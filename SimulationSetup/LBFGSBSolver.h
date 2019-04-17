#ifndef LBFGSB_H
#define LBFGSB_H
#include <iostream>
#include <vector>
#include <set>
#include <Eigen/LU>
#include <Eigen/Sparse>

#include "LineSearch.h"
#include "../MeshConnectivity.h"

template<typename TProblem>
class LbfgsbSolver
{
public:
    using Scalar = double;
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VariableTVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    
public:
    void setOptions(double eps, int maxItr, std::string path)
    {
        TOL = eps;
        MAXITR = maxItr;
        dataSavingPath = path;
    }
    

protected:
    // workspace matrices
    MatrixType W, M;
    double theta;
    int DIM;
    int m_historySize = 5;
    VariableTVector reductupperBound;
    VariableTVector reductlowerBound;
    
    double TOL;
    int MAXITR;
    std::string dataSavingPath;
    

    /**
     * @brief sort pairs (k,v) according v ascending
     * @details [long description]
     *
     * @param v [description]
     * @return [description]
     */
    std::vector<int> sort_indexes(const std::vector< std::pair<int, Scalar> > &v) {
        std::vector<int> idx(v.size());
        for (size_t i = 0; i != idx.size(); ++i)
            idx[i] = v[i].first;
        sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
            return v[i1].second < v[i2].second;
        });
        return idx;
    }
    /**
     * @brief Algorithm CP: Computation of the generalized Cauchy point
     * @details PAGE 8
     *
     * @param c [description]
     */
    void getGeneralizedCauchyPoint(const std::shared_ptr<TProblem> problem, const VariableTVector &x, const VariableTVector &g, VariableTVector &x_cauchy, VariableTVector &c) {
        const int DIM = x.rows();
        // Given x,l,u,g, and B = \theta I-WMW
        // {all t_i} = { (idx,value), ... }
        // TODO: use "std::set" ?
        std::vector<std::pair<int, Scalar> > SetOfT;
        // the feasible set is implicitly given by "SetOfT - {t_i==0}"
        VariableTVector d = -g;
        
        // n operations
        for (int j = 0; j < DIM; j++) {
            if (g(j) == 0) {
                SetOfT.push_back(std::make_pair(j, std::numeric_limits<double>::max()));
            } else {
                double tmp = 0;
                if (g(j) < 0) {
                    tmp = (x(j) - reductupperBound(j)) / g(j);
                } else {
                    tmp = (x(j) - reductlowerBound(j)) / g(j);
                }
                SetOfT.push_back(std::make_pair(j, tmp));
                if (tmp == 0) d(j) = 0;
            }
        }
        // sortedindices [1,0,2] means the minimal element is on the 1-st entry
        std::vector<int> sortedIndices = sort_indexes(SetOfT);
        x_cauchy = x;
        // Initialize
        // p :=     W^Scalar*p
        VariableTVector p = (W.transpose() * d);                     // (2mn operations)
        // c :=     0
        c = VariableTVector::Zero(W.cols());
        // f' :=    g^Scalar*d = -d^Td
        Scalar f_prime = -d.dot(d);                         // (n operations)
        // f'' :=   \theta*d^Scalar*d-d^Scalar*W*M*W^Scalar*d = -\theta*f' - p^Scalar*M*p
        Scalar f_doubleprime = (Scalar)(-1.0 * theta) * f_prime - p.dot(M * p); // (O(m^2) operations)
        f_doubleprime = std::max<Scalar>(std::numeric_limits<Scalar>::epsilon(), f_doubleprime);
        Scalar f_dp_orig = f_doubleprime;
        // \delta t_min :=  -f'/f''
        Scalar dt_min = -f_prime / f_doubleprime;
        // t_old :=     0
        Scalar t_old = 0;
        // b :=     argmin {t_i , t_i >0}
        int i = 0;
        for (int j = 0; j < DIM; j++) {
            i = j;
            if (SetOfT[sortedIndices[j]].second > 0)
                break;
        }
        int b = sortedIndices[i];
        // see below
        // t                    :=  min{t_i : i in F}
        Scalar t = SetOfT[b].second;
        // \delta Scalar             :=  t - 0
        Scalar dt = t ;
        // examination of subsequent segments
        while ((dt_min >= dt) && (i < DIM)) {
            if (d(b) > 0)
                x_cauchy(b) = reductupperBound(b);
            else if (d(b) < 0)
                x_cauchy(b) = reductlowerBound(b);
            // z_b = x_p^{cp} - x_b
            Scalar zb = x_cauchy(b) - x(b);
            // c   :=  c +\delta t*p
            c += dt * p;
            // cache
            VariableTVector wbt = W.row(b);
            f_prime += dt * f_doubleprime + (Scalar) g(b) * g(b) + (Scalar) theta * g(b) * zb - (Scalar) g(b) *
            wbt.transpose() * (M * c);
            f_doubleprime += (Scalar) - 1.0 * theta * g(b) * g(b)
            - (Scalar) 2.0 * (g(b) * (wbt.dot(M * p)))
            - (Scalar) g(b) * g(b) * wbt.transpose() * (M * wbt);
            f_doubleprime = std::max<Scalar>(std::numeric_limits<Scalar>::epsilon() * f_dp_orig, f_doubleprime);
            p += g(b) * wbt.transpose();
            d(b) = 0;
            dt_min = -f_prime / f_doubleprime;
            t_old = t;
            ++i;
            if (i < DIM) {
                b = sortedIndices[i];
                t = SetOfT[b].second;
                dt = t - t_old;
            }
        }
        dt_min = std::max<Scalar>(dt_min, (Scalar)0.0);
        t_old += dt_min;
#pragma omp parallel for
        for (int ii = i; ii < x_cauchy.rows(); ii++) {
            x_cauchy(sortedIndices[ii]) = x(sortedIndices[ii]) + t_old * d(sortedIndices[ii]);
        }
        c += dt_min * p;
    }
    /**
     * @brief find alpha* = max {a : a <= 1 and  l_i-xc_i <= a*d_i <= u_i-xc_i}
     * @details [long description]
     *
     * @param FreeVariables [description]
     * @return [description]
     */
    Scalar findAlpha(const std::shared_ptr<TProblem> problem, VariableTVector &x_cp, VariableTVector &du, std::vector<int> &FreeVariables) {
        Scalar alphastar = 1;
        const unsigned int n = FreeVariables.size();
        assert(du.rows() == n);
        for (unsigned int i = 0; i < n; i++) {
            if (du(i) > 0) {
                alphastar = std::min<Scalar>(alphastar, (reductupperBound(FreeVariables[i]) - x_cp(FreeVariables[i])) / du(i));
            } else {
                alphastar = std::min<Scalar>(alphastar, (reductlowerBound(FreeVariables[i]) - x_cp(FreeVariables[i])) / du(i));
            }
        }
        return alphastar;
    }
    /**
     * @brief solving unbounded probelm
     * @details [long description]
     *
     * @param SubspaceMin [description]
     */
    void SubspaceMinimization(const std::shared_ptr<TProblem> problem, Eigen::VectorXd &x_cauchy, Eigen::VectorXd &x, VariableTVector &c, Eigen::VectorXd &g,
                              VariableTVector &SubspaceMin) {
        Scalar theta_inverse = 1 / theta;
        std::vector<int> FreeVariablesIndex;
        
        for (int i = 0; i < x_cauchy.rows(); i++)
        {
            if ((x_cauchy(i) != reductupperBound(i)) && (x_cauchy(i) != reductlowerBound(i))) {
                FreeVariablesIndex.push_back(i);
            }
        }
        const int FreeVarCount = FreeVariablesIndex.size();
        MatrixType WZ = MatrixType::Zero(W.cols(), FreeVarCount);
        for (int i = 0; i < FreeVarCount; i++)
            WZ.col(i) = W.row(FreeVariablesIndex[i]);
        Eigen::VectorXd rr = (g + theta * (x_cauchy - x) - W * (M * c));
        // r=r(FreeVariables);
        MatrixType r = MatrixType::Zero(FreeVarCount, 1);
        for (int i = 0; i < FreeVarCount; i++)
            r.row(i) = rr.row(FreeVariablesIndex[i]);
        // STEP 2: "v = w^T*Z*r" and STEP 3: "v = M*v"
        VariableTVector v = M * (WZ * r);
        // STEP 4: N = 1/theta*W^T*Z*(W^T*Z)^T
        MatrixType N = theta_inverse * WZ * WZ.transpose();
        // N = I - MN
        N = MatrixType::Identity(N.rows(), N.rows()) - M * N;
        // STEP: 5
        // v = N^{-1}*v
        if (v.size() > 0)
            v = N.lu().solve(v);
        // STEP: 6
        // HERE IS A MISTAKE IN THE ORIGINAL PAPER!
        VariableTVector du = -theta_inverse * r - theta_inverse * theta_inverse * WZ.transpose() * v;
        // STEP: 7
        Scalar alpha_star = findAlpha(problem, x_cauchy, du, FreeVariablesIndex);
        // STEP: 8
        VariableTVector dStar = alpha_star * du;
        SubspaceMin = x_cauchy;
        for (int i = 0; i < FreeVarCount; i++) {
            SubspaceMin(FreeVariablesIndex[i]) = SubspaceMin(FreeVariablesIndex[i]) + dStar(i);
        }
    }
public:
    void setHistorySize(const int hs) { m_historySize = hs; }
  
    
    void minimize(std::shared_ptr<TProblem> problem, Eigen::VectorXd &x0, Eigen::MatrixXd &pos0)
    {
        if(!problem->isValid(x0))
            std::cerr << "start with invalid x0" << std::endl;
        reductupperBound = problem -> projM * problem->upperBound();
        reductlowerBound = problem -> projM * problem->lowerBound();
        DIM = x0.rows();
        theta = 1.0;
        W = MatrixType::Zero(DIM, 0);
        M = MatrixType::Zero(0, 0);
        MatrixType yHistory = MatrixType::Zero(DIM, 0);
        MatrixType sHistory = MatrixType::Zero(DIM, 0);
        
        
        Eigen::VectorXd fullx0, fullx;
        fullx0 = problem->getFullVariables(x0);
        Eigen::MatrixXd newPos = pos0;
        
        Eigen::VectorXd x = x0, g = x0;
        std::cout<<"Optimization begins!!"<<std::endl;
        problem->projectBack(x, newPos);
        
        Scalar f = problem->value(x, newPos);
        problem->gradient(x, newPos, g);
        
        double fmin = f;
        problem->save(x, newPos, dataSavingPath, true);
        // conv. crit.
        auto noConvergence =
        [&](Eigen::VectorXd &x, Eigen::VectorXd &g)->bool
        {
            return (((x - g).cwiseMax(reductlowerBound).cwiseMin(reductupperBound) - x).template lpNorm<Eigen::Infinity>() >= TOL);
        };
    
        int itr = 0;
        while (noConvergence(x, g) && itr < MAXITR)
        {
            Scalar f_old = f;
            Eigen::VectorXd x_old = x;
            Eigen::VectorXd g_old = g;
            // STEP 2: compute the cauchy point
            Eigen::VectorXd CauchyPoint = VariableTVector::Zero(DIM);
            VariableTVector c = VariableTVector::Zero(W.cols());
            getGeneralizedCauchyPoint(problem, x, g, CauchyPoint, c);
            std::cout<<"Compute the Cauchy point finished!"<<std::endl;
            // STEP 3: compute a search direction d_k by the primal method for the sub-problem
            Eigen::VectorXd SubspaceMin;
            SubspaceMinimization(problem, CauchyPoint, x, c, g, SubspaceMin);
             std::cout<<"compute a search direction d_k by the primal method for the sub-problem finished!"<<std::endl;
            // STEP 4: perform linesearch and STEP 5: compute gradient
            Scalar rate = 1.0;
            if(itr < m_historySize)
            {
                rate = std::min(x.norm() / sqrt(x.size()) * 1e-3, 1e-5/(SubspaceMin-x).norm() * x.norm());
            }
            bool is_success = lineSearch(problem, x, newPos, SubspaceMin-x,  rate);
            if(!is_success)
            {
                std::cout<<"L-BFGS-B failed, using GD instead!"<<std::endl;
                rate = std::min(x.norm() / sqrt(x.size()) * 1e-3, 1e-5/(SubspaceMin-x).norm() * x.norm());
                is_success = lineSearch(problem, x, newPos, -g, rate);
                if(!is_success)
                {
                    std::cout<<"Line Search failed to find a reasonable step"<<std::endl;
                    return;
                }
                x = x - rate * g;
            }
            // update current guess and function information
            else
            {
                x = x - rate*(x-SubspaceMin);
            }
            fullx = problem->getFullVariables(x);
            problem->projectBack(x, newPos);
            
            f = problem->value(x, newPos);
            problem->gradient(x, newPos, g);
            // prepare for next iteration
            Eigen::VectorXd newY = g - g_old;
            Eigen::VectorXd newS = x - x_old;
            // STEP 6:
            Scalar test = newS.dot(newY);
            test = (test < 0) ? -1.0 * test : test;
            if (test > TOL * newY.squaredNorm())
            {
                if (yHistory.cols() < m_historySize)
                {
                    yHistory.conservativeResize(DIM, yHistory.cols() + 1);
                    sHistory.conservativeResize(DIM, sHistory.cols() + 1);
                }
                else
                {
                    yHistory.leftCols(m_historySize - 1) = yHistory.rightCols(m_historySize - 1).eval();
                    sHistory.leftCols(m_historySize - 1) = sHistory.rightCols(m_historySize - 1).eval();
                }
                yHistory.rightCols(1) = newY;
                sHistory.rightCols(1) = newS;
                // STEP 7:
                theta = (Scalar)(newY.transpose() * newY) / (newY.transpose() * newS);
                W = MatrixType::Zero(yHistory.rows(), yHistory.cols() + sHistory.cols());
                W << yHistory, (theta * sHistory);
                MatrixType A = sHistory.transpose() * yHistory;
                MatrixType L = A.template triangularView<Eigen::StrictlyLower>();
                MatrixType MM(A.rows() + L.rows(), A.rows() + L.cols());
                MatrixType D = -1 * A.diagonal().asDiagonal();
                MM << D, L.transpose(), L, ((sHistory.transpose() * sHistory) * theta);
                M = MM.inverse();
            }
            
            Eigen::VectorXd fullg = problem->projM.transpose() * g;
            
            Eigen::VectorXd curL, curS;
            problem -> convertParams2LAndS(x, curL, curS);
            
            std::cout<<std::endl<< "iter: "<<itr<< ", \t Rate: "<<rate<< ", \t f = " <<  f<<", \t ||g||_inf = "<< g.lpNorm<Eigen::Infinity>()<<", \t ||Dir||_Inf = "<<(SubspaceMin-x_old).lpNorm<Eigen::Infinity>()<<std::endl;
            std::cout<<"||g_L||_inf = "<<fullg.segment(0, curL.rows()).lpNorm<Eigen::Infinity>()<<", \t ||g_s||_inf = "<<fullg.segment(curL.rows(), curS.rows()).lpNorm<Eigen::Infinity>()<<std::endl;
            std::cout<<"abar change "<<(fullx-fullx0).segment(0, curL.rows()).lpNorm<Eigen::Infinity>()<<", \t bbar changes: "<<(fullx-fullx0).segment(curL.rows(), curS.rows()).lpNorm<Eigen::Infinity>()<<std::endl;
            std::cout <<"Shape Difference = "<<problem->computeDifference(newPos)<<", \t Abar Penalty: "<<problem->computeAbarSmoothness(curL)<<", \t Bbar Penalty: "<<problem->computeBbarSmoothness(curL, curS, newPos)<<std::endl<<std::endl;
            itr++;
            if (fabs(f_old - f) < TOL)
            {
                // successive function values too similar
                break;
            }
            
            if(itr % 10 == 0)
            {
                double f = problem->value(x, newPos);
                if(f < fmin)
                {
                    problem->save(x, newPos, dataSavingPath, false);
                    fmin = f;
                }
            }
            
        }
        x0 = x;
        fullx0 = fullx;
        pos0 = newPos;
        f = problem->value(x, pos0);
        if(f < fmin)
        {
            problem->save(x, newPos, dataSavingPath, false);
        }
    }
};

#endif
