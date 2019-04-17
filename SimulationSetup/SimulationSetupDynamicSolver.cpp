#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include <igl/edge_lengths.h>
#include <iomanip>
#include <fstream>
#include <Eigen/CholmodSupport>

#include "SimulationSetupDynamicSolver.h"
#include "../ElasticShell.h"
#include "LBFGSBSolver.h"

#ifndef M_TOL
#define M_TOL 1e-8
#endif

#ifndef M_MIN
#define M_MIN 1e-15
#endif

#ifndef M_MAX
#define M_MAX 1e15
#endif

#ifndef MAX_ITR
#define MAX_ITR 1e4
#endif

bool SimulationSetupDynamicSolver::loadParams()
{
//        testValueAndGradient();
//        return true;
    //    testLBFGS();
    //    return true;
    auto op = std::shared_ptr<SensitiveAnalysis>();
    
    if(selectedDynamicType == "ABbarPos")
    {
        op = std::make_shared<SensitiveAnalysisABbarPos>();
    }
    
    std::ifstream infile(abarPath);
    //    std::ifstream infile("/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build/release/error.dat");
    if(!infile)
        return false;
    infile >> thickness;
    infile >> abarCoef;
    infile >> bbarCoef;
    infile >> smoothCoef;
    
    op->initialization(initialPos, targetPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
    int vertNum, paramsNum, faceNum;
    infile >> vertNum >> paramsNum;
    Eigen::VectorXd x(vertNum), curL, curS;
    
    if(selectedDynamicType == "ABbarPos")
    {
        faceNum = paramsNum / 4;
        curL.resize(3 * faceNum);
        curS.resize(faceNum);
    }
    double d;
    x.setZero();
    for(int i=0;i<vertNum; i++)
    {
        infile >> d;
        x(i) = d;
    }
    double d1, d2, d3;
    for(int i = 0; i < abars.size(); i++)
    {
        infile >> d1 >> d2 >> d3;
        Eigen::Matrix2d L;
        L << d1,0,
        d2,d3;
        abars[i] = L*L.transpose();
        
        curL(3*i) = d1;
        curL(3*i + 1) = d2;
        curL(3*i + 2) = d3;
        
    }
    
    if(selectedDynamicType == "ABbarPos")
    {
        for(int i = 0; i < bbars.size(); i++)
        {
            infile >> d1;
            bbars[i] = d1 * abars[i];
            
            curS(i) = d1 * op->_areaList(i);
            
        }
    }
    
    int nverts = vertNum / 3;
    Eigen::MatrixXd curPos(nverts,3);
    Eigen::VectorXd derivative;
    
    for(int i=0; i<nverts; i++)
    {
        curPos.row(i) = x.segment<3>(3*i);
    }
    
    Eigen::MatrixXd L(mesh.nFaces(), 3);
    igl::edge_lengths(curPos, mesh.faces(), L);
    
    
  
    
    
    Eigen::VectorXd sol(curL.size() + curS.size());
    sol.segment(0, curL.size()) = curL;
    sol.segment(curL.size(), curS.size()) = curS;
    auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisABbarPos>(op);
    if(op1->isValid(sol))
    {
        std::cout<<"Solution is feasible"<<std::endl;
    }
    else
    {
        
        std::cout<<"Solution is not feasible"<<std::endl;
        std::cout<<(sol - op1->lowerBound()).minCoeff()<<" "<<(sol - op1->upperBound()).maxCoeff()<<std::endl;
    }
    
    double differenceValue = op->computeDifference(curPos);
    double smoothnessAbarValue = op->computeAbarSmoothness(curL);
    //    double smoothnessPositionValue = op->computePositionSmoothness(curPos);
    Eigen::SparseMatrix<double> projM;
    
    Eigen::VectorXd params(curL.size() + curS.size());
    params.segment(0, curL.size()) = curL;
    params.segment(curL.size(), curS.size()) = curS;
    
    op->projM.resize(paramsNum, paramsNum);
    projM.setIdentity();
    
    Eigen::VectorXd grad;
    op->gradient(params, curPos, grad);
    std::cout<<grad.norm()<<std::endl;
    
    std::cout<<"Abar Penalty: "<<smoothnessAbarValue<<" Penalty Coefficient: "<<abarCoef<<std::endl;
    if(selectedDynamicType == "ABbarPos")
    {
        auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisABbarPos>(op);
        double smoothnessBbarValue = op1->computeBbarSmoothness(curL, curS, curPos);
        std::cout<<"Bbar Penalty: "<<smoothnessBbarValue<<" Penalty Coefficient: "<<bbarCoef<<std::endl;
    }
    //    std::cout<<"Delta q Penalty: "<<smoothnessPositionValue<<" Smootheness Coefficient: "<<smoothCoef<<std::endl;
    std::cout<<"Shape Changed Value: "<<differenceValue<<" Thickness: "<<thickness<<std::endl;
    
    MidedgeAverageFormulation sff;
    
    
    projM.resize(3*nverts, 3*nverts);
    projM.setIdentity();
    
    std::vector<Eigen::Triplet<double> > hessianTriplet;
    
    double energy = elasticEnergy(mesh, curPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, &hessianTriplet);
    
    Eigen::SparseMatrix<double> H(3*nverts, 3*nverts), resH;
    H.setFromTriplets(hessianTriplet.begin(), hessianTriplet.end());
    
    resH = projM * H * projM.transpose();
    
    Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > solver;
    
    Eigen::SparseMatrix<double> I(resH.rows(), resH.cols());
    I.setIdentity();
    
    resH += 1e-8 * I;
    
    solver.compute(resH);
    
    if(solver.info() == Eigen::ComputationInfo::Success )
    {
        std::cout<<"Rearch the local minimun at current state"<<std::endl;
    }
    else
    {
        std::cout<<"Saddle point"<<std::endl;
    }
    
    std::cout<<"Energy with current positions: "<<energy<<" Derivative with current positions: "<<(projM * derivative).norm()<<std::endl;
    
    energy = elasticEnergy(mesh, targetPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, NULL);
    std::cout<<"Energy with target positions: "<<energy<<" Derivative with target positions: "<<(projM * derivative).norm()<<std::endl;
    targetPosAfterFirstStep = curPos;
    
    std::cout<<curS.norm()<<std::endl;
    std::cout<<curL.norm()<<std::endl;
    
    return true;
}

void SimulationSetupDynamicSolver::buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
    int nfaces = mesh.nFaces();
    abars.resize(nfaces);
    bbars.resize(nfaces);
    
    initialAbars.resize(nfaces);
    
    lameAlpha = YoungsModulus * PoissonsRatio / (1.0 - PoissonsRatio * PoissonsRatio);
    lameBeta = YoungsModulus / 2.0 / (1.0 + PoissonsRatio);
    
    for (int i = 0; i < nfaces; i++)
    {
        bbars[i].setZero();
        initialAbars[i] = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
    }
    if(!loadParams() || _is_overwrite || _is_continue)
    {
        findFirstFundamentalForms(sff);
        bool flag = loadParams();
        if(flag)
            std::cout<<"Finished building rest fundamental forms!"<<std::endl;
    }
    
}



void SimulationSetupDynamicSolver::findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
    ////////////////////////////// Initialization ////////////////////////////////////////////
    auto op = std::shared_ptr<SensitiveAnalysis>();
    
    if(selectedDynamicType == "ABbarPos")
    {
        op = std::make_shared<SensitiveAnalysisABbarPos>();
    }
    op->initialization(initialPos, targetPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
    op->setPenalty(abarCoef, bbarCoef, smoothCoef);
    
    bool is_exist = false;
    int nfaces = mesh.nFaces();
    Eigen::VectorXd L(3*nfaces);
    Eigen::VectorXd S(nfaces);
    Eigen::MatrixXd pos;
    
    if(_is_continue)
    {
        is_exist = loadParams();
    }
    
    if(!is_exist)
    {
        std::cout<<"1"<<std::endl;
        for(int i=0;i<nfaces;i++)
        {
            //        Eigen::Matrix2d a = firstFundamentalForm(mesh, ellipsoid, i, NULL, NULL);
            Eigen::Matrix2d a = firstFundamentalForm(mesh, targetPos, i, NULL, NULL);
            L(3*i) = sqrt(a(0,0));
            L(3*i+1) = a(0,1)/L(3*i);
            L(3*i+2) = sqrt(a.determinant())/L(3*i);
            

    
       
          if(selectedDynamicType == "ABbarPos")
          {
              Eigen::VectorXd extraDOFs(0);
              Eigen::Matrix2d b = sff.secondFundamentalForm(mesh, targetPos, extraDOFs, i, NULL, NULL);
              S(i) = 0.5 * (a.inverse() * b).trace()  * op->_areaList(i);
              
            }
        }
        Eigen::VectorXd params(L.size() + S.size());
        params.segment(0, L.size()) = L;
        params.segment(L.size(), S.size()) = S;
        std::cout<<"2"<<std::endl;
        pos = targetPos;
        op->save(params, pos, abarPath, false);
    }
    else
    {
        std::cout<<"2"<<std::endl;
        for(int i=0;i<nfaces;i++)
        {
            //        Eigen::Matrix2d a = firstFundamentalForm(mesh, ellipsoid, i, NULL, NULL);
            Eigen::Matrix2d a = abars[i];
            L(3*i) = sqrt(a(0,0));
            L(3*i+1) = a(0,1)/L(3*i);
            L(3*i+2) = sqrt(a.determinant())/L(3*i);
            
            if(selectedDynamicType == "ABbarPos")
            {
                S.resize(nfaces);
                Eigen::Matrix2d b = bbars[i];
                S(i) = 0.5 * (a.inverse() * b).trace() * op->_areaList(i);
            }
        }
        pos = targetPosAfterFirstStep;
    }
    
    // First Stage: cut off L to satisfy the fibrication constraints
    std::set<int> cutAbar;
    Eigen::VectorXd infBound(4 * nfaces);  // 3*F = L.size, F = S.size;
    infBound.setConstant(std::numeric_limits<double>::infinity());
    int itr = 0;
    Eigen::VectorXd old_L = L;
    while(cutoffL(L, cutAbar))
    {
        std::cout<<"Outer iteration: "<<itr<<std::endl;
        LbfgsbSolver<SensitiveAnalysisABbarPos> solver;
        solver.setOptions(M_TOL, int (0.1 * MAX_ITR), abarPath);
        Eigen::VectorXd x(4*nfaces);
        x.segment(0, L.size()) = L;
        x.segment(L.size(), S.size()) = S;
        auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisABbarPos>(op);
        op1->setProjM(cutAbar);
        op1->updateFixedVariables(x);
        op1->setLowerBound(-infBound);
        op1->setUpperBound(infBound);
        for(int i = 1; i< 10; i++)
        {
            double rate =  0.1 * i;
            Eigen::VectorXd mid_L = (1 - rate) * old_L + rate * L;
            Eigen::VectorXd mid_x(4*nfaces);
            mid_x.segment(0, L.size()) = L;
            mid_x.segment(L.size(), S.size()) = S;
            
            op1->projectBack(op1->projM * mid_x, pos);
        }
        
        op1->_lambdaAbar = 1e-2;
        
        Eigen::VectorXd reductX = op1->projM * x;
        std::cout<<"L-BGFS solver"<<std::endl;
        solver.minimize(op1, reductX, pos);
        x = op1->getFullVariables(reductX);
        L = x.segment(0, L.size());
        S = x.segment(L.size(), S.size());
        itr ++;
        old_L = L;
        if(itr == 10)
            break;
    }
    
    
    LbfgsbSolver<SensitiveAnalysisABbarPos> solver;
    solver.setOptions(M_TOL, MAX_ITR, abarPath);
    Eigen::VectorXd x(4*nfaces);
    auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisABbarPos>(op);
    op1->projM.resize(x.size(), x.size());
    op1->projM.setIdentity();
    op1->_lambdaAbar = 1;
    for(int i=0;i<nfaces;i++)
    {
        op1->_atarget[i] = L2abar(L.segment(3*i, 3));
        op1->_starget(i) = S(i);
    }
    op1->updateBoxConstraint();
    for(int i=0;i<nfaces;i++)
    {
        double len = op->m_upperBound(3*nfaces + i) - op->m_lowerBound(3*nfaces + i);
        if(S(i) < op->m_lowerBound(3*nfaces + i))
            S(i) = op->m_lowerBound(3*nfaces + i) + 1e-2 * len;
        if(S(i) > op->m_upperBound(3*nfaces + i))
            S(i) = op->m_upperBound(3* nfaces + i) - 1e-2 * len;
    }
    x.segment(0, L.size()) = L;
    x.segment(L.size(), S.size()) = S;
    
    solver.minimize(op1, x, pos);
    
    
}

bool SimulationSetupDynamicSolver::cutoffL(Eigen::VectorXd &curL, std::set<int> &cutAbar)
{
    int nfaces = curL.size() / 3;
    Eigen::VectorXd valueList(nfaces);
    std::vector<Eigen::Matrix2d> abarList(nfaces);
    for(int i=0;i<nfaces;i++)
    {
        Eigen::Matrix2d T, abar;
        abar = L2abar(curL.segment(3*i, 3));
        T.col(0) = ( initialPos.row(mesh.faceVertex(i, 1)) - initialPos.row(mesh.faceVertex(i, 0)) ).segment(0, 2);
        T.col(1) = ( initialPos.row(mesh.faceVertex(i, 2)) - initialPos.row(mesh.faceVertex(i, 0)) ).segment(0, 2);
        abar = (T.transpose()).inverse() * abar * T.inverse();
        abarList[i] = abar;
        
        valueList(i) = abar.trace() / 2.0;
    }
    
    double min = valueList.minCoeff();
    double max = valueList.maxCoeff();
    
    if(min >= 0.6 * max)
        return false;
    
    
    int intervalNum = 10;
    double interval = ( max - 5.0/3 * min ) / (intervalNum * 1.0);  // maximum shrinking rate is 0.7 => min / max >= 0.7*2 * 11.0 / 9.0 (11.0 / 9.0 is for bbar soften)
    int bestCount = 0;
    double bestDivision = max;
    for(int i = 0; i<intervalNum;i++)
    {
        int count = 0;
        double end = max - i * interval;
        double begin = 0.6 * end;
        for(int j=0;j<nfaces;j++)
        {
            if(valueList(j) >= begin && valueList(j) <= end)
            {
                count ++;
            }
        }
        if(count > bestCount)
        {
            bestCount = count;
            bestDivision = end;
        }
    }
   
    cutAbar.clear();
    for(int i=0;i<nfaces;i++)
    {
        if(valueList(i) < 0.6 * bestDivision)
        {
            valueList(i) = 0.6 * bestDivision;
            cutAbar.insert(i);
        }
        else if(valueList(i) > bestDivision)
        {
            valueList(i) = bestDivision;
            cutAbar.insert(i);
        }
        
    }
    
    for(auto it=cutAbar.begin();it!=cutAbar.end();it++)
    {
        int fid = *it;
        abarList[fid] = abarList[fid] / (abarList[fid].trace() * 0.5) * valueList(fid);
        
        Eigen::Matrix2d T;
        T.col(0) = ( initialPos.row(mesh.faceVertex(fid, 1)) - initialPos.row(mesh.faceVertex(fid, 0)) ).segment(0, 2);
        T.col(1) = ( initialPos.row(mesh.faceVertex(fid, 2)) - initialPos.row(mesh.faceVertex(fid, 0)) ).segment(0, 2);
        Eigen::Matrix2d abar = T.transpose() * abarList[fid] * T;
        
        Eigen::Vector3d L = abar2L(abar);
        
        if(curL(3*fid) >= 0)
        {
            curL.segment(3*fid, 3) = L;
        }
        else
        {
            curL.segment(3*fid, 3) = -L;
        }
    }
    return true;
}

void SimulationSetupDynamicSolver::testValueAndGradient()
{
    SensitiveAnalysisABbarPos op;
    std::map<int, double> clampedDOFs;
    clampedDOFs.clear();
    op.initialization(initialPos, targetPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
    op.setPenalty(abarCoef, bbarCoef, smoothCoef);
    
    op.test();
}


bool SimulationSetupDynamicSolver::testLineSearch(Eigen::VectorXd x, Eigen::VectorXd dir, double &rate)
{
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = M_MIN;
    double beta = M_MAX;
    
    Eigen::VectorXd grad;
    double orig = testFunc(x, &grad);
    double deriv = dir.dot(grad);
    
    while (true)
    {
        Eigen::VectorXd newdE;
        Eigen::VectorXd newX = x + rate*dir;
        double newenergy = testFunc(newX, &newdE);
        
        std::cout << "Trying rate = " << rate << ", energy now " << newenergy << std::endl;
        
        if (std::isnan(newenergy) || newenergy > orig + rate*deriv*c1)
        {
            //            std::cout<<"Voilate the first Wolfe Condition"<<std::endl;
            beta = rate;
            rate = 0.5*(alpha + beta);
            if (beta - alpha < 1e-15)
            {
                rate = 1e-15;
                std::cout<<"Line Search failed, finished with Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                //                pos = newPos;
                return false;
            }
        }
        else if (newdE.dot(dir) < c2*deriv)
        {
            //            std::cout<<"Voilate the second Wolfe Condition"<<std::endl;
            alpha = rate;
            if (beta == M_MAX)
            {
                rate = 2 * alpha;
            }
            else
            {
                rate = 0.5*(alpha + beta);
            }
            
            if (beta - alpha < 1e-10)
            {
                if(newenergy > orig)
                {
                    std::cout<<"Line Search failed with beta - alph < 1e-10 without decreasing energy, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                    return false;
                }
                else
                {
                    std::cout<<"Line Search succeed with beta - alph < 1e-10 without sufficient decrease, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                    return false;
                }
            }
        }
        else
        {
            std::cout<<"Line Search Finished with Rate = "<<rate<<std::endl;
            return true;
        }
    }
}




void SimulationSetupDynamicSolver::testProjectBackSim()
{
//    SensitiveAnalysisABbarPos op;
//    op.initialization(initialPos, targetPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
//    op.setPenalty(abarCoef, bbarCoef, smoothCoef);
//
//
//    Eigen::VectorXd initialL(3*mesh.nFaces()), tarL(3*mesh.nFaces());
//    std::vector<Eigen::Matrix2d> as(mesh.nFaces());
//    for(int i =0; i<mesh.nFaces();i++)
//    {
//        Eigen::Matrix2d a = firstFundamentalForm(mesh, targetPos, i, NULL, NULL);
//        as[i] = a;
//        tarL(3*i) = sqrt(a(0,0));
//        tarL(3*i+1) = a(0,1) / tarL(3*i);
//        tarL(3*i+2) = sqrt(a.determinant()) / tarL(3*i);
//    }
//
//    Eigen::MatrixXd pos = targetPos;
//    Eigen::MatrixXd epsPos = Eigen::MatrixXd::Random(pos.rows(), pos.cols());
//    epsPos.row(0).setZero();
//    epsPos.row(1).setZero();
//    epsPos.row(2).setZero();
//    epsPos.row(3).setZero();
//    epsPos.normalized();
//
//    double eps = 1e-8;
//    MidedgeAverageFormulation sff;
//    Eigen::SparseMatrix<double> H0(3*targetPos.rows(), 3*targetPos.rows()), H1(3*targetPos.rows(), 3*targetPos.rows());
//    std::vector<Eigen::Triplet<double> > hessian;
//    Eigen::VectorXd edgeEOFs(0);
//    double  energy = 0;
//
//    for(int i=10;i<=10;i++)
//    {
//        double a = 0.1 * i;
//        pos = a * initialPos + (1-a) * targetPos;
//        Eigen::MatrixXd pos0 = pos;
//        Eigen::VectorXd S(0);
//        op.projectBack(tarL, S, pos0);
//
//
//        energy = elasticEnergy(mesh, pos0, edgeEOFs, lameAlpha, lameBeta, thickness, as, bbars, sff, NULL, &hessian);
//        H0.setFromTriplets(hessian.begin(), hessian.end());
//        Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
//        solver.compute(H0);
//        if(solver.info() == Eigen::ComputationInfo::Success)
//        {
//            std::cout<<"Local minimal"<<std::endl;
//        }
//        else
//        {
//            std::cout<<"saddle point"<<std::endl;
//        }
//
//        igl::writeOBJ("simed_saddle.obj", pos0, mesh.faces());
//        Eigen::MatrixXd pos1 = pos + eps * epsPos;
//        op.projectBack(tarL, S, pos1);
//
//        hessian.clear();
//        energy = elasticEnergy(mesh, pos1, edgeEOFs, lameAlpha, lameBeta, thickness, as, bbars, sff, NULL, &hessian);
//        H1.setFromTriplets(hessian.begin(), hessian.end());
//
//        solver.compute(H1);
//        if(solver.info() == Eigen::ComputationInfo::Success)
//        {
//            std::cout<<"Local minimal"<<std::endl;
//        }
//        else
//        {
//            std::cout<<"saddle point"<<std::endl;
//        }
//
//        igl::writeOBJ("simed_perturbed_saddle.obj", pos1, mesh.faces());
//        std::cout<<std::endl<<"The difference between two optimal solution is: "<<std::endl;
//        std::cout<<(pos1 - pos0).norm()<<std::endl<<std::endl;
//    }
    
}

double SimulationSetupDynamicSolver::testFunc(Eigen::VectorXd x, Eigen::VectorXd *grad)
{
    double fx = 0.0;
    Eigen::VectorXd df(x.size());
    for(int i = 0; i < x.size(); i += 2)
    {
        double t1 = 1.0 - x(i);
        double t2 = 10 * (x(i + 1) - x(i) * x(i));
        df(i + 1) = 20 * t2;
        df(i)     = -2.0 * (x(i) * df(i + 1) + t1);
        fx += t1 * t1 + t2 * t2;
    }
    *grad = df;
    return fx;
}

void SimulationSetupDynamicSolver::testLBFGS()
{
    int DIM = 20;
    int m = 10;
    Eigen::VectorXd grad(DIM), q(DIM), grad_old(DIM), s(DIM), y(DIM), z(DIM), x(DIM), x_old(DIM);
    Eigen::MatrixXd sVector = Eigen::MatrixXd ::Zero(DIM, m);
    Eigen::MatrixXd yVector = Eigen::MatrixXd::Zero(DIM, m);
    Eigen::Matrix<double, Eigen::Dynamic, 1> alpha = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(m);
    
    x.setZero();
    double f = testFunc(x, &grad);
    
    x_old = x;
    grad_old = grad;
    int iter = 0;
    
    
    while(true)
    {
        
        ///////////////////////////////////////////// L-BFGS /////////////////////////////////////////////////////////////////
        double H0k = 1;
        q = -grad;
        const int k = std::min<double>(m, iter);
        // for i = k − 1, k − 2, . . . , k − m§
        for (int i = k - 1; i >= 0; i--)
        {
            if(k < m)
                break;
            // alpha_i <- rho_i*s_i^T*q
            double rho = 1.0 / (sVector.col(i)).dot(yVector.col(i));
            //            std::cout<<rho<<std::endl;
            alpha(i) = rho * sVector.col(i).dot(q);
            // q <- q - alpha_i*y_i
            q = q - alpha(i) * yVector.col(i);
        }
        std::cout<<"q norm: "<<q.norm()<<std::endl;
        // z <- H_k^0*q
        // update the scaling factor
        if(k >= m)
        {
            H0k = yVector.col(m-1).dot(sVector.col(m-1)) / yVector.col(m-1).dot(yVector.col(m-1));
        }
        std::cout<<"H0K: "<<H0k<<std::endl;
        z = H0k * q;
        //for i k − m, k − m + 1, . . . , k − 1
        for (int i = 0; i <= k-1; i++)
        {
            if(k < m)
                break;
            // beta <- rho_i * y_i^T * r
            double rho = 1.0 / (sVector.col(i)).dot(yVector.col(i));
            double beta = rho * yVector.col(i).dot(z);
            // z <- z + s_i * ( alpha_i - beta)
            z = z + sVector.col(i) * (alpha(i) - beta);
        }
        std::cout<<"z: "<<z.norm()<<std::endl;
        double rate = 1.0;
        if(k<m)
        {
            rate = 1.0;
            bool isSuccess = testLineSearch(x, -grad, rate);
            x = x - rate * grad;
            f = testFunc(x, &grad);
            
            y = grad - grad_old;
            s = x - x_old;
            grad_old = grad;
            x_old = x;
        }
        else
        {
            rate = 1.0;
            bool isSuccess = testLineSearch(x, z, rate);
            if(isSuccess)
            {
                std::cout<<"L-BFGS succeeded!"<<std::endl;
                x = x + rate * z;
                f = testFunc(x, &grad);
                
                y = grad - grad_old;
                s = x - x_old;
                grad_old = grad;
                x_old = x;
            }
            else
            {
                rate = 1.0;
                bool isSuccess = testLineSearch(x, -grad, rate);
                x = x - rate * grad;
                f = testFunc(x, &grad);
                
                y = grad - grad_old;
                s = x - x_old;
                grad_old = grad;
                x_old = x;
            }
        }
        // update the history
        if (iter < m)
        {
            sVector.col(iter) = s;
            yVector.col(iter) = y;
        }
        else
        {
            sVector.leftCols(m - 1) = sVector.rightCols(m - 1).eval();
            sVector.rightCols(1) = s;
            yVector.leftCols(m - 1) = yVector.rightCols(m - 1).eval();
            yVector.rightCols(1) = y;
        }
        std::cout<<"s.dot(y) = "<<s.dot(y)<<std::endl;
        
        std::cout<<std::endl<< "iter: "<<iter<< ", Rate: "<<rate<< ", f = " <<  f <<", ||g||_inf "<< grad.lpNorm<Eigen::Infinity>()<<std::endl;
        
        if(iter == MAX_ITR)
        {
            std::cout<<"Maximun iteration reached!"<<std::endl;
            break;
        }
        if(grad.norm() <= M_TOL * std::max<double>(static_cast<double>(1.0), x.norm()))
        {
            std::cout<<"Force norm is less than "<<M_TOL<<std::endl;
            break;
        }
        if(s.template lpNorm<Eigen::Infinity>() <= M_TOL)
        {
            std::cout<<"x update is less than "<<M_TOL<<std::endl;
            break;
        }
        if(y.template lpNorm<Eigen::Infinity>() <= M_TOL)
        {
            std::cout<<"gradient update is less than "<<M_TOL<<std::endl;
            break;
        }
        iter++;
    }
    std::cout<<x<<std::endl;
}
