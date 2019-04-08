#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>
#include <igl/edge_lengths.h>
#include <iomanip>
#include <fstream>
#include <Eigen/CholmodSupport>

#include "SimulationSetupDynamicSolver.h"
#include "../ElasticShell.h"

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

bool SimulationSetupDynamicSolver::loadAbars()
{
//    testValueAndGradient();
//    return true;
    
    auto op = std::shared_ptr<SensitiveAnalysis>();
    
    if(selectedDynamicType == "AbarBbar")
    {
        op = std::make_shared<SensitiveAnalysisAbarBbar>();
    }
    else if (selectedDynamicType == "AbarPos")
    {
        op = std::make_shared<SensitiveAnalysisAbarPos>();
    }
    else if(selectedDynamicType == "ABbarPos")
    {
        op = std::make_shared<SensitiveAnalysisABbarPos>();
    }
    
    std::ifstream infile(abarPath);
//    std::ifstream infile("/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build/release/error.dat");
    if(!infile)
        return false;
    infile >> thickness;
    infile >> abarCoef;
    if(selectedDynamicType != "AbarPos")
        infile >> bbarCoef;
    infile >> smoothCoef;
    int vertNum, faceNum;
    infile >> vertNum >> faceNum;
    Eigen::VectorXd x(vertNum), curL, curS;
    // vertNum = 3*nverts; faceNum = 3 * nfaces;
    if(selectedDynamicType == "AbarBbar")
    {
        faceNum = faceNum / 2;
        curL.resize(faceNum);
        curS.resize(faceNum);
    }
    else if(selectedDynamicType == "ABbarPos")
    {
        faceNum = faceNum * 3 / 4;
        curL.resize(faceNum);
        curS.resize(faceNum / 3);
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
    if(selectedDynamicType == "AbarBbar")
    {
        for(int i = 0; i < bbars.size(); i++)
        {
            infile >> d1 >> d2 >> d3;
            bbars[i] = d1 * abars[i];
            
            curS(3*i) = d1;
            curS(3*i + 1) = d2;
            curS(3*i + 2) = d3;
            
        }
    }
    else if(selectedDynamicType == "ABbarPos")
    {
        for(int i = 0; i < bbars.size(); i++)
        {
            infile >> d1;
            bbars[i] = d1 * abars[i];
            
            curS(i) = d1;
            
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
    
    
    op->initialization(initialPos, targetPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
    
    double differenceValue = op->computeDifference(curPos);
    double smoothnessAbarValue = op->computeAbarSmoothness(curL);
//    double smoothnessPositionValue = op->computePositionSmoothness(curPos);
    Eigen::SparseMatrix<double> projM;
    
    
    Eigen::VectorXd grad;
    op->gradient(curL, curS, curPos, grad);
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
    if(!loadAbars())
    {
        findFirstFundamentalForms(sff);
        bool flag = loadAbars();
        if(flag)
            std::cout<<"Finished building rest fundamental forms!"<<std::endl;
    }
    
}

bool SimulationSetupDynamicSolver::lineSearch(std::shared_ptr<SensitiveAnalysis> op, Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd &pos, Eigen::VectorXd dir, double &rate)
{
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = M_MIN;
    double beta = M_MAX;
    
    Eigen::VectorXd grad;
    double orig = op->value(L, S, pos);
    op->gradient(L, S, pos, grad);
    double deriv = dir.dot(grad);
    
    while (true)
    {
        Eigen::VectorXd newdE;
        Eigen::VectorXd newL = L + rate*dir.segment(0, L.size());
        Eigen::MatrixXd newPos = pos;
        Eigen::VectorXd newS = S + rate*dir.segment(L.size(), S.size());
        op->projectBack(newL, newS, newPos);
        double newenergy = op->value(newL, newS, newPos);
        op->gradient(newL, newS, newPos, newdE);
        
        std::cout << "Trying rate = " << rate << ", energy now " << newenergy<<", L update "<<(newL-L).norm()<<", S update "<<(newS - S).norm() << std::endl;
        
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
                    std::cout<<"Line Search failed with beta - alph < 1e-10, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                    return false;
                }
                else
                {
                    std::cout<<"Line Search succeed with beta - alph < 1e-10, Rate = "<<rate<<" alpha: "<<alpha<<" beta: "<<beta<<std::endl;
                    return true;
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

void SimulationSetupDynamicSolver::findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
//    using namespace cppoptlib;
    auto op = std::shared_ptr<SensitiveAnalysis>();
//    srand((unsigned)time(NULL));
//    for(int i=0;i<initialPos.rows();i++)
//    {
//        initialPos(i,2) = (1e-6*rand())/RAND_MAX;
//    }
    
    if(selectedDynamicType == "AbarBbar")
    {
        op = std::make_shared<SensitiveAnalysisAbarBbar>();
    }
    else if (selectedDynamicType == "AbarPos")
    {
        op = std::make_shared<SensitiveAnalysisAbarPos>();
    }
    else if(selectedDynamicType == "ABbarPos")
    {
        op = std::make_shared<SensitiveAnalysisABbarPos>();
    }
    
    
    op->initialization(initialPos, targetPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
    op->setPenalty(abarCoef, bbarCoef, smoothCoef);
    
    Eigen::MatrixXd ellipsoid;
    Eigen::MatrixXi F;
//    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/shellbenchmarks/benchmarks/TestModels/middleCoarse/ellipsoid/ellipsoid_geometry.obj", ellipsoid, F);
    int nfaces = mesh.nFaces();
    Eigen::VectorXd L(3*nfaces);
    Eigen::VectorXd S(0);
    for(int i=0;i<nfaces;i++)
    {
//        Eigen::Matrix2d a = firstFundamentalForm(mesh, ellipsoid, i, NULL, NULL);
        Eigen::Matrix2d a = firstFundamentalForm(mesh, targetPos, i, NULL, NULL);
        L(3*i) = sqrt(a(0,0));
        L(3*i+1) = a(0,1)/L(3*i);
        L(3*i+2) = sqrt(a.determinant())/L(3*i);
        
        if(selectedDynamicType == "AbarBbar")
        {
            S.resize(3*nfaces);
            Eigen::VectorXd extraDOFs(0);
            Eigen::Matrix2d b = sff.secondFundamentalForm(mesh, targetPos, extraDOFs, i, NULL, NULL);
            auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisAbarBbar>(op);
            S.segment(3*i, 3) = op1->convertBar2s(a, b);
        }
        else if(selectedDynamicType == "ABbarPos")
        {
            S.resize(nfaces);
            Eigen::VectorXd extraDOFs(0);
            Eigen::Matrix2d b = sff.secondFundamentalForm(mesh, targetPos, extraDOFs, i, NULL, NULL);
            auto op1 = std::dynamic_pointer_cast<SensitiveAnalysisABbarPos>(op);
            S(i) = 0.5 * (a.inverse() * b).trace();
        }
    }
    if(selectedDynamicType == "AbarBbar")
    {
        double result = 0;
        for(int i=0;i<nfaces;i++)
        {
            std::cout<<S(3*i+1)<<std::endl;
            std::cout<<S(3*i+2)<<std::endl;
            result += S(3*i+1) * S(3*i+1) + S(3*i+2) * S(3*i+2);
    }
    std::cout<<result<<std::endl;
    }
    std::cout<<S.size()<<std::endl;
    Eigen::MatrixXd pos = targetPos;
    op->projectBack(L, S, pos);
    igl::writeOBJ("test.obj", pos, mesh.faces());
    
    int m = 10;
    int DIM = L.rows() + S.rows();
    Eigen::VectorXd grad(DIM), q(DIM), grad_old(DIM), s(DIM), y(DIM), z(DIM);
    Eigen::MatrixXd sVector = Eigen::MatrixXd ::Zero(DIM, m);
    Eigen::MatrixXd yVector = Eigen::MatrixXd::Zero(DIM, m);
    Eigen::Matrix<double, Eigen::Dynamic, 1> alpha = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(m);
    op->gradient(L, S, pos, grad);
    Eigen::VectorXd L_old = L;
    Eigen::VectorXd S_old = S;
    grad_old = grad;
    double fmin = op->value(L, S, pos);
    
    int iter = 0;
    double gradNorm = 0;
    std::cout<<"Simulation start (L-BFGS)!! Initial function value is: "<<fmin<<std::endl;
    while(true)
    {
        
//        const double relativeEpsilon = static_cast<double>(1e-6) * std::max<double>(static_cast<double>(1.0), L.norm());
//
//        if (grad.norm() < relativeEpsilon)
//            break;
//
///////////////////////////////////////////// L-BFGS /////////////////////////////////////////////////////////////////
        double H0k = 1;
        q = grad;
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
        std::cout<<"2: "<<q.norm()<<std::endl;
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
            // r <- r + s_i * ( alpha_i - beta)
            z = z + sVector.col(i) * (alpha(i) - beta);
        }
        std::cout<<"z: "<<z.norm()<<std::endl;
        double rate = std::min(L.norm() / sqrt(L.size()) * 1e-2, 1.0/grad.norm());
        if(k<m)
        {
            rate = std::min(L.norm() / sqrt(L.size()) * 1e-2, 1.0/grad.norm());
            bool isSuccess = lineSearch(op, L, S, pos, -grad, rate);
            L = L - rate * grad.segment(0,L.size());
            S = S - rate * grad.segment(L.size(), S.size());
            op->projectBack(L, S, pos);
            op->gradient(L, S, pos, grad);
            
            y = grad - grad_old;
            s.segment(0, L.size()) = L - L_old;
            s.segment(L.size(), S.size()) = S - S_old;
            grad_old = grad;
            L_old = L;
            S_old = S;
        }
        else
        {
            bool isSuccess = lineSearch(op, L, S, pos, -z, rate);
            if(isSuccess)
            {
                std::cout<<"L-BFGS succeeded!"<<std::endl;
                L = L - rate * z.segment(0,L.size());
                S = S - rate * z.segment(L.size(), S.size());
                op->projectBack(L, S, pos);
                op->gradient(L, S, pos, grad);
                
                y = grad - grad_old;
                s.segment(0, L.size()) = L - L_old;
                s.segment(L.size(), S.size()) = S - S_old;
                grad_old = grad;
                L_old = L;
                S_old = S;
            }
            else
            {
                std::cout<<"L-BFGS failed, using GD instead!!"<<std::endl;
                rate = std::min(L.norm() / sqrt(L.size()) * 1e-2, 1.0/grad.norm());
                bool isSuccess = lineSearch(op, L, S, pos, -grad, rate);
                L = L - rate * grad.segment(0,L.size());
                S = S - rate * grad.segment(L.size(), S.size());
                op->projectBack(L, S, pos);
                op->gradient(L, S, pos, grad);
                
                y = grad - grad_old;
                s.segment(0, L.size()) = L - L_old;
                s.segment(L.size(), S.size()) = S - S_old;
                grad_old = grad;
                L_old = L;
                S_old = S;
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
        
///////////////////////////////////////////// G-N /////////////////////////////////////////////////////////////////
//        double rate = std::min(L.norm() / sqrt(L.size()) * 1e-2, 1.0/grad.norm());
//        Eigen::SparseMatrix<double> H = (grad * grad.transpose()).sparseView();
//        Eigen::SparseMatrix<double> I(H.rows(), H.cols());
//        I.setIdentity();
//        H += 1e-6 * I;
//        Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
//        solver.compute(H);
//        Eigen::VectorXd dir = solver.solve(-grad);
//        bool isSuccess = lineSearch(op, L, pos, dir, rate);
//
//        L = L + rate * dir;
//        grad_old = grad;
//        projectBackOp(sff, L, pos);
//        op.gradient(L, pos, grad);
//
        
/////////////////////////////////////////////// GD /////////////////////////////////////////////////////////////////
//        // find steplength
//        double rate = std::min(L.norm() / sqrt(L.size()) * 1e-2, 1.0/grad.norm());
//        bool isSuccess = lineSearch(op, L, pos, -grad, rate);
//
//        L = L - rate * grad;
//        grad_old = grad;
//        projectBackOp(sff, L, pos);
//
//        op.gradient(L, pos, grad);
//
//        L_old = L;
        
        std::cout << "iter: "<<iter<< ", f = " <<  op->value(L, S, pos) << ", ||g||_inf "<<grad.template lpNorm<Eigen::Infinity>()<<", abar change "<<s.segment(0, L.size()).lpNorm<Eigen::Infinity>()<<", bbar changes: "<<s.segment(L.size(), S.size()).lpNorm<Eigen::Infinity>()<<", gradient change "<<y.lpNorm<Eigen::Infinity>()<< ", Rate: "<<rate<<std::endl<<std::endl;
        
        iter++;
        gradNorm = grad.template lpNorm<Eigen::Infinity>();
        if(isnan(gradNorm) || isnan(s.dot(y)))
        {
            std::cout<<"Something wrong happened!!"<<std::endl;
            std::cout<<s.norm()<<std::endl;
            std::cout<<y.norm()<<std::endl;
            std::ofstream outfile("error.dat", std::ios::trunc);
            int nverts = targetPos.rows();
            int nfaces = mesh.nFaces();
            
            outfile<<thickness<<"\n";
            outfile<<abarCoef<<"\n";
            if(selectedDynamicType != "AbarPos")
                 outfile<<bbarCoef<<"\n";
            outfile<<smoothCoef<<"\n";
            outfile<<3*nverts<<"\n";
            outfile<<3*nfaces<<"\n";
            
            std::cout<<3*nverts + 3*nfaces<<std::endl;
            
            for(int i=0;i<nverts;i++)
            {
                outfile<<std::setprecision(16)<<pos(i, 0)<<"\n";
                outfile<<std::setprecision(16)<<pos(i, 1)<<"\n";
                outfile<<std::setprecision(16)<<pos(i, 2)<<"\n";
            }
            
            for(int i=0;i<3*nfaces;i++)
            {
                outfile<<std::setprecision(16)<<L(i)<<"\n";
            }
            outfile<<std::setprecision(16)<<L(L.size()-1);
            outfile.close();
            igl::writeOBJ("resampled_error.obj", pos, mesh.faces());
        }
        if(iter == MAX_ITR)
        {
            std::cout<<"Maximun iteration reached!"<<std::endl;
            break;
        }
        if(gradNorm <= M_TOL * std::max<double>(static_cast<double>(1.0), L.norm()))
        {
            std::cout<<"Force norm is less than "<<M_TOL<<std::endl;
            break;
        }
        if(s.template lpNorm<Eigen::Infinity>() <= M_TOL)
        {
            std::cout<<"Abar update is less than "<<M_TOL<<std::endl;
            break;
        }
        if(y.template lpNorm<Eigen::Infinity>() <= M_TOL)
        {
            std::cout<<"gradient update is less than "<<M_TOL<<std::endl;
            break;
        }
        if(iter % 10 == 0)
        {
            double f = op->value(L, S, pos);
            if(f < fmin)
            {
                saveAbars(L, S, pos);
                fmin = f;
            }
        }
    }
    
    double f = op->value(L, S, pos);
    if(f < fmin)
    {
        saveAbars(L, S, pos);
        fmin = f;
    }
    
}

void SimulationSetupDynamicSolver::saveAbars(Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd pos)
{
    std::cout<<"Saving Abar path: "<<abarPath<<std::endl;
    std::ofstream outfile(abarPath, std::ios::trunc);
    int nverts = targetPos.rows();
    int nfaces = mesh.nFaces();
    
    outfile<<thickness<<"\n";
    outfile<<abarCoef<<"\n";
    if(selectedDynamicType != "AbarPos")
        outfile<<bbarCoef<<"\n";
    outfile<<smoothCoef<<"\n";
    outfile<<3*nverts<<"\n";
    
    int numEOFs = L.size() + S.size();
    outfile<<numEOFs<<"\n";
    
    std::cout<<3*nverts + 3*nfaces<<std::endl;
    
    for(int i=0;i<nverts;i++)
    {
        outfile<<std::setprecision(16)<<pos(i, 0)<<"\n";
        outfile<<std::setprecision(16)<<pos(i, 1)<<"\n";
        outfile<<std::setprecision(16)<<pos(i, 2)<<"\n";
    }
    
    for(int i=0;i<3*nfaces;i++)
    {
        outfile<<std::setprecision(16)<<L(i)<<"\n";
    }
    if(selectedDynamicType != "AbarPos")
    {
        for(int i=0;i<S.size();i++)
        {
            outfile<<std::setprecision(16)<<S(i)<<"\n";
        }
    }
    outfile.close();
    
    int startIdx, endIdx, expCoef;
    std::string subString = "";
    std::string resampledPath = abarPath;
    
    startIdx = resampledPath.rfind("/");
    endIdx = resampledPath.find("_");
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx + 1,resampledPath.begin() + endIdx, "resampled");
    
    // thickness
    if(thickness == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(thickness));
    startIdx = resampledPath.rfind("T");
    endIdx = resampledPath.rfind("A");
    subString = "";
    if(thickness > 0)
        subString = "T_1e" + std::to_string(expCoef);
    else
        subString = "T_0";
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    //Abar penalty
    if(abarCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(abarCoef));
    
    startIdx = resampledPath.rfind("A");
    endIdx = resampledPath.rfind("B");
    subString = "";
    if(abarCoef > 0)
        subString = "A_1e" + std::to_string(expCoef);
    else
        subString = "A_0";
    resampledPath= resampledPath .replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    // Bbar penalty
    if(bbarCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(bbarCoef));
    startIdx = resampledPath.rfind("B");
    endIdx = resampledPath.rfind("S");
    subString = "";
    if(bbarCoef > 0)
        subString = "B_1e" + std::to_string(expCoef);
    else
        subString = "B_0";
    resampledPath= resampledPath .replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    
    // smoothness
    if(smoothCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(smoothCoef));
    
    startIdx = resampledPath.rfind("S");
    endIdx = resampledPath.rfind(".");
    subString = "";
    if(smoothCoef > 0)
        subString = "S_1e" + std::to_string(expCoef);
    else
        subString = "S_0";
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx, subString);
    
    startIdx = resampledPath.rfind(".");
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.end(), ".obj");
    std::cout<<"Current abar loading path is: "<<resampledPath<<std::endl;
    igl::writeOBJ("resampled.obj", pos, mesh.faces());
    igl::writeOBJ(resampledPath, pos, mesh.faces());
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


void SimulationSetupDynamicSolver::testProjectBackSim()
{
    SensitiveAnalysisAbarPos op;
    op.initialization(initialPos, targetPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
    op.setPenalty(abarCoef, bbarCoef, smoothCoef);
    
   
    Eigen::VectorXd initialL(3*mesh.nFaces()), tarL(3*mesh.nFaces());
    std::vector<Eigen::Matrix2d> as(mesh.nFaces());
    for(int i =0; i<mesh.nFaces();i++)
    {
        Eigen::Matrix2d a = firstFundamentalForm(mesh, targetPos, i, NULL, NULL);
        as[i] = a;
        tarL(3*i) = sqrt(a(0,0));
        tarL(3*i+1) = a(0,1) / tarL(3*i);
        tarL(3*i+2) = sqrt(a.determinant()) / tarL(3*i);
    }
    
    Eigen::MatrixXd pos = targetPos;
    Eigen::MatrixXd epsPos = Eigen::MatrixXd::Random(pos.rows(), pos.cols());
    epsPos.row(0).setZero();
    epsPos.row(1).setZero();
    epsPos.row(2).setZero();
    epsPos.row(3).setZero();
    epsPos.normalized();
    
    double eps = 1e-8;
    MidedgeAverageFormulation sff;
    Eigen::SparseMatrix<double> H0(3*targetPos.rows(), 3*targetPos.rows()), H1(3*targetPos.rows(), 3*targetPos.rows());
    std::vector<Eigen::Triplet<double> > hessian;
    Eigen::VectorXd edgeEOFs(0);
    double  energy = 0;
    
    for(int i=10;i<=10;i++)
    {
        double a = 0.1 * i;
        pos = a * initialPos + (1-a) * targetPos;
        Eigen::MatrixXd pos0 = pos;
        Eigen::VectorXd S(0);
        op.projectBack(tarL, S, pos0);
        
        
        energy = elasticEnergy(mesh, pos0, edgeEOFs, lameAlpha, lameBeta, thickness, as, bbars, sff, NULL, &hessian);
        H0.setFromTriplets(hessian.begin(), hessian.end());
        Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        solver.compute(H0);
        if(solver.info() == Eigen::ComputationInfo::Success)
        {
            std::cout<<"Local minimal"<<std::endl;
        }
        else
        {
            std::cout<<"saddle point"<<std::endl;
        }
        
        igl::writeOBJ("simed_saddle.obj", pos0, mesh.faces());
        Eigen::MatrixXd pos1 = pos + eps * epsPos;
        op.projectBack(tarL, S, pos1);
        
        hessian.clear();
        energy = elasticEnergy(mesh, pos1, edgeEOFs, lameAlpha, lameBeta, thickness, as, bbars, sff, NULL, &hessian);
        H1.setFromTriplets(hessian.begin(), hessian.end());
        
        solver.compute(H1);
        if(solver.info() == Eigen::ComputationInfo::Success)
        {
            std::cout<<"Local minimal"<<std::endl;
        }
        else
        {
            std::cout<<"saddle point"<<std::endl;
        }
        
        igl::writeOBJ("simed_perturbed_saddle.obj", pos1, mesh.faces());
        std::cout<<std::endl<<"The difference between two optimal solution is: "<<std::endl;
        std::cout<<(pos1 - pos0).norm()<<std::endl<<std::endl;
    }
    
}
