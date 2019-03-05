#include <igl/writeOBJ.h>
#include <igl/edge_lengths.h>
#include <iomanip>
#include <Eigen/CholmodSupport>

#include "SimulationSetupDynamicSolver.h"
#include "../cppoptlib/solver/lbfgssolver.h"
#include "../cppoptlib/linesearch/morethuente.h"
#include "../ElasticShell.h"

#ifndef M_TOL
#define M_TOL 1e-8
#endif

#ifndef MAX_ITR
#define MAX_ITR 1e5
#endif

bool SimulationSetupDynamicSolver::loadAbars()
{
    std::ifstream infile(abarPath);
    if(!infile)
        return false;
    infile >> thickness;
    infile >> penaltyCoef;
    infile >> smoothCoef;
    int vertNum, faceNum;
    infile >> vertNum >> faceNum;
    
    double d;
    Eigen::VectorXd x(vertNum + faceNum);
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
        
        x(vertNum + 3*i) = d1;
        x(vertNum + 3*i + 1) = d2;
        x(vertNum + 3*i + 2) = d3;
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
    
    convertedProblem op;
    op.initialization(initialPos, targetPos, mesh, lameAlpha, lameBeta, thickness);
    
    double differenceValue = op.computeDifference(curPos);
    double smoothnessAbarValue = op.computeAbarSmoothness(x.segment(3*nverts, 3*mesh.nFaces()));
    double smoothnessPositionValue = op.computePositionSmoothness(curPos);
    
    std::cout<<"Abar Penalty: "<<smoothnessAbarValue<<" Penalty Coefficient: "<<penaltyCoef<<std::endl;
    std::cout<<"Delta q Penalty: "<<smoothnessPositionValue<<" Smootheness Coefficient: "<<smoothCoef<<std::endl;
    std::cout<<"Shape Changed Value: "<<differenceValue<<" Thickness: "<<thickness<<std::endl;
    
    MidedgeAverageFormulation sff;
    
    std::vector<Eigen::Triplet<double> > hessianTriplet;
    
    double energy = elasticEnergy(mesh, curPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, &hessianTriplet);
    
    Eigen::SparseMatrix<double> resH(3*nverts, 3*nverts);
    Eigen::SparseMatrix<double> H(3*nverts, 3*nverts);
    H.setFromTriplets(hessianTriplet.begin(), hessianTriplet.end());
    
    resH.setIdentity();
    resH = 1e-8 * resH + H;
    Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > solver(resH);
    
    if(solver.info() == Eigen::ComputationInfo::Success )
    {
        std::cout<<"Rearch the local minimun at current state"<<std::endl;
    }
    else
    {
        std::cout<<"Saddle point"<<std::endl;
    }
    
    std::cout<<"Energy with current positions: "<<energy<<" Derivative with current positions: "<<derivative.norm()<<std::endl;
    
    energy = elasticEnergy(mesh, targetPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, NULL);
    std::cout<<"Energy with target positions: "<<energy<<" Derivative with target positions: "<<derivative.norm()<<std::endl;
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
        findFirstFundamentalForms(sff);
}

double SimulationSetupDynamicSolver::lineSearch(convertedProblem op, Eigen::VectorXd L, Eigen::MatrixXd &pos, Eigen::VectorXd dir)
{
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 1e-8;
    double beta = 1e8;
    
    MidedgeAverageFormulation sff;
    Eigen::VectorXd grad;
    double orig = op.value(L, pos);
    double rate = std::min(L.norm() / sqrt(L.size()) * 1e-2, 1.0/dir.norm());
    op.gradient(L, pos, grad);
    double deriv = dir.dot(grad);
    
    while (true)
    {
        Eigen::VectorXd newdE;
        Eigen::VectorXd newL = L + rate*dir;
        Eigen::MatrixXd newPos = pos;
        projectBackOp(sff, newL, newPos);
        double newenergy = op.value(newL, newPos);
        op.gradient(newL, newPos, newdE);
        
//        std::cout << "Trying rate = " << rate << ", energy now " << newenergy << std::endl;
        
        if (std::isnan(newenergy) || newenergy > orig + rate*deriv*c1)
        {
//            std::cout<<"Voilate the first Wolfe Condition"<<std::endl;
            beta = rate;
            rate = 0.5*(alpha + beta);
            if (beta - alpha < 1e-8 || rate <= 1e-8)
            {
//                std::cout<<"Line Search failed, finished with Rate = "<<rate<<std::endl;
                pos = newPos;
                return rate;
            }
        }
        else if (newdE.dot(dir) < c2*deriv)
        {
//            std::cout<<"Voilate the second Wolfe Condition"<<std::endl;
            alpha = rate;
            if (beta == 1e8)
            {
                rate = 2 * alpha;
            }
            else
            {
                rate = 0.5*(alpha + beta);
            }
            
            if (beta - alpha < 1e-8 || rate <= 1e-8)
            {
//                std::cout<<"Line Search Finished with Rate = "<<rate<<std::endl;
                pos = newPos;
                return rate;
            }
        }
        else
        {
//            std::cout<<"Line Search Finished with Rate = "<<rate<<std::endl;
            pos = newPos;
            return rate;
        }
    }
}

void SimulationSetupDynamicSolver::findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
        testValueAndGradient();
        return;
    using namespace cppoptlib;
    convertedProblem op;
    srand((unsigned)time(NULL));
    for(int i=0;i<initialPos.rows();i++)
    {
        initialPos(i,2) = (1e-6*rand())/RAND_MAX;
    }
    op.initialization(initialPos, targetPos, mesh, lameAlpha, lameBeta, thickness);
    op.setPenalty(penaltyCoef, smoothCoef);
    
    int nfaces = mesh.nFaces();
    Eigen::VectorXd L(3*nfaces);
    
    for(int i=0;i<nfaces;i++)
    {
        Eigen::Matrix2d a = firstFundamentalForm(mesh, targetPos, i, NULL, NULL);
        L(3*i) = sqrt(a(0,0));
        L(3*i+1) = a(0,1)/L(3*i);
        L(3*i+2) = sqrt(a.determinant())/L(3*i);
    }
    Eigen::MatrixXd pos = targetPos;
    projectBackOp(sff, L, pos);
    igl::writeOBJ("test.obj", pos, mesh.faces());
    
    int m = 10;
    int DIM = L.rows();
    Eigen::VectorXd grad(DIM), q(DIM), grad_old(DIM), s(DIM), y(DIM), z(DIM);
    Eigen::MatrixXd sVector = Eigen::MatrixXd ::Zero(DIM, m);
    Eigen::MatrixXd yVector = Eigen::MatrixXd::Zero(DIM, m);
    Eigen::Matrix<double, Eigen::Dynamic, 1> alpha = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(m);
    op.gradient(L, pos, grad);
    Eigen::VectorXd L_old = L;
    
    int iter = 0;
    double gradNorm = 0;
    std::cout<<"Simulation start (L-BFGS)!!"<<std::endl;
    while(true)
    {
        
//        const double relativeEpsilon = static_cast<double>(1e-6) * std::max<double>(static_cast<double>(1.0), L.norm());
//
//        if (grad.norm() < relativeEpsilon)
//            break;
//
//        double H0k = 1;
//        q = grad;
//        const int k = std::min<double>(m, iter);
//        std::cout<<k<<std::endl;
//        std::cout<<q.norm()<<std::endl;
//        // for i = k − 1, k − 2, . . . , k − m§
//        for (int i = k - 1; i >= 0; i--)
//        {
//            if(k < m)
//                break;
//            // alpha_i <- rho_i*s_i^T*q
//            double rho = 1.0 / (sVector.col(i)).dot(yVector.col(i));
//            alpha(i) = rho * sVector.col(i).dot(q);
//            // q <- q - alpha_i*y_i
//            q = q - alpha(i) * yVector.col(i);
//        }
//        std::cout<<"2: "<<q.norm()<<std::endl;
//        // r <- H_k^0*q
//        // update the scaling factor
//        if(k >= m)
//        {
//            H0k = yVector.col(m-1).dot(sVector.col(m-1)) / yVector.col(m-1).dot(yVector.col(m-1));
//        }
//        std::cout<<"H0K: "<<H0k<<std::endl;
//        z = -H0k * q;
//        //for i k − m, k − m + 1, . . . , k − 1
//        for (int i = 0; i <= k-1; i++)
//        {
//            if(k < m)
//                break;
//            // beta <- rho_i * y_i^T * r
//            double rho = 1.0 / (sVector.col(i)).dot(yVector.col(i));
//            double beta = rho * yVector.col(i).dot(z);
//            // r <- r + s_i * ( alpha_i - beta)
//            z = z + sVector.col(i) * (alpha(i) - beta);
//        }
//        std::cout<<"z: "<<z.norm()<<std::endl;
//        double rate = lineSearch(op, L, pos, z);
//        L = L - rate * grad;
//        s = L - L_old;
//        y = grad - grad_old;
//
//        grad_old = grad;
//        op.gradient(L, pos, grad);
//        L_old = L;
//
//        // update the history
//        if (iter < m)
//        {
//            sVector.col(iter) = s;
//            yVector.col(iter) = y;
//        }
//        else
//        {
//            sVector.leftCols(m - 1) = sVector.rightCols(m - 1).eval();
//            sVector.rightCols(1) = s;
//            yVector.leftCols(m - 1) = yVector.rightCols(m - 1).eval();
//            yVector.rightCols(1) = y;
//        }
        
        // find steplength
        double rate = lineSearch(op, L, pos, -grad);
        L = L - rate * grad;
        grad_old = grad;
        op.gradient(L, pos, grad);

        L_old = L;
        std::cout <<std::endl << "iter: "<<iter<< ", f = " <<  op.value(L, pos) << ", ||g||_inf "<<grad.template lpNorm<Eigen::Infinity>()<< ", Rate: "<<rate<<std::endl;
        
        iter++;
        gradNorm = grad.template lpNorm<Eigen::Infinity>();
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
        if(iter % 1000 == 0)
            saveAbars(L, pos);
    }
    
    saveAbars(L, pos);
    
}

void SimulationSetupDynamicSolver::projectBackOp(const SecondFundamentalFormDiscretization &sff, Eigen::VectorXd &L, Eigen::MatrixXd &pos)
{
    double reg = 1e-6;
    int nverts = pos.rows();
    int nfaces = mesh.nFaces();
    
    std::vector<Eigen::Matrix2d> curAbars(nfaces);
    std::vector<Eigen::Matrix2d> bbars(nfaces);
    
    for(int i = 0; i < nfaces; i++)
    {
        curAbars[i] << L(3*i),0,
        L(3*i+1), L(3*i+2);
        curAbars[i] = curAbars[i] * curAbars[i].transpose();
        
        bbars[i].setZero();
    }
    
    Eigen::VectorXd edgeEOFs(0);
    double forceResidual = std::numeric_limits<double>::infinity();
    double updateMag = std::numeric_limits<double>::infinity();
    
    Eigen::SparseMatrix<double> projM;
    
    std::vector<Eigen::Triplet<double>> proj;
    int row = 0;
    for(int i=4;i<initialPos.rows();i++)
    {
        for(int j=0;j<3;j++)
        {
            proj.push_back(Eigen::Triplet<double>(row, 3*i+j, 1.0));
            row ++;
        }
    }
    projM.resize(3*initialPos.rows()-12, 3*initialPos.rows());
    projM.setFromTriplets(proj.begin(), proj.end());
    
    for(int i=0; i< MAX_ITR; i++)
    {
        while (true)
        {
            Eigen::VectorXd derivative;
            std::vector<Eigen::Triplet<double> > hessian;
            double energy = elasticEnergy(mesh, pos, edgeEOFs, lameAlpha, lameBeta, thickness, curAbars, bbars, sff, &derivative, &hessian);
            
            Eigen::SparseMatrix<double> H(3 * nverts, 3 * nverts);
            H.setFromTriplets(hessian.begin(), hessian.end());
            
            
            Eigen::VectorXd force = -derivative;
            
            Eigen::VectorXd reducedForce = projM * force;
            Eigen::SparseMatrix<double> reducedH = projM * H * projM.transpose();
            Eigen::SparseMatrix<double> I(reducedH.rows(), reducedH.cols());
            I.setIdentity();
            reducedH += reg * I;
            
            Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
            
            solver.compute(reducedH);
            while(solver.info() != Eigen::Success)
            {
                reg *=2;
                reducedH += reg * I;
                solver.compute(reducedH);
            }
            
            Eigen::VectorXd descentDir = solver.solve(reducedForce);
            Eigen::VectorXd fullDir = projM.transpose() * descentDir;
            
            Eigen::MatrixXd newPos = pos;
            for (int i = 0; i < nverts; i++)
            {
                newPos.row(i) += fullDir.segment<3>(3 * i);
            }
            
            double newenergy = elasticEnergy(mesh, newPos, edgeEOFs, lameAlpha, lameBeta, thickness, curAbars, bbars, sff, &derivative, NULL);
            force = -derivative;
            
            reducedForce = projM * force;
            forceResidual = reducedForce.norm();
            updateMag =  (projM.transpose() * descentDir).norm();
            if (newenergy <= energy)
            {
                //                std::cout << "Old energy: " << energy << " new energy: " << newenergy << " force residual " << forceResidual << " pos change " << updateMag << std::endl;
                pos = newPos;
                reg = std::max(1e-6, reg * 0.5);
                break;
            }
            else
            {
                reg *= 2.0;
                //                std::cout << "Old energy: " << energy << " new energy: " << newenergy << " lambda now: " << reg << std::endl;
            }
        }
        if(forceResidual < M_TOL || updateMag < M_TOL)
            break;
    }
    std::cout<<forceResidual<<std::endl;
}

void SimulationSetupDynamicSolver::saveAbars(Eigen::VectorXd L, Eigen::MatrixXd pos)
{
    std::cout<<"Saving Abar path: "<<abarPath<<std::endl;
    std::ofstream outfile(abarPath, std::ios::trunc);
    int nverts = targetPos.rows();
    int nfaces = mesh.nFaces();
    
    outfile<<thickness<<"\n";
    outfile<<penaltyCoef<<"\n";
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
    endIdx = resampledPath.rfind("P");
    subString = "";
    if(thickness > 0)
        subString = "T_1e" + std::to_string(expCoef);
    else
        subString = "T_0";
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    // penalty
    if(penaltyCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(penaltyCoef));
    
    startIdx = resampledPath.rfind("P");
    endIdx = resampledPath.rfind("S");
    subString = "";
    if(penaltyCoef > 0)
        subString = "P_1e" + std::to_string(expCoef);
    else
        subString = "P_0";
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
    using namespace cppoptlib;
    convertedProblem op;
    srand((unsigned)time(NULL));
    for(int i=0;i<initialPos.rows();i++)
    {
        initialPos(i,2) = (1e-6*rand())/RAND_MAX;
    }
    op.initialization(initialPos, targetPos, mesh, lameAlpha, lameBeta, thickness);
    op.setPenalty(penaltyCoef, smoothCoef);
    
    int nfaces = mesh.nFaces();
    Eigen::VectorXd L(3*nfaces);
    
    for(int i=0;i<nfaces;i++)
    {
        L(3*i) = sqrt(initialAbars[i](0,0));
        L(3*i+1) = initialAbars[i](0,1)/L(3*i);
        L(3*i+2) = sqrt(initialAbars[i].determinant())/L(3*i);
    }
    
    Eigen::VectorXd initialL(3*mesh.nFaces()), tarL(3*mesh.nFaces());
    for(int i =0; i<mesh.nFaces();i++)
    {
        Eigen::Matrix2d abar = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
        initialL(3*i) = sqrt(abar(0,0));
        initialL(3*i+1) = abar(0,1) / initialL(3*i);
        initialL(3*i+2) = sqrt(abar.determinant()) / initialL(3*i);
        
        Eigen::Matrix2d a = firstFundamentalForm(mesh, targetPos, i, NULL, NULL);
        tarL(3*i) = sqrt(a(0,0));
        tarL(3*i+1) = a(0,1) / tarL(3*i);
        tarL(3*i+2) = sqrt(a.determinant()) / tarL(3*i);
        
    }
    std::cout<<"dE/dq validation: "<<std::endl;
    op.testDifferenceGrad(initialPos);
//    std::cout<<std::endl<<"Position Smoothness Term validation"<<std::endl;
//    op.testPositionSmoothnessGrad(initialPos);
//
//    std::cout<<std::endl<<"Abar Smoothness Term validation"<<std::endl;
//    op.testAbarSmoothnessGrad(tarL);
    std::cout<<std::endl<<"dF/da validation"<<std::endl;
    op.testAbarDerivative(initialL, initialPos);
    std::cout<<std::endl<<"dq/da validation: "<<std::endl;
    op.testDerivativeQ2Abr(tarL, initialPos);
    
    //    double f = op(L);
    //    Eigen::VectorXd gradf;
    //    op.gradient(L, gradf);
    //
    //    Eigen::VectorXd epsVec = Eigen::VectorXd::Random(3*nfaces);
    //    epsVec = epsVec / epsVec.norm();
    //    MidedgeAverageFormulation sff;
    //    for(int i=5; i <= 12; i++ )
    //    {
    //        double eps = pow(10,-i);
    //        Eigen::VectorXd epsL = L + eps * epsVec;
    //        Eigen::MatrixXd pos = op._curPos;
    //        projectBackOp(sff, epsL, pos);
    //        op._curPos = pos;
    //        double epsF = op(L + eps * epsVec);
    //        std::cout<<"EPS is "<<eps<<std::endl;
    //        std::cout<< "Gradient is "<<gradf.dot(epsVec)<< " Finite difference is "<<(epsF - f)/eps<<std::endl;
    //        std::cout<< "The error is "<<abs(gradf.dot(epsVec) - (epsF - f)/eps)<<std::endl;
    //        igl::writeOBJ("test" + std::to_string(i)+".obj", pos, mesh.faces());
    //    }
}
