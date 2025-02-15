#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <Eigen/CholmodSupport>
#include <chrono>

#include "../ElasticShell.h"
#include "../GeometryDerivatives.h"
#include "../GeometryTools.h"
#include "SensitiveAnalysisAbarPos.h"

double SensitiveAnalysisAbarPos::value(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos)
{
    //    return computeDifference(curPos) + _lambda * computeAbarSmoothness(curL) + _mu * computePositionSmoothness(curPos);
    return computeDifference(curPos) + _lambdaAbar * computeAbarSmoothness(curL);
}

void SensitiveAnalysisAbarPos::gradient(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd &grad)
{
    Eigen::VectorXd gradDiff, gradAbarSmooth, gradPosSmooth, gradPos;
    Eigen::SparseMatrix<double> hessianL(3*curPos.rows(), curL.size());
    
    std::vector<Eigen::Triplet<double> > hessian;
    computeAbarDerivative(curL, curPos, &hessian);
    hessianL.setFromTriplets(hessian.begin(), hessian.end());
    
    //    computeDerivativeQ2Abr(curPos, curL, _convertedGrad);
    computeAbarSmoothnessGrad(curL, gradAbarSmooth);
    //    std::cout<<gradAbarSmooth.norm()<<std::endl;
    computeDifferenceGrad(curPos, gradDiff);
    //    computePositionSmoothnessGrad(curPos, gradPosSmooth);
    //    grad.transpose = (gradDiff + _mu * gradPosSmooth).transpose() *_convertedGrad + _lambda * gradAbarSmooth.transpose();
    //    grad.transpose() = gradDiff.transpose() * _convertedGrad;
    
    gradPos = gradDiff;
    //    std::cout<<gradPos.norm()<<std::endl;
    
    Eigen::VectorXd convertedGrad;
    computeConvertedGrad(curPos, curL, gradPos, convertedGrad);
    //    std::cout<<convertedGrad<<std::endl;
    grad.transpose() = convertedGrad.transpose() * hessianL + _lambdaAbar * gradAbarSmooth.transpose();
    
    
}

void SensitiveAnalysisAbarPos::computeConvertedGrad(Eigen::MatrixXd curPos, Eigen::VectorXd curL, Eigen::VectorXd grad, Eigen::VectorXd &convertedGrad)
{
    int nverts = curPos.rows();
    int nfaces =  _mesh.nFaces();
    std::vector<Eigen::Matrix2d> abars(nfaces);
    std::vector<Eigen::Matrix2d> bbars(nfaces);
    for(int i=0;i<nfaces;i++)
    {
        abars[i] << curL(3*i), 0,
        curL(3*i+1),curL(3*i+2);
        abars[i] = abars[i] * abars[i].transpose();
        
        bbars[i].setZero();
    }
    std::vector<Eigen::Triplet<double> > hessian;
    Eigen::VectorXd edgeDOFS(0);
    Eigen::VectorXd f;
    MidedgeAverageFormulation sff;
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &f, &hessian);
    //    std::cout<<f.norm()<<std::endl;
    Eigen::SparseMatrix<double> hessianQ(3*nverts, 3*nverts), hessianL(3*nverts, 3*nfaces);
    hessianQ.setFromTriplets(hessian.begin(), hessian.end());
    
    Eigen::SparseMatrix<double> reductQ = _projM * hessianQ * _projM.transpose();
    
    Eigen::SparseMatrix<double> I(reductQ.rows(), reductQ.cols());
    I.setIdentity();
    double reg = pow(10, -8);
    reductQ += reg*I;
    
    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(reductQ);
    
    Eigen::VectorXd reductGrad = solver.solve(-_projM * grad);
    
    convertedGrad = _projM.transpose() * reductGrad;
    
}

void SensitiveAnalysisAbarPos::computeDerivativeQ2Abr(Eigen::MatrixXd curPos, Eigen::VectorXd curL, Eigen::SparseMatrix<double> &derivative)
{
    int nverts = curPos.rows();
    int nfaces =  _mesh.nFaces();
    std::vector<Eigen::Matrix2d> abars(nfaces);
    std::vector<Eigen::Matrix2d> bbars(nfaces);
    for(int i=0;i<nfaces;i++)
    {
        abars[i] << curL(3*i), 0,
        curL(3*i+1),curL(3*i+2);
        abars[i] = abars[i] * abars[i].transpose();
        
        bbars[i].setZero();
    }
    std::vector<Eigen::Triplet<double> > hessian;
    Eigen::VectorXd edgeDOFS(0);
    Eigen::VectorXd f;
    MidedgeAverageFormulation sff;
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &f, &hessian);
    //    std::cout<<f.norm()<<std::endl;
    Eigen::SparseMatrix<double> hessianQ(3*nverts, 3*nverts), hessianL(3*nverts, 3*nfaces);
    hessianQ.setFromTriplets(hessian.begin(), hessian.end());
    
    Eigen::SparseMatrix<double> reductQ = _projM * hessianQ * _projM.transpose();
    //    std::cout<<hessianQ.norm()<<std::endl;
    
    // Compute the hessian w.r.t. abars.
    hessian.clear();
    computeAbarDerivative(curL, curPos, &hessian);
    hessianL.setFromTriplets(hessian.begin(), hessian.end());
    
    Eigen::SparseMatrix<double> reductL = _projM * hessianL;
    
    Eigen::SparseMatrix<double> I(reductQ.rows(), reductQ.cols());
    I.setIdentity();
    double reg = pow(10, -8);
    reductQ += reg*I;
    
    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(reductQ);
    
    Eigen::MatrixXd dev = solver.solve(-reductL.toDense());
    derivative = _projM.transpose() * dev.sparseView();
    
    //
    ////    std::cout<<hessianL.norm()<<std::endl;
    //    Eigen::SparseMatrix<double> I(3*nverts, 3*nverts);
    //    I.setIdentity();
    //    double reg = pow(10, -8);
    //    hessianQ += reg*I;
    //
    //    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    //    solver.compute(hessianQ);
    //
    //    Eigen::MatrixXd dev = solver.solve(-hessianL.toDense());
    //    derivative = dev.sparseView();
    
    //    std::cout<<derivative.norm()<<std::endl;
    
}

void SensitiveAnalysisAbarPos::computeAbarDerivative(Eigen::VectorXd curL, Eigen::MatrixXd curPos, std::vector<Eigen::Triplet<double> > *grad)
{
    int nfaces =  _mesh.nFaces();
    
    for(int i=0; i<nfaces; i++)
    {
        Eigen::Matrix2d L;
        L << curL(3*i), 0,
        curL(3*i + 1), curL(3*i + 2);
        
        Eigen::Matrix2d abar = L * L.transpose();
        
        Eigen::Matrix2d abarinv = abar.inverse();
        double dA = 0.5 * sqrt(abar.determinant());
        
        Eigen::Matrix<double, 4, 3> abarinvderiv;
        Eigen::Vector3d abarsqrtdetderiv(3);
        
        computeInvMatDeriv(L, abarinvderiv);
        computeSqrtDetDerv(L, abarsqrtdetderiv);
        
        Eigen::Matrix2d a;
        Eigen::Matrix<double, 4, 9> aderiv;
        
        a = firstFundamentalForm(_mesh, curPos, i, &aderiv, NULL);
        
        std::cout.precision(9);
        std::cout<<std::scientific;
        
        if(abs(a(0,1) - a(1,0)) > 1e-8)
        {
            std::cerr<<"Error with asymmetric a"<<std::endl;
            std::cout<<a<<std::endl;
        }
        
        if(abs(abar(0,1) - abar(1,0)) > 1e-8)
        {
            std::cerr<<"Error with asymmetric a"<<std::endl;
            std::cout<<abar<<std::endl;
        }
        
        double coeff = _thickness / 4.0;
        Eigen::Matrix2d M = abarinv * (a - abar);
        Eigen::Matrix<double, 1, 9> traceMderiv;
        traceMderiv.setZero();
        
        traceMderiv += abarinv(0,0) * aderiv.row(0).transpose();
        traceMderiv += abarinv(1,0) * aderiv.row(1).transpose();
        traceMderiv += abarinv(0,1) * aderiv.row(2).transpose();
        traceMderiv += abarinv(1,1) * aderiv.row(3).transpose();
        
        for(int k = 0; k < 3; k++)
        {
            Eigen::Matrix<double, 1, 9> result;
            Eigen::Matrix2d MderivL, MderivAinv, Mainvderiv;
            MderivL << abarinvderiv(0,k), abarinvderiv(1,k),
            abarinvderiv(2,k), abarinvderiv(3,k);
            MderivL = MderivL * a;
            
            result.setZero();
            Eigen::Matrix<double, 1, 4> abarinvderivVeck;
            abarinvderivVeck << abarinvderiv(0,k), abarinvderiv(2,k), abarinvderiv(1,k), abarinvderiv(3, k);
            result += coeff*dA * _lameAlpha * (MderivL.trace() * traceMderiv + M.trace() * abarinvderivVeck * aderiv);
            /*
             trace(A*B) = A11*B11 + A12*B21 + A21*B12 + A22*B22 = (A11, A21, A12, A22) * (B11, B12, B21, B22)
             */
            
            
            
            MderivAinv =  MderivL * abarinv;
            Eigen::Matrix<double, 1, 4> MderivAinvVec, MainvderivVec;
            MderivAinvVec << MderivAinv(0,0), MderivAinv(1,0), MderivAinv(0,1), MderivAinv(1,1);
            
            Eigen::Matrix2d tmpAbarinv;
            tmpAbarinv << abarinvderiv(0,k), abarinvderiv(1,k),
            abarinvderiv(2,k), abarinvderiv(3,k);
            Mainvderiv = M * tmpAbarinv;
            MainvderivVec <<  Mainvderiv(0,0),  Mainvderiv(1,0),  Mainvderiv(0,1),  Mainvderiv(1,1);
            
            
            result += coeff*dA* 2.0 * _lameBeta * (MderivAinvVec * aderiv +  MainvderivVec * aderiv);
            
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(0,0) * aderiv.row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(1,0) * aderiv.row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(0,1) * aderiv.row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(1,1) * aderiv.row(3).transpose();
            Eigen::Matrix2d Mainv = M*abarinv;
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(0, 0) * aderiv.row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(1, 0) * aderiv.row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(0, 1) * aderiv.row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(1, 1) * aderiv.row(3).transpose();
            
            for (int j = 0; j < 3; j++)
            {
                for(int r = 0; r < 3; r++)
                {
                    grad->push_back(Eigen::Triplet<double>(3 * _mesh.faceVertex(i,j) + r, 3*i+k, result(3*j + r)));
                }
            }
            
            
        }
        
        
        // Bending term
        coeff = _thickness * _thickness * _thickness / 12.0;
        
        Eigen::Matrix2d b;
        Eigen::MatrixXd bderiv;
        
        Eigen::VectorXd edgeDOFS(0);
        MidedgeAverageFormulation sff;
        
        
        b = sff.secondFundamentalForm(_mesh, curPos, edgeDOFS, i, &bderiv, NULL);
        
        Eigen::Matrix2d bbar;
        bbar.setZero();
        
        M = abarinv * (b - bbar);
        
        Eigen::Matrix<double, 1, 18> traceMderivb;
        traceMderivb.setZero();
        
        traceMderivb += abarinv(0,0) * bderiv.row(0).transpose();
        traceMderivb += abarinv(1,0) * bderiv.row(1).transpose();
        traceMderivb += abarinv(0,1) * bderiv.row(2).transpose();
        traceMderivb += abarinv(1,1) * bderiv.row(3).transpose();
        
        for(int k = 0; k < 3; k++)
        {
            Eigen::Matrix<double, 1, 18> result;
            Eigen::Matrix2d MderivL, MderivAinv, Mainvderiv;
            MderivL << abarinvderiv(0,k), abarinvderiv(1,k),
            abarinvderiv(2,k), abarinvderiv(3,k);
            MderivL = MderivL * b;
            
            result.setZero();
            Eigen::Matrix<double, 1, 4> abarinvderivVeck;
            abarinvderivVeck << abarinvderiv(0,k), abarinvderiv(2,k), abarinvderiv(1,k), abarinvderiv(3, k);
            result += coeff*dA * _lameAlpha * (MderivL.trace() * traceMderivb + M.trace() * abarinvderivVeck * bderiv);
            /*
             trace(A*B) = A11*B11 + A12*B21 + A21*B12 + A22*B22 = (A11, A21, A12, A22) * (B11, B12, B21, B22)
             */
            
            
            
            MderivAinv =  MderivL * abarinv;
            Eigen::Matrix<double, 1, 4> MderivAinvVec, MainvderivVec;
            MderivAinvVec << MderivAinv(0,0), MderivAinv(1,0), MderivAinv(0,1), MderivAinv(1,1);
            
            Eigen::Matrix2d tmpAbarinv;
            tmpAbarinv << abarinvderiv(0,k), abarinvderiv(1,k),
            abarinvderiv(2,k), abarinvderiv(3,k);
            Mainvderiv = M * tmpAbarinv;
            MainvderivVec <<  Mainvderiv(0,0),  Mainvderiv(1,0),  Mainvderiv(0,1),  Mainvderiv(1,1);
            
            result += coeff*dA* 2.0 * _lameBeta * (MderivAinvVec * bderiv +  MainvderivVec * bderiv);
            
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(0,0) * bderiv.row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(1,0) * bderiv.row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(0,1) * bderiv.row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * _lameAlpha * M.trace() * abarinv(1,1) * bderiv.row(3).transpose();
            Eigen::Matrix2d Mainv = M*abarinv;
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(0, 0) * bderiv.row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(1, 0) * bderiv.row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(0, 1) * bderiv.row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * _lameBeta * Mainv(1, 1) * bderiv.row(3).transpose();
            
            for (int j = 0; j < 3; j++)
            {
                for(int r = 0; r < 3; r++)
                {
                    grad->push_back(Eigen::Triplet<double>(3 * _mesh.faceVertex(i,j) + r, 3*i+k, result(3*j + r)));
                }
                
                int oppidx = _mesh.vertexOppositeFaceEdge(i, j);
                if(oppidx != -1)
                {
                    for(int r = 0; r < 3; r++)
                    {
                        grad->push_back(Eigen::Triplet<double>(3 * oppidx + r, 3*i+k, result(9 + 3*j + r)));
                    }
                }
            }
            
        }
    }
}

double SensitiveAnalysisAbarPos::computePositionSmoothness(Eigen::MatrixXd curPos)
// \int ||df||^2 = - \int f * lap(f) = - f^T * L * f (discretized)
{
    Eigen::Matrix3d resultMat;
    resultMat = -(curPos - _tarPos).transpose() * _laplacianMat * (curPos - _tarPos);
    
    double E = 0.5 * (resultMat(0,0) + resultMat(1,1) + resultMat(2,2));
    return E;
}


void SensitiveAnalysisAbarPos::computePositionSmoothnessGrad(Eigen::MatrixXd curPos, Eigen::VectorXd &grad)
{
    grad.resize(3*curPos.rows());
    grad.setZero();
    for(int i=0;i<3;i++)
    {
        grad += -selectedCoord[i] * _laplacianMat * (curPos - _tarPos).block(0, i, _tarPos.rows(), 1);
    }
    grad = _projM.transpose() * _projM * grad;
}


void SensitiveAnalysisAbarPos::test()
{
    int nfaces = _mesh.nFaces();
    Eigen::VectorXd initialL(3*_mesh.nFaces()), tarL(3*_mesh.nFaces());
    for(int i =0; i<nfaces;i++)
    {
        Eigen::Matrix2d abar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);
        initialL(3*i) = sqrt(abar(0,0));
        initialL(3*i+1) = abar(0,1) / initialL(3*i);
        initialL(3*i+2) = sqrt(abar.determinant()) / initialL(3*i);
        
        Eigen::Matrix2d a = firstFundamentalForm(_mesh, _tarPos, i, NULL, NULL);
        tarL(3*i) = sqrt(a(0,0));
        tarL(3*i+1) = a(0,1) / tarL(3*i);
        tarL(3*i+2) = sqrt(a.determinant()) / tarL(3*i);
        
    }
    std::cout<<"dE/dq validation: "<<std::endl;
    testDifferenceGrad(_initialPos);
    //    std::cout<<std::endl<<"Position Smoothness Term validation"<<std::endl;
    //    op.testPositionSmoothnessGrad(initialPos);
    //
    //    std::cout<<std::endl<<"Abar Smoothness Term validation"<<std::endl;
    //    op.testAbarSmoothnessGrad(tarL);
    std::cout<<std::endl<<"dF/da validation"<<std::endl;
    testAbarDerivative(initialL, _initialPos);
    std::cout<<std::endl<<"dE/da validation: "<<std::endl;
    Eigen::VectorXd S(0);
    testValueGrad(tarL, S, _tarPos);
}

void SensitiveAnalysisAbarPos::testDifferenceGrad(Eigen::MatrixXd curPos)
{
    int nverts = curPos.rows();
    
    srand((unsigned)time(NULL));
    
    Eigen::VectorXd epsVec(3*nverts);
    
    for(int i=0; i<nverts; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*i + k) = epsValue;
        }
    }
    epsVec = _projM.transpose() * _projM *epsVec;
    epsVec = epsVec.normalized();
    
    double value = computeDifference(curPos);
    Eigen::VectorXd J;
    computeDifferenceGrad(curPos, J);
    
    for(int i=3; i< 13; i++)
    {
        Eigen::MatrixXd epsPos = curPos;
        double eps =  pow(10,-i);
        for(int j=0; j<nverts; j++)
        {
            epsPos.row(j) += eps*epsVec.segment(3*j, 3).transpose();
        }
        double value1 = computeDifference(epsPos);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<J.dot(epsVec)<< " Finite difference is "<<(value1 - value)/eps<<std::endl;
        std::cout<< "The error is "<<abs(J.dot(epsVec) - (value1 - value)/eps)<<std::endl;
    }
    
}

void SensitiveAnalysisAbarPos::testAbarSmoothnessGrad(Eigen::VectorXd curL)
{
    int nfaces = _mesh.nFaces();
    
    srand((unsigned)time(NULL));
    
    Eigen::VectorXd epsVec(3*nfaces);
    
    for(int i=0; i<nfaces; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*i + k) = epsValue;
        }
    }
    
    epsVec = epsVec.normalized();
    
    double value = computeAbarSmoothness(curL);
    Eigen::VectorXd J;
    computeAbarSmoothnessGrad(curL, J);
    
    for(int i=3; i< 13; i++)
    {
        double eps =  pow(10,-i);
        double value1 = computeAbarSmoothness(curL + eps*epsVec);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<J.dot(epsVec)<< " Finite difference is "<<(value1 - value)/eps<<std::endl;
        std::cout<< "The error is "<<abs(J.dot(epsVec) - (value1 - value)/eps)<<std::endl;
    }
}

void SensitiveAnalysisAbarPos::testPositionSmoothnessGrad(Eigen::MatrixXd curPos)
{
    int nverts = curPos.rows();
    
    srand((unsigned)time(NULL));
    
    Eigen::VectorXd epsVec(3*nverts);
    
    for(int i=0; i<nverts; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*i + k) = epsValue;
        }
    }
    
    epsVec = epsVec.normalized();
    
    double value = computePositionSmoothness(curPos);
    Eigen::VectorXd J;
    computePositionSmoothnessGrad(curPos, J);
    
    
    for(int i=3; i< 13; i++)
    {
        Eigen::MatrixXd epsPos = curPos;
        double eps =  pow(10,-i);
        for(int j=0; j<nverts; j++)
        {
            epsPos.row(j) += eps*epsVec.segment(3*j, 3).transpose();
        }
        double value1 = computePositionSmoothness(epsPos);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<J.dot(epsVec)<< " Finite difference is "<<(value1 - value)/eps<<std::endl;
        std::cout<< "The error is "<<abs(J.dot(epsVec) - (value1 - value)/eps)<<std::endl;
    }
}

void SensitiveAnalysisAbarPos::testAbarDerivative(Eigen::VectorXd curL, Eigen::MatrixXd curPos)
{
    int nfaces = _mesh.nFaces();
    int nverts = curPos.rows();
    srand((unsigned)time(NULL));
    
    Eigen::VectorXd epsVec(3*nfaces);
    
    for(int i=0; i<nfaces; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*i + k) = epsValue;
        }
    }
    
    epsVec = epsVec.normalized();
    
    std::vector<Eigen::Triplet<double> > T;
    Eigen::SparseMatrix<double> J(3*nverts, 3*nfaces);
    Eigen::VectorXd F;
    Eigen::VectorXd edgeDOFS(0);
    MidedgeAverageFormulation sff;
    std::vector<Eigen::Matrix2d> abars(nfaces);
    std::vector<Eigen::Matrix2d> bbars(nfaces);
    for(int i=0;i<nfaces;i++)
    {
        abars[i] << curL(3*i), 0,
        curL(3*i+1),curL(3*i+2);
        abars[i] = abars[i] * abars[i].transpose();
        
        bbars[i].setZero();
    }
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &F, NULL);
    computeAbarDerivative(curL, curPos, &T);
    J.setFromTriplets(T.begin(), T.end());
    
    for(int i=3; i< 13; i++)
    {
        double eps =  pow(10,-i);
        Eigen::VectorXd epsL = curL + eps*epsVec;
        for(int i=0;i<nfaces;i++)
        {
            abars[i] << epsL(3*i), 0,
            epsL(3*i+1),epsL(3*i+2);
            abars[i] = abars[i] * abars[i].transpose();
            
        }
        Eigen::VectorXd epsF;
        elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &epsF, NULL);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<(J*epsVec).norm()<< " Finite difference is "<<(epsF - F).norm()/eps<<std::endl;
        std::cout<< "The error is "<<(J*epsVec - (epsF - F)/eps).norm()<<std::endl;
    }
}


void SensitiveAnalysisAbarPos::testDerivativeQ2Abr(Eigen::VectorXd curL, Eigen::MatrixXd curPos)
{
    Eigen::VectorXd S(0);
    projectBack(curL, S, curPos);
    Eigen::MatrixXd updatedPos = curPos;
    Eigen::SparseMatrix<double> J;
    computeDerivativeQ2Abr(curPos, curL, J);
    int nfaces = _mesh.nFaces();
    Eigen::VectorXd epsVec = Eigen::VectorXd::Random(3*nfaces);
    epsVec.normalized();
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    for(int i=4; i< 20; i++)
    {
        double eps = pow(4,-i);
        Eigen::VectorXd updatedL = curL + eps*epsVec;
        Eigen::MatrixXd tarPos = curPos;
        Eigen::VectorXd epsPos = eps * J * epsVec;
        for(int i=0; i<tarPos.rows(); i++)
        {
            tarPos.row(i) += epsPos.segment(3*i, 3);
        }
        //        updatedPos = tarPos;
        projectBack(updatedL, S, updatedPos);
        igl::writeOBJ("projected_"+std::to_string(i)+".obj", updatedPos, _mesh.faces());
        //        GeometryTools::rigidMotionTransformation(updatedPos, tarPos, _mesh, R, t);
        //        for(int i=0; i<updatedPos.rows();i++)
        //        {
        //            updatedPos.row(i).transpose() = R * updatedPos.row(i).transpose() + t;
        //        }
        //        std::cout<<R.transpose() * R<<std::endl;
        //        std::cout<<t<<std::endl;
        //        std::cout<<(updatedPos - tarPos).norm()/eps<<std::endl;
        Eigen::MatrixXd errorPos = curPos;
        for(int i=0;i<errorPos.rows();i++)
        {
            errorPos.row(i) = (J*epsVec).segment(3*i, 3);
        }
        
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<(J*epsVec).norm()<< " Finite difference is "<<(updatedPos - curPos).norm()/eps<<std::endl;
        std::cout<< "The error is "<<(errorPos - (updatedPos - curPos)/eps).norm()<<std::endl;
        igl::writeOBJ("iterate_pos"+std::to_string(i) +".obj", tarPos, _mesh.faces());
        
        updatedPos = curPos;
    }
    
}

void SensitiveAnalysisAbarPos::projectBack(Eigen::VectorXd L, Eigen::VectorXd &curS, Eigen::MatrixXd &pos)
{
    double reg = 1e-6;
    int nverts = pos.rows();
    int nfaces = _mesh.nFaces();
    
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
    MidedgeAverageFormulation sff;
    
    int itr;
    
    for(itr=0; itr < 1e4; itr++)
    {
        while (true)
        {
            Eigen::VectorXd derivative;
            std::vector<Eigen::Triplet<double> > hessian;
            double energy = elasticEnergy(_mesh, pos, edgeEOFs, _lameAlpha, _lameBeta, _thickness, curAbars, bbars, sff, &derivative, &hessian);
            
            Eigen::SparseMatrix<double> H(3 * nverts, 3 * nverts);
            H.setFromTriplets(hessian.begin(), hessian.end());
            
            
            Eigen::VectorXd force = -derivative;
            
            Eigen::VectorXd reducedForce = _projM * force;
            Eigen::SparseMatrix<double> reducedH = _projM * H * _projM.transpose();
            Eigen::SparseMatrix<double> I(reducedH.rows(), reducedH.cols());
            I.setIdentity();
            reducedH += reg * I;
            
            Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > solver;
            
            solver.compute(reducedH);
            while(solver.info() != Eigen::ComputationInfo::Success)
            {
                reg *=2;
                reducedH += reg * I;
                solver.compute(reducedH);
            }
            
            Eigen::VectorXd descentDir = solver.solve(reducedForce);
            Eigen::VectorXd fullDir = _projM.transpose() * descentDir;
            
            Eigen::MatrixXd newPos = pos;
            for (int i = 0; i < nverts; i++)
            {
                newPos.row(i) += fullDir.segment<3>(3 * i);
            }
            
            double newenergy = elasticEnergy(_mesh, newPos, edgeEOFs, _lameAlpha, _lameBeta, _thickness, curAbars, bbars, sff, &derivative, NULL);
            force = -derivative;
            
            reducedForce = _projM * force;
            forceResidual = reducedForce.norm();
            updateMag =  (_projM.transpose() * descentDir).norm();
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
        if(forceResidual < 1e-10 || updateMag < 1e-10)
            break;
    }
    std::cout<<"Finishwd with "<<"Force = "<<forceResidual<<" pos change = "<<updateMag<<" iteration time = "<<itr<<std::endl;
    
    if(_projM.rows() == _projM.cols())
    {
        Eigen::Matrix3d R;
        Eigen::Vector3d t;
        
        GeometryTools::rigidMotionTransformation(pos, _tarPos, _mesh, R, t);
        
        Eigen::MatrixXd ones(1, pos.rows());
        ones.setOnes();
        
        pos.transpose() = R * pos.transpose() + t * ones;
    }
    
    
}

