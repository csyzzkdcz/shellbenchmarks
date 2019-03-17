#include <Eigen/SPQRSupport>
#include "SensitiveAnalysisAbarBbar.h"
#include "../ElasticShell.h"

double SensitiveAnalysisAbarBbar::value(const Eigen::VectorXd &curL, const Eigen::MatrixXd curS)
{
    return 0.5 * (curS.transpose() * A * curS)(0, 0);
}

void SensitiveAnalysisAbarBbar::gradient(const Eigen::VectorXd curL, const Eigen::MatrixXd curS, Eigen::VectorXd &grad)
{
    Eigen::SparseMatrix<double> gradF2abar, gradF2s;
    std::vector<Eigen::Triplet<double>> T;
    computeAbarDerivative(curL, curS, &T);
    gradF2abar.resize(3*_tarPos.rows(), 3*_mesh.nFaces());
    gradF2abar.setFromTriplets(T.begin(), T.end());
    
    convertEquilibrium2MatrixForm(curL, NULL, gradF2s);
    Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
    solver.compute(gradF2s.transpose());
    
    Eigen::VectorXd g = solver.solve(-A.transpose() * curS);
    grad.transpose() = g.transpose() * gradF2abar;
    
}

void SensitiveAnalysisAbarBbar::projectBack(Eigen::VectorXd curL, Eigen::MatrixXd &curS)
{
    int nfaces = _mesh.nFaces();
    int nverts = _tarPos.rows();
    Eigen::VectorXd c, rhs(3*nfaces + 3*nverts);
    Eigen::SparseMatrix<double> W;
    convertEquilibrium2MatrixForm(curL, &c, W);
    rhs.setZero();
    rhs.segment(3*nfaces, 3*nverts) = c;
    
    Eigen::MatrixXd M(3*nfaces + 3*nverts, 3*nfaces + 3*nverts);
    M.setZero();
    M.block(0, 0, 3*nfaces, 3*nfaces) = A.toDense();
    M.block(0, 3*nfaces, 3*nfaces, 3*nverts) = W.transpose().toDense();
    M.block(3*nfaces, 0, 3*nverts, 3*nfaces) = W.toDense();
    
    
    Eigen::SparseMatrix<double> sparseM = M.sparseView();
    Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
    solver.compute(sparseM);
    Eigen::VectorXd sol = solver.solve(rhs);
    curS = sol.segment(0, 3*nfaces);
}

void SensitiveAnalysisAbarBbar::convertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd *C, Eigen::SparseMatrix<double> &W)
{
    int nfaces =  _mesh.nFaces();
    Eigen::VectorXd strechCoef(3*_tarPos.rows());
    Eigen::VectorXd bendingCoef(3*_tarPos.rows());
    std::vector<Eigen::Triplet<double>> T;
    
    for(int i=0;i<nfaces;i++)
    {
        Eigen::Matrix2d L;
        L << curL(3*i), 0,
        curL(3*i + 1), curL(3*i + 2);
        
        Eigen::Matrix2d abar = L * L.transpose();
        double dA = 0.5 * sqrt(abar.determinant());
        Eigen::Matrix2d abarinv = abar.inverse();
        Eigen::Matrix2d zeros;
        zeros.setZero();
        
        
        
        Eigen::Matrix<double, 1, 9> strechDeriv;
        Eigen::MatrixXd bendingDeriv;
        Eigen::Matrix<double, 1, 18> w1, w2, w3;
        
        MidedgeAverageFormulation sff;
        if(C)
        {
            double f = stretchingEnergy(_mesh, _tarPos, _lameAlpha, _lameBeta, _thickness, abar, i, &strechDeriv, NULL);
            
            for(int j=0;j<3;j++)
            {
                strechCoef.segment<3>(3*_mesh.faceVertex(i, j)) += strechDeriv.segment<3>(3*j);
            }
            
            f = bendingEnergy(_mesh, _tarPos, edgeDOFs, _lameAlpha, _lameBeta, _thickness, abar, zeros, i, sff, &bendingDeriv, NULL);
            
            for (int j = 0; j < 3; j++)
            {
                bendingCoef.segment<3>(3 * _mesh.faceVertex(i, j)).transpose() += bendingDeriv.block<1,3>(0, 3 * j);
                int oppidx = _mesh.vertexOppositeFaceEdge(i, j);
                if(oppidx != -1)
                    bendingCoef.segment<3>(3 * oppidx).transpose() += bendingDeriv.block<1,3>(0, 9 + 3 * j);
            }
        }
        double coef = _thickness*_thickness*_thickness / 12.0;
        
        w1.setZero();
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(0, 0) * bderivs[i].row(0);
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(1, 0) * bderivs[i].row(1);
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(0, 1) * bderivs[i].row(2);
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(1, 1) * bderivs[i].row(3);
        
        w2.setZero();
        Eigen::Matrix2d mat2;
        mat2 << 1,0, 0, -1;
        mat2 = mat2 * abarinv;
        w2 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat2(0, 0) * bderivs[i].row(0);
        w2 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat2(1, 0) * bderivs[i].row(1);
        w2 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat2(0, 1) * bderivs[i].row(2);
        w2 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat2(1, 1) * bderivs[i].row(3);
        
        w3.setZero();
        Eigen::Matrix2d mat3;
        mat3 << 0 ,1, 1, 0;
        mat3 = mat3 * abarinv;
        w3 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat3(0, 0) * bderivs[i].row(0);
        w3 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat3(1, 0) * bderivs[i].row(1);
        w3 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat3(0, 1) * bderivs[i].row(2);
        w3 += 2 * coef * dA * (_lameAlpha + _lameBeta) * mat3(1, 1) * bderivs[i].row(3);
        
        for(int j=0;j<3;j++)
        {
            for(int k=0; k<3;k++)
            {
                T.push_back(Eigen::Triplet<double>(3*_mesh.faceVertex(i, j) + k, 3*i, w1(0, 3*j + k)));
                T.push_back(Eigen::Triplet<double>(3*_mesh.faceVertex(i, j) + k, 3*i + 1, w2(0, 3*j + k)));
                T.push_back(Eigen::Triplet<double>(3*_mesh.faceVertex(i, j) + k, 3*i + 2, w3(0, 3*j + k)));
            }
            
            int oppidx = _mesh.vertexOppositeFaceEdge(i, j);
            if(oppidx != -1)
            {
                for(int k=0;k<3;k++)
                {
                    T.push_back(Eigen::Triplet<double>(3*oppidx + k, 3*i, w1(0, 3*j + k)));
                    T.push_back(Eigen::Triplet<double>(3*oppidx + k, 3*i + 1, w2(0, 3*j + k)));
                    T.push_back(Eigen::Triplet<double>(3*oppidx + k, 3*i + 2, w3(0, 3*j + k)));
                }
            }
        }
    }
    
    W.resize(3*_tarPos.rows(), 3*nfaces);
    W.setZero();
    W.setFromTriplets(T.begin(), T.end());
    if(C)
        *C = strechCoef + bendingCoef;
}

void SensitiveAnalysisAbarBbar::computeAbarDerivative(Eigen::VectorXd curL, Eigen::MatrixXd curS, std::vector<Eigen::Triplet<double> > *grad)
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
        
        a = firstFundamentalForm(_mesh, _tarPos, i, &aderiv, NULL);
        
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
        
        
        b = sff.secondFundamentalForm(_mesh, _tarPos, edgeDOFS, i, &bderiv, NULL);
        
        Eigen::Matrix2d bbar, mat1, mat2, mat3;
        
        mat1.setIdentity();
        mat2 << 1,0,0, -1;
        mat3 << 0,1,1,0;
        
        bbar = abar * (curS(3*i,0) * mat1 + curS(3*i+1, 0) * mat2 + curS(3*i+2, 0) * mat3);
        
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
