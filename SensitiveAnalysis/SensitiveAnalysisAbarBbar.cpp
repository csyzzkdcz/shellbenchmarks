#include <Eigen/SPQRSupport>
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#include "SensitiveAnalysisAbarBbar.h"
#include "../ElasticShell.h"
#include "../GeometryTools.h"

double SensitiveAnalysisAbarBbar::value(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos)
{
//    std::cout<<curS.rows()<<" "<<curS.cols()<<" "<<A.rows()<<" "<<A.cols()<<std::endl;
    return 0.5 * (curS.transpose() * A * curS)(0, 0) + computeDifference(curPos) +_lambdaAbar * computeAbarSmoothness(curL);
}

void SensitiveAnalysisAbarBbar::gradient(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd &grad)
{
    grad.resize(curL.rows() + curS.rows());
    Eigen::VectorXd gradDiff, gradAbarSmooth, gradPosSmooth, gradPos;
    Eigen::SparseMatrix<double> gradF2L(3*curPos.rows(), curL.size()), gradF2S(3*curPos.rows(), curS.size());
    
    std::vector<Eigen::Triplet<double> > hessian;
    computeAbarDerivative(curL, curS, curPos, &hessian);
    gradF2L.setFromTriplets(hessian.begin(), hessian.end());
    computeAbarSmoothnessGrad(curL, gradAbarSmooth);
    computeDifferenceGrad(curPos, gradDiff);
    gradPos = gradDiff;
    Eigen::VectorXd convertedGrad;
    computeConvertedGrad(curL, curS, curPos, gradPos, convertedGrad);
    Eigen::VectorXd gradAbar(curL.rows()), gradBbar(curS.rows());
    gradAbar.transpose() = convertedGrad.transpose() * gradF2L + _lambdaAbar * gradAbarSmooth.transpose();
    
    computeBbarDerivative(curL, curS, curPos, &gradF2S);
    gradBbar.transpose() = convertedGrad.transpose() * gradF2S + (A*curS).transpose();
    
    grad.segment(0, curL.size()) = gradAbar;
    grad.segment(curL.size(), curS.size()) = gradBbar;
    
}

void SensitiveAnalysisAbarBbar::projectBack(Eigen::VectorXd curL, Eigen::VectorXd &curS, Eigen::MatrixXd &curPos)
{
//    projectS(curL, curS, curPos);
//    std::cout<<"Project S Done!!"<<std::endl;
    projectPos(curL, curS, curPos);
//    std::cout<<"Project position Done!!"<<std::endl;
}

void SensitiveAnalysisAbarBbar::projectS(Eigen::VectorXd curL, Eigen::VectorXd &curS, Eigen::MatrixXd curPos)
{
    int nfaces = _mesh.nFaces();
    int nverts = _tarPos.rows();
    Eigen::VectorXd c, rhs(3*nfaces + 3*nverts);
    Eigen::SparseMatrix<double> W;
    convertEquilibrium2MatrixForm(curL, curS, curPos, &c, W);
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
    Eigen::VectorXd oldS = curS;
    curS = sol.segment(0, 3*nfaces);
    std::cout<<"bbar changes: "<<(curS - oldS).norm()<<std::endl;
}

void SensitiveAnalysisAbarBbar::projectPos(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd &curPos)
{
    double reg = 1e-6;
    int nverts = curPos.rows();
    int nfaces = _mesh.nFaces();
    
    std::vector<Eigen::Matrix2d> curAbars(nfaces);
    std::vector<Eigen::Matrix2d> bbars(nfaces);
    
    
    for(int i = 0; i < nfaces; i++)
    {
        curAbars[i] << curL(3*i),0,
        curL(3*i+1), curL(3*i+2);
        curAbars[i] = curAbars[i] * curAbars[i].transpose();
        
        bbars[i] = converts2Bbar(curAbars[i], curS.segment(3*i, 3));
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
            double energy = elasticEnergy(_mesh, curPos, edgeEOFs, _lameAlpha, _lameBeta, _thickness, curAbars, bbars, sff, &derivative, &hessian);
            
            Eigen::SparseMatrix<double> H(3 * nverts, 3 * nverts);
            H.setFromTriplets(hessian.begin(), hessian.end());
            
            
            Eigen::VectorXd force = -derivative;
            
            Eigen::VectorXd reducedForce = force;
            Eigen::SparseMatrix<double> reducedH = H;
            Eigen::SparseMatrix<double> I(reducedH.rows(), reducedH.cols());
            I.setIdentity();
            reducedH += reg * I;
            
            Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
            
            solver.compute(reducedH);
            while(solver.info() != Eigen::ComputationInfo::Success)
            {
                reg *=2;
                reducedH += reg * I;
                solver.compute(reducedH);
            }
            
            Eigen::VectorXd descentDir = solver.solve(reducedForce);
            Eigen::VectorXd fullDir =descentDir;
            
            Eigen::MatrixXd newPos = curPos;
            for (int i = 0; i < nverts; i++)
            {
                newPos.row(i) += fullDir.segment<3>(3 * i);
            }
            
            double newenergy = elasticEnergy(_mesh, newPos, edgeEOFs, _lameAlpha, _lameBeta, _thickness, curAbars, bbars, sff, &derivative, NULL);
            force = -derivative;
            
            reducedForce = force;
            forceResidual = reducedForce.norm();
            updateMag =  (descentDir).norm();
            if (newenergy <= energy)
            {
//                std::cout << "Old energy: " << energy << " new energy: " << newenergy << " force residual " << forceResidual << " pos change " << updateMag << std::endl;
                curPos = newPos;
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
    
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    
    GeometryTools::rigidMotionTransformation(curPos, _tarPos, _mesh, R, t);
    
    Eigen::MatrixXd ones(1, curPos.rows());
    ones.setOnes();
    
    curPos.transpose() = R * curPos.transpose() + t * ones;

}


void SensitiveAnalysisAbarBbar::computeConvertedGrad(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd grad, Eigen::VectorXd &convertedGrad)
// convertedGrad = -Hessian(q)^-1 * grad
{
    int nverts = _tarPos.rows();
    int nfaces =  _mesh.nFaces();
    std::vector<Eigen::Matrix2d> abars(nfaces);
    std::vector<Eigen::Matrix2d> bbars(nfaces);
    
    for(int i=0;i<nfaces;i++)
    {
        abars[i] << curL(3*i), 0,
        curL(3*i+1),curL(3*i+2);
        abars[i] = abars[i] * abars[i].transpose();
        
        bbars[i] = converts2Bbar(abars[i], curS.segment(3*i, 3));
    }
    std::vector<Eigen::Triplet<double> > hessian;
    Eigen::VectorXd edgeDOFS(0);
    Eigen::VectorXd f;
    MidedgeAverageFormulation sff;
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &f, &hessian);
    Eigen::SparseMatrix<double> hessianQ(3*nverts, 3*nverts);
    hessianQ.setFromTriplets(hessian.begin(), hessian.end());
    
    Eigen::SparseMatrix<double> reductQ = hessianQ;
    
    Eigen::SparseMatrix<double> I(reductQ.rows(), reductQ.cols());
    I.setIdentity();
    double reg = pow(10, -8);
    reductQ += reg*I;
    
    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(reductQ);
    Eigen::VectorXd reductGrad = solver.solve(-grad);
    convertedGrad = reductGrad;
}

void SensitiveAnalysisAbarBbar::convertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd *C, Eigen::SparseMatrix<double> &W)
{
    int nfaces =  _mesh.nFaces();
    Eigen::VectorXd strechCoef(3*curPos.rows());
    Eigen::VectorXd bendingCoef(3*curPos.rows());
    
    strechCoef.setZero();
    bendingCoef.setZero();
    
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
        Eigen::MatrixXd bendingDeriv(1, 18);
        Eigen::Matrix<double, 1, 18> w1, w2, w3, w4;
        
        Eigen::MatrixXd bderivs;
        
        MidedgeAverageFormulation sff;
        if(C)
        {
            double f = stretchingEnergy(_mesh, curPos, _lameAlpha, _lameBeta, _thickness, abar, i, &strechDeriv, NULL);

            for(int j=0;j<3;j++)
            {
                strechCoef.segment<3>(3*_mesh.faceVertex(i, j)) += strechDeriv.segment<3>(3*j);
            }

            f = bendingEnergy(_mesh, curPos, edgeDOFs, _lameAlpha, _lameBeta, _thickness, abar, zeros, i, sff, &bendingDeriv, NULL);

            for (int j = 0; j < 3; j++)
            {
                bendingCoef.segment<3>(3 * _mesh.faceVertex(i, j)).transpose() += bendingDeriv.block<1,3>(0, 3 * j);
                int oppidx = _mesh.vertexOppositeFaceEdge(i, j);
                if(oppidx != -1)
                    bendingCoef.segment<3>(3 * oppidx).transpose() += bendingDeriv.block<1,3>(0, 9 + 3 * j);
            }
        }
        double coef = _thickness*_thickness*_thickness / 12.0;
        
        Eigen::VectorXd edgeDOFS(0);
        Eigen::Matrix2d b = sff.secondFundamentalForm(_mesh, curPos, edgeDOFS, i, &bderivs, NULL);

        Eigen::Matrix2d mat2, mat3, mat4;
        mat2 << 0,0,1,0;
        mat3 << 1,0,0,-1;
        mat4 << 0,1,0,0;
        
        w1.setZero();
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(0, 0) * bderivs.row(0);
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(1, 0) * bderivs.row(1);
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(0, 1) * bderivs.row(2);
        w1 += 2 * coef * dA * (_lameAlpha + _lameBeta) * abarinv(1, 1) * bderivs.row(3);

        w2.setZero();
        mat2 = mat2 * abarinv;
        w2 += 2 * coef * dA * _lameBeta * mat2(0, 0) * bderivs.row(0);
        w2 += 2 * coef * dA * _lameBeta * mat2(1, 0) * bderivs.row(1);
        w2 += 2 * coef * dA * _lameBeta * mat2(0, 1) * bderivs.row(2);
        w2 += 2 * coef * dA * _lameBeta * mat2(1, 1) * bderivs.row(3);

        w3.setZero();
        mat3 = mat3 * abarinv;
        w3 += 2 * coef * dA * _lameBeta * mat3(0, 0) * bderivs.row(0);
        w3 += 2 * coef * dA * _lameBeta * mat3(1, 0) * bderivs.row(1);
        w3 += 2 * coef * dA * _lameBeta * mat3(0, 1) * bderivs.row(2);
        w3 += 2 * coef * dA * _lameBeta * mat3(1, 1) * bderivs.row(3);
        
        w4.setZero();
        mat4 = mat4 * abarinv;
        w4 += 2 * coef * dA * _lameBeta * mat4(0, 0) * bderivs.row(0);
        w4 += 2 * coef * dA * _lameBeta * mat4(1, 0) * bderivs.row(1);
        w4 += 2 * coef * dA * _lameBeta * mat4(0, 1) * bderivs.row(2);
        w4 += 2 * coef * dA * _lameBeta * mat4(1, 1) * bderivs.row(3);
        
        w2 += abar(1,1) / abar(0,0) * w4;
        w3 += 2 * abar(0,1) / abar(0,0) * w4;
        
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
                    T.push_back(Eigen::Triplet<double>(3*oppidx + k, 3*i, w1(0, 9 + 3*j + k)));
                    T.push_back(Eigen::Triplet<double>(3*oppidx + k, 3*i + 1, w2(0, 9 + 3*j + k)));
                    T.push_back(Eigen::Triplet<double>(3*oppidx + k, 3*i + 2, w3(0, 9 + 3*j + k)));
                }
            }
        }
    }

    W.resize(3*curPos.rows(), 3*nfaces);
    W.setZero();
    W.setFromTriplets(T.begin(), T.end());
    if(C)
        *C = strechCoef + bendingCoef;
}

void SensitiveAnalysisAbarBbar::computeAbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos,std::vector<Eigen::Triplet<double> > *grad)
{
    int nfaces =  _mesh.nFaces();
    Eigen::Matrix2d mat1, mat2, mat3, mat4;
    mat1.setIdentity();
    mat2 << 0,0,1,0;
    mat3 << 1,0,0,-1;
    mat4 << 0,1,0,0;

    for(int i=0; i<nfaces; i++)
    {
        Eigen::Matrix2d L;
        L << curL(3*i), 0,
        curL(3*i + 1), curL(3*i + 2);

        Eigen::Matrix2d abar = L * L.transpose();

        Eigen::Matrix2d abarinv = abar.inverse();
        double dA = 0.5 * sqrt(abar.determinant());

        Eigen::Matrix<double, 4, 3> abarinvderiv, abarinvbbarderiv;
        Eigen::Vector3d abarsqrtdetderiv(3);

        computeInvMatDeriv(L, abarinvderiv);
        computeSqrtDetDerv(L, abarsqrtdetderiv);
        computeAbarinvBbar2s(L, curS.segment(3*i, 3), abarinvbbarderiv);

        Eigen::Matrix2d a;
        Eigen::Matrix<double, 4, 9> aderiv;

        a = firstFundamentalForm(_mesh, curPos, i, &aderiv, NULL);

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
        
        double s1,s2,s3,s4;
        s1 = curS(3*i,0);
        s2 = curS(3*i+1,0);
        s3 = curS(3*i+2,0);
        s4 = abar(1,1) / abar(0,0) * s2 + 2 * abar(0,1) / abar(0,0) * s3;
        Eigen::Matrix2d bbar = abar * (curS(3*i,0) * mat1 + curS(3*i+1, 0) * mat2 + curS(3*i+2, 0) * mat3 + s4 * mat4);
        
        
        M = abarinv * (b-bbar);

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
            
            Eigen::Matrix2d Mderivabarinvbar;
            Mderivabarinvbar << abarinvbbarderiv(0,k), abarinvbbarderiv(1,k),
            abarinvbbarderiv(2,k), abarinvbbarderiv(3,k);
            
            MderivL = MderivL - Mderivabarinvbar;
//            std::cout<<bbar.norm()<<" "<<Mderivabarinvbar.norm()<<std::endl;

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
        
        // The term coming from bbar
        
//        for(int k = 0; k < 3; k++)
//        {
//            Eigen::Matrix<double, 1, 18> result;
//            Eigen::Matrix2d MderivL, MderivAinv, Mainvderiv;
//            MderivL << abarinvderiv(0,k), abarinvderiv(1,k),
//            abarinvderiv(2,k), abarinvderiv(3,k);
//            MderivL = MderivL * b;
//
//            result.setZero();
//
//            Eigen::Matrix2d tmpAbarinv, abarInvBbarAbarinvderiv, tmpABbarinv, abarInvBbarDerivAbarinv;
//            Eigen::Matrix<double, 1, 4> abarInvBbarAbarinvderivVec, abarInvBbarDerivAbarinvVec;
//
//            tmpAbarinv << abarinvderiv(0,k), abarinvderiv(1,k),
//            abarinvderiv(2,k), abarinvderiv(3,k);
//            abarInvBbarAbarinvderiv = abarinv * bbar * tmpAbarinv;
//            abarInvBbarAbarinvderivVec <<  abarInvBbarAbarinvderiv(0,0),  abarInvBbarAbarinvderiv(1,0),  abarInvBbarAbarinvderiv(0,1),  abarInvBbarAbarinvderiv(1,1);
//
//            tmpABbarinv << abarInvBbarAbarinvderivVec(0,k), abarInvBbarAbarinvderivVec(1,k),
//            abarInvBbarAbarinvderivVec(2,k), abarInvBbarAbarinvderivVec(3,k);
//            abarInvBbarDerivAbarinv = tmpABbarinv * abarinv;
//            abarInvBbarDerivAbarinvVec << abarInvBbarDerivAbarinv(0,0), abarInvBbarDerivAbarinv(1,0), abarInvBbarDerivAbarinv(0,1), abarInvBbarDerivAbarinv(1,1);
//
//            result += coeff*dA * _lameBeta * 2 * ( abarInvBbarAbarinvderivVec * bderiv + abarInvBbarDerivAbarinvVec * bderiv );
//
//
//
//            for (int j = 0; j < 3; j++)
//            {
//                for(int r = 0; r < 3; r++)
//                {
//                    grad->push_back(Eigen::Triplet<double>(3 * _mesh.faceVertex(i,j) + r, 3*i+k, -result(3*j + r)));
//                }
//
//                int oppidx = _mesh.vertexOppositeFaceEdge(i, j);
//                if(oppidx != -1)
//                {
//                    for(int r = 0; r < 3; r++)
//                    {
//                        grad->push_back(Eigen::Triplet<double>(3 * oppidx + r, 3*i+k, -result(9 + 3*j + r)));
//                    }
//                }
//            }
//
//
//        }
        
    }
}

void SensitiveAnalysisAbarBbar::computeBbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::SparseMatrix<double> *grad)
{
    convertEquilibrium2MatrixForm(curL, curS, curPos, NULL, *grad);
    *grad = - *grad;
}

Eigen::Vector3d SensitiveAnalysisAbarBbar::convertBar2s(Eigen::Matrix2d abar, Eigen::Matrix2d bbar)
{
    Eigen::Vector3d s;
    Eigen::Matrix2d mat;
    
    mat = abar.inverse() * bbar;
    s(0) = 0.5 * (mat(0,0) + mat(1,1));
    s(1) = mat(1,0);
    s(2) = 0.5 * (mat(0,0) - mat(1,1));
    
    return s;
}

Eigen::Matrix2d SensitiveAnalysisAbarBbar::converts2Bbar(Eigen::Matrix2d abar, Eigen::Vector3d s)
{
    Eigen::Matrix2d mat1, mat2, mat3, mat4;
    mat1.setIdentity();
    mat2 << 0,0,1,0;
    mat3 << 1,0,0,-1;
    mat4 << 0,1,0,0;
    double s4 = abar(1,1) / abar(0,0) * s(1) + 2 * abar(0,1) / abar(0,0) * s(2);
    Eigen::Matrix2d bbar = abar * (s(0) * mat1 + s(1) * mat2 + s(2) * mat3 + s4 * mat4);
    return bbar;
}

void SensitiveAnalysisAbarBbar::test()
{
    int nfaces = _mesh.nFaces();
    Eigen::VectorXd initialL(3*_mesh.nFaces()), tarL(3*_mesh.nFaces());
    Eigen::VectorXd initialS(3*_mesh.nFaces()), tarS(3*_mesh.nFaces());
    MidedgeAverageFormulation sff;
    
    for(int i =0; i<nfaces;i++)
    {
        Eigen::Matrix2d abar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);
        Eigen::Matrix2d bbar = sff.secondFundamentalForm(_mesh, _initialPos, edgeDOFs, i, NULL, NULL);
        initialL(3*i) = sqrt(abar(0,0));
        initialL(3*i+1) = abar(0,1) / initialL(3*i);
        initialL(3*i+2) = sqrt(abar.determinant()) / initialL(3*i);
        
        Eigen::Vector3d s = convertBar2s(abar, bbar);
        initialS.block(3*i, 0, 3, 1) = s;
        
        Eigen::Matrix2d a = firstFundamentalForm(_mesh, _tarPos, i, NULL, NULL);
        Eigen::Matrix2d b = sff.secondFundamentalForm(_mesh, _tarPos, edgeDOFs, i, NULL, NULL);
        testBbarAndSConvertion(a, b);
        tarL(3*i) = sqrt(a(0,0));
        tarL(3*i+1) = a(0,1) / tarL(3*i);
        tarL(3*i+2) = sqrt(a.determinant()) / tarL(3*i);
        
        tarS.segment(3*i, 3) = convertBar2s(a, b);
        
    }
//    std::cout<<std::endl<<"d(a^-1b)/da validation: "<<std::endl;
//    testAbarinvBbar2s(tarL, tarS, _tarPos);
    std::cout<<std::endl<<"dF/da validation"<<std::endl;
    testAbarDerivative(tarL, tarS, _initialPos);
    std::cout<<std::endl<<"Equlibrium validation: "<<std::endl;
    testConvertEquilibrium2MatrixForm(tarL, tarS, _initialPos);
    
    std::cout<<std::endl<<"dE/da validation"<<std::endl;
    Eigen::MatrixXd curPos = _tarPos;
    testValueGrad(tarL, tarS, curPos);
}

void SensitiveAnalysisAbarBbar::testAbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos)
{
    int nfaces = _mesh.nFaces();
    int nverts = _tarPos.rows();
    srand((unsigned)time(NULL));
    
    Eigen::VectorXd epsVec = Eigen::VectorXd::Random(3*nfaces);
    
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

        bbars[i] = converts2Bbar(abars[i], curS.segment(3*i, 3));
    }
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &F, NULL);
    computeAbarDerivative(curL, curS, curPos, &T);
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
            bbars[i] = converts2Bbar(abars[i], curS.segment(3*i, 3));
        }
        Eigen::VectorXd epsF;
        elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &epsF, NULL);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<(J*epsVec).norm()<< " Finite difference is "<<(epsF - F).norm()/eps<<std::endl;
        std::cout<< "The error is "<<(J*epsVec - (epsF - F)/eps).norm()<<std::endl;
    }
}

void SensitiveAnalysisAbarBbar::testConvertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos)
{
    int nfaces = _mesh.nFaces();
    std::vector<Eigen::Matrix2d> abars(nfaces);
    std::vector<Eigen::Matrix2d> bbars(nfaces);
    Eigen::VectorXd F, C;
    Eigen::SparseMatrix<double> W;
    MidedgeAverageFormulation sff;
    for(int i=0;i<nfaces;i++)
    {
        abars[i] << curL(3*i), 0,
        curL(3*i+1),curL(3*i+2);
        abars[i] = abars[i] * abars[i].transpose();
        bbars[i] = converts2Bbar(abars[i], curS.segment(3*i,3));
    }
    elasticEnergy(_mesh, curPos, edgeDOFs, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &F, NULL);
    
    convertEquilibrium2MatrixForm(curL, curS, curPos, &C, W);
    
    std::cout<<"F= "<<F.norm()<<" F - (C-W*S) = "<<(F-C + W*curS).norm()<<std::endl;
}


void SensitiveAnalysisAbarBbar::testBbarAndSConvertion(Eigen::Matrix2d abar, Eigen::Matrix2d bbar)
{
    Eigen::Vector3d s = convertBar2s(abar, bbar);
    Eigen::Matrix2d bbar_new = converts2Bbar(abar, s);
    double error = (bbar - bbar_new).norm();
    if(error > 1e-16)
        std::cout<<error<<std::endl;
}

Eigen::Matrix2d SensitiveAnalysisAbarBbar::computeAbarinvBbar(Eigen::Vector3d L, Eigen::Vector3d s)
{
    Eigen::Matrix2d mat1, mat2, mat3, mat4;
    mat1.setIdentity();
    mat2 << 0,0,1,0;
    mat3 << 1,0,0,-1;
    mat4 << 0,1,0,0;
    
    Eigen::Matrix2d abar;
    abar<< L(0), 0, L(1), L(2);
    abar = abar * abar.transpose();
    double s4 = abar(1,1) / abar(0,0) * s(1) + 2 * abar(0,1) / abar(0,0) * s(2);
    return s(0) * mat1 + s(1) * mat2 + s(2) * mat3 + s4 * mat4;
}

//void SensitiveAnalysisAbarBbar::testAbarinvBbar2s(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos)
//{
//    int i =  int(random() / RAND_MAX * curL.size());
//    Eigen::Vector3d L = curL.segment(3*i, 3);
//    Eigen::Matrix2d L_mat;
//    L_mat << L(0), 0, L(1), L(2);
//    Eigen::Vector3d s = curS.segment(3*i, 3);
//    Eigen::Matrix<double, 4, 3> dAbarinvBbar;
//
//    Eigen::Matrix2d abar, bbar, M;
//    abar = L_mat * L_mat.transpose();
//    bbar = converts2Bbar(abar, s);
//    M = computeAbarinvBbar(L, s);
//    Eigen::Vector4d f;
//    f << M(0,0), M(0,1), M(1,0), M(1,1);
//    s(1) = 0;
//    computeAbarinvBbar2s(L_mat, s, dAbarinvBbar);
//
//    std::cout<<s<<std::endl;
//    std::cout<<L<<std::endl;
//    std::cout<<dAbarinvBbar<<std::endl;
//
//    Eigen::Vector3d epsVec = Eigen::Vector3d::Random(3);
//
//    for(int i=3; i< 13; i++)
//    {
//        double eps =  pow(10,-i);
//        Eigen::Vector3d epsL = L + eps*epsVec;
//        L_mat << epsL(0), 0, epsL(1), epsL(2);
//        abar = L_mat * L_mat.transpose();
//        Eigen::Matrix2d epsM = computeAbarinvBbar(epsL, s);
//        Eigen::Vector4d epsf;
//        epsf << epsM(0,0), epsM(0,1), epsM(1,0), epsM(1,1);
//        std::cout<<"EPS is "<<eps<<std::endl;
//        std::cout<< "Gradient is "<<(dAbarinvBbar*epsVec).norm()<< " Finite difference is "<<(epsf - f).norm()/eps<<std::endl;
//        std::cout<< "The error is "<<(dAbarinvBbar*epsVec - (epsf - f)/eps).norm()<<std::endl;
//    }
//
//}

Eigen::VectorXd SensitiveAnalysisAbarBbar::getInplaneForce()
{
    int nfaces = _mesh.nFaces();
    Eigen::VectorXd initialL(3*_mesh.nFaces()), tarL(3*_mesh.nFaces());
    Eigen::MatrixXd initialS(3*_mesh.nFaces(), 1), tarS(3*_mesh.nFaces(), 1);
    MidedgeAverageFormulation sff;
    
    for(int i =0; i<nfaces;i++)
    {
        Eigen::Matrix2d abar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);
        Eigen::Matrix2d bbar = sff.secondFundamentalForm(_mesh, _initialPos, edgeDOFs, i, NULL, NULL);
        initialL(3*i) = sqrt(abar(0,0));
        initialL(3*i+1) = abar(0,1) / initialL(3*i);
        initialL(3*i+2) = sqrt(abar.determinant()) / initialL(3*i);
        
        Eigen::Vector3d s = convertBar2s(abar, bbar);
        initialS.block(3*i, 0, 3, 1) = s;
        
        Eigen::Matrix2d a = firstFundamentalForm(_mesh, _tarPos, i, NULL, NULL);
        Eigen::Matrix2d b = sff.secondFundamentalForm(_mesh, _tarPos, edgeDOFs, i, NULL, NULL);
        tarL(3*i) = sqrt(a(0,0));
        tarL(3*i+1) = a(0,1) / tarL(3*i);
        tarL(3*i+2) = sqrt(a.determinant()) / tarL(3*i);
        
        s = convertBar2s(a, b);
        tarS.block(3*i, 0, 3, 1) = s;
        
    }
    Eigen::VectorXd c;
    Eigen::SparseMatrix<double> W;
    Eigen::VectorXd epsVec = Eigen::VectorXd::Random(tarL.size());
    epsVec.normalized();
    Eigen::VectorXd L = tarL + 1e-5 * epsVec;
    convertEquilibrium2MatrixForm(L, tarS, _tarPos, &c, W);
    int nverts = _tarPos.rows();
    Eigen::VectorXd rhs(3*nfaces + 3*nverts);
    rhs.setZero();
    rhs.segment(3*nfaces, 3*nverts) = c;
    
    Eigen::MatrixXd M(3*nfaces + 3*nverts, 3*nfaces + 3*nverts);
    M.setZero();
    M.block(0, 0, 3*nfaces, 3*nfaces) = A.toDense();
    M.block(0, 3*nfaces, 3*nfaces, 3*nverts) = W.transpose().toDense();
    M.block(3*nfaces, 0, 3*nverts, 3*nfaces) = W.toDense();
    
    
    Eigen::SparseMatrix<double> sparseM = M.sparseView();
    Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
//    solver.compute(sparseM);
//    Eigen::VectorXd sol = solver.solve(rhs);
//    Eigen::VectorXd  s = sol.segment(0, 3*nfaces);
    solver.compute(W);
    Eigen::VectorXd  s = solver.solve(c);
    return c - W * s;
}
