#include "IfoptSolver.h"
#include "../ElasticShell.h"
#include <iomanip>
#include <iostream>
#include <fstream>

 #ifndef MAX_VALUE
 #define MAX_VALUE 1.0e20
 #endif

 using namespace ifopt;

int iter = 0;

void computeInvMatDeriv(Eigen::Matrix2d A, Eigen::Matrix<double, 4, 3> &dA)
{
    double x,y,z;
    x = A(0,0);
    y = A(1,0);
    z = A(1,1);
    Eigen::Matrix2d M, dA1, dA2, dA3;
    M << y*y+z*z,-x*y,
    -x*y,x*x;
    std::vector<Eigen::Matrix2d> C(3);
    C[0]<<0,-y,
    -y,2*x;
    C[1]<<2*y,-x,
    -x,0;
    C[2]<<2*z,0,
    0,0;

    dA1=1/(x*x*z*z)*C[0] - 2/(x*x*x*z*z)*M;
    dA2=1/(x*x*z*z)*C[1];
    dA3=1/(x*x*z*z)*C[2] - 2/(x*x*z*z*z)*M;

    dA.row(0) << dA1(0,0), dA2(0,0), dA3(0,0);
    dA.row(1) << dA1(0,1), dA2(0,1), dA3(0,1);
    dA.row(2) << dA1(1,0), dA2(1,0), dA3(1,0);
    dA.row(3) << dA1(1,1), dA2(1,1), dA3(1,1);
}

void computeSqrtDetDerv(Eigen::Matrix2d A, Eigen::Vector3d & diffSqrtDet)
{
    int sign = 1;
    if(A.determinant()<0)
        sign = -1;
    diffSqrtDet(0) = sign*A(1,1);
    diffSqrtDet(1) = 0;
    diffSqrtDet(2) = sign*A(0,0);
}

void optConstraint::convertVariable2ABbarsPos(Eigen::VectorXd x, std::vector<Eigen::Matrix2d> &abars, std::vector<Eigen::Matrix2d> &bbars, Eigen::MatrixXd &curPos) const
{
    int nfaces =  _mesh.nFaces();
    int nverts = x.size() / 3 - nfaces;

    abars.resize(nfaces);
    bbars.resize(nfaces);
    curPos.resize(nverts, 3);

    for(int i=0; i<nverts; i++)
    {
        curPos.row(i) = x.segment<3>(3*i);
    }

    for(int i=0; i< nfaces; i++)
    {
        abars[i] << x(3*i + 3*nverts), 0,
        x(3*i + 3*nverts + 1), x(3*i + 3*nverts + 2);
        abars[i] = abars[i]*abars[i].transpose();
        bbars[i].setZero();
    }
}

Eigen::VectorXd optConstraint::GetValues() const
{
    Eigen::VectorXd x = GetVariables()->GetComponent("var_set")->GetValues();
    return getValues(x);
}

Eigen::VectorXd optConstraint::getValues(Eigen::VectorXd x) const
{
    std::vector<Eigen::Matrix2d> abars;
    std::vector<Eigen::Matrix2d> bbars;
    Eigen::MatrixXd curPos;
    
    convertVariable2ABbarsPos(x, abars, bbars, curPos);
    
    Eigen::VectorXd derivative;
    Eigen::VectorXd edgeDOFS(0);
    MidedgeAverageFormulation sff;
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &derivative, NULL);
//    std::cout<<"The norm of the constrains is "<<derivative.norm()<<std::endl;
    return derivative;
}

void optConstraint::FillJacobianBlock(std::string var_set, Jacobian &jac_block) const
{
    Eigen::VectorXd x = GetVariables()->GetComponent("var_set")->GetValues();
    fillJacobianBlock(x, jac_block);
}

void optConstraint::fillJacobianBlock(Eigen::VectorXd x, Jacobian &jac_block) const
{
    std::vector<Eigen::Matrix2d> abars;
    std::vector<Eigen::Matrix2d> bbars;
    Eigen::MatrixXd curPos;
    
    convertVariable2ABbarsPos(x, abars, bbars, curPos);
    
    std::vector<Eigen::Triplet<double> > hessian;
    Eigen::VectorXd edgeDOFS(0);
    MidedgeAverageFormulation sff;
//    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, NULL, &hessian);
    
    int nfaces =  _mesh.nFaces();
    int nverts = x.size() / 3 - nfaces;

    for(int i=0; i<nfaces; i++)
    {
        Eigen::Matrix2d L;
        L << x(3*i + 3*nverts), 0,
        x(3*i + 3*nverts + 1), x(3*i + 3*nverts + 2);

        Eigen::Matrix2d abarinv = abars[i].inverse();
        double dA = 0.5 * sqrt(abars[i].determinant());

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

        if(abs(abars[i](0,1) - abars[i](1,0)) > 1e-8)
        {
            std::cerr<<"Error with asymmetric a"<<std::endl;
            std::cout<<abars[i]<<std::endl;
        }

        double coeff = _thickness / 4.0;
        Eigen::Matrix2d M = abarinv * (a - abars[i]);
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
                    hessian.push_back(Eigen::Triplet<double>(3 * _mesh.faceVertex(i,j) + r, 3*(i+nverts) + k, result(3*j + r)));
                }
            }


        }


        // Bending term
        coeff = _thickness * _thickness * _thickness / 12.0;

        Eigen::Matrix2d b;
        Eigen::MatrixXd bderiv;

        b = sff.secondFundamentalForm(_mesh, curPos, edgeDOFS, i, &bderiv, NULL);

        M = abarinv * (b - bbars[i]);

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
                    hessian.push_back(Eigen::Triplet<double>(3 * _mesh.faceVertex(i,j) + r, 3*(i+nverts)+k, result(3*j + r)));
                }

                int oppidx = _mesh.vertexOppositeFaceEdge(i, j);
                if(oppidx != -1)
                {
                    for(int r = 0; r < 3; r++)
                    {
                        hessian.push_back(Eigen::Triplet<double>(3 * oppidx + r, 3*(i+nverts)+k, result(9 + 3*j + r)));
                    }
                }
            }

        }
    }
    
//    Jacobian jacEntire(3*nverts, 3*(nverts+nfaces));
//    jacEntire.setFromTriplets(hessian.begin(), hessian.end());
    Eigen::VectorXi boundary = _mesh.getBoundaryLoop();
    jac_block.resize(3*nverts, 3*(nverts+nfaces));
    jac_block.setZero();
    jac_block.setFromTriplets(hessian.begin(), hessian.end());
    // std::vector<Eigen::Triplet<double>> nonzeroCols;
    // for(int i=0; i<boundary.size(); i++)
    // {
    //     for(int k=0; k<3; k++)
    //         nonzeroCols.push_back(Eigen::Triplet<double>(3*boundary(i) + k, 3*boundary(i) + k, 1));
    // }
    // for(int i=0; i<nfaces; i++)
    // {
    //     for(int k=0; k<3; k++)
    //         nonzeroCols.push_back(Eigen::Triplet<double>(3*i + 3*nverts + k, 3*i + 3*nverts + k, 1));
    // }
    // Jacobian nonzeroColMat;
    // nonzeroColMat.resize(3*(nverts+nfaces), 3*(nverts+nfaces));
    // nonzeroColMat.setFromTriplets(nonzeroCols.begin(), nonzeroCols.end());
    // jac_block = jac_block*nonzeroColMat;
}


void optConstraint::testValueJacobian(Eigen::VectorXd x)
{
    int nfaces = _mesh.nFaces();
    int nverts =  x.size() / 3 - nfaces;
    
    Eigen::VectorXi boundary = _mesh.getBoundaryLoop();
    
    Eigen::VectorXd epsVec(x.size());
    srand((unsigned)time(NULL));
    // for(int i=0; i<boundary.size(); i++)
    // {
    //     for(int k=0; k<3; k++)
    //     {
    //         double epsValue =  random();
    //         epsVec(3*boundary(i) + k) = epsValue;
    //     }
    // }

    for(int i=0; i<nverts; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*i + k) = epsValue;
        }
    }
    
    for(int i=0; i<nfaces; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*nverts + 3*i + k) = epsValue;
        }
    }
    
    epsVec = epsVec.normalized();
    
    Eigen::VectorXd value = getValues(x);
    Jacobian J;
    fillJacobianBlock(x, J);
    
    int selectedIdx = rand() % x.size();
    
    for(int i=3; i< 13; i++)
    {
        double eps =  pow(10,-i);
        Eigen::VectorXd x1 = x + eps * epsVec;
        Eigen::VectorXd value1 = getValues(x1);
        std::cout<<"EPS is "<<eps<<" Selected index is "<<selectedIdx<<std::endl;
        std::cout<<"Finite difference is "<< (value1(i) - value(i))/eps << " gradient is "<<(J*epsVec)(i)<<std::endl;
        std::cout<<"Error is "<<abs( (value1(i) - value(i))/eps - (J*epsVec)(i) )<<std::endl<<std::endl;;
        std::cout<<"The total norm of Finite difference is "<< ( (value1 - value)/eps ).norm()<<" gradient is " <<(J*epsVec).norm()<<std::endl;
        std::cout<<"The norm of difference is "<<( (value1 - value)/eps - J*epsVec ).norm()<<std::endl<<std::endl;
    }
}



// optCost Class
double optCost::GetCost() const
{
    VectorXd  x = GetVariables()->GetComponent("var_set")->GetValues();
    return getCost(x);

}

double optCost::getCost(Eigen::VectorXd x) const
{
    iter ++;
    if(iter % 1000 == 0)    // Saving x every 1000 evaluation
    {
        int nfaces = _mesh.nFaces();
        int nverts = x.size()/3 - nfaces;
        std::ofstream outfile("L_list.dat", std::ios::trunc);
        
        outfile<<1e-4<<"\n";
        outfile<<_lambda<<"\n";
        outfile<<_mu<<"\n";
        outfile<<3*nverts<<"\n";
        outfile<<3*nfaces<<"\n";
        
        for(int i=0;i<3*nverts;i++)
        {
            outfile<<std::setprecision(16)<<x(i)<<"\n";
        }
        
        for(int i=0;i<3*nfaces;i++)
        {
            outfile<<std::setprecision(16)<<x(3*nverts + i)<<"\n";
        }
        outfile<<std::setprecision(16)<<x(x.size()-1);
        outfile.close();
    }
    double E = 0;

//    E = getDifference(x);

    E += _lambda * getPenalty(x);
    
//    E += _mu * getSmoothness(x);

    return E; 
    
}


void optCost::FillJacobianBlock(std::string var_set, Jacobian &jac) const
{
    VectorXd  x = GetVariables()->GetComponent("var_set")->GetValues();
    fillJacobianBlock(x, jac);
}

void optCost::fillJacobianBlock(Eigen::VectorXd x, Jacobian &jac) const
{
    std::vector<Eigen::Triplet<double>> J;
    
    int nfaces =  _mesh.nFaces();
    int nverts = x.size() / 3 - nfaces;
    
    Eigen::Matrix<double, 2, 3> w;
    w << 1, 0, 1,
    -1, 1, 0;
    
    std::vector<Eigen::Matrix2d> Lderivs(3);
    
    Lderivs[0] << 1,0,
    0,0;
    
    Lderivs[1] << 0,0,
    1,0;
    
    Lderivs[2] << 0,0,
    0,1;
    
    
    for(int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix2d L;
        
        L << x(3*i + 3*nverts), 0,
        x(3*i+1 + 3*nverts), x(3*i+2 + 3*nverts);
        for(int j = 0; j < 3; j++)
        {
            int oppVerIdx =  _mesh.vertexOppositeFaceEdgeIndex(i, j);
            int oppFace = _mesh.faceOppositeVertex(i, j);
            if (oppFace != -1)
            {
                Eigen::Matrix2d Lj;
                Lj << x(3*oppFace + 3*nverts), 0,
                x(3*oppFace+1 + 3*nverts), x(3*oppFace+2 + 3*nverts);
                
                Eigen::Matrix2d R, oppR, abar, abarj;
                
                R.col(0) = ( _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                R.col(1) = ( _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                oppR.col(0) = ( _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                oppR.col(1) = ( _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                
                abar = (R.inverse()).transpose() * L * L.transpose() * R.inverse();
                
                abarj = (oppR.inverse()).transpose() * Lj * Lj.transpose() * oppR.inverse();
                
                Eigen::Matrix2d initialAbar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);
                
                
                
                for(int k = 0; k < 3; k++)
                {
                    Eigen::Matrix2d abarderiv, abarderivj;
                    abarderiv = (R.inverse()).transpose() *( Lderivs[k] * L.transpose() + L * Lderivs[k].transpose() ) * R.inverse();
                    abarderivj = (oppR.inverse()).transpose() * (Lderivs[k] * Lj.transpose() + Lj * Lderivs[k].transpose() ) * oppR.inverse();
                    
                    double result;
                    double area = ( _areaList(i) + _areaList(oppFace) ) / 3.0;
                    
                    result = 2.0 / _regionArea * _lambda * ( abarderiv * (abar - abarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                    
                    J.push_back(Eigen::Triplet<double>(0, 3 * i + k + 3*nverts, result));
                    
                    result = - 2.0 / _regionArea * _lambda * ( abarderivj * (abar - abarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                    
                    J.push_back(Eigen::Triplet<double>(0, 3 * oppFace + k + 3*nverts, result));
                }
                
            }
            
        }
    }
    // Eigen::VectorXi boundary =  _mesh.getBoundaryLoop();
    
    // for(int i=0; i<boundary.size(); i++)
    // {
    //     for(int k=0; k<3; k++)
    //     {
    //         double result = x(3*boundary(i)+k) - _tarPos(boundary(i), k);
    //         J.push_back(Eigen::Triplet<double>(0, 3 * boundary(i)+k, result));
    //     }
    // }

//    
//    for(int i=0; i<_tarPos.rows(); i++)
//    {
//        for(int k=0; k<3; k++)
//        {
//            double result = x(3*i + k) - _tarPos(i, k);
//            J.push_back(Eigen::Triplet<double>(0, 3 * i + k, result));
//        }
//    }
 
//    jac.resize(1, 3*(nverts + nfaces));
//    jac.setFromTriplets(J.begin(), J.end());
//
//
//    Jacobian smoothJ = -(selectedX.transpose()*L*selectedX*(x-_tarPosVec)).transpose().sparseView();
//    jac += _mu * smoothJ;
//
//    smoothJ = -(selectedY.transpose()*L*selectedY*(x-_tarPosVec)).transpose().sparseView();
//    jac += _mu * smoothJ;
//
//    smoothJ = -(selectedZ.transpose()*L*selectedZ*(x-_tarPosVec)).transpose().sparseView();
//    jac += _mu * smoothJ;
  
}


double optCost::getPenalty(Eigen::VectorXd x) const
{
    int nfaces =  _mesh.nFaces();
    int nverts = x.size() / 3 - nfaces;
    
    double E = 0;
    Eigen::Matrix<double, 2, 3> w;
    w << 1, 0, 1,
    -1, 1, 0;
    
    std::vector<Eigen::Matrix2d> Lderivs(3);
    
    Lderivs[0] << 1,0,
    0,0;
    
    Lderivs[1] << 0,0,
    1,0;
    
    Lderivs[2] << 0,0,
    0,1;
    
    
    for(int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix2d L;
        
        L << x(3*i + 3*nverts), 0,
        x(3*i+1 + 3*nverts), x(3*i+2 + 3*nverts);
        for(int j = 0; j < 3; j++)
        {
            int oppVerIdx =  _mesh.vertexOppositeFaceEdgeIndex(i, j);
            int oppFace = _mesh.faceOppositeVertex(i, j);
            if (oppFace != -1)
            {
                Eigen::Matrix2d Lj;
                Lj << x(3*oppFace + 3*nverts), 0,
                x(3*oppFace+1 + 3*nverts), x(3*oppFace+2 + 3*nverts);
                
                // Compute the tranfermation matrix M
                
                Eigen::Vector3d oppNormal, curNormal, oppEdge;
                
                oppNormal = faceNormal(_mesh, _initialPos, oppFace, oppVerIdx, NULL, NULL);
                curNormal = faceNormal(_mesh, _initialPos, i, j, NULL, NULL);
                
                oppNormal = oppNormal/oppNormal.norm();
                curNormal = curNormal/curNormal.norm();
                
                oppEdge = _initialPos.row(_mesh.faceVertex(i, (j + 1)%3)) - _initialPos.row(_mesh.faceVertex(i, (j + 2)%3));
                oppEdge = oppEdge/oppEdge.norm();
                
                Eigen::Matrix3d A, A1, T;
                A.col(0) = oppEdge;
                A.col(1) = curNormal;
                A.col(2) = oppEdge.cross(curNormal);
                
                A1.col(0) = oppEdge;
                A1.col(1) = oppNormal;
                A1.col(2) = oppEdge.cross(oppNormal);
                
                T = A1*A.inverse();
                
                Eigen::Matrix2d R, oppR, abar, abarj;
                
                R.col(0) = ( _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                R.col(1) = ( _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                oppR.col(0) = ( _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                oppR.col(1) = ( _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                
                abar = (R.transpose()).inverse() * L * L.transpose() * R.inverse();
                
                abarj = (oppR.transpose()).inverse() * Lj * Lj.transpose() * oppR.inverse();
                
                Eigen::Matrix2d initialAbar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);
                
                double area = ( _areaList(i) + _areaList(oppFace) ) / 3.0;
                
                E +=  1.0 / _regionArea * ( (abar - abarj) * (abar - abarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                
                
            }
            
        }
    }

    return E;
}

double optCost::getDifference(Eigen::VectorXd x) const
{
    double E = 0;
    Eigen::VectorXi boundary =  _mesh.getBoundaryLoop();
    
    // for(int i=0; i<boundary.size(); i++)
    // {
    //     E += 0.5 * (x.segment<3>(3*boundary(i)).transpose() - _tarPos.row(boundary(i))) * (x.segment<3>(3*boundary(i)).transpose() - _tarPos.row(boundary(i))).transpose();
    // }

    for(int i=0; i<_tarPos.rows(); i++)
    {
        E += 0.5 * (x.segment<3>(3*i).transpose() - _tarPos.row(i)) * (x.segment<3>(3*i).transpose() - _tarPos.row(i)).transpose(); 
    }

    return E;
}

double optCost::getSmoothness(Eigen::VectorXd x) const
{
    double E = 0;
    
    x = x - _tarPosVec;
    E += - 0.5 * (selectedX * x).transpose() * L * (selectedX * x);
    E += - 0.5 * (selectedY * x).transpose() * L * (selectedY * x);
    E += - 0.5 * (selectedZ * x).transpose() * L * (selectedZ * x);
    
    return E;
    
}

Eigen::SparseMatrix<double> optCost::computeSelectMatrix(int nVerts, int nFaces, int index)
{
    Eigen::SparseMatrix<double> M(nVerts, 3*(nVerts + nFaces));
    std::vector<Eigen::Triplet<double> > triplet;
    triplet.clear();
    for(int i=0;i<nVerts;i++)
    {
        triplet.push_back(Eigen::Triplet<double>(i, 3*i+index, 1));
    }
    M.setFromTriplets(triplet.begin(), triplet.end());
    return M;
}

void optCost::testCostJacobian(Eigen::VectorXd x)
{
    int nfaces = _mesh.nFaces();
    int nverts =  x.size() / 3 - nfaces;
    
    Eigen::VectorXi boundary = _mesh.getBoundaryLoop();
    
    Eigen::VectorXd epsVec(x.size());
    srand((unsigned)time(NULL));
    // for(int i=0; i<boundary.size(); i++)
    // {
    //     for(int k=0; k<3; k++)
    //     {
    //         double epsValue =  random();
    //         epsVec(3*boundary(i) + k) = epsValue;
    //     }
    // }

    for(int i=0; i<nverts; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*i + k) = epsValue;
        }
    }
    
    for(int i=0; i<nfaces; i++)
    {
        for(int k=0; k<3; k++)
        {
            double epsValue =  random();
            epsVec(3*nverts + 3*i + k) = epsValue;
        }
    }
    
    epsVec = epsVec.normalized();
    
    double value = getCost(x);
    Jacobian J;
    fillJacobianBlock(x, J);
    
    for(int i=3; i< 13; i++)
    {
        double eps =  pow(10,-i);
        Eigen::VectorXd x1 = x + eps * epsVec;
        double value1 = getCost(x1);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<(J*epsVec)(0,0)<< " Finite difference is "<<(value1 - value)/eps<<std::endl;
        std::cout<< "The error is "<<abs((J*epsVec)(0,0) - (value1 - value)/eps)<<std::endl;
    }
}
