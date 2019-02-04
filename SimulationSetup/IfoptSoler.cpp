 #include "IfoptSolver.h"
#include "../ElasticShell.h"

 #ifndef MAX_VALUE
 #define MAX_VALUE 1.0e20
 #endif

 using namespace ifopt;

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

    std::vector<Eigen::Matrix2d> abars;
    std::vector<Eigen::Matrix2d> bbars;
    Eigen::MatrixXd curPos;

    convertVariable2ABbarsPos(x, abars, bbars, curPos);

    Eigen::VectorXd derivative;
    Eigen::VectorXd edgeDOFS(0);
    MidedgeAverageFormulation sff;
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, &derivative, NULL);

    return derivative;
}

void optConstraint::FillJacobianBlock(std::string var_set, Jacobian &jac_block) const
{
    Eigen::VectorXd x = GetVariables()->GetComponent("var_set")->GetValues();

    std::vector<Eigen::Matrix2d> abars;
    std::vector<Eigen::Matrix2d> bbars;
    Eigen::MatrixXd curPos;

    convertVariable2ABbarsPos(x, abars, bbars, curPos);

    std::vector<Eigen::Triplet<double> > hessian;
    Eigen::VectorXd edgeDOFS(0);
    MidedgeAverageFormulation sff;
    elasticEnergy(_mesh, curPos, edgeDOFS, _lameAlpha, _lameBeta, _thickness, abars, bbars, sff, NULL, &hessian);

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
    
    Jacobian jacEntire(3*nverts, 3*(nverts+nfaces));
    jacEntire.setFromTriplets(hessian.begin(), hessian.end());
    Eigen::VectorXi boundary = _mesh.getBoundaryLoop();
    jac_block.resize(3*nverts, 3*(nverts+nfaces));
    jac_block.setZero();
    
    for(int i=0; i<boundary.size(); i++)
    {
        jac_block.row(3*boundary(i)) = jacEntire.row(3*boundary(i));
        jac_block.row(3*boundary(i) + 1) = jacEntire.row(3*boundary(i) + 1);
        jac_block.row(3*boundary(i) + 2) = jacEntire.row(3*boundary(i) + 2);
    }

}

double optCost::GetCost() const
{
    VectorXd  x = GetVariables()->GetComponent("var_set")->GetValues();

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

                Eigen::Matrix<double, 3, 2> R;
                Eigen::Matrix<double, 3, 2> oppR;

                R.col(0) = _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0));
                R.col(1) = _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0));

                oppR.col(0) = _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0));
                oppR.col(1) = _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0));

                Eigen::Matrix2d M = (oppR.transpose() * oppR).inverse() * oppR.transpose() * T * R;

                //                std::cout<<L * L.transpose()<<std::endl<<std::endl<<std::endl;
                //                std::cout<<R.transpose()*R<<std::endl<<std::endl;
                //                std::cout<<M.transpose() * Lj * Lj.transpose() * M<<std::endl<<std::endl;
                //                std::cout<<M<<std::endl<<std::endl;
                //
                //                std::cout<<T * R<<std::endl<<std::endl;;
                //                std::cout<<M(0,0) * oppR.col(0) + M(1,0) * oppR.col(1)<<std::endl<<std::endl;
                //                std::cout<<M(0,1) * oppR.col(0) + M(1,1) * oppR.col(1)<<std::endl<<std::endl;

                Eigen::Matrix2d initialAbar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);


                E += 1.0/2.0 * ( (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M) * (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M).transpose() ).trace() / sqrt(initialAbar.determinant()); // devided by det(A0) to make it scalar irrelavent

                //               double value =  1.0/2.0 * (w.col(j).transpose() * L * L.transpose() * w.col(j) - w.col(oppVerIdx).transpose() * Lj * Lj.transpose() * w.col(oppVerIdx)) * (w.col(j).transpose() * L * L.transpose() * w.col(j) - w.col(oppVerIdx).transpose() * Lj * Lj.transpose() * w.col(oppVerIdx));
                //
                //               E += value / sqrt(initialAbars[i].determinant());


            }

        }
    }

    E = _lambda * E;

    Eigen::VectorXi boundary =  _mesh.getBoundaryLoop();

    for(int i=0; i<boundary.size(); i++)
    {
        E += 0.5 * (x.segment<3>(3*boundary(i)).transpose() - _tarPos.row(boundary(i))) * (x.segment<3>(3*boundary(i)).transpose() - _tarPos.row(boundary(i))).transpose();
    }

    return E;

}

void optCost::FillJacobianBlock(std::string var_set, Jacobian &jac) const
{
    VectorXd  x = GetVariables()->GetComponent("var_set")->GetValues();
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

                Eigen::Matrix<double, 3, 2> R;
                Eigen::Matrix<double, 3, 2> oppR;

                R.col(0) = _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0));
                R.col(1) = _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0));

                oppR.col(0) = _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0));
                oppR.col(1) = _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0));

                Eigen::Matrix2d M = (oppR.transpose() * oppR).inverse() * oppR.transpose() * T * R;


                Eigen::Matrix2d initialAbar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);


                for(int k = 0; k < 3; k++)
                {
                    Eigen::Matrix2d abarderiv, abarderivj;
                    abarderiv = Lderivs[k] * L.transpose() + L * Lderivs[k].transpose();
                    abarderivj = Lderivs[k] * Lj.transpose() + Lj * Lderivs[k].transpose();

                    double result;

                    result = _lambda * ( abarderiv * (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M).transpose() ).trace() / sqrt(initialAbar.determinant());

                    J.push_back(Eigen::Triplet<double>(0, 3 * i + k + 3*nverts, result));

                    result = - _lambda * ( M.transpose() * abarderivj * M * (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M).transpose() ).trace() / sqrt(initialAbar.determinant());

                    J.push_back(Eigen::Triplet<double>(0, 3 * oppFace + k + 3*nverts, result));
                }


            }

        }
    }
    Eigen::VectorXi boundary =  _mesh.getBoundaryLoop();

    for(int i=0; i<boundary.size(); i++)
    {
        for(int k=0; k<3; k++)
        {
            double result = x(3*boundary(i)+k) - _tarPos(boundary(i), k);
            J.push_back(Eigen::Triplet<double>(0, 3 * boundary(i), result));
        }
    }
    
    jac.resize(1, 3*(nverts + nfaces));
    jac.setFromTriplets(J.begin(), J.end());
}
