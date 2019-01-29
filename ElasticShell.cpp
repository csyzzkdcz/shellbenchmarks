#include "ElasticShell.h"
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/Sparse>
#include "GeometryDerivatives.h"
#include <random>
#include <iostream>
#include <vector>
#include "MeshConnectivity.h"


double stretchingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,    
    double lameAlpha, double lameBeta, double thickness,
    const Eigen::Matrix2d &abar,
    int face,
    Eigen::Matrix<double, 1, 9> *derivative, // F(face, i)
    Eigen::Matrix<double, 9, 9> *hessian)
{
    double coeff = thickness / 4.0;
    Eigen::Matrix2d abarinv = abar.inverse();
    Eigen::Matrix<double, 4, 9> aderiv;
    std::vector<Eigen::Matrix<double, 9, 9> > ahess;
    Eigen::Matrix2d a = firstFundamentalForm(mesh, curPos, face, (derivative || hessian) ? &aderiv : NULL, hessian ? &ahess : NULL);
    Eigen::Matrix2d M = abarinv * (a - abar);
    double dA = 0.5 * sqrt(abar.determinant());

    double StVK = 0.5 * lameAlpha * M.trace() * M.trace() + lameBeta * (M*M).trace();
    double result = coeff * dA * StVK;

    if (derivative)
    {
        derivative->setZero();
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(0,0) * aderiv.row(0).transpose();
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(1,0) * aderiv.row(1).transpose();
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(0,1) * aderiv.row(2).transpose();
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(1,1) * aderiv.row(3).transpose();
        Eigen::Matrix2d Mainv = M*abarinv;
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(0, 0) * aderiv.row(0).transpose();
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(1, 0) * aderiv.row(1).transpose();
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(0, 1) * aderiv.row(2).transpose();
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(1, 1) * aderiv.row(3).transpose();
    }

    if (hessian)
    {
        hessian->setZero();
        Eigen::Matrix<double, 1, 9> inner = abarinv(0,0) * aderiv.row(0).transpose();
        inner += abarinv(1,0) * aderiv.row(1).transpose();
        inner += abarinv(0,1) * aderiv.row(2).transpose();
        inner += abarinv(1,1) * aderiv.row(3).transpose();
        *hessian += coeff*dA * lameAlpha * inner.transpose() * inner;
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(0,0) * ahess[0];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(1,0) * ahess[1];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(0,1) * ahess[2];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(1,1) * ahess[3];        
        Eigen::Matrix<double, 1, 9> inner00 = abarinv(0, 0) * aderiv.row(0) + abarinv(0, 1) * aderiv.row(2);        
        Eigen::Matrix<double, 1, 9> inner01 = abarinv(0, 0) * aderiv.row(1) + abarinv(0, 1) * aderiv.row(3);
        Eigen::Matrix<double, 1, 9> inner10 = abarinv(1, 0) * aderiv.row(0) + abarinv(1, 1) * aderiv.row(2);
        Eigen::Matrix<double, 1, 9> inner11 = abarinv(1, 0) * aderiv.row(1) + abarinv(1, 1) * aderiv.row(3);
        *hessian += coeff * dA * 2.0 * lameBeta * inner00.transpose() * inner00;
        *hessian += coeff * dA * 2.0 * lameBeta * inner01.transpose() * inner10;
        *hessian += coeff * dA * 2.0 * lameBeta * inner10.transpose() * inner01;
        *hessian += coeff * dA * 2.0 * lameBeta * inner11.transpose() * inner11;
        Eigen::Matrix2d Mainv = M*abarinv;
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(0, 0) * ahess[0];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(1, 0) * ahess[1];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(0, 1) * ahess[2];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(1, 1) * ahess[3];
    }

    return result;
}

double bendingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const Eigen::Matrix2d &abar, const Eigen::Matrix2d &bbar,
    int face,
    const SecondFundamentalFormDiscretization &sff,
    Eigen::MatrixXd *derivative, // F(face, i), then the three vertices opposite F(face,i), then the extra DOFs on oppositeEdge(face,i)
    Eigen::MatrixXd *hessian)
{
    double coeff = thickness*thickness*thickness / 12.0;
    int nedgedofs = sff.numExtraDOFs();
    Eigen::Matrix2d abarinv = abar.inverse();
    Eigen::MatrixXd bderiv(4, 18 + 3*nedgedofs);
    std::vector<Eigen::MatrixXd > bhess;
    Eigen::Matrix2d b = sff.secondFundamentalForm(mesh, curPos, extraDOFs, face, (derivative || hessian) ? &bderiv : NULL, hessian ? &bhess : NULL);
    Eigen::Matrix2d M = abarinv * (b - bbar);
    double dA = 0.5 * sqrt(abar.determinant());

    double StVK = 0.5 * lameAlpha * M.trace() * M.trace() + lameBeta * (M*M).trace();
    double result = coeff * dA * StVK;

    if (derivative)
    {
        derivative->setZero();
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(0,0) * bderiv.row(0);
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(1,0) * bderiv.row(1);
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(0,1) * bderiv.row(2);
        *derivative += coeff*dA * lameAlpha * M.trace() * abarinv(1,1) * bderiv.row(3);
        Eigen::Matrix2d Mainv = M*abarinv;
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(0, 0) * bderiv.row(0);
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(1, 0) * bderiv.row(1);
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(0, 1) * bderiv.row(2);
        *derivative += coeff*dA* 2.0 * lameBeta * Mainv(1, 1) * bderiv.row(3);
    }

    if (hessian)
    {
        hessian->setZero();
        Eigen::MatrixXd inner = abarinv(0,0) * bderiv.row(0);
        inner += abarinv(1,0) * bderiv.row(1);
        inner += abarinv(0,1) * bderiv.row(2);
        inner += abarinv(1,1) * bderiv.row(3);
        *hessian += coeff*dA * lameAlpha * inner.transpose() * inner;
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(0,0) * bhess[0];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(1,0) * bhess[1];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(0,1) * bhess[2];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(1,1) * bhess[3];        
        Eigen::MatrixXd inner00 = abarinv(0, 0) * bderiv.row(0) + abarinv(0, 1) * bderiv.row(2);        
        Eigen::MatrixXd inner01 = abarinv(0, 0) * bderiv.row(1) + abarinv(0, 1) * bderiv.row(3);
        Eigen::MatrixXd inner10 = abarinv(1, 0) * bderiv.row(0) + abarinv(1, 1) * bderiv.row(2);
        Eigen::MatrixXd inner11 = abarinv(1, 0) * bderiv.row(1) + abarinv(1, 1) * bderiv.row(3);
        *hessian += coeff * dA * 2.0 * lameBeta * inner00.transpose() * inner00;
        *hessian += coeff * dA * 2.0 * lameBeta * inner01.transpose() * inner10;
        *hessian += coeff * dA * 2.0 * lameBeta * inner10.transpose() * inner01;
        *hessian += coeff * dA * 2.0 * lameBeta * inner11.transpose() * inner11;
        Eigen::Matrix2d Mainv = M*abarinv;
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(0, 0) * bhess[0];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(1, 0) * bhess[1];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(0, 1) * bhess[2];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(1, 1) * bhess[3];
    }

    return result;
}

double elasticEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d> &abars, 
    const std::vector<Eigen::Matrix2d> &bbars,
    const SecondFundamentalFormDiscretization &sff,
    Eigen::VectorXd *derivative, // positions, then thetas
    std::vector<Eigen::Triplet<double> > *hessian)
{
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nverts = curPos.rows();

    if (derivative)
    {
        derivative->resize(3 * nverts + extraDOFs.size() * nedges);
        derivative->setZero();
    }
    if (hessian)
    {
        hessian->clear();
    }

    double result = 0;
    
    // stretching terms
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix<double, 1, 9> deriv;
        Eigen::Matrix<double, 9, 9> hess;
        result += stretchingEnergy(mesh, curPos, lameAlpha, lameBeta, thickness, abars[i], i, derivative ? &deriv : NULL, hessian ? &hess : NULL);
        if (derivative)
        {
            for (int j = 0; j < 3; j++)
                derivative->segment<3>(3 * mesh.faceVertex(i, j)) += deriv.segment<3>(3 * j);
        }
        if (hessian)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * mesh.faceVertex(i, k) + m, hess(3 * j + l, 3 * k + m)));
                        }
                    }
                }
            }
        }
    }
    
    
    // bending terms
    int nedgedofs = sff.numExtraDOFs();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::MatrixXd deriv(1, 18 + 3 * nedgedofs);
        Eigen::MatrixXd hess(18 + 3 * nedgedofs, 18 + 3 * nedgedofs);
        result += bendingEnergy(mesh, curPos, extraDOFs, lameAlpha, lameBeta, thickness, abars[i], bbars[i], i, sff, derivative ? &deriv : NULL, hessian ? &hess : NULL);
        if (derivative)
        {
            for (int j = 0; j < 3; j++)
            {
                derivative->segment<3>(3 * mesh.faceVertex(i, j)).transpose() += deriv.block<1,3>(0, 3 * j);
                int oppidx = mesh.vertexOppositeFaceEdge(i, j);
                if(oppidx != -1)
                    derivative->segment<3>(3 * oppidx).transpose() += deriv.block<1,3>(0, 9 + 3 * j);
                for (int k = 0; k < nedgedofs; k++)
                {
                    (*derivative)[3 * nverts + nedgedofs * mesh.faceEdge(i, j) + k] += deriv(0, 18 + nedgedofs *j + k);
                }
            }
        }
        if (hessian)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * mesh.faceVertex(i, k) + m, hess(3 * j + l, 3 * k + m)));
                            int oppidxk = mesh.vertexOppositeFaceEdge(i, k);
                            if(oppidxk != -1)
                                hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * oppidxk + m, hess(3 * j + l, 9 + 3 * k + m)));
                            int oppidxj = mesh.vertexOppositeFaceEdge(i, j);
                            if(oppidxj != -1)
                                hessian->push_back(Eigen::Triplet<double>(3 * oppidxj + l, 3 * mesh.faceVertex(i, k) + m, hess(9 + 3 * j + l, 3 * k + m)));
                            if(oppidxj != -1 && oppidxk != -1)
                                hessian->push_back(Eigen::Triplet<double>(3 * oppidxj + l, 3 * oppidxk + m, hess(9 + 3 * j + l, 9 + 3 * k + m)));
                        }
                        for (int m = 0; m < nedgedofs; m++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * mesh.faceVertex(i, j) + l, 3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, hess(3 * j + l, 18 + nedgedofs*k + m)));
                            hessian->push_back(Eigen::Triplet<double>(3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, 3 * mesh.faceVertex(i, j) + l, hess(18 + nedgedofs*k + m, 3 * j + l)));
                            int oppidxj = mesh.vertexOppositeFaceEdge(i, j);
                            if (oppidxj != -1)
                            {
                                hessian->push_back(Eigen::Triplet<double>(3 * oppidxj + l, 3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, hess(9 + 3 * j + l, 18 + nedgedofs * k + m)));
                                hessian->push_back(Eigen::Triplet<double>(3 * nverts + nedgedofs * mesh.faceEdge(i, k) + m, 3 * oppidxj + l, hess(18 + nedgedofs * k + m, 9 + 3 * j + l)));
                            }
                        }
                    }
                    for (int m = 0; m < nedgedofs; m++)
                    {
                        for (int n = 0; n < nedgedofs; n++)
                        {
                            hessian->push_back(Eigen::Triplet<double>(3 * nverts + nedgedofs * mesh.faceEdge(i, j) + m, 3 * nverts + nedgedofs * mesh.faceEdge(i, k) + n, hess(18 + nedgedofs * j + m, 18 + nedgedofs * k + n)));
                        }
                    }
                }
            }
        }
    }
    return result;
}




void testStretchingEnergy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    double eps = 1e-6;
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int ntests = 100;

    double lameAlpha = 1.0;
    double lameBeta = 2.0;
    double thickness = 3.0;
    Eigen::Matrix2d abar;
    abar.setIdentity();

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, nfaces);

    std::cout << "Testing " << ntests << " random faces" << std::endl;

    for(int i=0; i<ntests; i++)
    {
        int face = uni(rng);
        Eigen::Matrix<double, 1, 9> deriv;
        Eigen::Matrix<double, 9, 9> hess;
        double E = stretchingEnergy(mesh, V, lameAlpha, lameBeta, thickness, abar, face, &deriv, &hess);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(mesh.faceVertex(face, j), k) += 1e-6;
                Eigen::Matrix<double, 1, 9> derivpert;
                double Epert = stretchingEnergy(mesh, Vpert, lameAlpha, lameBeta, thickness, abar, face, &derivpert, NULL);
                double findiff = (Epert-E)/1e-6;
                double exact = deriv(0,3*j+k);
                std::cout << "q" << j << "[" << k <<"]: " << exact << " / " << findiff << std::endl;

                Eigen::Matrix<double, 1, 9> findiffhess = (derivpert-deriv)/1e-6;
                std::cout << " hess" << hess.col(3*j+k).transpose() << " / " << findiffhess << std::endl;;                
            }            
        }
    }
}

void testBendingEnergy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const SecondFundamentalFormDiscretization &sff)
{
    double eps = 1e-6;
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nedgedofs = sff.numExtraDOFs();
    Eigen::VectorXd edgeDOFs(nedgedofs*nedges);
    edgeDOFs.setConstant(0.5);
    int ntests = 100;

    double lameAlpha = 1.0;
    double lameBeta = 1.0;
    double thickness = 1.0;
    Eigen::Matrix2d bbar;
    bbar.setZero();

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, nfaces);

    std::cout << "Testing " << ntests << " random faces" << std::endl;
    for(int i=0; i<ntests; i++)
    {
        int face = uni(rng);
        Eigen::Matrix2d abar = firstFundamentalForm(mesh, V, face, NULL, NULL);
        Eigen::MatrixXd deriv(1, 18 + 3 * nedgedofs);
        Eigen::MatrixXd hess(18 + 3 * nedgedofs, 18 + 3 * nedgedofs);
        double E = bendingEnergy(mesh, V, edgeDOFs, lameAlpha, lameBeta, thickness, abar, bbar, face, sff, &deriv, &hess);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(mesh.faceVertex(face, j), k) += 1e-6;
                Eigen::MatrixXd derivpert(1, 18 + 3 * nedgedofs);
                double Epert = bendingEnergy(mesh, Vpert, edgeDOFs, lameAlpha, lameBeta, thickness, abar, bbar, face, sff, &derivpert, NULL);
                double findiff = (Epert-E)/1e-6;
                double exact = deriv(0,3*j+k);
                std::cout << "q" << j << "[" << k <<"]: " << exact << " / " << findiff << std::endl;

                Eigen::MatrixXd findiffhess = (derivpert-deriv)/1e-6;
                std::cout << " hess" << hess.col(3*j+k).transpose() << " / " << findiffhess << std::endl;;                
            }
            int edge = mesh.faceEdge(face, j);
            int ofaceidx = 0;
            if(mesh.edgeFace(edge, ofaceidx) == face)
                ofaceidx = 1;
            if(mesh.edgeFace(edge, ofaceidx) == -1)
                continue;
            int pidx = mesh.edgeOppositeVertex(edge, ofaceidx);
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(pidx, k) += 1e-6;
                Eigen::MatrixXd derivpert(1, 18 + 3 * nedgedofs);
                double Epert = bendingEnergy(mesh, Vpert, edgeDOFs, lameAlpha, lameBeta, thickness, abar, bbar, face, sff, &derivpert, NULL);
                double findiff = (Epert-E)/1e-6;
                double exact = deriv(0,9 + 3*j+k);
                std::cout << "p" << j << "[" << k <<"]: " << exact << " / " << findiff << std::endl;

                Eigen::MatrixXd findiffhess = (derivpert-deriv)/1e-6;

                std::cout << " hess" << hess.col(9 + 3*j+k).transpose() << " / " << findiffhess << std::endl;;                
            }
            Eigen::VectorXd edgepert(edgeDOFs);
            for (int k = 0; k < nedgedofs; k++)
            {
                edgepert[nedgedofs*edge + k] += 1e-6;
                Eigen::MatrixXd derivpert(1, 18 + 3*nedgedofs);
                double Epert = bendingEnergy(mesh, V, edgepert, lameAlpha, lameBeta, thickness, abar, bbar, face, sff, &derivpert, NULL);
                double findiff = (Epert - E) / 1e-6;
                double exact = deriv(0, 18 + nedgedofs *j + k);
                std::cout << "theta[" << j << "]: " << exact << " / " << findiff << std::endl;

                Eigen::MatrixXd findiffhess = (derivpert - deriv) / 1e-6;

                std::cout << " hess" << hess.col(18 + nedgedofs * j + k).transpose() << " / " << findiffhess << std::endl;;
            }
        }
    }
}

void testElasticEnergy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const SecondFundamentalFormDiscretization &sff)
{
    double eps = 1e-6;
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nedgedofs = sff.numExtraDOFs();
    int nverts = V.rows();
    Eigen::VectorXd extraDOFs(nedgedofs * nedges);
    extraDOFs.setConstant(0.3);
    int ntests = 100;

    double lameAlpha = 1.0;
    double lameBeta = 2.0;
    double thickness = 3.0;
    std::vector<Eigen::Matrix2d> bbars;
    Eigen::Matrix2d bbar;
    bbar.setZero();
    for (int i = 0; i < nfaces; i++)
        bbars.push_back(bbar);
    Eigen::Matrix2d abar;
    abar.setIdentity();
    std::vector<Eigen::Matrix2d> abars;
    for (int i = 0; i < nfaces; i++)
    {
        abars.push_back(abar);
    }

    Eigen::VectorXd deriv;
    std::vector<Eigen::Triplet<double> > hess;
    double E = elasticEnergy(mesh, V, extraDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &deriv, &hess);
    Eigen::SparseMatrix<double> H(3 * nverts + nedgedofs * nedges, 3 * nverts + nedgedofs * nedges);
    H.setFromTriplets(hess.begin(), hess.end());
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, 3*nverts);
    std::uniform_int_distribution<int> unie(0, nedgedofs*nedges);

    std::cout << "Testing " << ntests << " random DOFs" << std::endl;

    for (int i = 0; i < ntests; i++)
    {
        int dof = uni(rng);
        std::cout << "vert dof " << dof << std::endl;
        Eigen::MatrixXd Vpert(V);
        Vpert(dof / 3, dof % 3) += eps;
        Eigen::VectorXd derivpert;
        double Epert = elasticEnergy(mesh, Vpert, extraDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivpert, NULL);
        double findiff = (Epert-E)/eps;
        double exact = deriv[dof];
        std::cout << "q[" << dof <<"]: " << exact << " / " << findiff << std::endl;

        std::cout << "hessian: ";
        for (int j = 0; j < 3 * nverts + nedgedofs*nedges; j++)
        {
            std::cout << H.coeff(dof, j) << "/" << (derivpert[j] - deriv[j]) / eps << " ";
        }
        std::cout << std::endl;

        dof = unie(rng);
        std::cout << "edge dof " << dof << std::endl;
        Eigen::VectorXd epert = extraDOFs;
        epert[dof] += eps;
        derivpert;
        Epert = elasticEnergy(mesh, V, epert, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivpert, NULL);
        findiff = (Epert-E)/eps;
        exact = deriv[3*nverts + dof];
        std::cout << "q[" << dof <<"]: " << exact << " / " << findiff << std::endl;
        std::cout << "hessian: ";
        for (int j = 0; j < 3 * nverts + nedgedofs*nedges; j++)
        {
            std::cout << H.coeff(3 * nverts + dof, j) << "/" << (derivpert[j] - deriv[j]) / eps << " ";
        }
        std::cout << std::endl;
    }
}

void testElasticEnergy(const MeshConnectivity &mesh, const Eigen::MatrixXd &V, const Eigen::VectorXd &thetas, double lameAlpha, double lameBeta, double thickness, const std::vector<Eigen::Matrix2d> &abars, const std::vector<Eigen::Matrix2d> &bbars, const SecondFundamentalFormDiscretization &sff)
{
    double eps = 1e-6;

    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int nverts = V.rows();        
    int ntests = 100;

    Eigen::VectorXd deriv;
    std::vector<Eigen::Triplet<double> > hess;
    double E = elasticEnergy(mesh, V, thetas, lameAlpha, lameBeta, thickness, abars, bbars, sff, &deriv, &hess);
    Eigen::SparseMatrix<double> H(3 * nverts + nedges, 3 * nverts + nedges);
    H.setFromTriplets(hess.begin(), hess.end());
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, 3*nverts);

    std::cout << "Testing " << ntests << " random DOFs" << std::endl;

    for (int i = 0; i < ntests; i++)
    {
        int dof = uni(rng);
        Eigen::MatrixXd Vpert(V);
        Vpert(dof / 3, dof % 3) += eps;
        Eigen::VectorXd derivpert;
        double Epert = elasticEnergy(mesh, Vpert, thetas, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivpert, NULL);
        double findiff = (Epert-E)/eps;
        double exact = deriv[dof];
        std::cout << "q[" << dof <<"]: " << exact << " / " << findiff << std::endl;

        std::cout << "hessian: ";
        for (int j = 0; j < 3 * nverts + nedges; j++)
        {
            std::cout << H.coeff(dof, j) << "/" << (derivpert[j] - deriv[j]) / eps << " ";
        }
        std::cout << std::endl;
    }
}
