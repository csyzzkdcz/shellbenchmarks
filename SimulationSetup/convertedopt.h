#ifndef CONVERTEDOPT_H
#define CONVERTEDOPT_H

#include<iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/barycenter.h>
#include "../cppoptlib/solver/lbfgssolver.h"
#include "../MeshConnectivity.h"
#include "../SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "../GeometryDerivatives.h"

class convertedProblem
{
public:
    double value(const Eigen::VectorXd &curL, const Eigen::MatrixXd curPos);
    void gradient(const  Eigen::VectorXd curL, const Eigen::MatrixXd curPos, Eigen::VectorXd &grad);
    void update(const Eigen::VectorXd &curL);
    
    Eigen::SparseMatrix<double> getConvertedGrad()
    {
        return _convertedGrad;
    }
   
    void setPenalty(double abarPanalty, double deltaqPenalty)
    {
        _lambda = abarPanalty;
        _mu = deltaqPenalty;
    }
    
public:
    void initialization(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, double lameAlpha, double lameBeta, double thickness)
    {
        _lambda = 0;
        _mu = 0;
        _initialPos = initialPos;
        igl::cotmatrix(initialPos,mesh.faces(),_laplacianMat);
        igl::doublearea(initialPos, mesh.faces(), _areaList);
        igl::barycenter(initialPos, mesh.faces(), _bcPos);
        computeMassMatrix(_massVec, mesh, tarPos);
        _areaList = _areaList/2.0;
        _regionArea = _areaList.sum();
        
        _tarPos = tarPos;
        _mesh = mesh;
        _lameAlpha = lameAlpha;
        _lameBeta = lameBeta;
        _thickness = thickness;
        
//        int nfaces = _mesh.nFaces();
//
//        Eigen::VectorXd initialL(3*nfaces);
//        for(int i =0; i<nfaces;i++)
//        {
//            Eigen::Matrix2d abar = firstFundamentalForm(_mesh, _initialPos, i, NULL, NULL);
//            initialL(3*i) = sqrt(abar(0,0));
//            initialL(3*i+1) = abar(0,1) / initialL(3*i);
//            initialL(3*i+2) = sqrt(abar.determinant()) / initialL(3*i);
//
//        }
        
        selectedCoord.resize(3);
        for(int i =0; i<3; i++)
        {
           selectedCoord[i] = computeSelectMatrix(initialPos.rows(), i);
        }
//        computeDerivativeQ2Abr(_initialPos, initialL, _convertedGrad);
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
        _projM.resize(3*initialPos.rows()-12, 3*initialPos.rows());
        _projM.setFromTriplets(proj.begin(), proj.end());
    }
    
    void testAbarSmoothnessGrad(Eigen::VectorXd curL);
    void testPositionSmoothnessGrad(Eigen::MatrixXd curPos);
    void testDifferenceGrad(Eigen::MatrixXd curPos);
    void testAbarDerivative(Eigen::VectorXd curL, Eigen::MatrixXd curPos);
    void testDerivativeQ2Abr(Eigen::VectorXd curL, Eigen::MatrixXd curPos);
    
public:
    void computeDerivativeQ2Abr(Eigen::MatrixXd curPos, Eigen::VectorXd curL, Eigen::SparseMatrix<double> &derivative); // aAbr = L*L^T, we consider the drivative for L
    void computeAbarDerivative(Eigen::MatrixXd curPos, Eigen::VectorXd curL, std::vector<Eigen::Triplet<double> > *grad);
    
    double computeAbarSmoothness(Eigen::VectorXd curL);
    double computePositionSmoothness(Eigen::MatrixXd curPos);
    double computeDifference(Eigen::MatrixXd curPos);
    
    void computeAbarSmoothnessGrad(Eigen::VectorXd curL, Eigen::VectorXd &grad);
    void computePositionSmoothnessGrad(Eigen::MatrixXd curPos, Eigen::VectorXd &grad);
    void computeDifferenceGrad(Eigen::MatrixXd curPos, Eigen::VectorXd &grad);
    
    
    
private:
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
    
    Eigen::SparseMatrix<double> computeSelectMatrix(int nVerts, int index)
    {
        Eigen::SparseMatrix<double> M(3*nVerts, nVerts);
        std::vector<Eigen::Triplet<double> > triplet;
        triplet.clear();
        for(int i=0;i<nVerts;i++)
        {
            triplet.push_back(Eigen::Triplet<double>(3*i+index, i, 1));
        }
        M.setFromTriplets(triplet.begin(), triplet.end());
        return M;
    }
    
    void computeMassMatrix( Eigen::VectorXd &massVec, MeshConnectivity mesh, Eigen::MatrixXd V);
    void projectBack(Eigen::VectorXd curL, Eigen::MatrixXd &curPos);
    
private:
    Eigen::MatrixXd _initialPos;
    Eigen::MatrixXd _tarPos;
    Eigen::MatrixXd _bcPos;
    
    Eigen::VectorXd _prevL;
    Eigen::VectorXd _massVec;
    Eigen::SparseMatrix<double> _laplacianMat;
    Eigen::SparseMatrix<double> _convertedGrad;
    Eigen::SparseMatrix<double> _projM; // We fixed 4 corner points
    
    Eigen::VectorXd _areaList;
  
    
    MeshConnectivity _mesh;
    
    double _lameAlpha;
    double _lameBeta;
    double _thickness;
    double _regionArea;
    
    double _lambda;
    double _mu;
    
    std::vector<Eigen::SparseMatrix<double>> selectedCoord;
    

};
    


#endif
