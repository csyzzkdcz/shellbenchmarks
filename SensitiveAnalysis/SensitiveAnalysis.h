#ifndef SENSITIVEANALYSIS_H
#define SENSITIVEANALYSIS_H

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

class SensitiveAnalysis
/*
 Do sensitive analysis for min E(X, Y, Z), s.t. F(X,Y, Z) = 0
 */
{
public:
    virtual double value(Eigen::VectorXd X, Eigen::VectorXd Y, Eigen::MatrixXd Z) = 0;
    virtual void gradient(Eigen::VectorXd X, Eigen::VectorXd Y, Eigen::MatrixXd Z, Eigen::VectorXd &grad) = 0;
    virtual void projectBack(Eigen::VectorXd X, Eigen::VectorXd &Y, Eigen::MatrixXd &Z) = 0; 
    
    virtual void initialization(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, std::map<int, double> clampedDOFs, double lameAlpha, double lameBeta, double thickness) = 0;
    
    void setPenalty(double abarPenalty, double bbarPenalty, double deltaqPenalty)
    {
        _lambdaAbar = abarPenalty;
        _lambdaBbar = bbarPenalty;
        _mu = deltaqPenalty;
    }
    
    double computeAbarSmoothness(Eigen::VectorXd curL);
    double computeDifference(Eigen::MatrixXd curPos);
    void computeDifferenceGrad(Eigen::MatrixXd curPos, Eigen::VectorXd &grad);
    void computeAbarSmoothnessGrad(Eigen::VectorXd curL, Eigen::VectorXd &grad);
    void testValueGrad(Eigen::VectorXd X, Eigen::VectorXd Y, Eigen::MatrixXd Z);
    
    virtual void test() = 0;
    

protected:
    void generalInitialization(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, std::map<int, double> clampedDOFs, double lameAlpha, double lameBeta, double thickness)
    {
        _lambdaAbar = 0;
        _lambdaBbar = 0;
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
    }
    
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

    
protected:
    double _lambdaAbar;
    double _lambdaBbar;
    double _mu;
    
    Eigen::MatrixXd _initialPos;
    Eigen::MatrixXd _tarPos;
    Eigen::MatrixXd _bcPos;
    
    Eigen::SparseMatrix<double> _laplacianMat;
    
    MeshConnectivity _mesh;
    
    double _lameAlpha;
    double _lameBeta;
    double _thickness;
    double _regionArea;
    
    Eigen::VectorXd _areaList;
    Eigen::VectorXd _massVec;

};
    


#endif
