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
        
        int nfaces = mesh.nFaces();
        _atarget.resize(nfaces);
        _starget.resize(nfaces);
        
        MidedgeAverageFormulation sff;
        Eigen::VectorXd extraDOFs(0);
        for(int i=0;i<nfaces;i++)
        {
            _atarget[i] = firstFundamentalForm(_mesh, _tarPos, i, NULL, NULL);
            Eigen::Matrix2d b = sff.secondFundamentalForm(_mesh, _tarPos, extraDOFs, i, NULL, NULL);
            _starget(i) = (_atarget[i].inverse() * b).trace() / 2.0;
//            _starget(i) = 0;
        }
        
        
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

    
public:
    double _lambdaAbar;
    double _lambdaBbar;
    double _mu;
    
    Eigen::MatrixXd _initialPos;
    Eigen::MatrixXd _tarPos;
    Eigen::MatrixXd _bcPos;
    
    Eigen::SparseMatrix<double> _laplacianMat;
    std::vector<Eigen::Matrix2d> _atarget;
    Eigen::VectorXd _starget;
    
    MeshConnectivity _mesh;
    
    double _lameAlpha;
    double _lameBeta;
    double _thickness;
    double _regionArea;
    
    Eigen::VectorXd _areaList;
    Eigen::VectorXd _massVec;
    

// Box Constraints
public:
    Eigen::VectorXd m_lowerBound;
    Eigen::VectorXd m_upperBound;
    
public:
    const Eigen::VectorXd &lowerBound() const { return m_lowerBound; }
    void setLowerBound(const Eigen::VectorXd &lb) { m_lowerBound = lb; }
    const Eigen::VectorXd &upperBound() const { return m_upperBound; }
    void setUpperBound(const Eigen::VectorXd &ub) { m_upperBound = ub; }
    
    void setBoxConstraint(Eigen::VectorXd  lb, Eigen::VectorXd  ub) {
        setLowerBound(lb);
        setUpperBound(ub);
    }
    
    bool isValid(const Eigen::VectorXd &x)
    {
        return ((x - m_lowerBound).array() >= 0.0).all() && ((x - m_upperBound).array() <= 0.0).all();
    }
    
    double maxStep(Eigen::VectorXd L, Eigen::VectorXd S, Eigen::VectorXd dir)
    {
        Eigen::VectorXd x(L.size() + S.size());
        x.segment(0, L.size()) = L;
        x.segment(L.size(), S.size()) = S;
        if(!isValid(x))
        {
            std::cout<<"infeasible!"<<std::endl;
        }
        double step = 1e10;
        for(int i=0;i<S.size();i++)
        {
            double li = m_lowerBound(L.size() + i);
            double ui = m_upperBound(L.size() + i);
            double di = dir(L.size() + i);
            double curStep = 1e10;
            if(di>0)
            {
                curStep = ( ui - S(i) ) / di;
            }
            if(di < 0)
            {
                curStep = ( li - S(i) ) / di;
            }
            
            if(curStep < step)
            {
                step = curStep;
            }
        }
        return step;
    }

};
    


#endif
