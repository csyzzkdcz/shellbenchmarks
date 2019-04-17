#ifndef SENSITIVEANALYSISABBARPOS_H
#define SENSITIVEANALYSISABBARPOS_H

#include <set>
#include "SensitiveAnalysis.h"

class SensitiveAnalysisABbarPos : public SensitiveAnalysis
/*
 abar = L * L^T
 bbar = s * abar
 */
{
public:
    double value(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos) override;
    void gradient(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd &grad) override;
    void projectBack(Eigen::VectorXd curL, Eigen::VectorXd &curS, Eigen::MatrixXd &curPos) override;
    void projectS(Eigen::VectorXd curL, Eigen::VectorXd &curS, Eigen::MatrixXd curPos);
    void projectPos(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd &curPos);
    
    void computeConvertedGrad(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd grad, Eigen::VectorXd &convertedGrad);
    
    void initialization(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, std::map<int, double> clampedDOFs, double lameAlpha, double lameBeta, double thickness) override
    {
        
        generalInitialization(initialPos, tarPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
        std::vector<Eigen::Triplet<double> >T, Ts;
        int nfaces = _mesh.nFaces();
        
        MidedgeAverageFormulation sff;
        
        edgeDOFs.resize(0);
        
        for(int i=0;i<nfaces; i++)
        {
            T.push_back(Eigen::Triplet<double>(3*i+1, 3*i+1, 1));
            T.push_back(Eigen::Triplet<double>(3*i+2, 3*i+2, 1));
            
        }
        
        A.resize(3*nfaces, 3*nfaces);
        A.setFromTriplets(T.begin(), T.end());
        
        updateBoxConstraint();
      
        
        T.clear();
        Ts.clear();
        lapFaceMat.resize(nfaces, nfaces);
        Ws.resize(nfaces, nfaces);
        for(int i=0;i<nfaces;i++)
        {
            double totalNeiborArea = 0;
            for(int j=0;j<3;j++)
            {
                int oppFace = _mesh.faceVertex(i,j);
                if(oppFace!=-1)
                {
                    totalNeiborArea += _areaList(oppFace);
                }
            }
            for(int j=0;j<3;j++)
            {
                int oppFace = _mesh.faceVertex(i,j);
                if(oppFace!=-1)
                {
                    T.push_back(Eigen::Triplet<double>(i, oppFace, _areaList(oppFace) / totalNeiborArea));
                }
            }
            
            T.push_back(Eigen::Triplet<double>(i,i, -1));
            Ts.push_back(Eigen::Triplet<double>(i,i, 1.0 / _areaList(i)));
        }
        lapFaceMat.setFromTriplets(T.begin(), T.end());
        Ws.setFromTriplets(Ts.begin(), Ts.end());
        
        projM.resize(4*nfaces, 4*nfaces);
        projM.setIdentity();
        
        fixedVariables.resize(4*nfaces);
        fixedVariables.setZero();
    }
    
public:
    // Test functions
    void testConvertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    void testAbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    void testBbarAndSConvertion(Eigen::Matrix2d abar, Eigen::Matrix2d bbar);
    void testBbarSmoothAndGrad(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    
    //    void testAbarinvBbar2s(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    
    void test() override;
    
    Eigen::VectorXd getInplaneForce();
    
    Eigen::Vector3d convertBar2s(Eigen::Matrix2d abar, Eigen::Matrix2d bbar);
    Eigen::Matrix2d converts2Bbar(Eigen::Matrix2d abar, Eigen::Vector3d s);
    
    void setProjM(std::set<int> fixedFlags);
    void updateFixedVariables(Eigen::VectorXd L, Eigen::VectorXd S);
    Eigen::VectorXd getFullVariables(Eigen::VectorXd reductVariables);
    
public:
    void convertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd *C, Eigen::SparseMatrix<double> &W);
    // convert F(abars, bbars) = 0 to W*s = C
    
    void computeAbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, std::vector<Eigen::Triplet<double> > *grad);
    
    void computeBbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::SparseMatrix<double> *grad);
    
    double computeBbarSmoothness(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    
    void computeBbarSmoothnessGrad(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd &grad);
    
    Eigen::Matrix2d computeAbarinvBbar(Eigen::Vector3d L, Eigen::Vector3d s);
    
    void computeAbarinvBbar2s(Eigen::Matrix2d L, Eigen::Vector3d s, Eigen::Matrix<double, 4, 3> &dAbarinvBbar)
    {
        Eigen::Matrix2d dA1,dA2,dA3,M;
        M << 0,1,0,0;
        double x,y,z;
        x = L(0,0);
        y = L(1,0);
        z = L(1,1);
        dA1 = -(2 * (y*y + z*z)*s(1) / (x*x*x) + 2 * y * s(2) / (x*x)) * M;
        dA2 = (2 * y / (x*x) * s(1) + 2 / x * s(2)) * M;
        dA3 = (2 * z / (x*x) * s(1)) * M;
        
        //        std::cout<<dA1<<std::endl<<std::endl;
        //        std::cout<<dA2<<std::endl<<std::endl;
        //        std::cout<<dA3<<std::endl<<std::endl;
        
        dAbarinvBbar.row(0) << dA1(0,0), dA2(0,0), dA3(0,0);
        dAbarinvBbar.row(1) << dA1(0,1), dA2(0,1), dA3(0,1);
        dAbarinvBbar.row(2) << dA1(1,0), dA2(1,0), dA3(1,0);
        dAbarinvBbar.row(3) << dA1(1,1), dA2(1,1), dA3(1,1);
    }
    
    
    
public:
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd edgeDOFs;
    Eigen::SparseMatrix<double> lapFaceMat;
    Eigen::SparseMatrix<double> Ws;
    Eigen::SparseMatrix<double> projM;
    
    Eigen::VectorXd fixedVariables;
    
// Box Constraint
    
public:
    const Eigen::VectorXd &lowerBound() const { return m_lowerBound; }
    void setLowerBound(const Eigen::VectorXd &lb) { m_lowerBound = lb; }
    const Eigen::VectorXd &upperBound() const { return m_upperBound; }
    void setUpperBound(const Eigen::VectorXd &ub) { m_upperBound = ub; }
    
    void setBoxConstraint(Eigen::VectorXd  lb, Eigen::VectorXd  ub) {
        setLowerBound(lb);
        setUpperBound(ub);
    }
    
    void updateBoxConstraint()
    {
        int nfaces = _mesh.nFaces();
        Eigen::VectorXd infBound(4 * nfaces);  // 3*F = L.size, F = S.size;
        infBound.setConstant(std::numeric_limits<double>::infinity());
        m_lowerBound = -infBound;
        m_upperBound = infBound;
        
        for(int i=0;i<nfaces;i++)
        {
            Eigen::Matrix2d abar = _atarget[i];
            Eigen::Matrix2d TMat;
            TMat.col(0) = ( _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
            TMat.col(1) = ( _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
            abar = (TMat.transpose()).inverse() * abar * TMat.inverse();
            
            m_lowerBound(3*nfaces + i) = -abar.trace() / (10 * _thickness) * _areaList(i);
            m_upperBound(3*nfaces + i) = abar.trace() / (10 * _thickness) * _areaList(i);
        }
    }
    
    bool isValid(const Eigen::VectorXd x)
    {
        Eigen::VectorXd fullx = getFullVariables(x);
        return ((fullx - m_lowerBound).array() >= 0.0).all() && ((fullx - m_upperBound).array() <= 0.0).all();
    }
    
    double maxStep(Eigen::VectorXd x, Eigen::VectorXd dir)
    {
        Eigen::VectorXd reductlowerBound, reductupperBound;
        reductlowerBound = projM * m_lowerBound;
        reductupperBound = projM * m_upperBound;
        
        if(!isValid(x))
        {
            std::cout<<"infeasible!"<<std::endl;
        }
        double step = 1e10;
        for(int i=0;i<x.size();i++)
        {
            double li = reductlowerBound(i);
            double ui = reductupperBound(i);
            double di = dir(i);
            double curStep = 1e10;
            if(di>0)
            {
                curStep = ( ui - x(i) ) / di;
            }
            if(di < 0)
            {
                curStep = ( li - x(i) ) / di;
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
