#ifndef SENSITIVEANALYSISABARBBAR_H
#define SENSITIVEANALYSISABARBBAR_H

#include "SensitiveAnalysis.h"

class SensitiveAnalysisAbarBbar : public SensitiveAnalysis
/*
 abar = L * L^T
 bbar = s * abar
 */
{
public:
    double value(const Eigen::VectorXd &curL, const Eigen::MatrixXd curS) override;
    void gradient(const  Eigen::VectorXd curL, const Eigen::MatrixXd curS, Eigen::VectorXd &grad) override;
    void projectBack(Eigen::VectorXd curL, Eigen::MatrixXd &curS) override;
    void initialization(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, double lameAlpha, double lameBeta, double thickness)
    {
        _initialPos = initialPos;
        _tarPos = tarPos;
        _mesh = mesh;
        _lameAlpha = lameAlpha;
        _lameBeta = lameBeta;
        _thickness = thickness;
        
        std::vector<Eigen::Triplet<double> >T;
        int nfaces = _mesh.nFaces();
        aList.resize(nfaces);
        bList.resize(nfaces);
        aderivs.resize(nfaces);
        bderivs.resize(nfaces);
        
        MidedgeAverageFormulation sff;
        
        edgeDOFs.resize(0);
        
        for(int i=0;i<nfaces; i++)
        {
            T.push_back(Eigen::Triplet<double>(3*i+1, 3*i+1, 1));
            T.push_back(Eigen::Triplet<double>(3*i+2, 3*i+2, 1));
            
            aList[i] = firstFundamentalForm(_mesh, _tarPos, i, &aderivs[i], NULL);
            bderivs[i].resize(4, 18);
            bList[i] = sff.secondFundamentalForm(_mesh, _tarPos, edgeDOFs, i, &bderivs[i], NULL);
        }
        
        A.resize(3*nfaces, 3*nfaces);
        A.setFromTriplets(T.begin(), T.end());
        
        igl::cotmatrix(initialPos,mesh.faces(),_laplacianMat);
        
    }
private:
    void convertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd *C, Eigen::SparseMatrix<double> &W);
    // convert F(abars, bbars) = 0 to W*s = C
    
    void computeAbarDerivative(Eigen::VectorXd curL, Eigen::MatrixXd curS, std::vector<Eigen::Triplet<double> > *grad);
    
private:
    Eigen::SparseMatrix<double> A;
    std::vector<Eigen::Matrix2d> aList;
    std::vector<Eigen::Matrix2d> bList;
    std::vector<Eigen::Matrix<double, 4, 9>> aderivs;
    std::vector<Eigen::MatrixXd> bderivs;
    Eigen::VectorXd edgeDOFs;

    
};

#endif
