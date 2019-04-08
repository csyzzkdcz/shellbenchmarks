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
    double value(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos) override;
    void gradient(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd &grad) override;
    void projectBack(Eigen::VectorXd curL, Eigen::VectorXd &curS, Eigen::MatrixXd &curPos) override;
    void projectS(Eigen::VectorXd curL, Eigen::VectorXd &curS, Eigen::MatrixXd curPos);
    void projectPos(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd &curPos);
    
    void computeConvertedGrad(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd grad, Eigen::VectorXd &convertedGrad);
    
    void initialization(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, std::map<int, double> clampedDOFs, double lameAlpha, double lameBeta, double thickness) override
    {
       
        generalInitialization(initialPos, tarPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
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

    }
    
public:
    // Test functions
    void testConvertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    void testAbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    void testBbarAndSConvertion(Eigen::Matrix2d abar, Eigen::Matrix2d bbar);
    
//    void testAbarinvBbar2s(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos);
    
    void test() override;
    
    Eigen::VectorXd getInplaneForce();
    
    Eigen::Vector3d convertBar2s(Eigen::Matrix2d abar, Eigen::Matrix2d bbar);
    Eigen::Matrix2d converts2Bbar(Eigen::Matrix2d abar, Eigen::Vector3d s);
    
public:
    void convertEquilibrium2MatrixForm(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd *C, Eigen::SparseMatrix<double> &W);
    // convert F(abars, bbars) = 0 to W*s = C

    void computeAbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, std::vector<Eigen::Triplet<double> > *grad);
    
    void computeBbarDerivative(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::SparseMatrix<double> *grad);
    
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
    
    

private:
    Eigen::SparseMatrix<double> A;
    std::vector<Eigen::Matrix2d> aList;
    std::vector<Eigen::Matrix2d> bList;
    std::vector<Eigen::Matrix<double, 4, 9>> aderivs;
    std::vector<Eigen::MatrixXd> bderivs;
    Eigen::VectorXd edgeDOFs;


};

#endif
