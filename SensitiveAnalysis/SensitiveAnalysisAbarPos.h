#ifndef SENSITIVEANALYSISABARPOS_H
#define SENSITIVEANALYSISABARPOS_H

#include "SensitiveAnalysis.h"

class SensitiveAnalysisAbarPos : public SensitiveAnalysis
{
public:
    double value(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos) override;
    void gradient(Eigen::VectorXd curL, Eigen::VectorXd curS, Eigen::MatrixXd curPos, Eigen::VectorXd &grad) override;
    void projectBack(Eigen::VectorXd curL, Eigen::VectorXd &curS, Eigen::MatrixXd &curPos) override;
    
    Eigen::SparseMatrix<double> getConvertedGrad()
    {
        return _convertedGrad;
    }
    
public:
    void initialization(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, std::map<int, double> clampedDOFs, double lameAlpha, double lameBeta, double thickness) override
    {
        generalInitialization(initialPos, tarPos, mesh, clampedDOFs, lameAlpha, lameBeta, thickness);
        selectedCoord.resize(3);
        for(int i =0; i<3; i++)
        {
            selectedCoord[i] = computeSelectMatrix(initialPos.rows(), i);
        }
        //        computeDerivativeQ2Abr(_initialPos, initialL, _convertedGrad);
        std::vector<Eigen::Triplet<double>> proj;
        std::cout<<clampedDOFs.size()<<std::endl;
        int freeDOFs = 3 * initialPos.rows() - clampedDOFs.size();
        
        int row = 0;
        for(int i=0;i<initialPos.rows();i++)
        {
            for(int j=0;j<3;j++)
            {
                if (clampedDOFs.find(3*i+j) != clampedDOFs.end())
                    continue;
                proj.push_back(Eigen::Triplet<double>(row, 3*i+j, 1.0));
                row ++;
            }
        }
        _projM.resize(freeDOFs, 3*initialPos.rows());
        _projM.setFromTriplets(proj.begin(), proj.end());
    }
    
    void testAbarSmoothnessGrad(Eigen::VectorXd curL);
    void testPositionSmoothnessGrad(Eigen::MatrixXd curPos);
    void testDifferenceGrad(Eigen::MatrixXd curPos);
    void testAbarDerivative(Eigen::VectorXd curL, Eigen::MatrixXd curPos);
    void testDerivativeQ2Abr(Eigen::VectorXd curL, Eigen::MatrixXd curPos);
    
    void test() override;
    
public:
    void computeDerivativeQ2Abr(Eigen::MatrixXd curPos, Eigen::VectorXd curL, Eigen::SparseMatrix<double> &derivative); // aAbr = L*L^T, we consider the drivative for L
    void computeAbarDerivative(Eigen::VectorXd curL, Eigen::MatrixXd curPos, std::vector<Eigen::Triplet<double> > *grad);
    
//    double computeAbarSmoothness(Eigen::VectorXd curL);
    double computePositionSmoothness(Eigen::MatrixXd curPos);
    
    
//    void computeAbarSmoothnessGrad(Eigen::VectorXd curL, Eigen::VectorXd &grad);
    void computePositionSmoothnessGrad(Eigen::MatrixXd curPos, Eigen::VectorXd &grad);
    
    
    void computeConvertedGrad(Eigen::MatrixXd curPos, Eigen::VectorXd curL, Eigen::VectorXd grad, Eigen::VectorXd &convertedGrad);
    
    
    
public:
    Eigen::SparseMatrix<double> _projM; // We fixed 4 corner points
    
private:
    Eigen::VectorXd _prevL;
    Eigen::SparseMatrix<double> _convertedGrad;
    
    
    
    std::vector<Eigen::SparseMatrix<double>> selectedCoord;
    
    
};



#endif

