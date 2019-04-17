#ifndef SIMULATIONSETUPDYNAMICSOLVER_H
#define SIMULATIONSETUPDYNAMICSOLVER_H

#include "SimulationSetup.h"
#include "../SensitiveAnalysis/SensitiveAnalysisABbarPos.h"
#include "../SensitiveAnalysis/SensitiveAnalysis.h"

class SimulationSetupDynamicSolver : public SimulationSetup
{
public:
    virtual void buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff) override;
    void findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff);
    bool loadParams() override;
    void saveParams();
    
    void testValueAndGradient();
    double testFunc(Eigen::VectorXd x, Eigen::VectorXd *grad);
    void testLBFGS();
    bool testLineSearch(Eigen::VectorXd x, Eigen::VectorXd dir, double &rate);
    
private:
    void testProjectBackSim();
    
    bool cutoffL(Eigen::VectorXd &curL, std::set<int> &cutAbar); // cutAbar stores the faces where abar is cut of\

private:
    std::vector<Eigen::Matrix<double, 4, 9> > aderivs;
    std::vector<Eigen::MatrixXd > bderivs;
    std::vector<Eigen::Matrix2d> atargets;            // The first fundamental form of the target shape
    std::vector<Eigen::Matrix2d> btargets;            // The first fundamental form of the target shape
    double lameAlpha;
    double lameBeta;

};



#endif


