#ifndef SIMULATIONSETUPDYNAMICSOLVER_H
#define SIMULATIONSETUPDYNAMICSOLVER_H

#include "SimulationSetup.h"
#include "../SensitiveAnalysis/SensitiveAnalysisAbarPos.h"
#include "../SensitiveAnalysis/SensitiveAnalysisABbarPos.h"
#include "../SensitiveAnalysis/SensitiveAnalysisAbarBbar.h"
#include "../SensitiveAnalysis/SensitiveAnalysis.h"

class SimulationSetupDynamicSolver : public SimulationSetup
{
public:
    virtual void buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff) override;
    void findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff);
    bool loadAbars() override;   // The path is given by abarPath + "L_list.dat".  abar = L*L^T
    void saveAbars(Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd pos, bool is_target);
    
    void testValueAndGradient();
    double testFunc(Eigen::VectorXd x, Eigen::VectorXd *grad);
    void testLBFGS();
    bool testLineSearch(Eigen::VectorXd x, Eigen::VectorXd dir, double &rate);
    
private:
    bool lineSearch(std::shared_ptr<SensitiveAnalysis> op, Eigen::VectorXd L, Eigen::VectorXd S, Eigen::MatrixXd &Pos, Eigen::VectorXd dir, double &rate);
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


