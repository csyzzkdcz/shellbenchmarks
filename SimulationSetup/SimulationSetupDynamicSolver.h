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
    void saveAbars(Eigen::VectorXd L, Eigen::MatrixXd pos);
    
    void testValueAndGradient();

private:
    bool lineSearch(std::shared_ptr<SensitiveAnalysis> op, Eigen::VectorXd L, Eigen::MatrixXd &Pos, Eigen::VectorXd dir, double &rate);
    void testProjectBackSim();

private:
    std::vector<Eigen::Matrix<double, 4, 9> > aderivs;
    std::vector<Eigen::MatrixXd > bderivs;
    std::vector<Eigen::Matrix2d> atargets;            // The first fundamental form of the target shape
    std::vector<Eigen::Matrix2d> btargets;            // The first fundamental form of the target shape
    double lameAlpha;
    double lameBeta;

};



#endif


