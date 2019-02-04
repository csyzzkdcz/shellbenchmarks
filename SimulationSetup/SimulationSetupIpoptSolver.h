
#ifndef SIMULATIONSETUPIPOPTSOLVER_H
#define SIMULATIONSETUPIPOPTSOLVER_H

#include "SimulationSetup.h"
#include "IfoptSolver.h"

class SimulationSetupIpoptSolver : public SimulationSetup
{
public:
    virtual void buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff) override;
    void findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff);
    
private:
    std::vector<Eigen::Matrix<double, 4, 9> > aderivs;
    std::vector<Eigen::MatrixXd > bderivs;
    std::vector<Eigen::Matrix2d> atargets;            // The first fundamental form of the target shape
    std::vector<Eigen::Matrix2d> btargets;            // The first fundamental form of the target shape
    double lameAlpha;
    double lameBeta;
    
};


#endif
