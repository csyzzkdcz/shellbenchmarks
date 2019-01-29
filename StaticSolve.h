#ifndef STATICSOLVE_H
#define STATICSOLVE_H
#include <Eigen/Core>

class SimulationSetup;
class SimulationSetupNormal;
struct SimulationState;
class SecondFundamentalFormDiscretization;

void takeOneStep(const SimulationSetup &setup, SimulationState &state, const SecondFundamentalFormDiscretization &sff, double &reg, double interp,
    int &funcEvals, // number of times function was evaluated
    double &forceResidual, // nodal force residual ||F|| at beginning of step
    double &updateMag
);

void leadingEigenvector(const SimulationSetup &setup, SimulationState &state, const SecondFundamentalFormDiscretization &sff, Eigen::VectorXd &result);

#endif