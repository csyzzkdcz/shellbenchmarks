#ifndef SIMULATIONSTATE_H
#define SIMULATIONSTATE_H

#include <Eigen/Core>

struct SimulationState
{
    Eigen::MatrixXd curPos;
    Eigen::VectorXd curEdgeDOFs;
};

#endif