#ifndef SIMULATIONSETUP_H
#define SIMULATIONSETUP_H

#include "../MeshConnectivity.h"
#include <set>
#include <vector>
#include <map>

class SecondFundamentalFormDiscretization;

struct PostprocessTest
{
    int vertex;
    int coord;
    double wimDisplacement;
};

class SimulationSetup
{
    public:
    virtual ~SimulationSetup() = default;
    
    public:
    // Core data structures
    MeshConnectivity mesh; // mesh combinatorics
    Eigen::MatrixXd initialPos; // mesh vertices of the rest state
    Eigen::MatrixXd targetPos; //   mesh vertices of the target state
    Eigen::MatrixXd targetPosAfterFirstStep; 
    Eigen::VectorXd initialEdgeDOFs;
    std::map<int, double> clampedDOFs;
    Eigen::MatrixXd externalForces; // same size as initialPos
    
    std::vector<Eigen::Matrix2d> initialAbars; // Abars on the rest shape
    double thickness;
    double YoungsModulus;
    double PoissonsRatio;
    double penaltyCoef;     // default value is 0
    double smoothCoef;      // default value is 0
    
    std::string abarPath;
    std::vector<PostprocessTest> tests;

    // Derived from the above
    std::vector<Eigen::Matrix2d> abars;
    std::vector<Eigen::Matrix2d> bbars;

    virtual void buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff) = 0;
    virtual bool loadAbars() = 0;   // The path is given by abarPath + "L_list.dat".  abar = L*L^T
};



#endif
