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
    double abarCoef;     // default value is 0
    double bbarCoef;     // default value is 0
    double smoothCoef;      // default value is 0
    
    std::string abarPath;
    std::string resamplingPath;
    std::vector<PostprocessTest> tests;

    // Derived from the above
    std::vector<Eigen::Matrix2d> abars;
    std::vector<Eigen::Matrix2d> bbars;

    virtual void buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff) = 0;
    virtual bool loadAbars() = 0;   // The path is given by abarPath + "L_list.dat".  abar = L*L^T
    
    void remeshProcessing(Eigen::MatrixXd remeshedPos, Eigen::MatrixXi remeshedFaces);
    
    std::string selectedDynamicType;
    
    bool _is_overwrite;
    bool _is_continue;
    
    // Some operation
    Eigen::Vector3d abar2L(Eigen::Matrix2d abar)
    {
        Eigen::Vector3d L;
        L(0) = sqrt(abar(0,0));
        L(1) = abar(0,1)/L(0);
        L(2) = sqrt(abar.determinant())/L(0);
        return L;
    }
    
    Eigen::Matrix2d L2abar(Eigen::Vector3d L)
    {
        Eigen::Matrix2d abar;
        abar << L(0), 0, L(1), L(2);
        abar = abar * abar.transpose();
        return abar;
    }
};



#endif
