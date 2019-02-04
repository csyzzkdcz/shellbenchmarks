 #include <ifopt/problem.h>
 #include <ifopt/ipopt_solver.h>
 #include <ifopt/test_vars_constr_cost.h>


#include "SimulationSetupIpoptSolver.h"
#include "../GeometryDerivatives.h"

void SimulationSetupIpoptSolver::buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
    int nfaces = mesh.nFaces();
    abars.resize(nfaces);
    bbars.resize(nfaces);
    
    initialAbars.resize(nfaces);
    
    lameAlpha = YoungsModulus * PoissonsRatio / (1.0 - PoissonsRatio * PoissonsRatio);
    lameBeta = YoungsModulus / 2.0 / (1.0 + PoissonsRatio);
    
    for (int i = 0; i < nfaces; i++)
    {
        bbars[i].setZero();
        initialAbars[i] = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
    }
    findFirstFundamentalForms(sff);
}

void SimulationSetupIpoptSolver::findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
    int nfaces = mesh.nFaces();
    int nverts = initialPos.rows();
    
    using namespace ifopt;
    Problem nlp;
    auto variable_set = std::make_shared<optVariables>(3*(nfaces+nverts), mesh, initialPos, "var_set");
    auto constraint_set = std::make_shared<optConstraint>(3*nverts, mesh, lameAlpha, lameBeta, thickness, "constraint");
    auto cost_term = std::make_shared<optCost>(initialPos, targetPos, mesh, penaltyCoef);
    
    nlp.AddVariableSet(variable_set);
    nlp.AddConstraintSet(constraint_set);
    nlp.AddCostSet(cost_term);
    nlp.PrintCurrent();
    
    IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "mumps");
    ipopt.SetOption("jacobian_approximation", "exact");
    ipopt.SetOption("max_cpu_time", 1e6);
    ipopt.SetOption("tol", 1e-10);
    ipopt.SetOption("print_level", 5);
    
    ipopt.Solve(nlp);
    
    Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
    
    abars.resize(nfaces);
    
    for(int i=0; i< nfaces; i++)
    {
        abars[i] << x(3*i + 3*nverts), 0,
        x(3*i + 3*nverts + 1), x(3*i + 3*nverts + 2);
        abars[i] = abars[i]*abars[i].transpose();
    }
    
    
}
