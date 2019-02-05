 #include <ifopt/problem.h>
 #include <ifopt/ipopt_solver.h>
 #include <ifopt/test_vars_constr_cost.h>
#include <igl/writeOBJ.h>
#include <iomanip>


#include "SimulationSetupIpoptSolver.h"
#include "../GeometryDerivatives.h"
#include "../ElasticShell.h"

bool SimulationSetupIpoptSolver::loadAbars()
{
    std::ifstream infile(abarPath);
    if(!infile)
        return false;
    int num;
    infile >> num;
    double d1, d2, d3;
    for(int i = 0; i < abars.size(); i++)
    {
        infile >> d1 >> d2 >> d3;
        Eigen::Matrix2d L;
        L << d1,0,
        d2,d3;
        abars[i] = L*L.transpose();
    }
    return true;
}

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
    if(!loadAbars())
        findFirstFundamentalForms(sff);
}

void SimulationSetupIpoptSolver::findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
    int nfaces = mesh.nFaces();
    int nverts = initialPos.rows();
    
    using namespace ifopt;
    Problem nlp;
    auto variable_set = std::make_shared<optVariables>(3*(nfaces+nverts), mesh, targetPos, "var_set");
    auto constraint_set = std::make_shared<optConstraint>(3*nverts, mesh, lameAlpha, lameBeta, thickness, "constraint");
    auto cost_term = std::make_shared<optCost>(initialPos, targetPos, mesh, penaltyCoef);
    
//    Eigen::VectorXd x1 = variable_set->GetValues();
//
//    constraint_set->testValueJacobian(x1);
//    cost_term->testCostJacobian(x1);
//    return;
    
    nlp.AddVariableSet(variable_set);
    nlp.AddConstraintSet(constraint_set);
    nlp.AddCostSet(cost_term);
    nlp.PrintCurrent();
    
    IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "mumps");
    ipopt.SetOption("jacobian_approximation", "exact");
    ipopt.SetOption("max_cpu_time", 1e6);
    ipopt.SetOption("tol", 1e-14);
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
    
    Eigen::MatrixXd curPos(nverts,3);
    Eigen::VectorXd derivative;
    
    
    for(int i=0; i<nverts; i++)
    {
        curPos.row(i) = x.segment<3>(3*i);
    }
    
    double energy = elasticEnergy(mesh, curPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, NULL);
    
    std::cout<<std::setiosflags(std::ios::fixed)<<std::setprecision(16)<<energy<<" "<<derivative.norm()<<std::endl;
    
    std::ofstream outfile(abarPath, std::ios::trunc);
    std::cout<<abarPath<<std::endl;
    outfile<<3*nfaces<<"\n";
    for(int i=0;i<3*nfaces;i++)
    {
        outfile<<std::setprecision(16)<<x(3*nverts + i)<<"\n";
    }
    outfile<<std::setprecision(16)<<x(x.size()-1);
    outfile.close();
    
    
}
