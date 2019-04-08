#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/test_vars_constr_cost.h>
#include <igl/writeOBJ.h>
#include <igl/edge_lengths.h>
#include <iomanip>
#include <fstream>

#include "SimulationSetupIpoptSolver.h"
#include "../GeometryDerivatives.h"
#include "../ElasticShell.h"

bool SimulationSetupIpoptSolver::loadAbars()
{
    std::ifstream infile(abarPath);
    if(!infile)
        return false;
    infile >> thickness;
    infile >> abarCoef;
    infile >> smoothCoef;
    int vertNum, faceNum;
    infile >> vertNum >> faceNum;

    double d;
    Eigen::VectorXd x(vertNum + faceNum);
    x.setZero();
    for(int i=0;i<vertNum; i++)
    {
        infile >> d;
        x(i) = d;
    }
    double d1, d2, d3;
    for(int i = 0; i < abars.size(); i++)
    {
        infile >> d1 >> d2 >> d3;
        Eigen::Matrix2d L;
        L << d1,0,
        d2,d3;
        abars[i] = L*L.transpose();

        x(vertNum + 3*i) = d1;
        x(vertNum + 3*i + 1) = d2;
        x(vertNum + 3*i + 2) = d3;
    }

    auto cost_term = std::make_shared<ifopt::optCost>(initialPos, targetPos, mesh, abarCoef);
    double differenceValue = cost_term->getDifference(x);
    double penaltyValue = cost_term->getPenalty(x);
    double smoothnessValue = cost_term->getSmoothness(x);

    std::cout<<"Penalty: "<<penaltyValue<<"Abar Penalty Coefficient: "<<abarCoef<<std::endl;
    std::cout<<"Smoothness: "<<smoothnessValue<<" Smootheness Coefficient: "<<smoothCoef<<std::endl;
    std::cout<<"Shape Changed Value: "<<differenceValue<<" Thickness: "<<thickness<<std::endl;

    int nverts = vertNum / 3;
    Eigen::MatrixXd curPos(nverts,3);
    Eigen::VectorXd derivative;
    
    for(int i=0; i<nverts; i++)
    {
        curPos.row(i) = x.segment<3>(3*i);
    }
    
    Eigen::MatrixXd L(mesh.nFaces(), 3);
    igl::edge_lengths(curPos, mesh.faces(), L);
    
    
    MidedgeAverageFormulation sff;
    
    std::vector<Eigen::Triplet<double> > hessianTriplet;
    
    double energy = elasticEnergy(mesh, curPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, &hessianTriplet);
    
    Eigen::SparseMatrix<double> resH(3*nverts, 3*nverts);
    Eigen::SparseMatrix<double> H(3*nverts, 3*nverts);
    H.setFromTriplets(hessianTriplet.begin(), hessianTriplet.end());
    
    resH.setIdentity();
    resH = 1e-8 * resH + H;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(resH);
    
    if(solver.info() == Eigen::ComputationInfo::Success )
    {
        std::cout<<"Rearch the local minimun at current state"<<std::endl;
    }
    else
    {
        std::cout<<"Saddle point"<<std::endl;
    }
    
    std::cout<<"Energy with current positions: "<<energy<<" Derivative with current positions: "<<derivative.norm()<<std::endl;

    energy = elasticEnergy(mesh, targetPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, NULL);
    std::cout<<"Energy with target positions: "<<energy<<" Derivative with target positions: "<<derivative.norm()<<std::endl;
    targetPosAfterFirstStep = curPos;
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
//    auto variable_set = std::make_shared<optVariables>(3*(nfaces+nverts), mesh, initialPos, initialPos, "var_set");
    auto variable_set = std::make_shared<optVariables>(3*(nfaces+nverts), mesh, targetPos, targetPos, "var_set");
    auto constraint_set = std::make_shared<optConstraint>(3*nverts, mesh, lameAlpha, lameBeta, thickness, "constraint");
    auto cost_term = std::make_shared<optCost>(initialPos, targetPos, mesh, abarCoef);
    cost_term->_mu = smoothCoef;
    Eigen::VectorXd x1 = variable_set->GetValues();
    Eigen::MatrixXd curPos(nverts,3);
    Eigen::VectorXd derivative;

//    constraint_set->testValueJacobian(x1);
//    cost_term->testCostJacobian(x1);
//    return;
    
    
    abars.resize(nfaces);
    
    for(int i=0; i< nfaces; i++)
    {
        abars[i] << x1(3*i + 3*nverts), 0,
        x1(3*i + 3*nverts + 1), x1(3*i + 3*nverts + 2);
        abars[i] = abars[i]*abars[i].transpose();
    }
    
    double energy = elasticEnergy(mesh, targetPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, NULL);
    
    double bendingTerm = 0;
    Eigen::VectorXd bendingDerivative(3*nverts);
    bendingDerivative.setZero();
    
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::VectorXd extraDOFs(0);
        
        Eigen::MatrixXd deriv(1, 18);
        bendingTerm += bendingEnergy(mesh, targetPos, extraDOFs, lameAlpha, lameBeta, thickness, abars[i], bbars[i], i, sff,  &deriv, NULL);
        for (int j = 0; j < 3; j++)
        {
            bendingDerivative.segment<3>(3 * mesh.faceVertex(i, j)).transpose() += deriv.block<1,3>(0, 3 * j);
            int oppidx = mesh.vertexOppositeFaceEdge(i, j);
            if(oppidx != -1)
                bendingDerivative.segment<3>(3 * oppidx).transpose() += deriv.block<1,3>(0, 9 + 3 * j);
        }
    }
    
    std::cout<<"The initial bending energy is: "<<bendingTerm<<" bending derivative is: "<<bendingDerivative.norm()<<std::endl;
    std::cout<<"The derivative before initialization is: "<<derivative.lpNorm<Eigen::Infinity>()<<std::endl;
    
    double penaltyValue = cost_term->getPenalty(x1);
    double differenceValue = cost_term->getDifference(x1);
    double smoothnessValue = cost_term->getSmoothness(x1);

    std::cout<<"Initial Penalty: "<<penaltyValue<<" Penalty Coefficient: "<<abarCoef<<std::endl;
    std::cout<<"Initial Smoothness: "<<smoothnessValue<<" Smoothness coefficient: "<<smoothCoef<<std::endl;
    std::cout<<"Initial Difference: "<<differenceValue<<" Thickness: "<<thickness<<std::endl;
    
    nlp.AddVariableSet(variable_set);
    nlp.AddConstraintSet(constraint_set);
    nlp.AddCostSet(cost_term);
    nlp.PrintCurrent();
    
    IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "mumps");
    // ipopt.SetOption("linear_solver", "ma97");
    ipopt.SetOption("jacobian_approximation", "exact");
    ipopt.SetOption("max_cpu_time", 1e6);
    ipopt.SetOption("tol", std::min(1e-8, 1e-2*std::min(derivative.lpNorm<Eigen::Infinity>(), bendingDerivative.lpNorm<Eigen::Infinity>())));
    ipopt.SetOption("print_level", 5);
    ipopt.SetOption("max_iter", int(2e4));
    
    
    std::cout.precision(9);
    std::cout<<std::scientific;
    
//    if(derivative.lpNorm<Eigen::Infinity>() < 1e-14)
//    {
//        std::cout<<"Derivative is almost 0: "<<derivative.lpNorm<Eigen::Infinity>()<<std::endl;
//    }
//    else
    {
        ipopt.Solve(nlp);
    }
    
    Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

    penaltyValue = cost_term->getPenalty(x);
    differenceValue = cost_term->getDifference(x);
    smoothnessValue = cost_term->getSmoothness(x);


    std::cout<<"Final Penalty: "<<penaltyValue<<" Penalty Coefficient: "<<abarCoef<<std::endl;
    std::cout<<"Final Smoothness: "<<smoothnessValue<<" Smoothness Coefficient: "<<smoothCoef<<std::endl;
    std::cout<<"Final Difference: "<<differenceValue<<" Thickness: "<<thickness<<std::endl;

    
    abars.resize(nfaces);
    
    for(int i=0; i< nfaces; i++)
    {
        abars[i] << x(3*i + 3*nverts), 0,
        x(3*i + 3*nverts + 1), x(3*i + 3*nverts + 2);
        abars[i] = abars[i]*abars[i].transpose();
    }
    
    
    for(int i=0; i<nverts; i++)
    {
        curPos.row(i) = x.segment<3>(3*i);
    }
    
    std::vector<Eigen::Triplet<double> > hessianTriplet;
    
    energy = elasticEnergy(mesh, curPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, &hessianTriplet);
    
    Eigen::SparseMatrix<double> H(3*nverts, 3*nverts);
    H.setFromTriplets(hessianTriplet.begin(), hessianTriplet.end());
    
    Eigen::SparseMatrix<double> resH(3*nverts, 3*nverts);
    resH.setIdentity();
    resH = 1e-8 * resH + H;
    
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(resH);
    
    if(solver.info() == Eigen::ComputationInfo::Success )
    {
        std::cout<<"Rearch the local minimun at current state"<<std::endl;
    }
    else
    {
        std::cout<<"Saddle point"<<std::endl;
    }
    
    std::cout<<"Energy with current position: "<<energy<<" Derivative with current positions: "<<derivative.norm()<<std::endl;
    
    energy = elasticEnergy(mesh, targetPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, NULL);
    std::cout<<"Energy with target positions: "<<energy<<" Derivative with target positions: "<<derivative.norm()<<std::endl;

    std::cout<<"Saving Abar path: "<<abarPath<<std::endl;
    std::ofstream outfile(abarPath, std::ios::trunc);
    
    outfile<<thickness<<"\n";
    outfile<<abarCoef<<"\n";
    outfile<<smoothCoef<<"\n";
    outfile<<3*nverts<<"\n";
    outfile<<3*nfaces<<"\n";
    
    std::cout<<3*nverts + 3*nfaces<<std::endl;
    std::cout<<x.size()<<std::endl;

    for(int i=0;i<3*nverts;i++)
    {
        outfile<<std::setprecision(16)<<x(i)<<"\n";
    }

    for(int i=0;i<3*nfaces;i++)
    {
        outfile<<std::setprecision(16)<<x(3*nverts + i)<<"\n";
    }
    outfile<<std::setprecision(16)<<x(x.size()-1);
    outfile.close();

    int startIdx, endIdx, expCoef;
    std::string subString = "";
    std::string resampledPath = abarPath;
    
    startIdx = resampledPath.rfind("/");
    endIdx = resampledPath.find("_");
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx + 1,resampledPath.begin() + endIdx, "resampled");
    
    // thickness
    if(thickness == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(thickness));
    startIdx = resampledPath.rfind("T");
    endIdx = resampledPath.rfind("P");
    subString = "";
    if(thickness > 0)
        subString = "T_1e" + std::to_string(expCoef);
    else
        subString = "T_0";
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    // penalty
    if(abarCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(abarCoef));
    
    startIdx = resampledPath.rfind("P");
    endIdx = resampledPath.rfind("S");
    subString = "";
    if(abarCoef > 0)
        subString = "P_1e" + std::to_string(expCoef);
    else
        subString = "P_0";
    resampledPath= resampledPath .replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    // smoothness
    if(smoothCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(smoothCoef));
    
    startIdx = resampledPath.rfind("S");
    endIdx = resampledPath.rfind(".");
    subString = "";
    if(smoothCoef > 0)
        subString = "S_1e" + std::to_string(expCoef);
    else
        subString = "S_0";
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx, subString);
    
    startIdx = resampledPath.rfind(".");
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.end(), ".obj");
    std::cout<<"Current abar loading path is: "<<resampledPath<<std::endl;
    igl::writeOBJ("resampled.obj", curPos, mesh.faces());
    igl::writeOBJ(resampledPath, curPos, mesh.faces());
    targetPosAfterFirstStep = curPos;
}
