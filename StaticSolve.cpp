#include "StaticSolve.h"
#include "SimulationSetup/SimulationSetup.h"
#include "SimulationSetup/SimulationSetupNormal.h"
#include "SimulationSetup/SimulationSetupAlglibSolver.h"
#include "SimulationState.h"
#include "ElasticShell.h"
#include "GeometryDerivatives.h"
#include <iostream>
#include <igl/cotmatrix.h>



void takeOneStep(const SimulationSetup &setup, SimulationState &state, const SecondFundamentalFormDiscretization &sff, double &reg, double interp,
    int &funcEvals,
    double &forceResidual,
    double &updateMag
)
{
    double lameAlpha = setup.YoungsModulus * setup.PoissonsRatio / (1.0 - setup.PoissonsRatio * setup.PoissonsRatio);
    double lameBeta = setup.YoungsModulus / 2.0 / (1.0 + setup.PoissonsRatio);
    
    int nverts = state.curPos.rows();
    int nedges = setup.mesh.nEdges();
    int nedgedofs = sff.numExtraDOFs();

    int constrainedDOFs = setup.clampedDOFs.size();

    int freeDOFs = 3 * nverts + nedgedofs * nedges - constrainedDOFs;
    std::vector<Eigen::Triplet<double> > proj;
    int row = 0;
    for (int i = 4; i < nverts; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (setup.clampedDOFs.find(3*i+j) != setup.clampedDOFs.end())
                continue;
            proj.push_back(Eigen::Triplet<double>(row, 3 * i + j, 1.0));
            row++;
        }
    }
    for (int i = 0; i < nedges; i++)
    {
        for (int j = 0; j < nedgedofs; j++)
        {
            proj.push_back(Eigen::Triplet<double>(row, 3 * nverts + nedgedofs * i + j, 1.0));
            row++;
        }
    }
    assert(row == freeDOFs);
    Eigen::SparseMatrix<double> projM(freeDOFs, 3 * nverts + nedgedofs * nedges);
    projM.setFromTriplets(proj.begin(), proj.end());

    // enforce constrained DOFs
//    for (auto &it : setup.clampedDOFs)
//    {
//        int vid = it.first / 3;
//        int coord = it.first % 3;
//        state.curPos(vid, coord) = interp * it.second + (1.0 - interp) * setup.initialPos(vid, coord);
//    }
    
    // combine the optimized abar with the initial abars
    std::vector<Eigen::Matrix2d> curAbars(setup.abars.size());
    for(int i = 0; i < setup.abars.size(); i++)
    {
        Eigen::Matrix2d abar = firstFundamentalForm(setup.mesh, setup.initialPos, i, NULL, NULL);
        curAbars[i] = interp * setup.abars[i] + (1 - interp) * abar;
    }

    while (true)
    {
        Eigen::VectorXd derivative;
        std::vector<Eigen::Triplet<double> > hessian;
        
//        double energy = elasticEnergy(setup.mesh, state.curPos, state.curEdgeDOFs, lameAlpha, lameBeta, setup.thickness, setup.abars, setup.bbars, sff, &derivative, &hessian);
        double energy = elasticEnergy(setup.mesh, state.curPos, state.curEdgeDOFs, lameAlpha, lameBeta, setup.thickness, curAbars, setup.bbars, sff, &derivative, &hessian);

        funcEvals++;
        // work of external forces
//        for (int i = 0; i < nverts; i++)
//        {
//            energy -= setup.externalForces.row(i).dot(state.curPos.row(i) - setup.initialPos.row(i));
//        }
        
        Eigen::SparseMatrix<double> H(3 * nverts + nedgedofs * nedges, 3 * nverts + nedgedofs * nedges);
        H.setFromTriplets(hessian.begin(), hessian.end());


        Eigen::VectorXd force = -derivative;
//        for (int i = 0; i < nverts; i++)
//            force.segment<3>(3 * i) += setup.externalForces.row(i).transpose();
//        std::cout<<force.norm()<<std::endl;
        Eigen::VectorXd reducedForce = projM * force;
        Eigen::SparseMatrix<double> reducedH = projM * H * projM.transpose();
        Eigen::SparseMatrix<double> I(freeDOFs, freeDOFs);
        I.setIdentity();
        reducedH += reg * I;                
        
        //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
        //cg.compute(reducedH);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(reducedH);
        while(solver.info() != Eigen::ComputationInfo::Success )
        {
            reg *= 2.0;
            std::cout << "The Hessian is not positive definite " << energy << " lambda now: " << reg << std::endl;
            reducedH += reg * I;
            solver.compute(reducedH);
        }
        Eigen::VectorXd descentDir = solver.solve(reducedForce);
        //std::cout << "Solver residual: " << (reducedH*descentDir - reducedForce).norm() << std::endl;
        Eigen::VectorXd fullDir = projM.transpose() * descentDir;


        Eigen::MatrixXd newPos = state.curPos;
        for (int i = 0; i < nverts; i++)
        {
            newPos.row(i) += fullDir.segment<3>(3 * i);
        }    
        Eigen::VectorXd newEdgeDofs = state.curEdgeDOFs + fullDir.segment(3 * nverts, nedgedofs * nedges);



//        double newenergy = elasticEnergy(setup.mesh, newPos, newEdgeDofs, lameAlpha, lameBeta, setup.thickness, setup.abars, setup.bbars, sff, &derivative, NULL);
        double newenergy = elasticEnergy(setup.mesh, newPos, state.curEdgeDOFs, lameAlpha, lameBeta, setup.thickness, curAbars, setup.bbars, sff, &derivative, &hessian);
        funcEvals++;
        force = -derivative;

//        for (int i = 0; i < nverts; i++)
//        {
//            newenergy -= setup.externalForces.row(i).dot(newPos.row(i) - setup.initialPos.row(i));
//            force.segment<3>(3 * i) += setup.externalForces.row(i).transpose();
//        }
        reducedForce = projM * force;
        forceResidual = reducedForce.norm();
        updateMag = fullDir.norm();
        if (newenergy <= energy)
        {
            std::cout << "Old energy: " << energy << " new energy: " << newenergy << " force residual " << forceResidual << " pos change " << fullDir.segment(0, 3 * nverts).norm() << " theta change " << fullDir.segment(3 * nverts, nedgedofs*nedges).norm() << std::endl;
            state.curPos = newPos;
            state.curEdgeDOFs = newEdgeDofs;
            reg = std::max(1e-6, reg * 0.5);
            break;
        }
        else
        {
            reg *= 2.0;
            std::cout << "Old energy: " << energy << " new energy: " << newenergy << " lambda now: " << reg << std::endl;
        }
    }           
}

void inertiaTensor(const SimulationSetup &setup, SimulationState &state, const SecondFundamentalFormDiscretization &sff, Eigen::SparseMatrix<double> &result)
{
    int nverts = state.curPos.rows();
    int nedges = setup.mesh.nEdges();
    int nedgedofs = sff.numExtraDOFs();
    std::vector<Eigen::Triplet<double> > coeffs;
    for (int i = 0; i < 3 * nverts; i++)
    {
        coeffs.push_back(Eigen::Triplet<double>(i, i, 1.0));
    }
    result.resize(3 * nverts + nedges*nedgedofs, 3 * nverts + nedges*nedgedofs);
    result.setFromTriplets(coeffs.begin(), coeffs.end());
}

void leadingEigenvector(const SimulationSetup &setup, SimulationState &state, const SecondFundamentalFormDiscretization &sff, Eigen::VectorXd &result)
{
    double lameAlpha = setup.YoungsModulus * setup.PoissonsRatio / (1.0 - setup.PoissonsRatio * setup.PoissonsRatio);
    double lameBeta = setup.YoungsModulus / 2.0 / (1.0 + setup.PoissonsRatio);
    int nverts = state.curPos.rows();
    int nedges = setup.mesh.nEdges();
    int nedgedofs = sff.numExtraDOFs();

    int constrainedDOFs = setup.clampedDOFs.size();

    int freeDOFs = 3 * nverts + nedgedofs * nedges - constrainedDOFs;
    std::vector<Eigen::Triplet<double> > proj;
    int row = 0;
//    for (int i = 0; i < nverts; i++)
//    {
//        for (int j = 0; j < 3; j++)
//        {
//            if (setup.clampedDOFs.find(3*i+j) != setup.clampedDOFs.end())
//                continue;
//            proj.push_back(Eigen::Triplet<double>(row, 3 * i + j, 1.0));
//            row++;
//        }
//    }
    for (int i = 0; i < nedges; i++)
    {
        for (int j = 0; j < nedgedofs; j++)
        {
            proj.push_back(Eigen::Triplet<double>(row, 3 * nverts + nedgedofs * i + j, 1.0));
            row++;
        }
    }
    assert(row == freeDOFs);
    Eigen::SparseMatrix<double> projM(freeDOFs, 3 * nverts + nedgedofs * nedges);
    projM.setFromTriplets(proj.begin(), proj.end());
    std::vector<Eigen::Triplet<double> > hessian;

    elasticEnergy(setup.mesh, state.curPos, state.curEdgeDOFs, lameAlpha, lameBeta, setup.thickness, setup.abars, setup.bbars, sff, NULL, &hessian);

    Eigen::SparseMatrix<double> H(3 * nverts + nedgedofs * nedges, 3 * nverts + nedgedofs * nedges);
    H.setFromTriplets(hessian.begin(), hessian.end());
    
    double reg = 1.0;

    Eigen::SparseMatrix<double> reducedH = projM * H * projM.transpose();
    Eigen::SparseMatrix<double> I(freeDOFs, freeDOFs);
    I.setIdentity();
    reducedH += reg * I;                

    Eigen::SparseMatrix<double> M;
    inertiaTensor(setup, state, sff, M);
    Eigen::SparseMatrix<double> reducedM = projM * M * projM.transpose();
    Eigen::VectorXd v(freeDOFs);
    v.setRandom();
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(reducedH);
    for (int i = 0; i < 100; i++)
    {
        Eigen::VectorXd newv = solver.solve(v);
        newv = reducedM * newv;
        double n = newv.transpose() * (reducedM * newv);
        v = newv / sqrt(n);
    }
    double eval = v.transpose() * reducedH * v - reg;
    std::cout << "Eigenvalue: " << eval << std::endl;
    Eigen::VectorXd fullv = projM.transpose() * v;
    result.resize(nverts);
    for (int i = 0; i < nverts; i++)
        result[i] = fullv.segment<3>(3 * i).norm();    
}
