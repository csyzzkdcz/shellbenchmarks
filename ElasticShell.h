#ifndef ELASTICSHELL_H
#define ELASTICSHELL_H

#include <Eigen/Core>
#include "MeshConnectivity.h"
#include <vector>
#include <Eigen/Sparse>
#include "SecondFundamentalForm/SecondFundamentalFormDiscretization.h"

double stretchingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    double lameAlpha, double lameBeta, double thickness,
    const Eigen::Matrix2d &abar,
    int face,
    Eigen::Matrix<double, 1, 9> *derivative, // F(face, i)
    Eigen::Matrix<double, 9, 9> *hessian);
    
double bendingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &edgeDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const Eigen::Matrix2d &abar, const Eigen::Matrix2d &bbar,
    int face,
    const SecondFundamentalFormDiscretization &sff,
    Eigen::MatrixXd *derivative, // F(face, i), then the three vertices opposite F(face,i), then the thetas on oppositeEdge(face,i)
    Eigen::MatrixXd *hessian);

double elasticEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &edgeDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const std::vector<Eigen::Matrix2d> &abars,
    const std::vector<Eigen::Matrix2d> &bbars,
    const SecondFundamentalFormDiscretization &sff,
    Eigen::VectorXd *derivative, // positions, then thetas
    std::vector<Eigen::Triplet<double> > *hessian);


void testStretchingEnergy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
void testBendingEnergy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const SecondFundamentalFormDiscretization &sff);
void testElasticEnergy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const SecondFundamentalFormDiscretization &sff);
void testElasticEnergy(const MeshConnectivity &mesh, const Eigen::MatrixXd &V, const Eigen::VectorXd &thetas, double lameAlpha, double lameBeta, double thickness, const std::vector<Eigen::Matrix2d> &abars, const std::vector<Eigen::Matrix2d> &bbars, const SecondFundamentalFormDiscretization &sff);


#endif
