#ifndef BUILDMODELS_H
#define BUILDMODELS_H

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <Eigen/Dense>

void compute_sphere(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        /*
         x = R*sin(u / R)
         y = v
         z = R*cos(u / R) - R
         */
        double R = 0.5;
        double u = Vo(i,0);
        double v = Vo(i,1);
        double z = R - R*R/sqrt(R*R+u*u+v*v);
        double x = (R-z)/R*u;
        double y = (R-z)/R*v;
        Vo(i,0) = x;
        Vo(i,1) = y;
        Vo(i,2) = z;
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind-1) + "/sphere_geometry.obj", Vo, Fo);
    
}

void compute_hypar(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        Vo(i,2) = 32*Vo(i,1) * Vo(i,0);
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind-1) + "/hypar_geometry.obj", Vo, Fo);
}

void compute_cylinder(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        /*
         x = R*sin(u / R)
         y = v
         z = R*cos(u / R) - R
         */
        double R = 0.5;
        double x = R*sin(Vo(i,0)/R);
        double z = R*cos(Vo(i,0)/R) - R;
        Vo(i,0) = x;
        Vo(i,2) = z;
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind-1) + "/cylinder_geometry.obj", Vo, Fo);
}

void compute_saddle(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        double x = Vo(i,0);
        double y = Vo(i,1);
        double z = x*x -y*y;
        Vo(i,2) = z;
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind-1) + "/saddle_geometry.obj", Vo, Fo);
}

void meshResample(std::string filepath, int targetResolution)
{
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi E;
//    Eigen::MatrixXd H;
//
//    // Triangulated interior
//    Eigen::MatrixXd V2;
//    Eigen::MatrixXi F2;
//    V.resize(8,1);
//    E.resize(8,1);
//
//    V << -1,-1, 1,-1, 1,1, -1, 1,
//
//    E << 0,1, 1,2, 2,3, 3,0,
//
//
//    // Triangulate the interior
//    igl::triangle::triangulate(V,E,H,"a0.005q",V2,F2);
    
    // Plot the generated mesh
    Eigen::MatrixXd Vcurr;
    Eigen::MatrixXi Fcurr;
    igl::readOBJ(filepath, Vcurr, Fcurr);
    while ( Fcurr.rows() < targetResolution * 2)
    {
        Eigen::MatrixXd Vtmp = Vcurr;
        Eigen::MatrixXi Ftmp = Fcurr;
        igl::upsample(Vtmp, Ftmp, Vcurr, Fcurr, 1);
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi J;

    igl::decimate(Vcurr, Fcurr, targetResolution, V, F, J);
    int ind = filepath.rfind(".");
    igl::writeOBJ(filepath.substr(0, ind-1) + "_resampled.obj", V, F);
    
}

#endif
