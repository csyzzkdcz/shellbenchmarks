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

void compute_expanded_rect(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        Vo.row(i) = 2 * Vo.row(i);
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/expandedRect_geometry.obj", Vo, Fo);
}

void compute_extended_rect(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        Vo(i, 0) = 2*Vo(i, 0);
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/extendedRect_geometry.obj", Vo, Fo);
}

Eigen::Matrix3d compute_mapping(Eigen::Matrix<double, 4, 2> rectCorners, Eigen::Matrix<double, 4, 2> trapeCorners)
{
    Eigen::Matrix<double, 8, 8> A;
    A.setZero();
    
    A.block(0,0,4, 2) = rectCorners;
    A.block(4, 3, 4, 2) = rectCorners;
    
    Eigen::Matrix<double, 4, 1> ones;
    ones.setConstant(1);
    
    A.block(0, 2, 4, 1) = ones;
    A.block(4, 5, 4, 1) = ones;
    
    for(int i = 0; i<4; i++)
    {
        A(i, 6) = -rectCorners(i, 0) * trapeCorners(i, 0);
        A(i, 7) = -rectCorners(i, 1) * trapeCorners(i, 0);
        A(4+i, 6) = -rectCorners(i, 0) * trapeCorners(i, 1);
        A(4+i, 7) = -rectCorners(i, 1) * trapeCorners(i, 1);
    }
    
//    std::cout<<A<<std::endl;
    
    Eigen::Matrix<double, 8, 1> b;
    for(int i=0; i<4; i++)
    {
        b(i, 0) = trapeCorners(i, 0);
        b(4+i, 0) = trapeCorners(i , 1);
    }
    
//    std::cout<<b<<std::endl;
//    std::cout<<A.inverse()<<std::endl;
    
    Eigen::Matrix<double, 8, 1> sol = A.inverse() * b;
    
//    std::cout<<sol<<std::endl;
    
    Eigen::Matrix3d transM;
    transM << sol(0, 0), sol(1, 0), sol(2, 0),
    sol(3, 0), sol(4, 0), sol(5, 0),
    sol(6, 0), sol(7, 0), 1;
    
//    std::cout<<transM<<std::endl;
    return transM;
}


void compute_trapezoid(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);

    Eigen::Matrix<double, 4, 2> rectCorners;
    Eigen::Matrix<double, 4, 2> trapeCorners;
    
    rectCorners << -0.5, -0.5,
                    -0.5, 0.5,
                    0.5, 0.5,
                    0.5, -0.5;
    
    trapeCorners << -1, -1,
                    -0.5, 0.5,
                    0.5, 0.5,
                    1, -1;

    Eigen::Matrix3d M = compute_mapping(rectCorners, trapeCorners);

    for(int i=0;i<Vo.rows();i++)
    {
        Eigen::Vector3d x;
        x << Vo(i, 0), Vo(i, 1), 1;
        x = M * x;
        Vo(i, 0) = x(0)/x(2);
        Vo(i, 1) = x(1)/x(2);
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/trapeZoid_geometry.obj", Vo, Fo);
}


void compute_circle(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    
    double R = 0.5;
    
    for(int i=0;i<Vo.rows();i++)
    {
        double x = Vo(i, 0);
        double y = Vo(i, 1);
        
        Vo(i, 0) = x * sqrt(1 - y*y / (2 * R * R));
        Vo(i, 1) = y * sqrt(1 - x*x / (2 * R * R));
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/circle_geometry.obj", Vo, Fo);
}




#endif
