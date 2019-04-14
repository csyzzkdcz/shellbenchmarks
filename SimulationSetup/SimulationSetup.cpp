#include <igl/boundary_loop.h>
#include <igl/barycenter.h>
#include <igl/triangle/triangulate.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include "SimulationSetup.h"


void SimulationSetup::remeshProcessing(Eigen::MatrixXd remeshedPos, Eigen::MatrixXi remeshedFaces)
{
    int nfaces = mesh.nFaces();
    Eigen::MatrixXd newTarPos(remeshedPos.rows(),3);
    newTarPos.setZero();
    std::vector<Eigen::Matrix2d> convertedAbars(nfaces);
    std::vector<Eigen::Matrix2d> newAbars;
    std::vector<Eigen::Matrix2d> newBbars;
    
    for(int i=0;i<nfaces;i++)
    {
        Eigen::Matrix2d T;
        T.col(0) = ( initialPos.row(mesh.faceVertex(i, 1)) - initialPos.row(mesh.faceVertex(i, 0)) ).segment(0, 2);
        T.col(1) = ( initialPos.row(mesh.faceVertex(i, 2)) - initialPos.row(mesh.faceVertex(i, 0)) ).segment(0, 2);
        convertedAbars[i] = (T.transpose()).inverse() * abars[i] * T.inverse();
    }
    
    newAbars.resize(remeshedFaces.rows());
    newBbars.resize(remeshedFaces.rows());
    
    Eigen::MatrixXd remeshedBC;
    
    igl::barycenter(remeshedPos, remeshedFaces, remeshedBC);
    
    MeshConnectivity newMesh(remeshedFaces);
    
    for(int i=0;i<remeshedFaces.rows();i++)
    {
        Eigen::Matrix2d T;
        T.col(0) = ( remeshedPos.row(newMesh.faceVertex(i, 1)) - remeshedPos.row(newMesh.faceVertex(i, 0)) ).segment(0, 2);
        T.col(1) = ( remeshedPos.row(newMesh.faceVertex(i, 2)) - remeshedPos.row(newMesh.faceVertex(i, 0)) ).segment(0, 2);
        int flag = -1;
        for(int j=0;j<nfaces;j++)
        {
            Eigen::Matrix3d A;
            Eigen::Vector3d b;
            b << remeshedBC(i,0), remeshedBC(i,1),1;
            for(int k=0;k<3;k++)
            {
                A.col(k) = initialPos.row(mesh.faceVertex(j, k)).transpose();
            }
            A.row(2).setOnes();
            Eigen::Vector3d sol = A.inverse() * b;
            if(sol(0) >=0 && sol(1)>=0 && sol(2)>=0)
            {
                flag = j;
                break;
            }
        }
        newAbars[i] = T.transpose() * convertedAbars[flag] * T;
        newBbars[i] = abars[flag].inverse() * bbars[flag] * newAbars[i];
    }
    
    
    for(int i=0;i<remeshedPos.rows();i++)
    {
        int flag = -1;
        for(int j=0;j<nfaces;j++)
        {
            Eigen::Matrix3d A;
            Eigen::Vector3d b;
            b << remeshedPos(i,0), remeshedPos(i,1),1;
            for(int k=0;k<3;k++)
            {
                A.col(k) = initialPos.row(mesh.faceVertex(j, k)).transpose();
            }
            A.row(2).setOnes();
            Eigen::Vector3d sol = A.inverse() * b;
            if(i == remeshedPos.rows() - 2)
            {
//                std::cout<<j<<" "<<initialPos.row(mesh.faceVertex(flag, 0))<<" "<<initialPos.row(mesh.faceVertex(flag, 1))<<" "<<initialPos.row(mesh.faceVertex(flag, 2))<<std::endl;
                std::cout<<sol.transpose()<<std::endl;
            }
            if(sol(0) >=-1e-3 && sol(1)>=-1e-3 && sol(2)>=-1e-3)
            {
                flag = j;
                for(int k=0;k<3;k++)
                    newTarPos.row(i) += sol(k) * targetPos.row(mesh.faceVertex(flag, k));
                break;
            }
        }
        if(flag == -1)
        {
            double min = std::numeric_limits<double>::infinity();
            for(int k=0;k<initialPos.rows();k++)
            {
                double dist = (remeshedPos.row(i) - initialPos.row(k)).norm();
                if(dist < min)
                {
                    min = dist;
                    flag = k;
                }
            }
            newTarPos.row(i) = targetPos.row(flag);
        }
        
    }
    std::cout<<newTarPos<<std::endl;
    abars.swap(newAbars);
    bbars.swap(newBbars);
    mesh = newMesh;
    initialPos = remeshedPos;
    targetPos = newTarPos;
    igl::writeOBJ("remeshedTarget.obj", newTarPos, remeshedFaces);
}
