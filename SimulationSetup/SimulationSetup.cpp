#include <igl/boundary_loop.h>
#include <igl/barycenter.h>
#include <igl/triangle/triangulate.h>
#include "SimulationSetup.h"


void SimulationSetup::remeshProcessing(Eigen::MatrixXd remeshedPos, Eigen::MatrixXi remeshedFaces)
{
    int nfaces = mesh.nFaces();
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
        int flag = 0;
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
    
    abars.swap(newAbars);
    bbars.swap(newBbars);
    mesh = newMesh;
    initialPos = remeshedPos;
    
}
