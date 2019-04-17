
#include <igl/doublearea.h>
#include <igl/readOBJ.h>
#include <iostream>
#include "GeometryTools.h"

using namespace GeometryTools;

void GeometryTools::rigidMotionTransformation(Eigen::MatrixXd pos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, Eigen::Matrix3d &R, Eigen::Vector3d &t)
{
    int nverts = tarPos.rows();
    int nfaces = mesh.nFaces();
    Eigen::VectorXd massVec;
    massVec.resize(nverts);
    massVec.setZero();
    
    Eigen::VectorXd areaList;
    igl::doublearea(tarPos, mesh.faces(), areaList);
    areaList = areaList / 2;
    
    for(int i=0; i < nfaces; i++)
    {
        double faceArea = areaList(i);
        for(int j=0; j<3; j++)
        {
            int vertIdx = mesh.faceVertex(i, j);
            massVec(vertIdx) += faceArea / 3;
        }
    }
    
    massVec = massVec / 3;
    massVec = massVec / massVec.maxCoeff();
    
    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(nverts, nverts);
    W.diagonal() = massVec;
    
    Eigen::Vector3d avePos, aveTarPos;
    avePos.setZero();
    aveTarPos.setZero();
    
    for(int i=0;i<nverts;i++)
    {
        avePos += massVec(i) * pos.row(i);
        aveTarPos += massVec(i) * tarPos.row(i);
    }
    
    avePos = avePos / massVec.sum();
    aveTarPos = aveTarPos / massVec.sum();
    
    Eigen::MatrixXd onesMat(nverts,1);
    onesMat.setOnes();
    
    pos = pos - onesMat * avePos.transpose();
    tarPos = tarPos - onesMat * aveTarPos.transpose();
    
    Eigen::MatrixXd S = pos.transpose() * W * tarPos;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(3);
    
    
    Eigen::BDCSVD<Eigen::MatrixXd> solver(S);
    solver.compute(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U = solver.matrixU();
    Eigen::MatrixXd V = solver.matrixV();
    
    Eigen::MatrixXd middleMat(S.rows(), S.cols());
    middleMat.setIdentity();
    middleMat(S.rows()-1,S.cols()-1) = (V*U.transpose()).determinant();
    
    R = V * middleMat * U.transpose();
    t = aveTarPos - R*avePos;
    

}


void GeometryTools::testRigigMotionTransformation()
{
    Eigen::MatrixXd pos,tarPos;
    Eigen::MatrixXi F;
    igl::readOBJ("../../benchmarks/TestModels/veryCoarse/sphere/sphere_geometry.obj", pos, F);
    MeshConnectivity mesh(F);
    tarPos = pos;
    
    Eigen::Quaterniond q(2, 0, 1, -3);
    q.normalize();
    
    Eigen::Matrix3d R = q.toRotationMatrix(); // convert a quaternion to a 3x3 rotation matrix
    
    std::cout << "Rotation matrix is: " << std::endl << R << std::endl;
    
    Eigen::Vector3d t = Eigen::Vector3d::Random();
    
    std::cout << "Transition is : "<<std::endl<<t<<std::endl;
    
    int nverts = pos.rows();
    for(int i=0;i<nverts;i++)
    {
        tarPos.row(i).transpose() = R*pos.row(i).transpose() + t;
    }
    
    rigidMotionTransformation(pos, tarPos, mesh, R, t);
    
    std::cout << "Rotation matrix is: " << std::endl << R << std::endl;;
    std::cout << "Transition is : "<<std::endl<<t<<std::endl;
    
}
