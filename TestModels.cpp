#include <vector>
#include <igl/barycenter.h>
#include <igl/doublearea.h>

#include "TestModels.h"
#include "BuildModels.h"
#include "MeshConnectivity.h"

void testTrapezoid(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
   int nfaces = F.rows();
   int nverts = V.rows();
   
   std::vector<Eigen::Matrix2d> oldAbars;
   std::vector<Eigen::Matrix2d> abars;
   std::vector<Eigen::Matrix2d> bbars;
   Eigen::VectorXd areaList;
   igl::doublearea(V, F, areaList);
   
   areaList = areaList / 2;
   
   double regionArea = areaList.sum();
   
   Eigen::MatrixXd BC;
   igl::barycenter(V, F, BC);
   
   MeshConnectivity mesh(F);
   
   Eigen::Matrix<double, 4, 2> rectCorners;
   Eigen::Matrix<double, 4, 2> trapeCorners;
   
   rectCorners << -0.5, -0.5,
   -0.5, 0.5,
   0.5, 0.5,
   0.5, -0.5;
   
   trapeCorners << -1, -0.5,
   -0.5, 0.5,
   0.5, 0.5,
   1, -0.5;
   
   
   Eigen::Matrix3d M = compute_mapping(rectCorners, trapeCorners);
   
   //                std::cout<<M<<std::endl;
   
   abars.resize(nfaces);
   bbars.resize(nfaces);
   for(int i = 0; i < nfaces; i++)
   {
       Eigen::Vector3d bcPos = BC.row(i);
       bcPos(2) = 1;
       Eigen::Matrix2d dr; // dr = (rx, ry)
       Eigen::Vector2d rx, ry;
       dr.setZero();
       for(int i = 0; i < 2; i++)
       {
           for(int k = 0; k < 2; k++ )
               dr(k, i) = ( M(k,i) * M.row(2).dot(bcPos) - M(2,i) * M.row(0).dot(bcPos) ) / ( M.row(2).dot(bcPos) * M.row(2).dot(bcPos) );
       }
       rx << M(0,0) * M.row(2).dot(bcPos) - M(2,0)*M.row(0).dot(bcPos),
       M(1,0) * M.row(2).dot(bcPos) - M(2,0)*M.row(1).dot(bcPos);
       
       rx = rx / ( M.row(2).dot(bcPos) * M.row(2).dot(bcPos) );
       
       ry << M(0,1) * M.row(2).dot(bcPos) - M(2,1)*M.row(0).dot(bcPos),
       M(1,1) * M.row(2).dot(bcPos) - M(2,1)*M.row(1).dot(bcPos);
       
       ry = ry / ( M.row(2).dot(bcPos) * M.row(2).dot(bcPos) );
       
       Eigen::Matrix2d abar;
       
       abar << rx.dot(rx), rx.dot(ry),
       ry.dot(rx), ry.dot(ry);
       
       Eigen::Matrix2d newAbar, T;
       T.col(0) = V.row(mesh.faceVertex(i, 1)).segment(0, 2) - V.row(mesh.faceVertex(i, 0)).segment(0, 2);
       T.col(1) = V.row(mesh.faceVertex(i, 2)).segment(0, 2) - V.row(mesh.faceVertex(i, 0)).segment(0, 2);
       newAbar = T.transpose() * abar * T;
       abars[i] = newAbar;
       bbars[i].setZero();
       oldAbars.push_back(abar);
       //                    std::cout<<"double(subs(a, [x,y], ["<<bcPos(0)<<","<<bcPos(1)<<"]"<<"))"<<std::endl<<std::endl;
       //                    std::cout<<abar<<std::endl<<std::endl;
   }
   // compute the penalty term
   double E = 0;
   for(int i=0; i< nfaces; i++)
   {
       for(int j = 0; j < 3; j++)
       {
           int oppFace = mesh.faceOppositeVertex(i, j);
           if (oppFace != -1)
           {
               double area = ( areaList(i) + areaList(oppFace) ) / 3;
               Eigen::Matrix2d finiteDifference =  ( oldAbars[i] - oldAbars[oppFace] ) / (BC.row(i) - BC.row(oppFace)).norm();
               E += ( finiteDifference * finiteDifference.transpose() ).trace() * area / regionArea;
           }
       }
   }
   std::cout<<E<<std::endl;
}

void testSphere(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
   int nfaces = F.rows();
   int nverts = V.rows();

   std::vector<Eigen::Matrix2d> oldAbars;
   std::vector<Eigen::Matrix2d> abars;
   std::vector<Eigen::Matrix2d> bbars;
   Eigen::VectorXd areaList;
   igl::doublearea(V, F, areaList);

   areaList = areaList / 2;

   double regionArea = areaList.sum();

   Eigen::MatrixXd BC;
   igl::barycenter(V, F, BC);

   MeshConnectivity mesh(F);

   abars.resize(nfaces);
   bbars.resize(nfaces);
   for(int i = 0; i < nfaces; i++)
   {
       double R = 0.5;
//        double z = 0.5*(R - R*R/sqrt(R*R+u*u+v*v));
//        double x = 0.5*((R-z)/R*u);
//        double y = 0.5*(R-z)/R*v;

       Eigen::Vector3d bcPos = BC.row(i);
       bcPos(2) = 1;
       double u = bcPos(0);
       double v = bcPos(1);
       
       Eigen::Vector3d ru, rv;
       ru <<  1/(2*sqrt(u*u + v*v + 1.0/4.0)) - u*u/(2*pow( (u*u + v*v + 1.0/4.0), 3.0/2.0)),
        -(u*v)/(2*pow( (u*u + v*v + 1.0/4.0), 3.0/2.0)),
       u/(4*pow( (u*u + v*v + 1.0/4.0), 3.0/2.0));

       rv << -(u*v)/(2*pow( (u*u + v*v + 1.0/4.0), 3.0/2.0)),
       1/(2*sqrt(u*u + v*v + 1.0/4.0)) - v*v/(2*pow( (u*u + v*v + 1.0/4.0), 3.0/2.0)),
       v/(4*pow( (u*u + v*v + 1.0/4.0), 3.0/2.0));

       ru = ru / 2.0;
       rv = rv / 2.0;
       
       Eigen::Matrix2d abar;

       abar << ru.dot(ru), ru.dot(rv),
       rv.dot(ru), rv.dot(rv);

       Eigen::Matrix2d newAbar, T;
       T.col(0) = V.row(mesh.faceVertex(i, 1)).segment(0, 2) - V.row(mesh.faceVertex(i, 0)).segment(0, 2);
       T.col(1) = V.row(mesh.faceVertex(i, 2)).segment(0, 2) - V.row(mesh.faceVertex(i, 0)).segment(0, 2);
       newAbar = T.transpose() * abar * T;
       abars[i] = newAbar;
       bbars[i].setZero();
       oldAbars.push_back(abar);
       
//       std::cout<<"double(subs(A, [u,v], ["<<bcPos(0)<<","<<bcPos(1)<<"]"<<"))"<<std::endl<<std::endl;
//       std::cout<<abar<<std::endl<<std::endl;
   }
   // compute the penalty term
   double E = 0;
   for(int i=0; i< nfaces; i++)
   {
       for(int j = 0; j < 3; j++)
       {
           int oppFace = mesh.faceOppositeVertex(i, j);
           if (oppFace != -1)
           {
               double area = ( areaList(i) + areaList(oppFace) ) / 3;
               Eigen::Matrix2d finiteDifference =  ( oldAbars[i] - oldAbars[oppFace] ) / (BC.row(i) - BC.row(oppFace)).norm();
               E += ( finiteDifference * finiteDifference.transpose() ).trace() * area / regionArea;
           }
       }
   }
   std::cout<<E<<std::endl;
}

void testHypar(Eigen::MatrixXd V, Eigen::MatrixXi F)
{

}

void testSaddle(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    
}
