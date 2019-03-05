#ifndef IFOPTSOLVER_H
#define IFOPTSOLVER_H
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/barycenter.h>


#include "../MeshConnectivity.h"
#include "../GeometryDerivatives.h"
#include "../SecondFundamentalForm/MidedgeAverageFormulation.h"

 #ifndef EPS_BOUND
 #define EPS_BOUND 1e-10
 #endif

 namespace ifopt
 {
     class optVariables: public VariableSet
     {
     public:
         optVariables(int num, MeshConnectivity mesh, Eigen::MatrixXd initialPos, Eigen::MatrixXd targetPos, const std::string& name) : VariableSet(num, name)
         {
             _x.resize(num);
             _x.setZero();
             int nfaces =  mesh.nFaces();
             int nverts = initialPos.rows();

             if(num != 3* (nfaces +  nverts))
             {
                 std::cout<<"Wrong variable size"<<std::endl;
             }
             else
             {
                 for (int i=0; i<nverts; i++)
                 {
                     _x(3*i) = targetPos(i, 0);
                     _x(3*i+1) = targetPos(i,1);
                     _x(3*i+2) = targetPos(i,2);
                 }

                 for ( int i=0; i < nfaces; i++)
                 {
                     Eigen::Matrix2d A = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
                     _x(3*i + 3*nverts) = sqrt(A(0,0));
                     _x(3*i + 3*nverts + 1) = A(0,1)/_x(3*i + 3*nverts);
                     _x(3*i + 3*nverts + 2) = sqrt(A.determinant()) / _x(3*i + 3*nverts);
                 }


             }
         }

         void SetVariables(const VectorXd& x) override
         {
             _x = x;
         }

         VectorXd GetValues() const override
         {
             return _x;
         }

         VecBound GetBounds() const override
         {
             VecBound bounds(GetRows());
             int num = GetRows()/3;
             for(int i=0;i<num;i++)
             {
                 bounds.at(3*i) = NoBound;
                 bounds.at(3*i+1) = NoBound;
                 bounds.at(3*i+2) = NoBound;
             }
             return bounds;
         }

     private:
         VectorXd _x;
     };

     class optConstraint: public ConstraintSet
     {
     public:
         optConstraint(int num, MeshConnectivity mesh, double lameAlpha, double lameBeta, double thickness, const std::string& name) : ConstraintSet(num, name)
         {
             _mesh = mesh;
             _lameAlpha = lameAlpha;
             _lameBeta = lameBeta;
             _thickness = thickness;
             
             int nfaces = mesh.nFaces();
             int nverts = num / 3;
             
             _boundaryMatrix.resize(3*(nverts + nfaces), 3*(nverts + nfaces));
             _boundaryMatrix.setZero();
             Eigen::VectorXi boundaryLoop = mesh.getBoundaryLoop();
             std::vector<Eigen::Triplet<double> > triplet;

             for(int i=0;i<boundaryLoop.size();i++)
             {
                 for(int j=0;j<3;j++)
                 {
                     triplet.push_back(Eigen::Triplet<double>(3*boundaryLoop(i) + j, 3*boundaryLoop(i) + j, 1));
                 }
             }

             for(int i=0; i<nfaces; i++)
             {
                 for(int j=0;j<3;j++)
                 {
                     triplet.push_back(Eigen::Triplet<double>(3*(nverts + i) + j, 3*(nverts + i) + j, 1));
                 }
             }

             _boundaryMatrix.setFromTriplets(triplet.begin(), triplet.end());

         }

         VectorXd GetValues() const override;
         
         VectorXd getValues(Eigen::VectorXd x) const;

         VecBound GetBounds() const override
         {
             VecBound b(GetRows());
             for(int i=0;i<GetRows();i++)
             {
                b.at(i) = Bounds(0.0,0.0);
             }
             return b;
         }

         void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override;
         
         void fillJacobianBlock(Eigen::VectorXd x, Jacobian& jac_block) const;
         
         void testValueJacobian(Eigen::VectorXd x);

     private:
         void convertVariable2ABbarsPos(Eigen::VectorXd x, std::vector<Eigen::Matrix2d> &abars, std::vector<Eigen::Matrix2d> &bbars, Eigen::MatrixXd &curPos) const;

     private:
         MeshConnectivity _mesh;
         double _lameAlpha;
         double _lameBeta;
         double _thickness;
         
         Eigen::SparseMatrix<double> _boundaryMatrix;

     };

     class optCost: public CostTerm
     {
     public:
         optCost(Eigen::MatrixXd initialPos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, double lambda) : optCost("cost_term1")
         {
             _tarPos = tarPos;
             _initialPos = initialPos;
             _mesh = mesh;
             _lambda = lambda;
//             _mu = lambda_2;
             igl::cotmatrix(initialPos,mesh.faces(),L);
             selectedX = computeSelectMatrix(tarPos.rows(), mesh.nFaces(), 0);
             selectedY = computeSelectMatrix(tarPos.rows(), mesh.nFaces(), 1);
             selectedZ = computeSelectMatrix(tarPos.rows(), mesh.nFaces(), 2);
             
             int nverts = tarPos.rows();
             int nfaces = mesh.nFaces();
             
             _tarPosVec.resize(3*(nverts + nfaces));
             _tarPosVec.setZero();
             
             for (int i=0; i<nverts; i++)
             {
                 _tarPosVec(3*i) = tarPos(i, 0);
                 _tarPosVec(3*i+1) = tarPos(i,1);
                 _tarPosVec(3*i+2) = tarPos(i,2);
             }
             igl::doublearea(initialPos, mesh.faces(), _areaList);
             igl::barycenter(initialPos, mesh.faces(), _bcPos);
             
             _areaList = _areaList / 2;
             
             _regionArea = _areaList.sum();
             
             computeMassMatrix(_massVec, mesh, tarPos);
             
             boundaryMatrix.resize(3*(nverts + nfaces), 3*(nverts + nfaces));
             boundaryMatrix.setZero();
             Eigen::VectorXi boundaryLoop = mesh.getBoundaryLoop();
             std::vector<Eigen::Triplet<double> > triplet;
             
             for(int i=0;i<boundaryLoop.size();i++)
             {
                 for(int j=0;j<3;j++)
                 {
                     triplet.push_back(Eigen::Triplet<double>(3*boundaryLoop(i) + j, 3*boundaryLoop(i) + j, 1));
                 }
             }
             
             for(int i=0; i<nfaces; i++)
             {
                 for(int j=0;j<3;j++)
                 {
                     triplet.push_back(Eigen::Triplet<double>(3*(nverts + i) + j, 3*(nverts + i) + j, 1));
                 }
             }

             boundaryMatrix.setFromTriplets(triplet.begin(), triplet.end());
             
                 
         }
         optCost(const std::string& name) : CostTerm(name){}

         double GetCost() const override;
         
         double getCost(Eigen::VectorXd x) const;
         
         double getPenalty(Eigen::VectorXd x) const;
         
         double getSmoothness(Eigen::VectorXd x) const;

         double getDifference(Eigen::VectorXd x) const;

         void FillJacobianBlock(std::string var_set, Jacobian &jac) const override;
         
         void fillJacobianBlock(Eigen::VectorXd x, Jacobian &jac) const;
         
         void testCostJacobian(Eigen::VectorXd x);
         
     private:
         Eigen::SparseMatrix<double> computeSelectMatrix(int nVerts, int nFaces, int index);
         void computeMassMatrix( Eigen::VectorXd &massVec, MeshConnectivity mesh, Eigen::MatrixXd V);

     private:
         Eigen::MatrixXd _tarPos;
         Eigen::MatrixXd _initialPos;
         Eigen::MatrixXd _bcPos;
         Eigen::VectorXd _areaList;
         MeshConnectivity _mesh;
         Eigen::VectorXd _tarPosVec;
         Eigen::VectorXd _massVec;
         
         Eigen::SparseMatrix<double> L;
         Eigen::SparseMatrix<double> selectedX;
         Eigen::SparseMatrix<double> selectedY;
         Eigen::SparseMatrix<double> selectedZ;
         Eigen::SparseMatrix<double> boundaryMatrix;
         
         
         double _regionArea;
         double _lambda;
     public:
         double  _mu;


     };

 }

#endif
