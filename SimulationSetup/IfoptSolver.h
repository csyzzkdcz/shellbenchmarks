#ifndef IFOPTSOLVER_H
#define IFOPTSOLVER_H
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <iostream>

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
         optVariables(int num, MeshConnectivity mesh, Eigen::MatrixXd initialPos, const std::string& name) : VariableSet(num, name)
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
                     _x(3*i) = initialPos(i, 0);
                     _x(3*i+1) = initialPos(i,1);
                     _x(3*i+2) = initialPos(i,2);
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

         }

         VectorXd GetValues() const override;
         
         VectorXd getValues(Eigen::VectorXd x) const;

         VecBound GetBounds() const override
         {
             VecBound b(GetRows());
             for(int i=0;i<GetRows();i++)
             {
                b.at(i) = Bounds(-1e-7,1e-7);
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
         }
         optCost(const std::string& name) : CostTerm(name){}

         double GetCost() const override;
         
         double getCost(Eigen::VectorXd x) const;
         
         void FillJacobianBlock(std::string var_set, Jacobian &jac) const override;
         
         void fillJacobianBlock(Eigen::VectorXd x, Jacobian &jac) const;
         
         void testCostJacobian(Eigen::VectorXd x);

     private:
         Eigen::MatrixXd _tarPos;
         Eigen::MatrixXd _initialPos;
         MeshConnectivity _mesh;
         double _lambda;


     };

 }

#endif
