#include <igl/boundary_loop.h>
#include "SensitiveAnalysis.h"

double SensitiveAnalysis::computeAbarSmoothness(Eigen::VectorXd curL)
{
    int nfaces =  _mesh.nFaces();
    
    double E = 0;
    Eigen::Matrix<double, 2, 3> w;
    w << 1, 0, 1,
    -1, 1, 0;
    
    std::vector<Eigen::Matrix2d> Lderivs(3);
    
    Lderivs[0] << 1,0,
    0,0;
    
    Lderivs[1] << 0,0,
    1,0;
    
    Lderivs[2] << 0,0,
    0,1;
    
    
    for(int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix2d L;
        
        L << curL(3*i), 0,
        curL(3*i+1), curL(3*i+2);
        for(int j = 0; j < 3; j++)
        {
            int oppVerIdx =  _mesh.vertexOppositeFaceEdgeIndex(i, j);
            int oppFace = _mesh.faceOppositeVertex(i, j);
            if (oppFace != -1)
            {
                Eigen::Matrix2d Lj;
                Lj << curL(3*oppFace), 0,
                curL(3*oppFace+1), curL(3*oppFace+2);
                
                // Compute the tranfermation matrix M
                
                Eigen::Vector3d oppNormal, curNormal, oppEdge;
                
                oppNormal = faceNormal(_mesh, _initialPos, oppFace, oppVerIdx, NULL, NULL);
                curNormal = faceNormal(_mesh, _initialPos, i, j, NULL, NULL);
                
                oppNormal = oppNormal/oppNormal.norm();
                curNormal = curNormal/curNormal.norm();
                
                oppEdge = _initialPos.row(_mesh.faceVertex(i, (j + 1)%3)) - _initialPos.row(_mesh.faceVertex(i, (j + 2)%3));
                oppEdge = oppEdge/oppEdge.norm();
                
                Eigen::Matrix3d A, A1, T;
                A.col(0) = oppEdge;
                A.col(1) = curNormal;
                A.col(2) = oppEdge.cross(curNormal);
                
                A1.col(0) = oppEdge;
                A1.col(1) = oppNormal;
                A1.col(2) = oppEdge.cross(oppNormal);
                
                T = A1*A.inverse();
                
                Eigen::Matrix2d R, oppR, abar, abarj;
                
                R.col(0) = ( _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                R.col(1) = ( _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                oppR.col(0) = ( _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                oppR.col(1) = ( _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                
                abar = (R.transpose()).inverse() * L * L.transpose() * R.inverse();
                
                abarj = (oppR.transpose()).inverse() * Lj * Lj.transpose() * oppR.inverse();
                
                double area = ( _areaList(i) + _areaList(oppFace) ) / 3.0;
                
                E +=  1.0 / _regionArea * ( (abar - abarj) * (abar - abarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                
                
            }
            
        }
    }
    
    return E;
}

void SensitiveAnalysis::computeAbarSmoothnessGrad(Eigen::VectorXd curL, Eigen::VectorXd &grad)
{
    int nfaces =  _mesh.nFaces();
    grad.resize(3*nfaces);
    grad.setZero();
    
    Eigen::Matrix<double, 2, 3> w;
    w << 1, 0, 1,
    -1, 1, 0;
    
    std::vector<Eigen::Matrix2d> Lderivs(3);
    
    Lderivs[0] << 1,0,
    0,0;
    
    Lderivs[1] << 0,0,
    1,0;
    
    Lderivs[2] << 0,0,
    0,1;
    
    
    for(int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix2d L;
        
        L << curL(3*i), 0,
        curL(3*i+1), curL(3*i+2);
        for(int j = 0; j < 3; j++)
        {
            int oppFace = _mesh.faceOppositeVertex(i, j);
            if (oppFace != -1)
            {
                Eigen::Matrix2d Lj;
                Lj << curL(3*oppFace), 0,
                curL(3*oppFace+1), curL(3*oppFace+2);
                
                Eigen::Matrix2d R, oppR, abar, abarj;
                
                R.col(0) = ( _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                R.col(1) = ( _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                oppR.col(0) = ( _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                oppR.col(1) = ( _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                
                abar = (R.inverse()).transpose() * L * L.transpose() * R.inverse();
                
                abarj = (oppR.inverse()).transpose() * Lj * Lj.transpose() * oppR.inverse();
                
                for(int k = 0; k < 3; k++)
                {
                    Eigen::Matrix2d abarderiv, abarderivj;
                    abarderiv = (R.inverse()).transpose() *( Lderivs[k] * L.transpose() + L * Lderivs[k].transpose() ) * R.inverse();
                    abarderivj = (oppR.inverse()).transpose() * (Lderivs[k] * Lj.transpose() + Lj * Lderivs[k].transpose() ) * oppR.inverse();
                    
                    double result;
                    double area = ( _areaList(i) + _areaList(oppFace) ) / 3.0;
                    
                    result = 2.0 / _regionArea * ( abarderiv * (abar - abarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                    
                    grad(3*i+k) += result;
                    
                    
                    result = - 2.0 / _regionArea * ( abarderivj * (abar - abarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                    
                    grad(3 * oppFace + k) += result;
                }
                
            }
            
        }
    }
}


void SensitiveAnalysis::testValueGrad(Eigen::VectorXd X, Eigen::VectorXd Y, Eigen::MatrixXd Z)
{
    projectBack(X, Y, Z);
    std::cout<<"Projection Done!!"<<std::endl;
    Eigen::VectorXd updatedY = Y;
    Eigen::MatrixXd updatedZ = Z;
    double f = value(X, Y, Z);
    Eigen::VectorXd df;
    gradient(X, Y, Z, df);
    Eigen::VectorXd epsVec = Eigen::VectorXd::Random(X.size() + Y.size());
    
    epsVec.normalized();
    for(int i = 6; i< 14; i++ )
    {
        double eps = pow(4, -i);
        Eigen::VectorXd updatedX = X + eps * epsVec.segment(0, X.size());
        Eigen::VectorXd updatedY = Y + eps * epsVec.segment(X.size(), Y.size());
        projectBack(updatedX, updatedY, updatedZ);
        double f1 = value(updatedX, updatedY, updatedZ);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<df.dot(epsVec)<< " Finite difference is "<<(f1 - f)/eps<<std::endl;
        std::cout<< "The error is "<<std::abs(df.dot(epsVec) - (f1 - f)/eps )<<std::endl;
        updatedY = Y;
        updatedZ = Z;
    }
    
}

void SensitiveAnalysis::computeMassMatrix(Eigen::VectorXd &massVec, MeshConnectivity mesh, Eigen::MatrixXd V)
{
    int nverts = V.rows();
    int nfaces = mesh.nFaces();
    massVec.resize(nverts);
    massVec.setZero();
    
    Eigen::VectorXd areaList;
    igl::doublearea(V, mesh.faces(), areaList);
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
    
    Eigen::VectorXi boundary;
    igl::boundary_loop(mesh.faces(), boundary);
    for(int i=0;i<boundary.size();i++)
    {
        int vidx = boundary(i);
        massVec(vidx) *= 1e-3;
    }
}


double SensitiveAnalysis::computeDifference(Eigen::MatrixXd curPos)
{
    double differece = 0;
    int nverts = curPos.rows();
    
    for(int i = 0; i < nverts; i++)
    {
        differece += 0.5 * _massVec(i) * (curPos.row(i) - _tarPos.row(i)).squaredNorm();
    }
    return differece;
}

void SensitiveAnalysis::computeDifferenceGrad(Eigen::MatrixXd curPos, Eigen::VectorXd &grad)
{
    int nverts = curPos.rows();
    grad.resize(3*nverts);
    grad.setZero();
    for(int i = 0; i<nverts;i++)
    {
        grad.segment(3*i, 3) = _massVec(i) * (curPos.row(i) - _tarPos.row(i));
    }
    grad = grad;
}
