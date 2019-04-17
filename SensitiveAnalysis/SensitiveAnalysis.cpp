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
                
                Eigen::Matrix2d R, oppR, deltaAbar, deltaAbarj;
                
                R.col(0) = ( _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                R.col(1) = ( _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                oppR.col(0) = ( _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                oppR.col(1) = ( _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                
                deltaAbar = (R.transpose()).inverse() * ( L * L.transpose() - _atarget[i] ) * R.inverse();
                
                deltaAbarj = (oppR.transpose()).inverse() * ( Lj * Lj.transpose() - _atarget[oppFace] ) * oppR.inverse();
                
                double area = ( _areaList(i) + _areaList(oppFace) ) / 2.0;
                
                E +=  1.0 / _regionArea * ( (deltaAbar - deltaAbarj) * (deltaAbar - deltaAbarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                
                
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
                
                Eigen::Matrix2d R, oppR, deltaAbar, deltaAbarj;
                
                R.col(0) = ( _initialPos.row(_mesh.faceVertex(i, 1)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                R.col(1) = ( _initialPos.row(_mesh.faceVertex(i, 2)) - _initialPos.row(_mesh.faceVertex(i, 0)) ).segment(0, 2);
                
                oppR.col(0) = ( _initialPos.row(_mesh.faceVertex(oppFace, 1)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                oppR.col(1) = ( _initialPos.row(_mesh.faceVertex(oppFace, 2)) - _initialPos.row(_mesh.faceVertex(oppFace, 0)) ).segment(0, 2);
                
                deltaAbar = (R.inverse()).transpose() * ( L * L.transpose() - _atarget[i] ) * R.inverse();
                
                deltaAbarj = (oppR.inverse()).transpose() * ( Lj * Lj.transpose() - _atarget[oppFace] ) * oppR.inverse();
                
                for(int k = 0; k < 3; k++)
                {
                    Eigen::Matrix2d abarderiv, abarderivj;
                    abarderiv = (R.inverse()).transpose() *( Lderivs[k] * L.transpose() + L * Lderivs[k].transpose() ) * R.inverse();
                    abarderivj = (oppR.inverse()).transpose() * (Lderivs[k] * Lj.transpose() + Lj * Lderivs[k].transpose() ) * oppR.inverse();
                    
                    double result;
                    double area = ( _areaList(i) + _areaList(oppFace) ) / 2.0;
                    
                    result = 2.0 / _regionArea * ( abarderiv * (deltaAbar - deltaAbarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                    
                    grad(3*i+k) += result;
                    
                    
                    result = - 2.0 / _regionArea * ( abarderivj * (deltaAbar - deltaAbarj).transpose() ).trace() * area / ( _bcPos.row(i) - _bcPos.row(oppFace) ).squaredNorm();
                    
                    grad(3 * oppFace + k) += result;
                }
                
            }
            
        }
    }
}


void SensitiveAnalysis::testValueGrad(Eigen::VectorXd params, Eigen::MatrixXd pos)
{
    projectBack(params, pos);
    std::cout<<"Projection Done!!"<<std::endl;
    double f = value(params, pos);
    Eigen::VectorXd df;
    gradient(params, pos, df);
    Eigen::VectorXd epsVec = Eigen::VectorXd::Random(params.size());
    
    epsVec.normalized();
    for(int i = 6; i< 20; i++ )
    {
        double eps = pow(4, -i);
        Eigen::VectorXd updatedParams = params + eps * params;
        Eigen::MatrixXd updatedPos = pos;
        projectBack(updatedParams, updatedPos);
        double f1 = value(updatedParams, updatedPos);
        std::cout<<"EPS is "<<eps<<std::endl;
        std::cout<< "Gradient is "<<df.dot(epsVec)<< " Finite difference is "<<(f1 - f)/eps<<std::endl;
        std::cout<< "The error is "<<std::abs(df.dot(epsVec) - (f1 - f)/eps )<<std::endl;
        updatedParams = params;
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
}


void SensitiveAnalysis::updateFixedVariables(Eigen::VectorXd variables)
{
    if(projM.rows() == projM.cols())
        return;
    
    fixedVariables = variables - projM.transpose() * projM * variables;
}

Eigen::VectorXd SensitiveAnalysis::getFullVariables(Eigen::VectorXd reductVariables)
{
    if(projM.rows() == projM.cols())
        return reductVariables;
    else
        return projM.transpose() * reductVariables + fixedVariables;
}


void SensitiveAnalysis::setProjM(std::set<int> fixedFlags)
{
    std::vector<Eigen::Triplet<double>> T;
    int nfaces = _mesh.nFaces();
    if(fixedFlags.size() == 0)
    {
        return;
    }
    int totalEOFs = projM.cols();
    int freeEOFs = totalEOFs - fixedFlags.size();
    projM.resize(freeEOFs, totalEOFs);
    
    int row = 0;
    for (int i = 0; i < nfaces; i++)
    {
        if (fixedFlags.find(i) != fixedFlags.end())
            continue;
        if(totalEOFs == 4 * nfaces)
        {
            for (int j = 0; j < 3; j++)
            {
                T.push_back(Eigen::Triplet<double>(row, 3 * i + j, 1.0));
                row++;
            }
        }
        else
        {
            T.push_back(Eigen::Triplet<double>(row, i, 1.0));
            row++;
        }
    }
    
    for(int i =0;i<nfaces;i++)
    {
        if (fixedFlags.find(i) != fixedFlags.end())
            continue;
        T.push_back(Eigen::Triplet<double>(row, totalEOFs - nfaces + i, 1.0));
        row++;
    }
    
    projM.setFromTriplets(T.begin(), T.end());
    
}

void SensitiveAnalysis::save(Eigen::VectorXd params, Eigen::MatrixXd pos, std::string path, bool is_initial)
{
    if(is_initial == false)
    {
        std::cout<<"Saving Abar path: "<<path<<std::endl;
        std::ofstream outfile(path, std::ios::trunc);
        int nverts = pos.rows();
        int nfaces = _mesh.nFaces();
        
        outfile<<_thickness<<"\n";
        outfile<<_lambdaAbar<<"\n";
        outfile<<_lambdaBbar<<"\n";
        outfile<<_mu<<"\n";
        outfile<<3*nverts<<"\n";
        
        Eigen::VectorXd L, S;
        convertParams2LAndS(params, L, S);
        
        int numEOFs = L.size() + S.size();
        outfile<<numEOFs<<"\n";
        
        //        std::cout<<3*nverts + 3*nfaces<<std::endl;
        
        for(int i=0;i<nverts;i++)
        {
            outfile<<std::setprecision(16)<<pos(i, 0)<<"\n";
            outfile<<std::setprecision(16)<<pos(i, 1)<<"\n";
            outfile<<std::setprecision(16)<<pos(i, 2)<<"\n";
        }
        
        for(int i=0;i<L.size();i++)
        {
            outfile<<std::setprecision(16)<<L(i)<<"\n";
        }
        
        for(int i=0;i<S.size();i++)
        {
            outfile<<std::setprecision(16)<<S(i)<<"\n";
        }
        
        outfile.close();
    }
    int startIdx, endIdx, expCoef;
    std::string subString = "";
    std::string resampledPath = path;
    
    startIdx = resampledPath.rfind("/");
    endIdx = resampledPath.find("_");
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx + 1,resampledPath.begin() + endIdx, "resampled");
    
    // thickness
    if(_thickness == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(_thickness));
    startIdx = resampledPath.rfind("T");
    endIdx = resampledPath.rfind("A");
    subString = "";
    if(_thickness > 0)
        subString = "T_1e" + std::to_string(expCoef);
    else
        subString = "T_0";
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    //Abar penalty
    if(_lambdaAbar == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(_lambdaAbar));
    
    startIdx = resampledPath.rfind("A");
    endIdx = resampledPath.rfind("B");
    subString = "";
    if(_lambdaAbar > 0)
        subString = "A_1e" + std::to_string(expCoef);
    else
        subString = "A_0";
    resampledPath= resampledPath .replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    // Bbar penalty
    if(_lambdaBbar == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(_lambdaBbar));
    startIdx = resampledPath.rfind("B");
    endIdx = resampledPath.rfind("S");
    subString = "";
    if(_lambdaAbar > 0)
        subString = "B_1e" + std::to_string(expCoef);
    else
        subString = "B_0";
    resampledPath= resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);
    
    
    // smoothness
    if(_mu == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(_mu));
    
    startIdx = resampledPath.rfind("S");
    endIdx = resampledPath.rfind(".");
    subString = "";
    if(_mu > 0)
        subString = "S_1e" + std::to_string(expCoef);
    else
        subString = "S_0";
    resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx, subString);
    
    startIdx = resampledPath.rfind(".");
    if(is_initial)
        resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.end(), "_target.obj");
    else
        resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.end(), ".obj");
    std::cout<<"Current abar loading path is: "<<resampledPath<<std::endl;
    igl::writeOBJ(resampledPath, pos, _mesh.faces());
}

