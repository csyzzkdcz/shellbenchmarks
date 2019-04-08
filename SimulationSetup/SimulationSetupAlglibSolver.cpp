#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include "SimulationSetupAlglibSolver.h"
#include "../GeometryDerivatives.h"
#include "../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../ElasticShell.h"

int itrTimes = 0;
typedef Eigen::Triplet<double> T;

void callback_func_coef(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df, void* ptr)
{
    SimulationSetupAlglibSolver* callback = (SimulationSetupAlglibSolver*) ptr;
    callback->getFunctionGradient(x,f,df);
}

bool SimulationSetupAlglibSolver::loadAbars()
{
    std::ifstream infile(abarPath);
    if(!infile)
        return false;
    int num;
    infile >> num;
    double d1, d2, d3;
    for(int i = 0; i < abars.size(); i++)
    {
        infile >> d1 >> d2 >> d3;
        Eigen::Matrix2d L;
        L << d1,0,
        d2,d3;
        abars[i] = L*L.transpose();
    }
    return true;
}

void SimulationSetupAlglibSolver::saveAbars()
{
    
}

void SimulationSetupAlglibSolver::buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
//    testFucntionGradient(sff);
//    return;
    int nfaces = mesh.nFaces();
    abars.resize(nfaces);
    bbars.resize(nfaces);
    
    initialAbars.resize(nfaces);

    lameAlpha = YoungsModulus * PoissonsRatio / (1.0 - PoissonsRatio * PoissonsRatio);
    lameBeta = YoungsModulus / 2.0 / (1.0 + PoissonsRatio);

    for (int i = 0; i < nfaces; i++)
    {
        bbars[i].setZero();
        initialAbars[i] = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
    }
    
//    testPenaltyTerm();
//    return;
    
    if(!loadAbars())
        findFirstFundamentalForms(sff);

}

void SimulationSetupAlglibSolver::findFirstFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
    int nfaces = mesh.nFaces();
    atargets.resize(nfaces);
    btargets.resize(nfaces);
    aderivs.resize(nfaces);
    bderivs.resize(nfaces);
    
    alglib::real_1d_array x;
    double epsg = 0;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0;
    alglib::ae_int_t maxits = 1e5;
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;
    x.setlength(3*nfaces);

    for (int i = 0; i < nfaces; i++)
    {
        atargets[i] = firstFundamentalForm(mesh, targetPos, i, &aderivs[i], NULL);
        btargets[i] = sff.secondFundamentalForm(mesh, targetPos, initialEdgeDOFs, i, &bderivs[i], NULL);
        abars[i] = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
        
        x[3*i] = sqrt(atargets[i](0,0));
        x[3*i+1] = atargets[i](0,1)/x[3*i];
        x[3*i+2] = sqrt(atargets[i].determinant())/x[3*i];

//        x[3*i] = sqrt(abars[i](0,0));
//        x[3*i+1] = abars[i](0,1)/x[3*i];
//        x[3*i+2] = sqrt(abars[i].determinant())/x[3*i];
    }
    
    alglib::minlbfgscreate(x.length(),x, state);
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    alglib::minlbfgssetstpmax(state, stpmax);
    alglib::minlbfgsoptimize(state, callback_func_coef, NULL, this);
    alglib::minlbfgsresults(state, x, rep);
    
    printf("%d\n", int(rep.terminationtype));
    
    //    Rep     -   optimization report:
    //    * Rep.TerminationType completetion code:
    //    * -8    internal integrity control  detected  infinite
    //    or NAN values in  function/gradient.  Abnormal
    //    termination signalled.
    //    * -7    gradient verification failed.
    //    See MinLBFGSSetGradientCheck() for more information.
    //        * -2    rounding errors prevent further improvement.
    //        X contains best point found.
    //        * -1    incorrect parameters were specified
    //        *  1    relative function improvement is no more than
    //        EpsF.
    //        *  2    relative step is no more than EpsX.
    //        *  4    gradient norm is no more than EpsG
    //        *  5    MaxIts steps was taken
    //        *  7    stopping conditions are too stringent,
    //        further improvement is impossible
    //        *  8    terminated by user who called minlbfgsrequesttermination().
    //        X contains point which was "current accepted" when
    //        termination request was submitted.
    //        * Rep.IterationsCount contains iterations count
    //        * NFEV countains number of function calculations
    
    for(int i = 0; i < nfaces; i++)
    {
        abars[i] << x[3*i], 0,
        x[3*i+1], x[3*i+2];
        abars[i] = abars[i] * abars[i].transpose();
    }
    
    std::ofstream outfile(abarPath,std::ios::trunc);
    outfile<<x.length()<<"\n";
    for(int i=0;i<x.length()-1;i++)
    {
        outfile<<std::setprecision(16)<<x[i]<<"\n";
    }
    outfile<<std::setprecision(16)<<x[x.length()-1];
    outfile.close();
    

}

void SimulationSetupAlglibSolver::getFunctionGradient(const alglib::real_1d_array &x, double &f, alglib::real_1d_array &df)
{
    
    itrTimes++;
    
    if(itrTimes%100 == 0)
    {
        std::ofstream outfile(abarPath,std::ios::trunc);
        outfile<<x.length()<<"\n";
        for(int i=0;i<x.length()-1;i++)
        {
            outfile<<std::setprecision(16)<<x[i]<<"\n";
        }
        outfile<<std::setprecision(16)<<x[x.length()-1];
        outfile.close();
    }

    int nfaces = mesh.nFaces();
    int nverts = targetPos.rows();
    
    Eigen::VectorXd dE(3*nverts);
    dE.setZero();
    
     // Compute the energy f = 1/2||dE||^2, where E is the SV-energy
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix2d abar;
        abar << x[3*i], 0,
                x[3*i+1], x[3*i+2];
        abar = abar * abar.transpose();
        
        Eigen::Matrix2d abarinv = abar.inverse();
        double dA =  0.5 * sqrt(abar.determinant());

        // streching term
        double coeff = thickness / 4.0;
        Eigen::Matrix2d M = abarinv * (atargets[i] - abar);
        Eigen::Matrix<double, 1, 9> strechderiv;
        strechderiv.setZero();
        strechderiv += coeff*dA * lameAlpha * M.trace() * abarinv(0,0) * aderivs[i].row(0).transpose();
        strechderiv += coeff*dA * lameAlpha * M.trace() * abarinv(1,0) * aderivs[i].row(1).transpose();
        strechderiv += coeff*dA * lameAlpha * M.trace() * abarinv(0,1) * aderivs[i].row(2).transpose();
        strechderiv += coeff*dA * lameAlpha * M.trace() * abarinv(1,1) * aderivs[i].row(3).transpose();
        Eigen::Matrix2d Mainv = M*abarinv;
        strechderiv += coeff*dA* 2.0 * lameBeta * Mainv(0, 0) * aderivs[i].row(0).transpose();
        strechderiv += coeff*dA* 2.0 * lameBeta * Mainv(1, 0) * aderivs[i].row(1).transpose();
        strechderiv += coeff*dA* 2.0 * lameBeta * Mainv(0, 1) * aderivs[i].row(2).transpose();
        strechderiv += coeff*dA* 2.0 * lameBeta * Mainv(1, 1) * aderivs[i].row(3).transpose();
        for (int j = 0; j < 3; j++)
        {
            dE.segment<3>(3 * mesh.faceVertex(i,j)) += strechderiv.segment<3>(3 * j);
        }
        
        // bending term
        coeff = thickness*thickness*thickness / 12.0;
        M = abarinv * (btargets[i] - bbars[i]);
        Eigen::Matrix<double, 1, 18> bendingderiv;
        bendingderiv.setZero();
        bendingderiv += coeff*dA * lameAlpha * M.trace() * abarinv(0,0) * bderivs[i].row(0).transpose();
        bendingderiv += coeff*dA * lameAlpha * M.trace() * abarinv(1,0) * bderivs[i].row(1).transpose();
        bendingderiv += coeff*dA * lameAlpha * M.trace() * abarinv(0,1) * bderivs[i].row(2).transpose();
        bendingderiv += coeff*dA * lameAlpha * M.trace() * abarinv(1,1) * bderivs[i].row(3).transpose();
        Mainv = M*abarinv;
        bendingderiv += coeff*dA* 2.0 * lameBeta * Mainv(0, 0) * bderivs[i].row(0).transpose();
        bendingderiv += coeff*dA* 2.0 * lameBeta * Mainv(1, 0) * bderivs[i].row(1).transpose();
        bendingderiv += coeff*dA* 2.0 * lameBeta * Mainv(0, 1) * bderivs[i].row(2).transpose();
        bendingderiv += coeff*dA* 2.0 * lameBeta * Mainv(1, 1) * bderivs[i].row(3).transpose();
        for (int j = 0; j < 3; j++)
        {
            dE.segment<3>(3 * mesh.faceVertex(i,j)) += bendingderiv.segment<3>(3 * j);
            int oppidx = mesh.vertexOppositeFaceEdge(i,j);
            if(oppidx != -1)
                dE.segment<3>(3 * oppidx) += bendingderiv.segment<3>(9 + 3 * j);
        }
    }

    f = 0.5 * dE.transpose() * dE;
    

    // Compute df
    Eigen::SparseMatrix<double> gradvec(dE.rows(), 3*nfaces);
    std::vector<T> tripletlist;
    gradvec.setZero();
    Eigen::Matrix2d I;
    I.setIdentity();

    for(int i=0; i < nfaces; i++)
    {
        Eigen::Matrix2d abar, L;
        L << x[3*i], 0,
             x[3*i+1], x[3*i+2];
        abar = L * L.transpose();
        Eigen::Matrix2d abarinv = abar.inverse();
        double dA = 0.5 * sqrt(abar.determinant());
        Eigen::Matrix<double, 4, 3> abarinvderiv;
        Eigen::Vector3d abarsqrtdetderiv(3);

        computeInvMatDeriv(L, abarinvderiv);
        computeSqrtDetDerv(L, abarsqrtdetderiv);

        // Streching term
        double coeff = thickness / 4.0;
        Eigen::Matrix2d M = abarinv * (atargets[i] - abar);
        Eigen::Matrix<double, 1, 9> traceMderiv;
        traceMderiv.setZero();

        traceMderiv += abarinv(0,0) * aderivs[i].row(0).transpose();
        traceMderiv += abarinv(1,0) * aderivs[i].row(1).transpose();
        traceMderiv += abarinv(0,1) * aderivs[i].row(2).transpose();
        traceMderiv += abarinv(1,1) * aderivs[i].row(3).transpose();

        for(int k = 0; k < 3; k++)
        {
            Eigen::Matrix<double, 1, 9> result;
            Eigen::Matrix2d MderivL, MderivAinv, Mainvderiv;
            MderivL << abarinvderiv(0,k), abarinvderiv(1,k),
                        abarinvderiv(2,k), abarinvderiv(3,k);
            MderivL = MderivL * atargets[i];
            
            result.setZero();
            Eigen::Matrix<double, 1, 4> abarinvderivVeck;
            abarinvderivVeck << abarinvderiv(0,k), abarinvderiv(2,k), abarinvderiv(1,k), abarinvderiv(3, k);
            result += coeff*dA * lameAlpha * (MderivL.trace() * traceMderiv + M.trace() * abarinvderivVeck * aderivs[i]);
            /*
            trace(A*B) = A11*B11 + A12*B21 + A21*B12 + A22*B22 = (A11, A21, A12, A22) * (B11, B12, B21, B22)
            */
            
            
            
            MderivAinv =  MderivL * abarinv;
            Eigen::Matrix<double, 1, 4> MderivAinvVec, MainvderivVec;
            MderivAinvVec << MderivAinv(0,0), MderivAinv(1,0), MderivAinv(0,1), MderivAinv(1,1);

            Eigen::Matrix2d tmpAbarinv;
            tmpAbarinv << abarinvderiv(0,k), abarinvderiv(1,k),
                        abarinvderiv(2,k), abarinvderiv(3,k);
            Mainvderiv = M * tmpAbarinv;
            MainvderivVec <<  Mainvderiv(0,0),  Mainvderiv(1,0),  Mainvderiv(0,1),  Mainvderiv(1,1);
            

            result += coeff*dA* 2.0 * lameBeta * (MderivAinvVec * aderivs[i] +  MainvderivVec * aderivs[i]);

            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(0,0) * aderivs[i].row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(1,0) * aderivs[i].row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(0,1) * aderivs[i].row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(1,1) * aderivs[i].row(3).transpose();
            Eigen::Matrix2d Mainv = M*abarinv;
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(0, 0) * aderivs[i].row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(1, 0) * aderivs[i].row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(0, 1) * aderivs[i].row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(1, 1) * aderivs[i].row(3).transpose();

            for (int j = 0; j < 3; j++)
            {
                for(int r = 0; r < 3; r++)
                {
                    tripletlist.push_back(T(3 * mesh.faceVertex(i,j) + r, 3*i+k, result(3*j + r)));
                }
            }


        }
        
        // Bending term
        coeff = thickness * thickness * thickness / 12.0;
        M = abarinv * (btargets[i] - bbars[i]);

        Eigen::Matrix<double, 1, 18> traceMderivb;
        traceMderivb.setZero();
        
        traceMderivb += abarinv(0,0) * bderivs[i].row(0).transpose();
        traceMderivb += abarinv(1,0) * bderivs[i].row(1).transpose();
        traceMderivb += abarinv(0,1) * bderivs[i].row(2).transpose();
        traceMderivb += abarinv(1,1) * bderivs[i].row(3).transpose();

        for(int k = 0; k < 3; k++)
        {
            Eigen::Matrix<double, 1, 18> result;
            Eigen::Matrix2d MderivL, MderivAinv, Mainvderiv;
            MderivL << abarinvderiv(0,k), abarinvderiv(1,k),
            abarinvderiv(2,k), abarinvderiv(3,k);
            MderivL = MderivL * btargets[i];
            
            result.setZero();
            Eigen::Matrix<double, 1, 4> abarinvderivVeck;
            abarinvderivVeck << abarinvderiv(0,k), abarinvderiv(2,k), abarinvderiv(1,k), abarinvderiv(3, k);
            result += coeff*dA * lameAlpha * (MderivL.trace() * traceMderivb + M.trace() * abarinvderivVeck * bderivs[i]);
            /*
             trace(A*B) = A11*B11 + A12*B21 + A21*B12 + A22*B22 = (A11, A21, A12, A22) * (B11, B12, B21, B22)
             */
            
            
            
            MderivAinv =  MderivL * abarinv;
            Eigen::Matrix<double, 1, 4> MderivAinvVec, MainvderivVec;
            MderivAinvVec << MderivAinv(0,0), MderivAinv(1,0), MderivAinv(0,1), MderivAinv(1,1);
            
            Eigen::Matrix2d tmpAbarinv;
            tmpAbarinv << abarinvderiv(0,k), abarinvderiv(1,k),
            abarinvderiv(2,k), abarinvderiv(3,k);
            Mainvderiv = M * tmpAbarinv;
            MainvderivVec <<  Mainvderiv(0,0),  Mainvderiv(1,0),  Mainvderiv(0,1),  Mainvderiv(1,1);

            result += coeff*dA* 2.0 * lameBeta * (MderivAinvVec * bderivs[i] +  MainvderivVec * bderivs[i]);

            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(0,0) * bderivs[i].row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(1,0) * bderivs[i].row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(0,1) * bderivs[i].row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * 0.5 * lameAlpha * M.trace() * abarinv(1,1) * bderivs[i].row(3).transpose();
            Eigen::Matrix2d Mainv = M*abarinv;
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(0, 0) * bderivs[i].row(0).transpose();
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(1, 0) * bderivs[i].row(1).transpose();
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(0, 1) * bderivs[i].row(2).transpose();
            result += coeff * abarsqrtdetderiv(k) * lameBeta * Mainv(1, 1) * bderivs[i].row(3).transpose();

            for (int j = 0; j < 3; j++)
            {
                for(int r = 0; r < 3; r++)
                {
                    tripletlist.push_back(T(3 * mesh.faceVertex(i,j) + r, 3*i+k, result(3*j + r)));
                }

                int oppidx = mesh.vertexOppositeFaceEdge(i, j);
                if(oppidx != -1)
                {
                    for(int r = 0; r < 3; r++)
                    {
                        tripletlist.push_back(T(3 * oppidx + r, 3*i+k, result(9 + 3*j + r)));
                    }
                }
            }

        }
    }
    gradvec.setFromTriplets(tripletlist.begin(), tripletlist.end());
    Eigen::VectorXd gradf = gradvec.transpose()*dE;
    
    std::cout<<"Iteration Times: "<<itrTimes<<std::endl;
    std::cout<<"Elastic Term: "<<f<<"\t"<<"||Elastic Gradient||: "<<gradf.norm()<<std::endl;
    
    Eigen::VectorXd penaltyDerivative;
    double penaltyTerm = computePenaltyTerm(x, penaltyDerivative);
    
    std::cout<<"Abar Penalty Coefficient: "<<abarCoef<<"\t"<<"Penalty Term: "<<penaltyTerm<<"\t"<<"||Penalty Gradient||: "<<penaltyDerivative.norm()<<std::endl;
    
    f += abarCoef * penaltyTerm;
    gradf += abarCoef * penaltyDerivative;
    
    for(int i=0;i<df.length();i++)
    {
        df[i] = gradf(i);
    }
    
    std::cout<<"Objective function: "<<f<<"\t"<<"||Gradient||: "<<gradf.norm()<<std::endl;
}

double SimulationSetupAlglibSolver::computePenaltyTerm(const alglib::real_1d_array &x, Eigen::VectorXd &derivative)
{
    int nfaces =  mesh.nFaces();
    if(3*nfaces != x.length())
        return -1;
    derivative.resize(3 * nfaces);
    derivative.setZero();
    
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
        
        L << x[3*i], 0,
        x[3*i+1], x[3*i+2];
        for(int j = 0; j < 3; j++)
        {
            int oppVerIdx =  mesh.vertexOppositeFaceEdgeIndex(i, j);
            int oppFace = mesh.faceOppositeVertex(i, j);
            if (oppFace != -1)
            {
                Eigen::Matrix2d Lj;
                Lj << x[3*oppFace], 0,
                    x[3*oppFace+1], x[3*oppFace+2];
                
                // Compute the tranfermation matrix M
                
                Eigen::Vector3d oppNormal, curNormal, oppEdge;
                
                oppNormal = faceNormal(mesh, initialPos, oppFace, oppVerIdx, NULL, NULL);
                curNormal = faceNormal(mesh, initialPos, i, j, NULL, NULL);
                
                oppNormal = oppNormal/oppNormal.norm();
                curNormal = curNormal/curNormal.norm();
                
                oppEdge = initialPos.row(mesh.faceVertex(i, (j + 1)%3)) - initialPos.row(mesh.faceVertex(i, (j + 2)%3));
                oppEdge = oppEdge/oppEdge.norm();
                
                Eigen::Matrix3d A, A1, T;
                A.col(0) = oppEdge;
                A.col(1) = curNormal;
                A.col(2) = oppEdge.cross(curNormal);
                
                A1.col(0) = oppEdge;
                A1.col(1) = oppNormal;
                A1.col(2) = oppEdge.cross(oppNormal);
                
                T = A1*A.inverse();
                
                Eigen::Matrix<double, 3, 2> R;
                Eigen::Matrix<double, 3, 2> oppR;
                
                R.col(0) = initialPos.row(mesh.faceVertex(i, 1)) - initialPos.row(mesh.faceVertex(i, 0));
                R.col(1) = initialPos.row(mesh.faceVertex(i, 2)) - initialPos.row(mesh.faceVertex(i, 0));
                
                oppR.col(0) = initialPos.row(mesh.faceVertex(oppFace, 1)) - initialPos.row(mesh.faceVertex(oppFace, 0));
                oppR.col(1) = initialPos.row(mesh.faceVertex(oppFace, 2)) - initialPos.row(mesh.faceVertex(oppFace, 0));
                
                Eigen::Matrix2d M = (oppR.transpose() * oppR).inverse() * oppR.transpose() * T * R;
                
//                std::cout<<L * L.transpose()<<std::endl<<std::endl<<std::endl;
//                std::cout<<R.transpose()*R<<std::endl<<std::endl;
//                std::cout<<M.transpose() * Lj * Lj.transpose() * M<<std::endl<<std::endl;
//                std::cout<<M<<std::endl<<std::endl;
//
//                std::cout<<T * R<<std::endl<<std::endl;;
//                std::cout<<M(0,0) * oppR.col(0) + M(1,0) * oppR.col(1)<<std::endl<<std::endl;
//                std::cout<<M(0,1) * oppR.col(0) + M(1,1) * oppR.col(1)<<std::endl<<std::endl;
                
               
                
                 E += 1.0/2.0 * ( (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M) * (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M).transpose() ).trace() / sqrt(initialAbars[i].determinant()); // devided by det(A0) to make it scalar irrelavent
                
//               double value =  1.0/2.0 * (w.col(j).transpose() * L * L.transpose() * w.col(j) - w.col(oppVerIdx).transpose() * Lj * Lj.transpose() * w.col(oppVerIdx)) * (w.col(j).transpose() * L * L.transpose() * w.col(j) - w.col(oppVerIdx).transpose() * Lj * Lj.transpose() * w.col(oppVerIdx));
//
//               E += value / sqrt(initialAbars[i].determinant());
                

                for(int k = 0; k < 3; k++)
                {
                    Eigen::Matrix2d abarderiv, abarderivj;
                    abarderiv = Lderivs[k] * L.transpose() + L * Lderivs[k].transpose();
                    abarderivj = Lderivs[k] * Lj.transpose() + Lj * Lderivs[k].transpose();
                    
                     derivative(3 * i + k) += ( abarderiv * (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M).transpose() ).trace() / sqrt(initialAbars[i].determinant());

                     derivative(3 * oppFace + k) += - ( M.transpose() * abarderivj * M * (L * L.transpose() - M.transpose() * Lj * Lj.transpose() * M).transpose() ).trace() / sqrt(initialAbars[i].determinant());

//                   value = (w.col(j).transpose() * L * L.transpose() * w.col(j) - w.col(oppVerIdx).transpose() * Lj * Lj.transpose() * w.col(oppVerIdx)) * (w.col(j).transpose() * abarderiv * w.col(j));
//
//                   derivative(3 * i + k) += value / sqrt(initialAbars[i].determinant());
//
//                   value = -(w.col(j).transpose() * L * L.transpose() * w.col(j) - w.col(oppVerIdx).transpose() * Lj * Lj.transpose() * w.col(oppVerIdx)) * (w.col(oppVerIdx).transpose() * abarderivj * w.col(oppVerIdx));
//
//                   derivative(3 * oppFace + k) += value / sqrt(initialAbars[i].determinant());
                }
            }
            
        }
    }
    
    return E;
}

void SimulationSetupAlglibSolver::computeInvMatDeriv(Eigen::Matrix2d A, Eigen::Matrix<double, 4, 3> &dA)
{
    double x,y,z;
    x = A(0,0);
    y = A(1,0);
    z = A(1,1);
    Eigen::Matrix2d M, dA1, dA2, dA3;
    M << y*y+z*z,-x*y,
    -x*y,x*x;
    std::vector<Eigen::Matrix2d> C(3);
    C[0]<<0,-y,
    -y,2*x;
    C[1]<<2*y,-x,
    -x,0;
    C[2]<<2*z,0,
    0,0;
    
    dA1=1/(x*x*z*z)*C[0] - 2/(x*x*x*z*z)*M;
    dA2=1/(x*x*z*z)*C[1];
    dA3=1/(x*x*z*z)*C[2] - 2/(x*x*z*z*z)*M;

    dA.row(0) << dA1(0,0), dA2(0,0), dA3(0,0);
    dA.row(1) << dA1(0,1), dA2(0,1), dA3(0,1);
    dA.row(2) << dA1(1,0), dA2(1,0), dA3(1,0);
    dA.row(3) << dA1(1,1), dA2(1,1), dA3(1,1);
}

void SimulationSetupAlglibSolver::computeSqrtDetDerv(Eigen::Matrix2d A, Eigen::Vector3d & diffSqrtDet)
{
    int sign = 1;
    if(A.determinant()<0)
        sign = -1;
    diffSqrtDet(0) = sign*A(1,1);
    diffSqrtDet(1) = 0;
    diffSqrtDet(2) = sign*A(0,0);
}

void SimulationSetupAlglibSolver::testFucntionGradient(const SecondFundamentalFormDiscretization &sff)
{
    std::string resMeshName = "../../benchmarks/TestModels/sphere/draped_rect_geometry.obj";
//    std::string resMeshName = "../../benchmarks/TestModels/test_square.obj";
    Eigen::MatrixXi F, F1;
    Eigen::MatrixXd V, V1;
    if (!igl::readOBJ(resMeshName, V, F))
        return;
    initialPos = V;
    mesh = MeshConnectivity(F);
    
    std::string tarMeshName = "../../benchmarks/TestModels/sphere/sphere_geometry.obj";
//    std::string tarMeshName = "../../benchmarks/TestModels/test_square_bended.obj";
    if (!igl::readOBJ(tarMeshName, V1, F1))
        return;
    targetPos = V1;
    
    thickness = 0.1;

    int nfaces = mesh.nFaces();
    abars.resize(nfaces);
    bbars.resize(nfaces);
    atargets.resize(nfaces);
    btargets.resize(nfaces);
    aderivs.resize(nfaces);
    bderivs.resize(nfaces);
    
    lameAlpha = YoungsModulus * PoissonsRatio / (1.0 - PoissonsRatio * PoissonsRatio);
    lameBeta = YoungsModulus / 2.0 / (1.0 + PoissonsRatio);
    
    
    for (int i = 0; i < nfaces; i++)
    {
        bbars[i].setZero();
    }
    
    double E, E1;
    Eigen::VectorXd dE(nfaces*3),dE2,dE3;
    alglib::real_1d_array alg_x,alg_x1,alg_df, alg_df1;
    Eigen::VectorXd x(nfaces*3);
    alg_x.setlength(nfaces*3);
    alg_df.setlength(nfaces*3);
    
    for (int i = 0; i < nfaces; i++)
    {
        atargets[i] = firstFundamentalForm(mesh, targetPos, i, &aderivs[i], NULL);
        btargets[i] = sff.secondFundamentalForm(mesh, targetPos, initialEdgeDOFs, i, &bderivs[i], NULL);
        abars[i] = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
        
        alg_x[3*i] = sqrt(abars[i](0,0));
        alg_x[3*i+1] = abars[i](0,1)/alg_x[3*i];
        alg_x[3*i+2] = sqrt(abars[i].determinant())/alg_x[3*i];
    }
    
    getFunctionGradient(alg_x, E, alg_df);
    for(int i=0;i<alg_df.length();i++)
    {
        dE(i) = alg_df[i];
    }
    
    Eigen::VectorXd derivative;
    
    double energy = elasticEnergy(mesh, targetPos, initialEdgeDOFs, lameAlpha, lameBeta, thickness, abars, bbars, sff, &derivative, NULL);
//    std::cout<<energy<<std::endl;
//    std::cout<<E<<" "<<derivative.transpose() * derivative * 0.5 <<std::endl;
//    std::cout<<targetPos.rows()*3<<" "<<derivative.rows()<<std::endl;
    
    srand((unsigned)time(NULL));
    int selected_i = rand()%(nfaces*3);
    Eigen::VectorXd eps(nfaces*3);
    eps.setZero();
    alg_x1 = alg_x;
    for(int k=3;k<16;k++)
    {
        eps(selected_i) = pow(10,-k);
        alg_x1[selected_i] += pow(10,-k);
        getFunctionGradient(alg_x1, E1, alg_df1);
        std::cout<<"Selected index is: "<<selected_i<<" Eplison is: "<<eps(selected_i)<<std::endl;
        std::cout<<"finite difference is: "<<(E1-E)/(pow(10,-k))<<std::endl;
        std::cout<<"gradient projection: "<<dE.dot(eps)/eps(selected_i)<<std::endl;
        std::cout<<"The different between the finite difference and gradient is: "<<std::abs((E1-E)/pow(10,-k) - dE.dot(eps)/eps(selected_i))<<std::endl;
        alg_x1 = alg_x;
        
    }
}

void SimulationSetupAlglibSolver::testPenaltyTerm()
{
    std::string testMeshName = "../../benchmarks/TestModels/sphere/sphere_geometry.obj";
//    std::string testMeshName = "../../benchmarks/TestModels/backup/test_square.obj";
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(testMeshName, V, F);
    MeshConnectivity testMesh(F);
    mesh = testMesh;
    std::vector<Eigen::Matrix2d> abarTest;
    initialPos = V;
    int nfaces = F.rows();
    abarTest.resize(nfaces);
    alglib::real_1d_array x;
    x.setlength(3*nfaces);
    
    for(int i=0; i<nfaces; i++)
    {
        abarTest[i] = firstFundamentalForm(testMesh, V, i, NULL, NULL);
        initialAbars[i] = abarTest[i];
        srand((unsigned)time(NULL));
//        x[3*i] = sqrt(abarTest[i](0,0)) + 1e-3 * rand()/RAND_MAX;
//        x[3*i+1] = abarTest[i](0,1) / x[3*i] + 1e-3 * rand()/RAND_MAX;
//        x[3*i+2] = sqrt(abarTest[i].determinant())/x[3*i] + 1e-3 * rand()/RAND_MAX;
        
        x[3*i] = sqrt(abarTest[i](0,0));
        x[3*i+1] = abarTest[i](0,1) / x[3*i];
        x[3*i+2] = sqrt(abarTest[i].determinant())/x[3*i];
    }
    
    Eigen::VectorXd derivative, derivative1;
    
    double E = computePenaltyTerm(x, derivative);
    std::cout<<"The energy is "<<E<<std::endl;
    std::cout<<"The norm of derivative is "<<derivative.norm()<<std::endl;
    
    srand((unsigned)time(NULL));
    int selected_i = rand()%(nfaces*3);
    Eigen::VectorXd eps(nfaces*3);
    eps.setZero();
    alglib::real_1d_array x1 = x;
    for(int k=1;k<16;k++)
    {
        eps(selected_i) = pow(10,-k);
        x1[selected_i] += pow(10,-k);
        double E1 = computePenaltyTerm(x1, derivative1);
        std::cout<<"Selected index is: "<<selected_i<<" Eplison is: "<<eps(selected_i)<<std::endl;
        std::cout<<"finite difference is: "<<(E1-E)/(pow(10,-k))<<std::endl;
        std::cout<<"gradient projection: "<<derivative.dot(eps)/eps(selected_i)<<std::endl;
        std::cout<<"The different between the finite difference and gradient is: "<<std::abs((E1-E)/pow(10,-k) - derivative.dot(eps)/eps(selected_i))<<std::endl;
        x1 = x;
        
    }
    
    std::cout<< E << std::endl;
}
