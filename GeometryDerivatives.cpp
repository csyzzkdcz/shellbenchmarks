#include "GeometryDerivatives.h"
#include "MeshConnectivity.h"
#include <iostream>
#include <random>
#include <Eigen/Geometry>

Eigen::Matrix3d crossMatrix(Eigen::Vector3d v)
{
    Eigen::Matrix3d ret;
    ret << 0, -v[2], v[1],
        v[2], 0, -v[0],
        -v[1], v[0], 0;
    return ret;
}

double angle(const Eigen::Vector3d &v, const Eigen::Vector3d &w, const Eigen::Vector3d &axis,
    Eigen::Matrix<double, 1, 6> *derivative, // v, w
    Eigen::Matrix<double, 6, 6> *hessian
)
{
    double theta = 2.0 * atan2((v.cross(w).dot(axis)), v.dot(w) + v.norm() * w.norm());

    if (derivative)
    {
        derivative->segment(0, 3) = -axis.cross(v) / v.squaredNorm();
        derivative->segment(3, 3) = axis.cross(w) / w.squaredNorm();
    }
    if (hessian)
    {
        hessian->setZero();
        hessian->block(0, 0, 3, 3) += 2.0 * (axis.cross(v))*v.transpose() / v.squaredNorm() / v.squaredNorm();
        hessian->block(3, 3, 3, 3) += -2.0 * (axis.cross(w))*w.transpose() / w.squaredNorm() / w.squaredNorm();
        hessian->block(0, 0, 3, 3) += -crossMatrix(axis) / v.squaredNorm();
        hessian->block(3, 3, 3, 3) += crossMatrix(axis) / w.squaredNorm();

        double sigma = 1.0;
        if (v.cross(w).dot(axis) < 0)
            sigma = -1.0;

        double vwnorm = v.cross(w).norm();
        if (vwnorm > 1e-8)
        {
            Eigen::Matrix3d da = sigma * (1.0 / vwnorm * Eigen::Matrix3d::Identity() - 1.0 / vwnorm / vwnorm / vwnorm * (v.cross(w)) * (v.cross(w)).transpose());
            hessian->block(0, 0, 3, 3) += crossMatrix(v) / v.squaredNorm() * da * -crossMatrix(w);
            hessian->block(3, 0, 3, 3) += crossMatrix(v) / v.squaredNorm() * da * crossMatrix(v);
            hessian->block(0, 3, 3, 3) += -crossMatrix(w) / w.squaredNorm() * da * -crossMatrix(w);
            hessian->block(3, 3, 3, 3) += -crossMatrix(w) / w.squaredNorm() * da * crossMatrix(v);
        }
    }

    return theta;
}

Eigen::Vector3d faceNormal(const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    int face, int startidx,
    Eigen::Matrix<double, 3, 9> *derivative,
    std::vector<Eigen::Matrix<double, 9, 9> > *hessian)
{
    if (derivative)
        derivative->setZero();

    if (hessian)
    {
        hessian->resize(3);
        for (int i = 0; i < 3; i++) (*hessian)[i].setZero();
    }

    int v0 = startidx % 3;
    int v1 = (startidx + 1) % 3;
    int v2 = (startidx + 2) % 3;
    Eigen::Vector3d qi0 = curPos.row(mesh.faceVertex(face, v0)).transpose();
    Eigen::Vector3d qi1 = curPos.row(mesh.faceVertex(face, v1)).transpose();
    Eigen::Vector3d qi2 = curPos.row(mesh.faceVertex(face, v2)).transpose();
    Eigen::Vector3d n = (qi1 - qi0).cross(qi2 - qi0);

    if (derivative)
    {
        derivative->block(0, 0, 3, 3) += crossMatrix(qi2 - qi1);
        derivative->block(0, 3, 3, 3) += crossMatrix(qi0 - qi2);
        derivative->block(0, 6, 3, 3) += crossMatrix(qi1 - qi0);
    }

    if (hessian)
    {
        for (int j = 0; j < 3; j++)
        {
            Eigen::Vector3d ej(0, 0, 0);
            ej[j] = 1.0;
            Eigen::Matrix3d ejc = crossMatrix(ej);
            (*hessian)[j].block(0, 3, 3, 3) -= ejc;
            (*hessian)[j].block(0, 6, 3, 3) += ejc;
            (*hessian)[j].block(3, 6, 3, 3) -= ejc;
            (*hessian)[j].block(3, 0, 3, 3) += ejc;
            (*hessian)[j].block(6, 0, 3, 3) -= ejc;
            (*hessian)[j].block(6, 3, 3, 3) += ejc;
        }
    }

    return n;
}

double triangleAltitude(const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    int face,
    int edgeidx,
    Eigen::Matrix<double, 1, 9> *derivative,
    Eigen::Matrix<double, 9, 9> *hessian)
{
    if (derivative)
        derivative->setZero();
    if (hessian)
        hessian->setZero();

    Eigen::Matrix<double, 3, 9> nderiv;
    std::vector<Eigen::Matrix<double, 9, 9> > nhess;
    Eigen::Vector3d n = faceNormal(mesh, curPos, face, edgeidx, (derivative || hessian ? &nderiv : NULL), hessian ? &nhess : NULL);

    int v2 = (edgeidx + 2) % 3;
    int v1 = (edgeidx + 1) % 3;
    Eigen::Vector3d q2 = curPos.row(mesh.faceVertex(face, v2)).transpose();
    Eigen::Vector3d q1 = curPos.row(mesh.faceVertex(face, v1)).transpose();

    Eigen::Vector3d e = q2 - q1;
    double nnorm = n.norm();
    double enorm = e.norm();
    double h = nnorm / enorm;

    if (derivative)
    {
        for (int i = 0; i < 3; i++)
        {
            *derivative += nderiv.row(i) * n[i] / nnorm / enorm;
        }
        derivative->block(0, 6, 1, 3) += -nnorm / enorm / enorm / enorm * e.transpose();
        derivative->block(0, 3, 1, 3) += nnorm / enorm / enorm / enorm * e.transpose();
    }

    if (hessian)
    {
        for (int i = 0; i < 3; i++)
        {
            *hessian += nhess[i] * n[i] / nnorm / enorm;
        }
        Eigen::Matrix3d P = Eigen::Matrix3d::Identity() / nnorm - n*n.transpose() / nnorm / nnorm / nnorm;
        *hessian += nderiv.transpose() * P * nderiv / enorm;
        hessian->block(6, 0, 3, 9) += -e * n.transpose() * nderiv / nnorm / enorm / enorm / enorm;
        hessian->block(3, 0, 3, 9) += e * n.transpose() * nderiv / nnorm / enorm / enorm / enorm;
        hessian->block(0, 6, 9, 3) += -nderiv.transpose() * n * e.transpose() / nnorm / enorm / enorm / enorm;
        hessian->block(0, 3, 9, 3) += nderiv.transpose() * n * e.transpose() / nnorm / enorm / enorm / enorm;
        hessian->block(6, 6, 3, 3) += -nnorm / enorm / enorm / enorm * Eigen::Matrix3d::Identity();
        hessian->block(6, 3, 3, 3) += nnorm / enorm / enorm / enorm * Eigen::Matrix3d::Identity();
        hessian->block(3, 6, 3, 3) += nnorm / enorm / enorm / enorm * Eigen::Matrix3d::Identity();
        hessian->block(3, 3, 3, 3) += -nnorm / enorm / enorm / enorm * Eigen::Matrix3d::Identity();
        Eigen::Matrix3d outer = e*e.transpose() * 3.0 * nnorm / enorm / enorm / enorm / enorm / enorm;
        hessian->block(6, 6, 3, 3) += outer;
        hessian->block(6, 3, 3, 3) += -outer;
        hessian->block(3, 6, 3, 3) += -outer;
        hessian->block(3, 3, 3, 3) += outer;
    }

    return h;
}

Eigen::Matrix2d firstFundamentalForm(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd curPos,
    int face,
    Eigen::Matrix<double, 4, 9> *derivative, // F(face, i)
    std::vector <Eigen::Matrix<double, 9, 9> > *hessian)
{
    Eigen::Vector3d q0 = curPos.row(mesh.faceVertex(face, 0));
    Eigen::Vector3d q1 = curPos.row(mesh.faceVertex(face, 1));
    Eigen::Vector3d q2 = curPos.row(mesh.faceVertex(face, 2));
    Eigen::Matrix2d result;
    result << (q1 - q0).dot(q1 - q0), (q1 - q0).dot(q2 - q0),
        (q2 - q0).dot(q1 - q0), (q2 - q0).dot(q2 - q0);

    if (derivative)
    {
        derivative->setZero();
        derivative->block<1, 3>(0, 3) += 2.0 * (q1 - q0).transpose();
        derivative->block<1, 3>(0, 0) -= 2.0 * (q1 - q0).transpose();
        derivative->block<1, 3>(1, 6) += (q1 - q0).transpose();
        derivative->block<1, 3>(1, 3) += (q2 - q0).transpose();
        derivative->block<1, 3>(1, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
        derivative->block<1, 3>(2, 6) += (q1 - q0).transpose();
        derivative->block<1, 3>(2, 3) += (q2 - q0).transpose();
        derivative->block<1, 3>(2, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
        derivative->block<1, 3>(3, 6) += 2.0 * (q2 - q0).transpose();
        derivative->block<1, 3>(3, 0) -= 2.0 * (q2 - q0).transpose();
    }

    if (hessian)
    {
        hessian->resize(4);
        for (int i = 0; i < 4; i++)
        {
            (*hessian)[i].setZero();
        }
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        (*hessian)[0].block<3, 3>(0, 0) += 2.0*I;
        (*hessian)[0].block<3, 3>(3, 3) += 2.0*I;
        (*hessian)[0].block<3, 3>(0, 3) -= 2.0*I;
        (*hessian)[0].block<3, 3>(3, 0) -= 2.0*I;

        (*hessian)[1].block<3, 3>(3, 6) += I;
        (*hessian)[1].block<3, 3>(6, 3) += I;
        (*hessian)[1].block<3, 3>(0, 3) -= I;
        (*hessian)[1].block<3, 3>(0, 6) -= I;
        (*hessian)[1].block<3, 3>(3, 0) -= I;
        (*hessian)[1].block<3, 3>(6, 0) -= I;
        (*hessian)[1].block<3, 3>(0, 0) += 2.0*I;

        (*hessian)[2].block<3, 3>(3, 6) += I;
        (*hessian)[2].block<3, 3>(6, 3) += I;
        (*hessian)[2].block<3, 3>(0, 3) -= I;
        (*hessian)[2].block<3, 3>(0, 6) -= I;
        (*hessian)[2].block<3, 3>(3, 0) -= I;
        (*hessian)[2].block<3, 3>(6, 0) -= I;
        (*hessian)[2].block<3, 3>(0, 0) += 2.0*I;

        (*hessian)[3].block<3, 3>(0, 0) += 2.0*I;
        (*hessian)[3].block<3, 3>(6, 6) += 2.0*I;
        (*hessian)[3].block<3, 3>(0, 6) -= 2.0*I;
        (*hessian)[3].block<3, 3>(6, 0) -= 2.0*I;
    }

    return result;
}


void testAngle()
{
    double eps = 1e-6;
    int ntests = 100;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::normal_distribution<double> dis;

    std::cout << "Testing " << ntests << " random vectors" << std::endl;

    for (int i = 0; i < ntests; i++)
    {
        Eigen::Vector3d v;
        Eigen::Vector3d w;
        for (int i = 0; i < 3; i++)
        {
            v[i] = dis(rng);
            w[i] = dis(rng);
        }
        double dir = dis(rng);
        dir = (dir < 0 ? -1.0 : 1.0);
        Eigen::Vector3d axis = dir * v.cross(w) / v.cross(w).norm();
        Eigen::Matrix<double, 1, 6> deriv;
        Eigen::Matrix<double, 6, 6> hess;
        double theta = angle(v, w, axis, &deriv, &hess);        
        for (int i = 0; i < 3; i++)
        {
            Eigen::Vector3d vpert = v;

            vpert[i] += eps;
            Eigen::Vector3d axispert = dir * vpert.cross(w) / vpert.cross(w).norm();
            Eigen::Matrix<double, 1, 6> derivpert;
            double newtheta = angle(vpert, w, axispert, &derivpert, NULL);
            double exact = deriv[i];
            double findiff = (newtheta - theta) / eps;
            std::cout << "v[" << i << "]: " << exact << " / " << findiff << std::endl;
            Eigen::Matrix<double, 1, 6> findiffhess = (derivpert - deriv) / eps;
            for (int j = 0; j < 6; j++)
            {
                std::cout << "hess[" << j << "]: " << hess(i, j) << " / " << findiffhess[j] << std::endl;
            }
        }
        for (int i = 0; i < 3; i++)
        {
            Eigen::Vector3d wpert = w;

            wpert[i] += eps;
            Eigen::Vector3d axispert = dir * v.cross(wpert) / v.cross(wpert).norm();
            Eigen::Matrix<double, 1, 6> derivpert;
            double newtheta = angle(v, wpert, axispert,  &derivpert, NULL);
            double exact = deriv[3+i];
            double findiff = (newtheta - theta) / eps;
            std::cout << "w[" << i << "]: " << exact << " / " << findiff << std::endl;
            Eigen::Matrix<double, 1, 6> findiffhess = (derivpert - deriv) / eps;
            for (int j = 0; j < 6; j++)
            {
                std::cout << "hess[" << j << "]: " << hess(3+i, j) << " / " << findiffhess[j] << std::endl;
            }
        }
    }   
}

void testFaceNormals(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    double eps = 1e-6;
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int ntests = 100;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, nfaces);
    std::uniform_int_distribution<int> idxdist(0, 2);

    std::cout << "Testing " << ntests << " random faces" << std::endl;

    for(int i=0; i<ntests; i++)
    {
        int face = uni(rng);
        int idx = idxdist(rng);
        Eigen::Matrix<double, 3, 9> deriv;
        std::vector<Eigen::Matrix<double, 9, 9> > hess;
        Eigen::Vector3d n = faceNormal(mesh, V, face, idx, &deriv, &hess);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(mesh.faceVertex(face, (idx + j)%3), k) += 1e-6;
                Eigen::Matrix<double, 3, 9> derivpert;
                Eigen::Vector3d npert = faceNormal(mesh, Vpert, face, idx, &derivpert, NULL);
                Eigen::Vector3d findiff = (npert-n)/1e-6;
                Eigen::Vector3d exact = deriv.col(3*j+k);
                std::cout << "q" << j << "[" << k <<"]: " << exact.transpose() << " / " << findiff.transpose() << std::endl;

                Eigen::Matrix<double, 3, 9> findiffhess = (derivpert-deriv)/1e-6;
                for(int l=0; l<3; l++)
                {
                    std::cout << " hess[" << l << "]: " << hess[l].col(3*j+k).transpose() << " / " << findiffhess.row(l) << std::endl;;
                }
            }            
        }
    }
}

void testAltitudes(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    double eps = 1e-6;
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int ntests = 100;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, nfaces);
    std::uniform_int_distribution<int> idxdist(0, 2);

    std::cout << "Testing " << ntests << " random altitudes" << std::endl;

    for(int i=0; i<ntests; i++)
    {
        int face = uni(rng);
        int idx = idxdist(rng);
        Eigen::Matrix<double, 1, 9> deriv;
        Eigen::Matrix<double, 9, 9> hess;
        double h = triangleAltitude(mesh, V, face, idx, &deriv, &hess);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(mesh.faceVertex(face, (idx + j)%3), k) += 1e-6;
                Eigen::Matrix<double, 1, 9> derivpert;
                double hpert = triangleAltitude(mesh, Vpert, face, idx, &derivpert, NULL);
                double findiff = (hpert-h)/1e-6;
                double exact = deriv[3 * j + k];
                std::cout << "q" << j << "[" << k <<"]: " << exact << " / " << findiff << std::endl;

                Eigen::Matrix<double, 1, 9> findiffhess = (derivpert-deriv)/1e-6;
                std::cout << " hess: " << hess.col(3*j+k).transpose() << " / " << findiffhess << std::endl;;                
            }            
        }
    }
}


void testFirstFundamentalForm(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    double eps = 1e-6;
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    int nedges = mesh.nEdges();
    int ntests = 100;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, nfaces);

    std::cout << "Testing " << ntests << " random faces" << std::endl;

    for(int i=0; i<ntests; i++)
    {
        int face = uni(rng);
        Eigen::Matrix<double, 4, 9> deriv;
        std::vector<Eigen::Matrix<double, 9, 9> > hess;
        Eigen::Matrix2d a = firstFundamentalForm(mesh, V, face, &deriv, &hess);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                Eigen::MatrixXd Vpert(V);
                Vpert(mesh.faceVertex(face, j), k) += 1e-6;
                Eigen::Matrix<double, 4, 9> derivpert;
                Eigen::Matrix2d apert = firstFundamentalForm(mesh, Vpert, face, &derivpert, NULL);
                Eigen::Matrix2d findiff = (apert-a)/1e-6;
                Eigen::Vector4d findiffvec;
                findiffvec << findiff(0,0), findiff(0,1), findiff(1,0), findiff(1,1);
                Eigen::Vector4d exact = deriv.col(3*j+k);
                std::cout << "q" << j << "[" << k <<"]: " << exact.transpose() << " / " << findiffvec.transpose() << std::endl;

                Eigen::Matrix<double, 4, 9> findiffhess = (derivpert-deriv)/1e-6;
                for(int l=0; l<4; l++)
                {
                    std::cout << " hess[" << l << "]: " << hess[l].col(3*j+k).transpose() << " / " << findiffhess.row(l) << std::endl;;
                }
            }            
        }
    }
}

