#include <igl/opengl/glfw/Viewer.h>
#include <random>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/decimate.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/upsample.h>
#include <igl/doublearea.h>
#include <igl/arap.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/file_dialog_open.h>
#include <igl/colormap.h>
#include <igl/edges.h>
#include <igl/lscm.h>
#include <imgui/imgui.h>
#include <igl/triangle/triangulate.h>
#include <igl/hessian_energy.h>
#include <igl/massmatrix.h>
#include <igl/barycenter.h>
#include <igl/boundary_loop.h>
#include <memory>
#include <iomanip>
#include <fstream>
//#include "SecondFundamentalForm/MidedgeAngleSinFormulation.h"
//#include "SecondFundamentalForm/MidedgeAngleTanFormulation.h"
#include "SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "MeshConnectivity.h"
#include "GeometryDerivatives.h"
//#include "BuildModels.h"
//#include "TestModels.h"
#include "ElasticShell.h"
#include "SimulationSetup/SimulationSetup.h"
#include "SimulationSetup/SimulationSetupNormal.h"
#include "SimulationSetup/SimulationSetupAlglibSolver.h"
#include "SimulationSetup/SimulationSetupIpoptSolver.h"
#include "SimulationSetup/SimulationSetupDynamicSolver.h"
#include "ParseWimFiles.h"
#include "StaticSolve.h"
#include "SimulationState.h"
#include "GeometryTools.h"
#include "SensitiveAnalysis/SensitiveAnalysisAbarBbar.h"
#include "SensitiveAnalysis/SensitiveAnalysisABbarPos.h"




std::unique_ptr<SimulationSetup> setup;
std::vector<Eigen::Matrix2d> oldAbars;
SimulationState curState;
Eigen::VectorXd evec;
Eigen::MatrixXd _initialPos;
MeshConnectivity _initialMesh;
int numSteps;
double tolerance;
bool isShowVerField = false;
bool isShowFaceColor = false;
bool isShowAbar =  false;
bool isStartFromMiddleShape = false;
bool isPlotForce = false;
bool isOverWrite = false;
bool isContinue = false;

enum Methods {Normal=0, alglibSolver, ipoptSolver, dynamicSolver};
static Methods methodType =  dynamicSolver;

enum faceColorTypes {optimal2target, optimal2plate, bbar, topLayer, bottomLayer};
static faceColorTypes faceColorType =  optimal2target;

enum paramatrizationTypes {ARAP = 0, Conformal};
static paramatrizationTypes paramatrizationType = Conformal;

enum dynamicSolverTypes {AbarPos = 0, AbarBbar, ABbarPos};
static dynamicSolverTypes dynamicSolverType = ABbarPos;

std::string selectedMethod = "dynamicSolver";
std::string selectedDynamicType = "ABbarPos";

std::string resShape = "";
std::string tarShape = "";
std::string curPath = "";
std::string selectedType = "sphere";

bool isFixedConer = false;

double thickness = 1e-4;
double abarCoef = 0;
double bbarCoef = 0;
double smoothnessCoef = 0;
int resampledResolution = 1000;
int expCoef;


void testTriangulation()
{
    Eigen::MatrixXd initialPos;
    Eigen::MatrixXi F;

    igl::readOBJ("/Users/chenzhen/UT/Research/Projects/shellbenchmarks/benchmarks/TestModels/coarse/Nefertiti/draped_rect_geometry.obj", initialPos, F);
    MeshConnectivity mesh(F);
    int nfaces = mesh.nFaces();
    std::vector<Eigen::Matrix2d> convertedAbars(nfaces);
    std::vector<Eigen::Matrix2d> newAbars;
    std::vector<Eigen::Matrix2d> newBbars;

    Eigen::VectorXi boundaryLoop;
    igl::boundary_loop(mesh.faces(), boundaryLoop);
    Eigen::MatrixXd boundaryVert(boundaryLoop.size(),2);
    Eigen::MatrixXi boundaryEdge(boundaryLoop.size(),2);

    for(int i=0;i<boundaryLoop.size();i++)
    {
        int vertid = boundaryLoop(i);
        boundaryVert.row(i) << initialPos(vertid,0), initialPos(vertid, 1);
    }
    for(int i=0;i<boundaryLoop.size()-1;i++)
    {
        int startid = boundaryLoop(i);
        int endid = boundaryLoop(i+1);
        boundaryEdge.row(i) << startid, endid;
    }
    boundaryEdge.row(boundaryLoop.size()-1) << boundaryLoop(boundaryLoop.size() - 1), boundaryLoop(0);

    std::cout<<boundaryLoop<<std::endl<<std::endl;
    std::cout<<std::endl<<boundaryEdge<<std::endl<<std::endl;
    std::cout<<std::endl<<boundaryVert<<std::endl<<std::endl;

    Eigen::MatrixXd remeshedPos;
    Eigen::MatrixXi remeshedFaces;
    Eigen::MatrixXd H;
    H.resize(0, 0);
    igl::triangle::triangulate(boundaryVert,boundaryEdge,H,"a0.0005q",remeshedPos,remeshedFaces);

    Eigen::MatrixXd V(remeshedPos.rows(),3);
    V.setZero();
    V.block(0, 0, remeshedPos.rows(), 2) = remeshedPos;
    igl::writeOBJ("remeshed.obj", V, remeshedFaces);

    std::ofstream outfile("boundary.txt", std::ios::trunc);

    outfile<<boundaryVert.rows()<<"\n";

    for(int i=0;i<boundaryVert.rows();i++)
    {
        outfile<<std::setprecision(16)<<boundaryVert(i, 0)<<" "<<boundaryVert(i,1)<<" "<<0<<"\n";
    }

}

void meshResampling(std::string filepath, int targetResolution)
{
    // Plot the generated mesh
    Eigen::MatrixXd Vcurr;
    Eigen::MatrixXi Fcurr;
    igl::readOBJ(filepath, Vcurr, Fcurr);
    Eigen::VectorXi bc;
    igl::boundary_loop(Fcurr, bc);
    for(int i=0;i<bc.size();i++)
    {
        std::cout<<Vcurr.row(bc(i))<<std::endl;
    }
    while ( Fcurr.rows() < targetResolution * 2)
    {
        Eigen::MatrixXd Vtmp = Vcurr;
        Eigen::MatrixXi Ftmp = Fcurr;
        igl::upsample(Vtmp, Ftmp, Vcurr, Fcurr, 1);
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi J;

    igl::decimate(Vcurr, Fcurr, targetResolution, V, F, J);
    int ind = filepath.rfind(".");

    typedef Eigen::SparseMatrix<double> SparseMat;

    //Read our mesh
    Eigen::MatrixXi E;
    igl::edges(F,E);

    //Constructing an exact function to smooth
    Eigen::VectorXd zexact = V.block(0,2,V.rows(),1).array()
    + 0.5*V.block(0,1,V.rows(),1).array()
    + V.block(0,1,V.rows(),1).array().pow(2)
    + V.block(0,2,V.rows(),1).array().pow(3);

    //Make the exact function noisy
    srand(5);
    const double s = 0.2*(zexact.maxCoeff() - zexact.minCoeff());
    Eigen::VectorXd znoisy = zexact + s*Eigen::VectorXd::Random(zexact.size());

    //Constructing the squared Laplacian and squared Hessian energy
    SparseMat L, M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    Eigen::SimplicialLDLT<SparseMat> solver(M);
    SparseMat MinvL = solver.solve(L);
    SparseMat QL = L.transpose()*MinvL;
    SparseMat QH;
    igl::hessian_energy(V, F, QH);
    const double al = 8e-4;
    Eigen::SimplicialLDLT<SparseMat> lapSolver(al*QL + (1.-al)*M);
    Eigen::VectorXd zl = lapSolver.solve(al*M*znoisy);
    const double ah = 5e-6;
    Eigen::SimplicialLDLT<SparseMat> hessSolver(ah*QH + (1.-ah)*M);
    Eigen::VectorXd zh = hessSolver.solve(ah*M*znoisy);
    igl::boundary_loop(F, bc);
//    for(int i=0;i<bc.size();i++)
//    {
//        V(bc(i), 0) = 0;
//    }
    igl::writeOBJ(filepath.substr(0, ind) + "_resampled.obj", V, F);

}


void computeSphere(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        /*
         x = R*sin(u / R)
         y = v
         z = R*cos(u / R) - R
         */
        double R = 0.5;
        double u = Vo(i,0);
        double v = Vo(i,1);
        double z = R - R*R/sqrt(R*R+u*u+v*v);
        double x = (R-z)/R*u;
        double y = (R-z)/R*v;
        Vo(i,0) = 0.5*x;
        Vo(i,1) = 0.5*y;
        Vo(i,2) = 0.5*z;
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/sphere_geometry.obj", Vo, Fo);

}

void computeEllipsoid(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        double R = 0.5;
        double u = Vo(i,0);
        double v = Vo(i,1);
        double z = R - R*R/sqrt(R*R+u*u+v*v);
        double x = (R-z)/R*u * 0.5;
        double y = (R-z)/R*v;
        Vo(i,0) = R*x;
        Vo(i,1) = R*y;
        Vo(i,2) = R*z;
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/ellipsoid_geometry.obj", Vo, Fo);

}


void computeCylinder(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        /*
         x = R*sin(u / R)
         y = v
         z = R*cos(u / R) - R
         */
        double R = 0.5;
        double x = R*sin(Vo(i,0)/R);
        double z = R*cos(Vo(i,0)/R) - R;
        Vo(i,0) = x;
        Vo(i,2) = z;
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/cylinder_geometry.obj", Vo, Fo);
}

void computeSaddle(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        double x = Vo(i,0);
        double y = Vo(i,1);
        double z = x*x -y*y;
        Vo(i,2) = z;
    }
    int ind = rectPath.rfind("/");
    igl::writeOBJ(rectPath.substr(0, ind) + "/saddle_geometry.obj", Vo, Fo);
}


void postProcess()
{
    for (PostprocessTest &t : setup->tests)
    {
        double disp = curState.curPos(t.vertex, t.coord) - setup->initialPos(t.vertex, t.coord);
        std::cout << t.vertex << '\t' << t.coord << '\t' << t.wimDisplacement << "\tvs\t" << disp << std::endl;
    }
}

void jitter(double magnitude)
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-magnitude, magnitude);
    for (int i = 0; i < curState.curPos.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            curState.curPos(i, j) += distribution(generator);
        }
    }
}

void reset()
{
    std::cout << std::endl << "Reset" << std::endl << std::endl;
    curState.curPos = setup->initialPos;
    curState.curEdgeDOFs = setup->initialEdgeDOFs;
    evec.resize(curState.curPos.rows());
    evec.setZero();

}

void restore()
{
    std::cout << std::endl << "Restore" << std::endl << std::endl;
    curState.curPos = _initialPos;
    curState.curEdgeDOFs = setup->initialEdgeDOFs;
    evec.resize(curState.curPos.rows());
    evec.setZero();
}

void setTarget()
{
    std::cout << std::endl << "Set Target" << std::endl << std::endl;
    curState.curPos = setup->targetPos;
    curState.curEdgeDOFs = setup->initialEdgeDOFs;
    numSteps = 1;
}

void setMiddle()
{
    std::cout << std::endl << "Set Middle" << std::endl << std::endl;
    curState.curPos = setup->targetPosAfterFirstStep;
    curState.curEdgeDOFs = setup->initialEdgeDOFs;
    numSteps = 1;

}

void repaint(igl::opengl::glfw::Viewer &viewer)
{
    viewer.data().clear();
    viewer.core.background_color<<1,1,1,1;
    if(curState.curPos.rows() == 0)
        return;
    viewer.data().set_mesh(curState.curPos, setup->mesh.faces());
    Eigen::MatrixXd colors(setup->initialPos.rows(), 3);

    colors.col(0).setConstant(1.0);
    colors.col(1).setConstant(1.0);
    colors.col(2).setConstant(0);

    viewer.data().set_colors(colors);
    viewer.data().line_width = 2;

    if(isPlotForce != 0)
    {
        MidedgeAverageFormulation sff;
        SensitiveAnalysisAbarBbar op;
        bool ok = parseWimFiles(resShape, tarShape, *setup, sff);
        setup->thickness = thickness;
        setup->abarCoef = abarCoef;
        setup->bbarCoef = bbarCoef;
        setup->smoothCoef = smoothnessCoef;
        double lameAlpha = setup->YoungsModulus * setup->PoissonsRatio / (1.0 - setup->PoissonsRatio * setup->PoissonsRatio);
        double lameBeta = setup->YoungsModulus / 2.0 / (1.0 + setup->PoissonsRatio);
        op.initialization(setup->initialPos, setup->targetPos, setup->mesh, setup->clampedDOFs, lameAlpha, lameBeta, setup->thickness);
        op.setPenalty(abarCoef, bbarCoef, smoothnessCoef);
        Eigen::VectorXd inplaneForce = op.getInplaneForce();
        Eigen::MatrixXd force(viewer.data().V.rows(), 3);
        for(int i=0;i<force.rows();i++)
        {
            force.row(i) = inplaneForce.segment(3*i, 3);
        }
        force = force / force.norm() * viewer.data().V.norm();
        const Eigen::RowVector3d red(1,0,0), black(0,0,0);
        Eigen::MatrixXd V = setup->targetPos;
        Eigen::MatrixXi F = setup->mesh.faces();
        Eigen::MatrixXd FN, VN;
        igl::per_face_normals(V, F, FN);
        igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, VN);
        VN = VN / VN.norm() * viewer.data().V.norm();
        viewer.data().add_edges(viewer.data().V,viewer.data().V+force, red);
        viewer.data().add_edges(viewer.data().V, viewer.data().V+VN, black);
        Eigen::VectorXd theta(force.rows());
        for(int i=0;i<force.rows();i++)
        {
            double cos = force.row(i).dot(VN.row(i)) / (force.row(i).norm() * VN.row(i).norm());
            theta(i) = acos(cos) / 3.1415926 * 180;
        }
        Eigen::VectorXi bl;
        igl::boundary_loop(F, bl);
        for(int i=0;i<bl.size();i++)
        {
            theta(bl(i)) = 90.0;
        }
        for(int i=0;i<force.rows();i++)
        {
            std::cout<<force(i)<<" "<<theta(i)<<std::endl;
        }
        std::cout<<theta.minCoeff()<<" "<<theta.maxCoeff()<<" "<<theta.sum()/ theta.size()<<std::endl;
    }
    else
    {

        if(isShowFaceColor != 0)
        {
            viewer.data().set_mesh(setup->initialPos, setup->mesh.faces());

            int nfaces = setup->mesh.faces().rows();
            igl::ColorMapType vizColor = igl::COLOR_MAP_TYPE_PARULA;
            Eigen::VectorXd Z(nfaces);
            Eigen::MatrixXd faceColors(nfaces, 3);

            double max = -1e10;
            double min = 1e10;

            for(int i=0;i<nfaces;i++)
            {
                if(faceColorType == 2)
                {
                    Eigen::Matrix2d abar, bbar;
                    abar = setup->abars[i];
                    bbar = setup->bbars[i];
                    Z(i) = (abar.inverse() * bbar).trace()/2;
                }
                else if(faceColorType == 3 || faceColorType == 4)
                {
                    Eigen::Matrix2d abar = setup->abars[i];
                    Eigen::Matrix2d T;
                    T.col(0) = ( setup->initialPos.row(setup->mesh.faceVertex(i, 1)) - setup->initialPos.row(setup->mesh.faceVertex(i, 0)) ).segment(0, 2);
                    T.col(1) = ( setup->initialPos.row(setup->mesh.faceVertex(i, 2)) - setup->initialPos.row(setup->mesh.faceVertex(i, 0)) ).segment(0, 2);
                    abar = (T.transpose()).inverse() * abar * T.inverse();
                    double s = (setup->abars[i].inverse() * setup->bbars[i]).trace() / 2.0;
//                    std::cout<<std::endl<<abar / abar.trace() * 2<<" "<<s<<std::endl;
//                    Eigen::Matrix2d bbar = s * abar;

                    if(faceColorType == 3)
                    {
                        s = abar.trace() / 2 - setup->thickness * s;
                    }
                    else
                    {
                        s = abar.trace() / 2 + setup->thickness * s;
                    }
                    Z(i) = s;
                }
                else
                {
                    Eigen::Matrix2d a, abar, A;
                    if(faceColorType == 0)      // from target to the optimal in first stage
                    {
                        a = setup->abars[i];
                        abar = firstFundamentalForm(setup->mesh, setup->targetPos, i, NULL, NULL);

                        A = abar.inverse()*a;
                        Eigen::EigenSolver<Eigen::MatrixXd> es(A);
                        double eigValue1 = es.eigenvalues()[0].real();

                        double eigValue2 = es.eigenvalues()[1].real();
                        Z(i) = std::max(eigValue1 , eigValue2);
                        if (Z(i) > max)
                            max = Z(i);
                        if(Z(i) < min)
                            min = Z(i);
                    }
                    else if(faceColorType == 1) // from plate to the optimal in first stage
                    {
                        Eigen::Matrix2d abar = setup->abars[i];
                        Eigen::Matrix2d T;
                        T.col(0) = ( setup->initialPos.row(setup->mesh.faceVertex(i, 1)) - setup->initialPos.row(setup->mesh.faceVertex(i, 0)) ).segment(0, 2);
                        T.col(1) = ( setup->initialPos.row(setup->mesh.faceVertex(i, 2)) - setup->initialPos.row(setup->mesh.faceVertex(i, 0)) ).segment(0, 2);
                        abar = (T.transpose()).inverse() * abar * T.inverse();

                        double s = abar.trace() / 2.0;

                        Z(i) = s;
                        if (Z(i) > max)
                            max = Z(i);
                        if (Z(i) < min)
                            min = Z(i);
                    }

                }
            }
            min = Z.minCoeff();
            max = Z.maxCoeff();
//            Z = Z / Z.lpNorm<Eigen::Infinity>();
            std::cout<<max<<" "<<min<<" "<<max - min <<std::endl;
            std::cout<<Z<<std::endl;
            igl::colormap(vizColor, Z, true, faceColors); // true here means libigl will automatically normalize Z, which may or may not be what you want.
            viewer.data().set_colors(faceColors);

        }

        if(isShowVerField)
        {
            Eigen::MatrixXd BC, Vec1, Vec2;
            viewer.data().set_mesh(setup->initialPos, setup->mesh.faces());
            igl::barycenter(viewer.data().V, setup->mesh.faces(), BC);

            int nfaces = setup->mesh.faces().rows();
            Vec1.resize(nfaces, 3);
            Vec2.resize(nfaces, 3);
            for(int i=0;i<nfaces;i++)
            {
                Eigen::Matrix2d a, abar, A;

                if(faceColorType == 0)  // from target to the optimal in first stage
                {
                    a = setup->abars[i];
                    abar = firstFundamentalForm(setup->mesh, setup->targetPos, i, NULL, NULL);
                }
                else if(faceColorType == 1) // from plate to the optimal in first stage
                {
                    a = setup->abars[i];
                    abar = firstFundamentalForm(setup->mesh, setup->initialPos, i, NULL, NULL);
                }
                A = abar.inverse()*a;
                //A = IU.inverse()*IM;
                Eigen::EigenSolver<Eigen::MatrixXd> es(A);
                double eigValue1 = es.eigenvalues()[0].real();
                double eigValue2 = es.eigenvalues()[1].real();

                int flag = 0;

                if(eigValue1 < eigValue2)
                {
                    flag = 1;
                }

                eigValue1 = es.eigenvalues()[flag].real();
                eigValue2 = es.eigenvalues()[1-flag].real();

                Eigen::VectorXd eigVec1 = es.eigenvectors().col(flag).real();

                Eigen::VectorXd eigVec2 = es.eigenvectors().col(1-flag).real();

                Vec1.row(i) = eigValue1*(eigVec1(0)*(setup->initialPos.row(setup->mesh.faces()(i,1))-setup->initialPos.row(setup->mesh.faces()(i,0))) + eigVec1(1)*(setup->initialPos.row(setup->mesh.faces()(i,2))-setup->initialPos.row(setup->mesh.faces()(i,0))));
                Vec2.row(i) = eigValue2*(eigVec2(0)*(setup->initialPos.row(setup->mesh.faces()(i,1))-setup->initialPos.row(setup->mesh.faces()(i,0))) + eigVec2(1)*(setup->initialPos.row(setup->mesh.faces()(i,2))-setup->initialPos.row(setup->mesh.faces()(i,0))));
                if(eigValue1 < eigValue2)
                    std::cout<<eigValue1<<" "<<eigValue2<<std::endl;
            }
            const Eigen::RowVector3d red(1,0,0), black(0,0,0);
            viewer.data().add_edges(BC,BC+Vec1/2, red);
            viewer.data().add_edges(BC,BC-Vec1/2, red);
            viewer.data().add_edges(BC,BC+Vec2/2, black);
            viewer.data().add_edges(BC,BC-Vec2/2, black);
        }
    }

}

void updateAbarPath(std::string modelPath, std::string modelType, std::string methodType, std::string dynamicSolveType, double thickness, double abarCoef, double bbarCoef, double smoothnessCoef, std::string &abarPath, bool isFixedCorners)
{
    std::string cornerStat = "Free";
    if(isFixedConer)
        cornerStat = "PinedCorners";
    abarPath = modelPath + "/" + methodType + "/" + dynamicSolveType + "/" + cornerStat + "/" + modelType + "_L_list_T_0_A_0_B_0_S_0.dat";

    int startIdx, endIdx;
    std::string subString = "";
    int expCoef;
    // thickness
    if(thickness == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(thickness));
    startIdx = abarPath.rfind("T");
    endIdx = abarPath.rfind("A");
    subString = "";
    if(thickness > 0)
        subString = "T_1e" + std::to_string(expCoef);
    else
        subString = "T_0";
    abarPath = abarPath.replace(abarPath.begin() + startIdx, abarPath.begin()+endIdx-1, subString);

    // Abar penalty
    if(abarCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(abarCoef));

    startIdx = abarPath.rfind("A");
    endIdx = abarPath.rfind("B");
    subString = "";
    if(abarCoef > 0)
        subString = "A_1e" + std::to_string(expCoef);
    else
        subString = "A_0";
    abarPath = abarPath.replace(abarPath.begin() + startIdx, abarPath.begin()+endIdx-1, subString);

    // bbar penalty
    if(bbarCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(bbarCoef));

    startIdx = abarPath.rfind("B");
    endIdx = abarPath.rfind("S");
    subString = "";
    if(bbarCoef > 0)
        subString = "B_1e" + std::to_string(expCoef);
    else
        subString = "B_0";
    abarPath = abarPath.replace(abarPath.begin() + startIdx, abarPath.begin()+endIdx-1, subString);

    // smoothness
    if(smoothnessCoef == 0)
        expCoef = 0;
    else
        expCoef = int(std::log10(smoothnessCoef));

    startIdx = abarPath.rfind("S");
    endIdx = abarPath.rfind(".");
    subString = "";
    if(smoothnessCoef > 0)
        subString = "S_1e" + std::to_string(expCoef);
    else
        subString = "S_0";
    abarPath = abarPath.replace(abarPath.begin() + startIdx, abarPath.begin()+endIdx, subString);
    std::cout<<"Current abar loading path is: "<<abarPath<<std::endl;
}

void conformalParametrization(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd &UV)
{
    Eigen::VectorXi bnd,b(2,1);
    igl::boundary_loop(F,bnd);
    b(0) = bnd(0);
    b(1) = bnd(round(bnd.size()/2));
    Eigen::MatrixXd bc(2,2);
    bc<<0,1,1,0;

    // LSCM parametrization
    if(igl::lscm(V,F,b,bc,UV))
    {
        std::cout<<"Parametrization Succeeded!"<<std::endl;
    }
    else
    {
        std::cout<<"Parametrization Failed!"<<std::endl;
        std::cout<<"Try harmonic Parametrization"<<std::endl;
        Eigen::MatrixXd bnd_uv;
        Eigen::MatrixXd initial_guess;
        igl::map_vertices_to_circle(V,bnd,bnd_uv);

        if(igl::harmonic(V,F,bnd,bnd_uv,1,UV))
        {
             std::cout<<"Parametrization Succeeded!"<<std::endl;
        }
    }

}

void conformalParametrizationARAP(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd &UV)
{
    // Compute the initial solution for ARAP (harmonic parametrization)
    Eigen::VectorXi bnd;
    igl::boundary_loop(F,bnd);
    Eigen::MatrixXd bnd_uv;
    Eigen::MatrixXd initial_guess;
    igl::map_vertices_to_circle(V,bnd,bnd_uv);

    conformalParametrization(V, F, initial_guess);

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;
    Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

    // Initialize ARAP
    arap_data.max_iter = 100;
    // 2 means that we're going to *solve* in 2d
    arap_precomputation(V,F,2,b,arap_data);


    // Solve arap using the harmonic map as initial guess
    UV = initial_guess;

    arap_solve(bc,arap_data,UV);

}

int main(int argc, char *argv[])
{
//    meshResampling("/Users/chenzhen/UT/Research/Projects/shellbenchmarks/benchmarks/TestModels/coarse/Nefertiti/Nefertiti_geometry.obj", 5000);
//    return 1;
//    testTriangulation();
//    return 1;
    //    Eigen::MatrixXd testV;
    //    Eigen::MatrixXi testF;
    //
    //    igl::readOBJ("../../benchmarks/TestModels/coarse/trapeZoid/draped_rect_geometry.obj", testV, testF);
    //    testSphere(testV, testF);
    //    igl::readOBJ("../../benchmarks/TestModels/fine/trapeZoid/draped_rect_geometry.obj", testV, testF);
    //    testSphere(testV, testF);
    //    return 0;
    //
    //    numSteps = 30;

    //    computeSphere("../../benchmarks/TestModels/middleCoarse/sphere/draped_disk_geometry.obj");
    ////    computeSaddle("../../benchmarks/TestModels/middleCoarse/saddle/draped_rect_geometry.obj");
    ////    computeCylinder("../../benchmarks/TestModels/middleCoarse/cylinder/draped_rect_geometry.obj");
    //    computeEllipsoid("../../benchmarks/TestModels/coarse/ellipsoid/draped_rect_geometry.obj");
    //    computeEllipsoid("../../benchmarks/TestModels/fine/ellipsoid/draped_rect_geometry.obj");
    //    return 0;

//

    tolerance = 1e-10;

    MidedgeAverageFormulation sff;

    setup = std::make_unique<SimulationSetupDynamicSolver>();

    igl::opengl::glfw::Viewer viewer;

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    //     Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        // Workspace
        if (ImGui::CollapsingHeader("Workspace", ImGuiTreeNodeFlags_DefaultOpen))
        {
            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;
            if (ImGui::Button("Load##Workspace", ImVec2((w-p)/2.f, 0)))
            {
                viewer.load_scene();
            }
            ImGui::SameLine(0, p);
            if (ImGui::Button("Save##Workspace", ImVec2((w-p)/2.f, 0)))
            {
                viewer.save_scene();
            }
        }

        // Mesh
        if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
        {
            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;
            if (ImGui::Button("Load##Mesh", ImVec2((w-p)/2.f, 0)))
            {
                viewer.data().clear();
                curPath = igl::file_dialog_open();
                Eigen::MatrixXi F;

                if (curPath.length() == 0)
                {
                    std::cout<<"Loading mesh failed"<<std::endl;
                }
                else
                {
                    int idx = curPath.rfind("/");
                    curPath = curPath.substr(0,idx);
                    idx = curPath.rfind("/");
                    selectedType = curPath.substr(idx+1, curPath.size()-1);
                    tarShape = curPath + "/" + selectedType ;
                    resShape = curPath + "/" + "draped_rect";
                    igl::readOBJ(tarShape + "_geometry.obj", curState.curPos, F);
                    setup->mesh = MeshConnectivity(F); // Just for repainting
                    std::cout<<curPath<<std::endl;
                    _initialPos = curState.curPos;
                    _initialMesh = setup->mesh;
                    //                    updateAbarPath(curPath, selectedType, selectedMethod, thickness, penaltyCoef, smoothnessCoef, setup->abarPath);
                }
                repaint(viewer);
            }
            ImGui::SameLine(0, p);
            if (ImGui::Button("Save##Mesh", ImVec2((w-p)/2.f, 0)))
            {
                viewer.open_dialog_save_mesh();
            }
        }

        if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Center object", ImVec2(-1, 0)))
            {
                viewer.core.align_camera_center(viewer.data().V, viewer.data().F);
            }
            if (ImGui::Button("Snap canonical view", ImVec2(-1, 0)))
            {
                viewer.snap_to_canonical_quaternion();
            }

            // Select rotation type
            int rotation_type = static_cast<int>(viewer.core.rotation_type);
            static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
            static bool orthographic = true;
            if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
            {
                using RT = igl::opengl::ViewerCore::RotationType;
                auto new_type = static_cast<RT>(rotation_type);
                if (new_type != viewer.core.rotation_type)
                {
                    if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
                    {
                        trackball_angle = viewer.core.trackball_angle;
                        orthographic = viewer.core.orthographic;
                        viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
                        viewer.core.orthographic = true;
                    }
                    else if (viewer.core.rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
                    {
                        viewer.core.trackball_angle = trackball_angle;
                        viewer.core.orthographic = orthographic;
                    }
                    viewer.core.set_rotation_type(new_type);
                }
            }

            // Orthographic view
            ImGui::Checkbox("Orthographic view", &(viewer.core.orthographic));
            //            ImGui::PopItemWidth();
        }
        // Parametrization
        if (ImGui::CollapsingHeader("Mesh Editor", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Combo("Parametrization Method", (int *)(&paramatrizationType), "ARAP\0Comformal\0\0"))
            {

            }
            if (ImGui::Button("Apply Parametrization", ImVec2(-1, 0)))
            {
                if (paramatrizationType == 0)   // ARAP
                {
                    std::cout<<"ARAP"<<std::endl;
                    Eigen::MatrixXd V, UV, V_P;
                    Eigen::MatrixXi F;
                    const std::string path = igl::file_dialog_open();
                    std::cout<<path<<std::endl;
                    igl::readOBJ(path, V, F);
                    V = V / V.lpNorm<Eigen::Infinity>();
                    igl::writeOBJ(path, V, F);
                    std::cout<<"Rescaling fnished!"<<std::endl;
                    V_P = V;
                    conformalParametrizationARAP(V, F, UV);
                    double norm_uv = UV.norm();
                    double norm_V = V_P.norm();
                    std::cout<<norm_uv<<" "<<norm_V<<std::endl;
                    for(int i = 0; i<UV.rows();i++)
                    {
                        V_P(i,0) = sqrt(2) * UV(i,0) * norm_V / norm_uv;
                        V_P(i,1) = sqrt(2) * UV(i,1) * norm_V / norm_uv;
                        V_P(i,2) = 0;
                    }
                    int idx =  path.rfind("/");
                    std::string planePath = path.substr(0, idx) + "/draped_rect_geometry.obj";
                    igl::writeOBJ(planePath, V_P, F);
                }
                if (paramatrizationType == 1)   // Comformal
                {
                    std::cout<<"conformal"<<std::endl;
                    Eigen::MatrixXd V, UV, V_P;
                    Eigen::MatrixXi F;
                    const std::string path = igl::file_dialog_open();
                    std::cout<<path<<std::endl;
                    igl::readOBJ(path, V, F);
                    V = V / V.lpNorm<Eigen::Infinity>();
                    igl::writeOBJ(path, V, F);
                    std::cout<<"Rescaling fnished!"<<std::endl;
                    V_P = V;
                    conformalParametrization(V, F, UV);
                    double norm_uv = UV.norm();
                    double norm_V = V_P.norm();
                    std::cout<<norm_uv<<" "<<norm_V<<std::endl;
                    for(int i = 0; i<UV.rows();i++)
                    {
                        V_P(i,0) = sqrt(2) * UV(i,0) * norm_V / norm_uv;
                        V_P(i,1) = sqrt(2) * UV(i,1) * norm_V / norm_uv;
                        V_P(i,2) = 0;
                    }
                    int idx =  path.rfind("/");
                    std::string planePath = path.substr(0, idx) + "/draped_rect_geometry.obj";
                    igl::writeOBJ(planePath, V_P, F);
                }
            }
            if (ImGui::InputInt("Resample Resolution", &resampledResolution))
            {

            }
            if (ImGui::Button("Mesh Resample", ImVec2(-1, 0)))
            {
                std::string path = igl::file_dialog_open();
                meshResampling(path, resampledResolution);
            }
            if (ImGui::Button("Quad to Tri", ImVec2(-1, 0)))
            {
                std::string path = igl::file_dialog_open();
                Eigen::MatrixXd V;
                Eigen::MatrixXi F;
                igl::readOBJ(path, V, F);
                assert(F.cols() == 4);

                int nfaces =  F.rows();
                Eigen::MatrixXi newF(2*nfaces, 3);

                for(int i=0;i<nfaces;i++)
                {
                    newF.row(2*i) << F(i, 0), F(i,1), F(i,2);
                    newF.row(2*i+1)<< F(i, 2), F(i,3), F(i,0);
                }
                int idx =  path.rfind("/");
                std::string planePath = path.substr(0, idx) + "/Tri_geometry.obj";
                igl::writeOBJ(planePath, V, newF);
            }
        }
        // Draw options
        if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if(ImGui::Checkbox("Vector Field", &isShowVerField))
            {
                repaint(viewer);
            }
            if(ImGui::Checkbox("Face Colors", &isShowFaceColor))
            {
                repaint(viewer);
            }
            if(ImGui::Checkbox("Plot Force", &isPlotForce))
            {
                repaint(viewer);
            }
            if(ImGui::Checkbox("Fixed Corners", &isFixedConer))
            {

            }
            ImGui::Checkbox("Show texture", &(viewer.data().show_texture));
            if (ImGui::Checkbox("Invert normals", &(viewer.data().invert_normals)))
            {
                viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
            }
            ImGui::Checkbox("Show overlay", &(viewer.data().show_overlay));
            ImGui::Checkbox("Show overlay depth", &(viewer.data().show_overlay_depth));
            ImGui::ColorEdit4("Background", viewer.core.background_color.data(),
                              ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
            ImGui::ColorEdit4("Line color", viewer.data().line_color.data(),
                              ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
            ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
            ImGui::DragFloat("Shininess", &(viewer.data().shininess), 0.05f, 0.0f, 100.0f);
            //ImGui::PopItemWidth();
        }

        // Overlays
        if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
            ImGui::Checkbox("Fill", &(viewer.data().show_faces));
//            ImGui::Checkbox("Show vertex labels", &(viewer.data().show_vertid));
//            ImGui::Checkbox("Show faces labels", &(viewer.data().show_faceid));
        }

        // Face color type
        if (ImGui::Combo("Face Color", (int *)(&faceColorType), "optimal2target\0optimal2plate\0bbar\0topLayer\0bottomLayer0\0"))
        {
            repaint(viewer);
        }
    };

    menu.callback_draw_custom_window = [&]()
    {
        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(0.0, 0.0), ImGuiSetCond_FirstUseEver);
        ImGui::Begin(
                     "Optimization", nullptr,
                     ImGuiWindowFlags_NoSavedSettings
                     );


        if (ImGui::Combo("Methods", (int *)(&methodType), "Normal\0alglibSolver\0ifOptSolver\0dynamicSolver\0\0"))
        {
            if (methodType == 0)
            {
                setup = std::make_unique<SimulationSetupNormal>();
            }
            else if (methodType == 1)
            {
                setup = std::make_unique<SimulationSetupAlglibSolver>();
                setup->abarCoef = abarCoef;
                setup->bbarCoef = bbarCoef;
                setup->smoothCoef = smoothnessCoef;
                setup->thickness = thickness;
                selectedMethod = "alglibSolver";
            }
            else if (methodType == 2)
            {
                setup = std::make_unique<SimulationSetupIpoptSolver>();
                setup->abarCoef = abarCoef;
                setup->bbarCoef = bbarCoef;
                setup->smoothCoef = smoothnessCoef;
                setup->thickness = thickness;
                selectedMethod = "ifOptSolver";
            }
            else if (methodType == 3)
            {
                setup = std::make_unique<SimulationSetupDynamicSolver>();
                setup->abarCoef = abarCoef;
                setup->bbarCoef = bbarCoef;
                setup->smoothCoef = smoothnessCoef; // for delta q
                setup->thickness = thickness;
                selectedMethod = "dynamicSolver";
            }
        }

        if (ImGui::Combo("Types", (int *)(&dynamicSolverType), "AbarPos\0AbarBbar\0ABbarPos\0\0"))
        {
            if(dynamicSolverType == 0)
            {
                selectedDynamicType = "AbarPos";
            }
            else if(dynamicSolverType == 1)
            {
                selectedDynamicType = "AbarBbar";
            }
            else if(dynamicSolverType == 2)
            {
                selectedDynamicType = "ABbarPos";
            }
        }

        if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen))
        {
            thickness = setup->thickness;
            setup->abarCoef = abarCoef;
            setup->bbarCoef = bbarCoef;
            smoothnessCoef = setup->smoothCoef;

            if (ImGui::InputDouble("Thickness", &thickness))
            {
                setup->thickness = thickness;
                if(thickness == 0)
                    expCoef = 0;
                else
                    expCoef = int(std::log10(thickness));
                if(setup->abarPath != "")
                {
                    int startIdx = setup->abarPath.rfind("T");
                    int endIdx = setup->abarPath.rfind("P");
                    std::string subString = "";
                    if(thickness > 0)
                        subString = "T_1e" + std::to_string(expCoef);
                    else
                        subString = "T_0";
                    setup->abarPath = setup->abarPath.replace(setup->abarPath.begin() + startIdx, setup->abarPath.begin()+endIdx-1, subString);
                    std::cout<<setup->abarPath<<std::endl;
                }

            }
            if (ImGui::InputDouble("Abar Penalty", &abarCoef))
            {
                setup->abarCoef = abarCoef;
                if(abarCoef == 0)
                    expCoef = 0;
                else
                    expCoef = int(std::log10(abarCoef));
                if(setup->abarPath != "")
                {
                    int startIdx = setup->abarPath.rfind("A");
                    int endIdx = setup->abarPath.rfind("B");
                    std::string subString = "";
                    if(abarCoef > 0)
                        subString = "A_1e" + std::to_string(expCoef);
                    else
                        subString = "A_0";
                    setup->abarPath = setup->abarPath.replace(setup->abarPath.begin() + startIdx, setup->abarPath.begin()+endIdx-1, subString);
                    std::cout<<setup->abarPath<<std::endl;
                }

            }

            if (ImGui::InputDouble("Bbar Penalty", &bbarCoef))
            {
                setup->bbarCoef = bbarCoef;
                if(bbarCoef == 0)
                    expCoef = 0;
                else
                    expCoef = int(std::log10(bbarCoef));
                if(setup->abarPath != "")
                {
                    int startIdx = setup->abarPath.rfind("B");
                    int endIdx = setup->abarPath.rfind("S");
                    std::string subString = "";
                    if(bbarCoef > 0)
                        subString = "B_1e" + std::to_string(expCoef);
                    else
                        subString = "B_0";
                    setup->abarPath = setup->abarPath.replace(setup->abarPath.begin() + startIdx, setup->abarPath.begin()+endIdx-1, subString);
                    std::cout<<setup->abarPath<<std::endl;
                }

            }
            if(ImGui::Checkbox("Overwrite Date", &isOverWrite))
            {

            }
            if(ImGui::Checkbox("Continue Unfinished OP", &isContinue))
            {

            }
            if (ImGui::InputDouble("Pos Smoothness Penalty", &smoothnessCoef))
            {
                setup->smoothCoef = smoothnessCoef;
                if(smoothnessCoef == 0)
                    expCoef = 0;
                else
                    expCoef = int(std::log10(smoothnessCoef));
                if(setup->abarPath != "")
                {
                    int startIdx = setup->abarPath.rfind("S");
                    int endIdx = setup->abarPath.rfind(".");
                    std::string subString = "";
                    if(smoothnessCoef > 0)
                        subString = "S_1e" + std::to_string(expCoef);
                    else
                        subString = "S_0";
                    setup->abarPath = setup->abarPath.replace(setup->abarPath.begin() + startIdx, setup->abarPath.begin()+endIdx, subString);
                    std::cout<<setup->abarPath<<std::endl;
                }

            }
        }

        if (ImGui::Button("Reset", ImVec2(-1, 0)))
        {
            reset();
            repaint(viewer);
        }
        if(ImGui::Button("Set Target", ImVec2(-1,0)))
        {
            setTarget();
            repaint(viewer);
        }

        if(ImGui::Button("Set Middle", ImVec2(-1,0)))
        {
            setMiddle();
            repaint(viewer);
        }

        if(ImGui::Button("Set target metric", ImVec2(-1, 0)))
        {
            updateAbarPath(curPath, selectedType, selectedMethod, selectedDynamicType, thickness, abarCoef, bbarCoef, smoothnessCoef, setup->abarPath, isFixedConer);
            bool ok = parseWimFiles(resShape, tarShape, *setup, sff);
            if (!ok)
            {
                std::cerr << "Couldn't load problem: " << std::endl;
                std::cerr << "Rest Shape: "<< resShape << std::endl;
                std::cerr << "Target Shape: "<< tarShape << std::endl;
                return -1;
            }
            
            if(isFixedConer)
            {
                for(int vid = 0; vid < 4; vid++)
                {
                    setup->clampedDOFs[3*vid + 0] = setup->targetPos(vid,0);
                    setup->clampedDOFs[3*vid + 1] = setup->targetPos(vid,0);
                    setup->clampedDOFs[3*vid + 2] = setup->targetPos(vid,0);
                }
            }
            else
            {
                setup->clampedDOFs.clear();
            }
            
            std::cout<<thickness<<std::endl;
            setup->thickness = thickness;
            setup->abarCoef = abarCoef;
            setup->bbarCoef = bbarCoef;
            setup->smoothCoef = smoothnessCoef;
            setup->selectedDynamicType = selectedDynamicType;
            setup->_is_overwrite = isOverWrite;
            setup->_is_continue = isContinue;
            int nfaces =  setup->mesh.nFaces();
            setup->abars.resize(nfaces);
            setup->bbars.resize(nfaces);
            for(int i=0;i<nfaces;i++)
            {
                Eigen::Matrix2d bbar;
                setup->abars[i] = firstFundamentalForm(setup->mesh, setup->targetPos, i, NULL, NULL);
                MidedgeAverageFormulation sff;
                Eigen::VectorXd vec(0);
                bbar = sff.secondFundamentalForm(setup->mesh, setup->targetPos, vec, i, NULL, NULL);
                setup->bbars[i] = ( setup->abars[i].inverse() * bbar ).trace() / 2.0 * setup->abars[i];
            }
            Eigen::MatrixXd V_temp;
            Eigen::MatrixXi F_temp;
            std::string path = setup->abarPath;
        }

        if (ImGui::Button("load and Compute", ImVec2(-1, 0)))
        {
            updateAbarPath(curPath, selectedType, selectedMethod, selectedDynamicType, thickness, abarCoef, bbarCoef, smoothnessCoef, setup->abarPath, isFixedConer);
            bool ok = parseWimFiles(resShape, tarShape, *setup, sff);
            if (!ok)
            {
                std::cerr << "Couldn't load problem: " << std::endl;
                std::cerr << "Rest Shape: "<< resShape << std::endl;
                std::cerr << "Target Shape: "<< tarShape << std::endl;
                return -1;
            }

            if(isFixedConer)
            {
                for(int vid = 0; vid < 4; vid++)
                {
                    setup->clampedDOFs[3*vid + 0] = setup->targetPos(vid,0);
                    setup->clampedDOFs[3*vid + 1] = setup->targetPos(vid,0);
                    setup->clampedDOFs[3*vid + 2] = setup->targetPos(vid,0);
                }
            }
            else
            {
                setup->clampedDOFs.clear();
            }

            std::cout<<thickness<<std::endl;
            setup->thickness = thickness;
            setup->abarCoef = abarCoef;
            setup->bbarCoef = bbarCoef;
            setup->smoothCoef = smoothnessCoef;
            setup->selectedDynamicType = selectedDynamicType;
            setup->_is_overwrite = isOverWrite;
            setup->_is_continue = isContinue;
            setup->buildRestFundamentalForms(sff);


            /*
             if(selectedType == "sphere")
             {
             Eigen::MatrixXd V;
             Eigen::MatrixXi F;
             igl::readOBJ("../../benchmarks/TestModels/coarse/trapeZoid/draped_rect_geometry.obj", V, F);

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



             //                std::cout<<M<<std::endl;

             abars.resize(nfaces);
             bbars.resize(nfaces);
             for(int i = 0; i < nfaces; i++)
             {
             double R = 0.5;
             //        double z = R - R*R/sqrt(R*R+u*u+v*v);
             //        double x = (R-z)/R*u;
             //        double y = (R-z)/R*v;

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
             setup->abars = abars;
             }
             */

        }
        if (ImGui::Button("Loading Remeshed plane", ImVec2(-1, 0)))
        {
            const std::string path = igl::file_dialog_open();
            Eigen::MatrixXd remeshedPos;
            Eigen::MatrixXi remeshedFaces;
            igl::readOBJ(path, remeshedPos, remeshedFaces);
//            igl::readOBJ(path, remeshedPos, remeshedPos);
            setup->remeshProcessing(remeshedPos, remeshedFaces);
            std::cout<<"Remesh finished!"<<std::endl;
        }


        if (ImGui::InputInt("Interpolation Steps", &numSteps))
        {

        }
        if (ImGui::InputDouble("Convergence Tolerance", &tolerance))
        {

        }
        if (ImGui::Button("Optimize Some Step", ImVec2(-1,0)))
        {
            //sff.testSecondFundamentalForm(setup->initialPos, setup->mesh.faces());
            double reg = 1e-10;
            int funcEvals = 0;
            double updateMag = std::numeric_limits<double>::infinity();
            double forceResidual = std::numeric_limits<double>::infinity();
//            if(numSteps > 1)
//            {
//                srand((unsigned)time(NULL));
//                curState.curPos = setup->initialPos;
//                for(int i=0;i<curState.curPos.rows();i++)
//                {
//                    curState.curPos(i,2) = (1e-4*rand())/RAND_MAX;
//                }
//            }
            setup->thickness = thickness;
            setup->abarCoef = abarCoef;
            setup->bbarCoef = bbarCoef;
            setup->smoothCoef = smoothnessCoef;
            for (int j = 1; j <= numSteps; j++)
            {
                double interp = double(j) / double(numSteps);
                std::cout<<interp<<std::endl;
                for (int i = 0; i < 10000; i++)
                {
                    takeOneStep(*setup, curState, sff, reg, interp, funcEvals, forceResidual, updateMag);
                    repaint(viewer);

                    if (forceResidual < tolerance || updateMag < tolerance)
                        break;
                }
            }
            std::cout << "Finished with " << funcEvals << " evaluations" << std::endl;
            Eigen::Matrix3d R;
            Eigen::Vector3d t;
            GeometryTools::rigidMotionTransformation(curState.curPos, setup->targetPos, setup->mesh, R, t);
            Eigen::MatrixXd ones(1, curState.curPos.rows());
            ones.setOnes();

            curState.curPos.transpose() = R * curState.curPos.transpose() + t * ones;
            /*
             if(selectedType == "trapeZoid")
             {
             std::cout<<"Saving Abar path: "<<setup->abarPath<<std::endl;
             std::ofstream outfile(setup->abarPath, std::ios::trunc);

             outfile<<thickness<<"\n";
             outfile<<penaltyCoef<<"\n";
             outfile<<smoothnessCoef<<"\n";

             int nverts = curState.curPos.rows();
             int nfaces = setup->mesh.nFaces();
             outfile<<3*nverts<<"\n";
             outfile<<3*nfaces<<"\n";

             for(int i=0;i<nverts;i++)
             {
             outfile<<std::setprecision(16)<<curState.curPos(i, 0)<<"\n";
             outfile<<std::setprecision(16)<<curState.curPos(i, 1)<<"\n";
             outfile<<std::setprecision(16)<<curState.curPos(i, 2)<<"\n";
             }

             for(int i=0;i<nfaces;i++)
             {
             double x = sqrt(setup->abars[i](0,0));
             double y = setup->abars[i](0,1)/x;
             double z = sqrt(setup->abars[i].determinant())/x;
             outfile<<std::setprecision(16)<<x<<"\n";
             outfile<<std::setprecision(16)<<y<<"\n";
             if(i != nfaces - 1)
             outfile<<std::setprecision(16)<<z<<"\n";
             else
             outfile<<std::setprecision(16)<<z;
             }
             outfile.close();

             int startIdx, endIdx, expCoef;
             std::string subString = "";
             std::string resampledPath = setup->abarPath;

             startIdx = resampledPath.rfind("/");
             endIdx = resampledPath.find("_");
             resampledPath = resampledPath.replace(resampledPath.begin() + startIdx + 1,resampledPath.begin() + endIdx, "resampled");

             // thickness
             if(thickness == 0)
             expCoef = 0;
             else
             expCoef = int(std::log10(thickness));
             startIdx = resampledPath.rfind("T");
             endIdx = resampledPath.rfind("P");
             subString = "";
             if(thickness > 0)
             subString = "T_1e" + std::to_string(expCoef);
             else
             subString = "T_0";
             resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);

             // penalty
             if(penaltyCoef == 0)
             expCoef = 0;
             else
             expCoef = int(std::log10(penaltyCoef));

             startIdx = resampledPath.rfind("P");
             endIdx = resampledPath.rfind("S");
             subString = "";
             if(penaltyCoef > 0)
             subString = "P_1e" + std::to_string(expCoef);
             else
             subString = "P_0";
             resampledPath= resampledPath .replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx - 1, subString);

             // smoothness
             if(smoothnessCoef == 0)
             expCoef = 0;
             else
             expCoef = int(std::log10(smoothnessCoef));

             startIdx = resampledPath.rfind("S");
             endIdx = resampledPath.rfind(".");
             subString = "";
             if(smoothnessCoef > 0)
             subString = "S_1e" + std::to_string(expCoef);
             else
             subString = "S_0";
             resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.begin() + endIdx, subString);

             startIdx = resampledPath.rfind(".");
             resampledPath = resampledPath.replace(resampledPath.begin() + startIdx,resampledPath.end(), ".obj");
             std::cout<<"Current abar loading path is: "<<resampledPath<<std::endl;
             igl::writeOBJ(resampledPath, curState.curPos, setup->mesh.faces());
             }
             */
        }
        if (ImGui::Button("Run Postprocessing", ImVec2(-1, 0)))
        {
            postProcess();
        }
        if (ImGui::Button("Leading Eigenvector", ImVec2(-1, 0)))
        {
            leadingEigenvector(*setup, curState, sff, evec);
            repaint(viewer);
        }
         ImGui::End();
    };



    viewer.data().set_face_based(false);
    repaint(viewer);
    viewer.launch();
}


