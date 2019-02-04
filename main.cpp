#include <igl/opengl/glfw/Viewer.h>
#include "MeshConnectivity.h"
#include <random>
#include "ElasticShell.h"
#include "SimulationSetup/SimulationSetup.h"
#include "SimulationSetup/SimulationSetupNormal.h"
#include "SimulationSetup/SimulationSetupFindAbar.h"
#include "SimulationSetup/SimulationSetupIpoptSolver.h"
#include "ParseWimFiles.h"
#include "StaticSolve.h"
#include "SimulationState.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/decimate.h>
#include <igl/upsample.h>
#include <imgui/imgui.h>
//#include <igl/triangle/triangulate.h>
//#include "SecondFundamentalForm/MidedgeAngleSinFormulation.h"
//#include "SecondFundamentalForm/MidedgeAngleTanFormulation.h"
#include "SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "GeometryDerivatives.h"
#include <memory>


std::unique_ptr<SimulationSetup> setup;
SimulationState curState;
Eigen::VectorXd evec;
int numSteps;
double tolerance;
bool isShowVerField = false;
bool isShowFaceColor = false;
bool isShowAbar =  false;
enum tarShapes { sphere=0, saddle, hypar, cylinder };
static tarShapes selected = sphere;
std::string resShape = "";
std::string tarShape = "";

void compute_sphere(std::string rectPath)
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
        Vo(i,0) = x;
        Vo(i,1) = y;
        Vo(i,2) = z;
    }
    igl::writeOBJ("../../benchmarks/TestModels/coarse/sphere/sphere_geometry.obj", Vo, Fo);
    
}

void compute_hypar(std::string rectPath)
{
    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;
    igl::readOBJ(rectPath, Vo, Fo);
    for(int i=0;i<Vo.rows();i++)
    {
        Vo(i,2) = 32*Vo(i,1) * Vo(i,0);
    }
    igl::writeOBJ("../../benchmarks/TestModels/coarse/hypar/hypar_geometry.obj", Vo, Fo);
}

void compute_cylinder(std::string rectPath)
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
    igl::writeOBJ("../../benchmarks/TestModels/coarse/cylinder/cylinder_geometry.obj", Vo, Fo);
}

void meshResample(int targetResolution)
{
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi E;
//    Eigen::MatrixXd H;
//
//    // Triangulated interior
//    Eigen::MatrixXd V2;
//    Eigen::MatrixXi F2;
//    V.resize(8,1);
//    E.resize(8,1);
//
//    V << -1,-1, 1,-1, 1,1, -1, 1,
//
//    E << 0,1, 1,2, 2,3, 3,0,
//
//
//    // Triangulate the interior
//    igl::triangle::triangulate(V,E,H,"a0.005q",V2,F2);
    
    // Plot the generated mesh
    Eigen::MatrixXd Vcurr;
    Eigen::MatrixXi Fcurr;
    igl::readOBJ("../../benchmarks/DrapedRect/3876_triangles/draped_rect_geometry.obj", Vcurr, Fcurr);
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
//    igl::writeOBJ("../../benchmarks/TestModels/resampled/draped_rect_geometry.obj",V2,F2);
    
    
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

void repaint(igl::opengl::glfw::Viewer &viewer)
{
    viewer.data().clear();
    viewer.data().set_mesh(curState.curPos, setup->mesh.faces());
    Eigen::MatrixXd colors(setup->initialPos.rows(), 3);
    //    igl::jet(evec, true, colors);
    //    for (auto &it : setup->clampedDOFs)
    //    {
    //        int vid = it.first / 3;
    //        colors(vid, 0) = 1.0;
    //        colors(vid, 1) = 0.0;
    //        colors(vid, 2) = 0.0;
    //    }
    //    viewer.data().set_colors(colors);
    
    colors.col(0).setConstant(1.0);
    colors.col(1).setConstant(1.0);
    colors.col(2).setConstant(0);
    //    for (int i=0;i<op->p_fixed_index.size();i++)
    //    {
    //        colors(op->p_fixed_index[i], 0) = 1.0;
    //        colors(op->p_fixed_index[i], 1) = 0.0;
    //        colors(op->p_fixed_index[i], 2) = 0.0;
    //    }
    viewer.data().set_colors(colors);
    viewer.data().line_width = 2;
    viewer.core.background_color<<1,1,1,1;
    
    if(isShowVerField)
    {
        Eigen::MatrixXd BC, Vec1, Vec2;
        
        if(isShowAbar)
        {
            viewer.data().set_mesh(setup->initialPos, setup->mesh.faces());
        }
        igl::barycenter(viewer.data().V, setup->mesh.faces(), BC);
        
        int nfaces = setup->mesh.faces().rows();
        Vec1.resize(nfaces, 3);
        Vec2.resize(nfaces, 3);
        for(int i=0;i<nfaces;i++)
        {
            Eigen::Matrix2d a, abar, A;
            
            if(isShowAbar)
            {
                a = setup->abars[i];
                abar = firstFundamentalForm(setup->mesh, setup->initialPos, i, NULL, NULL);
        
                
                A = abar.inverse()*a;
                //A = IU.inverse()*IM;
                Eigen::EigenSolver<Eigen::MatrixXd> es(A);
                double eigValue1 = es.eigenvalues()[0].real();
                Eigen::VectorXd eigVec1 = es.eigenvectors().col(0).real();
                
                double eigValue2 = es.eigenvalues()[1].real();
                Eigen::VectorXd eigVec2 = es.eigenvectors().col(1).real();
                
                Vec1.row(i) = eigValue1*(eigVec1(0)*(setup->initialPos.row(setup->mesh.faces()(i,1))-setup->initialPos.row(setup->mesh.faces()(i,0))) + eigVec1(1)*(setup->initialPos.row(setup->mesh.faces()(i,2))-setup->initialPos.row(setup->mesh.faces()(i,0))));
                Vec2.row(i) = eigValue2*(eigVec2(0)*(setup->initialPos.row(setup->mesh.faces()(i,1))-setup->initialPos.row(setup->mesh.faces()(i,0))) + eigVec2(1)*(setup->initialPos.row(setup->mesh.faces()(i,2))-setup->initialPos.row(setup->mesh.faces()(i,0))));
                
            }
            else
            {
                a = firstFundamentalForm(setup->mesh, curState.curPos, i, NULL, NULL);
                abar = firstFundamentalForm(setup->mesh, setup->initialPos, i, NULL, NULL);
                
                A = abar.inverse()*a;
                //A = IU.inverse()*IM;
                Eigen::EigenSolver<Eigen::MatrixXd> es(A);
                double eigValue1 = es.eigenvalues()[0].real();
                Eigen::VectorXd eigVec1 = es.eigenvectors().col(0).real();
                
                double eigValue2 = es.eigenvalues()[1].real();
                Eigen::VectorXd eigVec2 = es.eigenvectors().col(1).real();
                
                Vec1.row(i) = eigValue1*(eigVec1(0)*(curState.curPos.row(setup->mesh.faces()(i,1))-curState.curPos.row(setup->mesh.faces()(i,0)))+ eigVec1(1)*(curState.curPos.row(setup->mesh.faces()(i,2))-curState.curPos.row(setup->mesh.faces()(i,0))));
                Vec2.row(i) = eigValue2*(eigVec2(0)*(curState.curPos.row(setup->mesh.faces()(i,1))-curState.curPos.row(setup->mesh.faces()(i,0))) + eigVec2(1)*(curState.curPos.row(setup->mesh.faces()(i,2))-curState.curPos.row(setup->mesh.faces()(i,0))));
                
            }
            
           
            
        }
        const Eigen::RowVector3d red(1,0,0), black(0,0,0);
        viewer.data().add_edges(BC,BC+Vec1, red);
        viewer.data().add_edges(BC,BC+Vec2, black);
        
    }
    
    if(isShowFaceColor)
    {
        //        load_L_list("../../benchmarks/TestModels/L_list_cylinder.dat");
        int nfaces = setup->mesh.faces().rows();
        igl::ColorMapType vizColor = igl::COLOR_MAP_TYPE_PARULA;
        Eigen::VectorXd Z(nfaces);
        Eigen::MatrixXd faceColors(nfaces, 3);
        
        if(isShowAbar)
        {
            viewer.data().set_mesh(setup->initialPos, setup->mesh.faces());
        }
        
        for(int i=0;i<nfaces;i++)
        {
            Eigen::Matrix2d a, abar, A;
            
            if(isShowAbar)
            {
//                Eigen::MatrixXd V;
//                Eigen::MatrixXd F;
//                igl::readOBJ(tarShape + "_geometry.obj", V, F);
//                a = firstFundamentalForm(setup->mesh, V, i, NULL, NULL);
                a = setup->abars[i];
                abar = firstFundamentalForm(setup->mesh, setup->initialPos, i, NULL, NULL);
            }
            else
            {
                a = firstFundamentalForm(setup->mesh, curState.curPos, i, NULL, NULL);
                abar = firstFundamentalForm(setup->mesh, setup->initialPos, i, NULL, NULL);
            }
            A = abar.inverse()*a;
            //A = IU.inverse()*IM;
            Eigen::EigenSolver<Eigen::MatrixXd> es(A);
            double eigValue1 = es.eigenvalues()[0].real();
            
            double eigValue2 = es.eigenvalues()[1].real();
            Z(i) = eigValue1 + eigValue2;
            
        }
        
        igl::colormap(vizColor, Z, true, faceColors); // true here means libigl will automatically normalize Z, which may or may not be what you want.
        viewer.data().set_colors(faceColors);
    }
}

int main(int argc, char *argv[])
{
//    compute_cylinder("../../benchmarks/TestModels/coarse/cylinder/draped_rect_geometry.obj");
//    compute_sphere("../../benchmarks/TestModels/coarse/sphere/draped_rect_geometry.obj");
//    compute_hypar("../../benchmarks/TestModels/coarse/hypar/draped_rect_geometry.obj");
//    return;
    numSteps = 30;
    tolerance = 1e-6;
    
    std::string selectedType = "sphere";
    
    MidedgeAverageFormulation sff;
    //MidedgeAngleTanFormulation sff;
    //MidedgeAngleSinFormulation sff;
    //std::string problem = "benchmarks/ClampedCylindricalShell/2496_triangles/cantilever_cylindrical_shell";
    //std::string problem = "C:/Users/evouga/Documents/wim/tests/trivial/trivial";
    //std::string problem = "D:/Dropbox/Dropbox/ShellBenchmarks/CantileverPlate/17640_triangles/cantilever_plate";
    //std::string problem = "benchmarks/CrushedCylinder/4900_triangles/crushed_cylinder";
    
    //    setup = std::make_unique<SimulationSetupFindAbar>();
    
    //    Eigen::MatrixXd Vt;
    //    Eigen::MatrixXi Ft;
    //    igl::readOBJ("../../benchmarks/TestModels/sphere_geometry.obj", Vt, Ft);
    //    curState.curPos = Vt;
    //jitter(1e-6);
    resShape = "../../benchmarks/TestModels/coarse/" + selectedType + "/draped_rect";
    tarShape = "../../benchmarks/TestModels/coarse/" + selectedType + "/" + selectedType;
    
//    meshResample(512);
//    compute_hypar(resShape + "_geometry.obj");
//    compute_cylinder(resShape + "_geometry.obj");
//    compute_sphere(resShape + "_geometry.obj");
    
    setup = std::make_unique<SimulationSetupIpoptSolver>();
    

    setup->penaltyCoef = 0;
    
    bool ok = parseWimFiles(resShape, tarShape, *setup, sff);
    if (!ok)
    {
        std::cerr << "Couldn't load problem: " << std::endl;
        std::cerr << "Rest Shape: "<< resShape << std::endl;
        std::cerr << "Target Shape: "<< tarShape << std::endl;
        return -1;
    }
    
    reset();
    igl::opengl::glfw::Viewer viewer;
    
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    //     Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        // menu.draw_viewer_menu();
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
//                viewer.open_dialog_load_mesh();
                std::string fname = igl::file_dialog_open();
                Eigen::MatrixXi F;

                if (fname.length() == 0)
                {
                    std::cout<<"Loading mesh failed"<<std::endl;
                }
                else
                {
                    igl::readOBJ(fname, curState.curPos, F);
                }
//                curState.curPos = viewer.data().F;
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
        //
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
            if(ImGui::Checkbox("Abars", &isShowAbar))
            {
                repaint(viewer);
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
            ImGui::Checkbox("Show vertex labels", &(viewer.data().show_vertid));
            ImGui::Checkbox("Show faces labels", &(viewer.data().show_faceid));
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
        
        
        if (ImGui::Combo("Target Shapes", (int *)(&selected), "sphere\0saddle\0hypar\0cylinder\0\0"))
        {
            if (selected == 0)
                selectedType = "sphere";
            else if (selected == 1)
                selectedType = "saddle";
            else if (selected == 2)
                selectedType = "hypar";
            else if (selected == 3)
                selectedType = "cylinder";
            
            resShape = "../../benchmarks/TestModels/coarse/" + selectedType + "/draped_rect";
            tarShape = "../../benchmarks/TestModels/coarse/" + selectedType + "/" + selectedType;
            
            bool ok = parseWimFiles(resShape, tarShape, *setup, sff);
            if (!ok)
            {
                std::cerr << "Couldn't load problem: " << std::endl;
                std::cerr << "Rest Shape: "<< resShape << std::endl;
                std::cerr << "Target Shape: "<< tarShape << std::endl;
                return -1;
            }
            
           // if (selected == 3 || selected == 2)
            {
                Eigen::MatrixXd Vt;
                Eigen::MatrixXi Ft;
                igl::readOBJ(tarShape + "_geometry.obj", Vt, Ft);
                curState.curPos = Vt;
                numSteps = 1;
            }
        }
        
        if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen))
        {
            double thickness = setup->thickness;
            double penaltyCoef = setup->penaltyCoef;
            if (ImGui::InputDouble("Thickness", &thickness))
            {
                setup->thickness = thickness;
            }
            if (ImGui::InputDouble("Penalty Coefficient", &penaltyCoef))
            {
                setup->penaltyCoef = penaltyCoef;
            }
        }
        
        if (ImGui::Button("Reset", ImVec2(-1, 0)))
        {
            reset();
            repaint(viewer);
        }
        if (ImGui::Button("Recompute Abars", ImVec2(-1, 0)))
        {
            parseWimFiles(resShape, tarShape, *setup, sff);
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
            double reg = 1e-6;
            int funcEvals = 0;
            double updateMag = std::numeric_limits<double>::infinity();
            double forceResidual = std::numeric_limits<double>::infinity();
            if(numSteps > 1)
            {
                srand((unsigned)time(NULL));
                for(int i=0;i<curState.curPos.rows();i++)
                {
                    curState.curPos(i,2) = (1e-6*rand())/RAND_MAX;
                }
            }
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
