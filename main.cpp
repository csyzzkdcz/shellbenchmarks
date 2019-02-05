#include <igl/opengl/glfw/Viewer.h>

#include <random>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/decimate.h>
#include <igl/upsample.h>
#include <imgui/imgui.h>
#include <memory>
//#include <igl/triangle/triangulate.h>
//#include "SecondFundamentalForm/MidedgeAngleSinFormulation.h"
//#include "SecondFundamentalForm/MidedgeAngleTanFormulation.h"
#include "SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "MeshConnectivity.h"
#include "GeometryDerivatives.h"
#include "BuildModels.h"
#include "ElasticShell.h"
#include "SimulationSetup/SimulationSetup.h"
#include "SimulationSetup/SimulationSetupNormal.h"
#include "SimulationSetup/SimulationSetupFindAbar.h"
#include "SimulationSetup/SimulationSetupIpoptSolver.h"
#include "ParseWimFiles.h"
#include "StaticSolve.h"
#include "SimulationState.h"



std::unique_ptr<SimulationSetup> setup;
SimulationState curState;
Eigen::VectorXd evec;
int numSteps;
double tolerance;
bool isShowVerField = false;
bool isShowFaceColor = false;
bool isShowAbar =  false;

enum Methods {Normal=0, alglibSolver, ipoptSolver};
static Methods methodType =  ipoptSolver;

std::string resShape = "";
std::string tarShape = "";
std::string curPath = "";


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
    numSteps = 1;
    tolerance = 1e-6;
    
    std::string selectedType = "sphere";
    
    MidedgeAverageFormulation sff;
    
    setup = std::make_unique<SimulationSetupIpoptSolver>();
    
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
                    std::cout<<curPath<<std::endl;
                    setup->mesh = MeshConnectivity(F);
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
        
        
        if (ImGui::Combo("Methods", (int *)(&methodType), "Normal\0alglibSolver\0ifOptSolver\0\0"))
        {
            if (methodType == 0)
            {
                setup = std::make_unique<SimulationSetupNormal>();
            }
            else if (methodType == 1)
            {
                setup = std::make_unique<SimulationSetupFindAbar>();
                setup->abarPath = curPath + "/alglibSolver" + "/" + selectedType + "_L_list_" +std::to_string(int(-log(setup->penaltyCoef))) + ".dat";
                std::cout<<setup->abarPath<<std::endl;
            }
            else if (methodType == 2)
            {
                setup = std::make_unique<SimulationSetupIpoptSolver>();
                setup->abarPath = curPath + "/ifOptSolver" + "/" + selectedType + "_L_list_0" +std::to_string(int(-log(setup->penaltyCoef)))+ ".dat";
                std::cout<<setup->abarPath<<std::endl;
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
        if (ImGui::Button("load and Compute", ImVec2(-1, 0)))
        {
            bool ok = parseWimFiles(resShape, tarShape, *setup, sff);
            if (!ok)
            {
                std::cerr << "Couldn't load problem: " << std::endl;
                std::cerr << "Rest Shape: "<< resShape << std::endl;
                std::cerr << "Target Shape: "<< tarShape << std::endl;
                return -1;
            }
            setup->buildRestFundamentalForms(sff);
            setup->penaltyCoef = 0;
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
