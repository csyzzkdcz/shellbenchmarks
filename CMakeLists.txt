cmake_minimum_required(VERSION 3.1)
project(shellbenchmark)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(LAPACK REQUIRED)


# Use C++14
set (CMAKE_CXX_STANDARD 14)

set(CMAKE_PREFIX_PATH "$ENV{HOME}/.local")

# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" ON)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

find_package(LIBIGL REQUIRED QUIET)
find_package(ifopt REQUIRED)
find_package(eigen3 REQUIRED)
find_package(MOSEK REQUIRED)

# Add your project files
file(GLOB SRCFILES *.cpp  *.h SecondFundamentalForm/*.cpp SimulationSetup/*.cpp SecondFundamentalForm/*.h SimulationSetup/*.h cppoptlib/*.cpp cppoptlib/*.h SensitiveAnalysis/*.h SensitiveAnalysis/*.cpp)

# find_package(ifopt REQUIRED)
message("Found IFOPT at ${ifopt_INCLUDE_DIRS} from $ENV(ifopt}")
include_directories( ${ifopt_INCLUDE_DIRS} )

if (LAPACK_FOUND)
  include_directories(${LAPACK_INCLUDE_DIRS})
endif()

find_package(SuiteSparse REQUIRED)
if (SUITESPARSE_FOUND)
  include_directories(${SUITESPARSE_INCLUDE_DIRS})
  message("found SuiteSparse at ${SUITESPARSE_INCLUDE_DIRS}")
endif()

find_package(Umfpack REQUIRED)
if (UMFPACK_FOUND)
  include_directories(${UMFPACK_INCLUDES})
  message("found umfpack at ${UMFPACK_INCLUDES}")
endif()


message("Found MOSEK at ${MOSEK_INCLUDE_DIR}")
include_directories(${MOSEK_INCLUDE_DIR})

message("Found Eigen at ${EIGEN3_INCLUDE_DIR}")
message("Eigen version is ${EIGEN3_VERSION}")

message(" MOSEK Libs are ${MOSEK_LIBRARIES} ${MOSEK_FUSION_LIBRARIES}")


# External libs
add_subdirectory(alglib)

add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::triangle)
target_link_libraries(${PROJECT_NAME}_bin alglib ifopt::ifopt_ipopt)
target_link_libraries(${PROJECT_NAME}_bin ${MOSEK_LIBRARIES} ${MOSEK_FUSION_LIBRARIES})
target_link_libraries(${PROJECT_NAME}_bin ${SUITESPARSE_LIBRARIES} ${UMFPACK_LIBRARIES})