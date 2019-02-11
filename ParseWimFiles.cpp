#include "ParseWimFiles.h"
#include <Eigen/Core>
#include <igl/readOBJ.h>
#include "SimulationSetup/SimulationSetup.h"
#include "SimulationSetup/SimulationSetupNormal.h"
#include <fstream>
#include "SecondFundamentalForm/SecondFundamentalFormDiscretization.h"

bool parseWimFiles(const std::string &prefixRes, const std::string &prefixTar, SimulationSetup &parsedSetup, const SecondFundamentalFormDiscretization &sff)
{
    std::string resMeshName = prefixRes + std::string("_geometry.obj");
    Eigen::MatrixXi F, F1;
    Eigen::MatrixXd V, V1;
    if (!igl::readOBJ(resMeshName, V, F))
        return false;
    parsedSetup.initialPos = V;
    parsedSetup.mesh = MeshConnectivity(F);

    std::string tarMeshName = prefixTar + std::string("_geometry.obj");
    if (!igl::readOBJ(tarMeshName, V1, F1))
        return false;
    parsedSetup.targetPos = V1;

    if (V1.rows()!=V.rows() || F1.rows()!= F.rows())
        return false;
    
    
    int nedgedofs = sff.numExtraDOFs();
    parsedSetup.initialEdgeDOFs.resize(nedgedofs * parsedSetup.mesh.nEdges());
    sff.initializeExtraDOFs(parsedSetup.initialEdgeDOFs, parsedSetup.mesh, parsedSetup.initialPos);    
    
    /*
    parsedSetup.clampedDOFs.clear();
    std::string clampedName = prefixRes + std::string("_clamped_vertices.dat");
    std::ifstream ifs(clampedName);
    if(!ifs)
        return false;
    int nclamped;
    ifs >> nclamped;
    char dummy;
    ifs >> dummy;
    if (!ifs)
        return false;
    ifs.ignore(std::numeric_limits<int>::max(), '\n');
    for (int i = 0; i < nclamped; i++)
    {
        std::string line;
        std::getline(ifs, line);
        std::stringstream ss(line);

        int vid;
        ss >> vid;
        if (!ss || vid < 0 || vid >= parsedSetup.initialPos.rows())
            return false;
        double x, y, z;
        ss >> x >> y >> z;
        if (!ss)
        {
            // no displacement specified, use original position
            x = V(vid, 0);
            y = V(vid, 1);
            z = V(vid, 2);
        }
        parsedSetup.clampedDOFs[3*vid + 0] = x;
        parsedSetup.clampedDOFs[3*vid + 1] = y;
        parsedSetup.clampedDOFs[3*vid + 2] = z;        
    }
    if (!ifs)
        return false;*/

    /*int nfaces = parsedSetup.mesh.nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        int nclamped = 0;
        for (int j = 0; j < 3; j++)
        {
            if (parsedSetup.clampedVertices.count(parsedSetup.mesh.faceVertex(i, j)))
                nclamped++;
        }
        if (nclamped == 2)
        {
            for (int j = 0; j < 3; j++)
            {
                if (!parsedSetup.clampedVertices.count(parsedSetup.mesh.faceVertex(i, j)))
                    std::cout << parsedSetup.mesh.faceVertex(i, j) << std::endl;
            }
        }
    }

    while (true);*/

    /*
    parsedSetup.externalForces.resize(parsedSetup.initialPos.rows(), 3);
    parsedSetup.externalForces.setZero();
    std::string loadName = prefix + std::string("_loaded_vertices.dat");
    std::ifstream lfs(loadName);
    if (!lfs)
        return false;
    int nloaded;
    lfs >> nloaded;
    lfs >> dummy;
    if (!lfs)
        return false;
    for (int i = 0; i < nloaded; i++)
    {
        int vidx;
        lfs >> vidx;
        if (!lfs || vidx < 0 || vidx >= parsedSetup.externalForces.rows())
            return false;
        for (int j = 0; j < 3; j++)
        {
            double val;
            lfs >> val;
            parsedSetup.externalForces(vidx, j) = val;
        }
    }
    if (!lfs)
        return false;
    */

    std::string matName = prefixRes + std::string("_material.dat");
    std::ifstream mfs(matName);
    if (!mfs)
        return false;
    double t;
//    mfs >> parsedSetup.thickness;
    mfs >> t;
    mfs >> parsedSetup.YoungsModulus;
    mfs >> parsedSetup.PoissonsRatio;
    if (!mfs)
        return false;
    /*
    parsedSetup.tests.clear();
    std::string testName = prefix + std::string("_postprocess_vertices.dat");
    std::ifstream tfs(testName);
    if (!tfs)
        return false;
    int numtests;
    tfs >> numtests;
    tfs >> dummy;
    if (!tfs)
        return false;

    tfs.ignore(std::numeric_limits<int>::max(), '\n');
    for (int i = 0; i < numtests; i++)
    {        
        std::string line;
        std::getline(tfs, line);
        std::stringstream ss(line);
        PostprocessTest test;

        ss >> test.vertex;
        ss >> test.coord;
        if (!ss)
            return false;
        ss >> test.wimDisplacement;
        if (!ss)
        {
            // no Wim data
            test.wimDisplacement = 0.0;
        }
        parsedSetup.tests.push_back(test);
    }
    */
    
    //parsedSetup.buildRestFundamentalForms(sff);

    return true;
}
