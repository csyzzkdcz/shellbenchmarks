#ifndef GEOMETRYTOOLS_H
#define GEOMETRYTOOLS_H

#include "MeshConnectivity.h"
#include <Eigen/Dense>

namespace GeometryTools
{
    void rigidMotionTransformation(Eigen::MatrixXd pos, Eigen::MatrixXd tarPos, MeshConnectivity mesh, Eigen::Matrix3d &R, Eigen::Vector3d &t);
    
    void testRigigMotionTransformation();
}


#endif
