import pymesh
import numpy as np
from math import *

N = 30;
vert = [];
for i in range(N):
    vert.append([0.5*cos(i/N *2*pi), 0.5*sin(i/N * 2*pi), 0.0]);
vertices = np.array(vert);
# vertices = np.array([[-0.5,-0.5, 0.0], [0.5, -0.5, 0.0], [0.5,0.5, 0.0], [-0.5,0.5,0.0]]);
tri = pymesh.triangle();
tri = pymesh.triangle();
tri.points = vertices;
tri.split_boundary = True;
tri.conforming_delaunay = True;
tri.max_area = 0.003;
tri.run()
mesh = tri.mesh
pymesh.save_mesh("draped_disk_geometry.obj", mesh, use_double =True)
# print(mesh.vertices)
# vert = []
# for vertex in mesh.vertices:
#     vertex = vertex * 2;
#     vert.append(vertex);
# vertices = np.array(vert);
# meshExpanded = pymesh.form_mesh(vertices, mesh.faces)
# pymesh.save_mesh("expanded_circle_geometry.obj", meshExpanded, use_double =True)