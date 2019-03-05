import pymesh
import numpy as np

vertices = np.array([[-0.5,-0.5, 0.0], [0.5, -0.5, 0.0], [0.5,0.5, 0.0], [-0.5,0.5,0.0]]);
tri = pymesh.triangle();
tri = pymesh.triangle();
tri.points = vertices;
tri.split_boundary = True;
tri.conforming_delaunay = True;
tri.max_area = 0.01;
tri.run()
mesh = tri.mesh
pymesh.save_mesh("draped_rect_geometry.obj", mesh, use_double =True)