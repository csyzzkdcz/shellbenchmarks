import pymesh
import numpy as np

path  = "../benchmarks/TestModels/coarse/bunnyHead/draped_rect_geometry_bilayer__sf.obj";
mesh = pymesh.load_mesh(path);
nfaces = mesh.num_faces;
nverts =  mesh.num_vertices;
vertices = mesh.vertices.tolist();
faces =  mesh.faces;
new_faces_0 = [];
new_vertices_0 = [];
new_faces_1 = [];
new_vertices_1 = [];
vertId_0 = 0;
vertId_1 = 0;
for face in faces:
    v = [];
    for vertid in face:
        v.append(vertid);
    if vertices[v[0]][2] == 0 and vertices[v[1]][2] == 0 and vertices[v[2]][2] == 0:
        new_face = [];
        for k in range(3):
            if vertices[v[k]] in new_vertices_0:
                new_face.append(new_vertices_0.index(vertices[v[k]]));
            else:
                new_face.append(vertId_0);
                new_vertices_0.append(vertices[v[k]]);
                vertId_0 = vertId_0 + 1;
        new_faces_0.append(new_face);
new_faces_0 = np.array(new_faces_0);
new_vertices_0 = np.array(new_vertices_0);

subpath = path[0: path.rfind("_")-1];
path_0 = subpath[0:subpath.rfind("_")]+ "_remeshed.obj";
mesh_0 = pymesh.form_mesh(new_vertices_0, new_faces_0);

pymesh.save_mesh(path_0, mesh_0, use_double=True);