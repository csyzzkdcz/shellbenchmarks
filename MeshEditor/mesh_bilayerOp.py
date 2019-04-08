import pymesh
import numpy as np

path  =  "../benchmarks/TestModels/fine/kitten/draped_rect_geometry.obj";
mesh = pymesh.load_mesh(path);
nfaces = mesh.num_faces;
nverts =  mesh.num_vertices;
vertices = mesh.vertices;
for vert in vertices:
    vert_new = [];
    vert_new.append(vert[0]);
    vert_new.append(vert[1]);
    vert_new.append(-1);
    vert_new = np.array([vert_new]);
    vertices = np.concatenate((vertices, vert_new));
faces = mesh.faces;
for face in faces:
    face_new = [];
    face_new.append( face[0] + nverts);
    face_new.append( face[2] + nverts);
    face_new.append( face[1] + nverts);
    face_new = np.array([face_new]);
    faces = np.concatenate((faces, face_new));

mesh_new = pymesh.form_mesh(vertices, faces)
path_new = path[0:path.rfind(".")]+"_bilayer.obj";
pymesh.save_mesh(path_new, mesh_new, use_double=True)