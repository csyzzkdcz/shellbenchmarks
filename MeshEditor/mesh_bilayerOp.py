import pymesh
import numpy as np

path  =  "../benchmarks/TestModels/coarse/bunnyHead/draped_rect_geometry.obj";
mesh = pymesh.load_mesh(path);
nfaces = mesh.num_faces;
nverts =  mesh.num_vertices;
vertices = mesh.vertices;
boundary_loops = mesh.boundary_loops;

for vert in vertices:
    vert_new = [];
    vert_new.append(vert[0]);
    vert_new.append(vert[1]);
    vert_new.append(-0.1);
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

for boundary_loop in boundary_loops:
    loop_size = boundary_loop.size;
    for i in range(loop_size):
        face_new = [boundary_loop[i], boundary_loop[i] + nverts, boundary_loop[(i+1)%loop_size] + nverts];
        face_new = np.array([face_new]);
        faces = np.concatenate((faces, face_new));

        face_new = [boundary_loop[(i+1)%loop_size] + nverts, boundary_loop[(i+1)%loop_size], boundary_loop[i]];
        face_new = np.array([face_new]);
        faces = np.concatenate((faces, face_new));


mesh_new = pymesh.form_mesh(vertices, faces)
path_new = path[0:path.rfind(".")]+"_bilayer.obj";
pymesh.save_mesh(path_new, mesh_new, use_double=True)