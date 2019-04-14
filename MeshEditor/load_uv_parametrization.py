# import pymesh
import numpy as np

path  = "../benchmarks/TestModels/fine/bunnyHead/BFF_parametrization.obj";

with open(path, "r") as myMesh:
    uv = [];
    faces = [];
    lines = myMesh.read();
    lines = lines.split("\n");
    for line in lines:
        # print(line)
        line = line.split();
        # print(line)
        if(len(line) == 0):
            break;
        if line[0] == "vt":
            uv.append([float(line[1]), float(line[2]), 0.0])
        if line[0] == "f":
            faceInfo = line[1:len(line)];
            face = [];
            for vertex in faceInfo:
                vertices = vertex.split("/");
                if(vertices[1] == '2448'):
                    print(faceInfo);
                    print(face)
                face.append(int(vertices[1]));
            faces.append(face);
    # uv = np.array(uv);
    # faces = np.array(faces);
    # print(faces);
    path_0 = path[0:path.rfind("/")]+ "/uv_1.obj";
    # mesh = pymesh.form_mesh(uv, faces);
    # print(mesh.faces);
    # pymesh.save_mesh(path_0, mesh, use_double=True);

    with open(path_0, "w") as wf:
        for v in uv:
            wf.write("{0:<3}{1:<16.8f}{2:<16.8f}{3:<16.8f}\n".format("v", v[0],v[1],v[2]));
        for face in faces:
            wf.write("{0:<3}{1:<5d}{2:<5d}{3:<5d}\n".format("f", face[0],face[1],face[2]));