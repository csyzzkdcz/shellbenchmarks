import pymesh
import numpy as np


with open("/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/boundary.txt","r") as f:
    rows = [];
    line1 = f.readline();
    line1.split("\n");
    print(line1)
    nverts = int(line1[0]);
    print(nverts)
    f = f.read();
    f = f.split("\n");
    i = 0;
    for line in f:
        if(line == ''):
            break;
        line = line.split(" ")
        column = [];
        print(line)
        column.append(line[0])
        column.append(line[1])
        column.append(line[2])
        rows.append(column)
        i = i+1;
    rows = [[-0.5,-0.5,0],[-0.5,0.5,0],[0.5,0.5,0],[0.5,-0.5,0]]
    vertices = np.array(rows);
    print(vertices)
    tri = pymesh.triangle();
    tri = pymesh.triangle();
    tri.points = vertices;
    tri.split_boundary = True;
    tri.conforming_delaunay = False;
    tri.max_area = 0.001;
    tri.run()
    mesh = tri.mesh
    pymesh.save_mesh("/Users/chenzhen/UT/Research/Projects/shellbenchmarks/draped_rect_geometry_remeshed.obj", mesh, use_double =True)