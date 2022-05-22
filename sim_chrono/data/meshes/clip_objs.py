import numpy as np
import matplotlib.pyplot as plt

w = 0.2
h = 0.4
gap = 0.05
cs_r = 0.01
f_r = 0.05

is_capped = True

n_cs = 8


skleton = np.array([
    [0.5*w, 0.5*gap, 0],

    [0.5*w, 0.5*h-f_r, 0],
    [0.5*w-(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), 0.5*h-(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), 0],
    [0.5*w-f_r, 0.5*h, 0],

    [-0.5*w+f_r, 0.5*h, 0],
    [-0.5*w+(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), 0.5*h-(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), 0],
    [-0.5*w, 0.5*h-f_r, 0],

    [-0.5*w, -0.5*h+f_r, 0],
    [-0.5*w+(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), -0.5*h+(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), 0],
    [-0.5*w+f_r, -0.5*h, 0],

    [0.5*w-f_r, -0.5*h, 0],
    [0.5*w-(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), -0.5*h+(np.sqrt(2.0)-1.0)*f_r/np.sqrt(2.0), 0],
    [0.5*w, -0.5*h+f_r, 0],

    [0.5*w, -0.5*gap, 0]
    ], dtype=float)

cs_normal = np.array([
    [[1., 0, 0], [0, 1., 0]], 
    [[1., 0, 0], [0, 1., 0]], 
    [[1./np.sqrt(2.), 1./np.sqrt(2.), 0], [-1./np.sqrt(2.), 1./np.sqrt(2.), 0]], 
    [[0, 1., 0], [-1., 0, 0]],
    [[0, 1., 0], [-1., 0, 0]],
    [[-1./np.sqrt(2.), 1./np.sqrt(2.), 0], [-1./np.sqrt(2.), -1./np.sqrt(2.), 0],],
    [[-1., 0, 0], [0, -1., 0]],
    [[-1., 0, 0], [0, -1., 0]],
    [[-1./np.sqrt(2.), -1./np.sqrt(2.), 0], [1./np.sqrt(2.), -1./np.sqrt(2.), 0]],
    [[0, -1., 0], [1., 0, 0]],
    [[0, -1., 0], [1., 0, 0]],
    [[1./np.sqrt(2.), -1./np.sqrt(2.), 0], [1./np.sqrt(2.), 1./np.sqrt(2.), 0]], 
    [[1., 0, 0], [0, 1., 0]],
    [[1., 0, 0], [0, 1., 0]]
])

def rotated(vec, axis, angle):
    return np.dot(axis, vec) * axis + np.cos(angle) * np.cross(np.cross(axis, vec), axis) + np.sin(angle) * np.cross(axis, vec)

# plt.plot(skleton[:, 0], skleton[:, 1])
# plt.axis('equal')
# plt.show()
n_layers = skleton.shape[0]
cs_nodes = np.zeros((n_layers, n_cs, 3))
for n_cntr in range(n_layers):
    for n_bound in range(n_cs):
        cs_nodes[n_cntr, n_bound, :] = skleton[n_cntr, :] + \
            cs_r * rotated(cs_normal[n_cntr, 0, :], cs_normal[n_cntr, 1, :], n_bound * 2.*np.pi/n_cs)

cs_nodes = np.reshape(cs_nodes, (n_cs*n_layers, 3))
# fig = plt.figure()
# ax = fig.add_subplot(111,projection='3d')
# ax.plot(cs_nodes[:, 0], cs_nodes[:, 1], cs_nodes[:, 2], "ko")
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# # ax.axis('equal')
# plt.show()

faces = np.zeros((2*n_cs*(n_layers-1), 3))
for cs in range(1, n_layers):
    for bd in range(1, n_cs):
        faces[2*(bd-1 + n_cs*(cs-1)), :] = np.array([(bd-1)+n_cs*(cs-1)+1, (bd-1)+n_cs*cs+1, bd+n_cs*(cs-1)+1])
        faces[2*(bd-1 + n_cs*(cs-1))+1, :] = np.array([(bd-1)+n_cs*cs+1, bd+n_cs*cs+1, bd+n_cs*(cs-1)+1])
    #     print(2*(bd-1 + n_cs*(cs-1)), 2*(bd-1 + n_cs*(cs-1))+1)
    
    faces[2*(n_cs-1 + n_cs*(cs-1)), :] = np.array([(n_cs-1)+n_cs*(cs-1)+1, (n_cs-1)+n_cs*cs+1, n_cs*(cs-1)+1])
    faces[2*(n_cs-1 + n_cs*(cs-1))+1, :] = np.array([(n_cs-1)+n_cs*cs+1, n_cs*cs+1, n_cs*(cs-1)+1])
    # print(2*(n_cs-1 + n_cs*(cs-1)), 2*(n_cs-1 + n_cs*(cs-1))+1)

if is_capped and faces.shape[0]==2*n_cs*(n_layers-1):
    id_init = cs_nodes.shape[0]+1
    id_end = cs_nodes.shape[0]+2
    cs_nodes = np.vstack([cs_nodes, skleton[0]])
    cs_nodes = np.vstack([cs_nodes, skleton[-1]])
    for bd in range(1, n_cs):
        faces = np.vstack([faces, np.array([bd+1, id_init, bd])])
        faces = np.vstack([faces, np.array([id_init-bd, id_end, id_init-bd-1])])
    faces = np.vstack([faces, np.array([1, id_init, n_cs])])
    faces = np.vstack([faces, np.array([id_init-n_cs, id_end, id_init-1])])


with open("sample.obj", "w") as ofile:
    for ni in cs_nodes:
        ofile.write("v " + str(ni[0]) + "  " + str(ni[1]) + " " + str(ni[2]) + "\n")
    ofile.write("s off \n")
    for nf in faces:
        ofile.write("f " + str(nf[0]) + "  " + str(nf[1]) + " " + str(nf[2]) + "\n")

ofile.close()
