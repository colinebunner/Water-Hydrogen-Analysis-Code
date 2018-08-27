import numpy as np

# Written by Colin Bunner (bunne043@umn.edu) some time in 2017

def unit_vector(v):
    # Normalize a vector
    v = np.array(v)
    return v/np.sqrt(np.dot(v,v))

def unit_vec_between(v1,v2):
    # Points from v1 to v2
    x = float(v2[0]) - float(v1[0])
    y = float(v2[1]) - float(v1[1])
    z = float(v2[2]) - float(v1[2])

    vec = unit_vector((x,y,z))

    return vec

def vec_angle(v1,v2):

    # Return the angle between two vectors (in degrees)
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    return (np.arccos(np.dot(v1_u,v2_u)))*(180./np.pi)

def pbc_distance(bead1_coords,bead2_coords,boxlen):

    xdist = bead1_coords[0]-bead2_coords[0]
    ydist = bead1_coords[1]-bead2_coords[1]
    zdist = bead1_coords[2]-bead2_coords[2]

    xdist = xdist - boxlen*int(round(xdist/boxlen))
    ydist = ydist - boxlen*int(round(ydist/boxlen))
    zdist = zdist - boxlen*int(round(zdist/boxlen))

    return np.sqrt(xdist**2 + ydist**2 + zdist**2)


