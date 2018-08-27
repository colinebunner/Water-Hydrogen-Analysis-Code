import numpy as np
import matplotlib.pyplot as plt
import pdb_reader 
import math_tools
import tip4p_energy as wpie

# Written by Colin Bunner some time in 2017

#### Cluster Class ############################################
# Contains the size of the cluster and the water objects that #
# form the cluster, including all of the detailed information #
# that comes with those objects                               #
###############################################################

class Cluster:

    # Initialize cluster object.
    # Cluster size: how many molecules are aggregated
    # Waters: list of water objects that comprise the cluster
    def __init__(self,cluster_size=0,waters=[]):
        self.size = cluster_size
        self.waters = waters
    
    # Getters
    def get_size(self):
        return self.size
    def get_waters(self):
        return waters

    # Setters
    def set_size(self,size):
        self.size = size
    def append_water(self,water):
        self.waters.append(water)

#############################
# H-bond Criteria Functions #
#############################

def is_hydrogenbond_ellipse(water_ref,water_obs,boxlen,hwidth,hheight,xc,yc,debug=False):

    #Input: A reference Water object and an observed Water object (in that order), as well as the
    #       boxlength. hwidth = half-width, hheight= half-height, xc = xcenter, yc = ycenter
    #
    # Please note that this analysis assumes the reference water is a potential hbond donor, while 
    # the observed water is a potential hbond acceptor

    # O-O distance
    oo_dist = math_tools.pbc_distance(water_ref.o_coords,water_obs.o_coords,boxlen)

    # Unit vector between the two oxygen atoms
    # Points from reference oxygen to observed oxygen
    oo_vec = math_tools.unit_vec_between(water_ref.o_coords,water_obs.o_coords)

    oh1_dist = math_tools.pbc_distance(water_ref.h1_coords,water_obs.o_coords,boxlen)
    oh2_dist = math_tools.pbc_distance(water_ref.h2_coords,water_obs.o_coords,boxlen)

    # Hydrogen bonding conditions from RADFs, using same ellipse approach as Wernet et al.
    if oh1_dist <= oh2_dist:
        angle = math_tools.vec_angle(water_ref.oh1_vec,oo_vec)
        ellipse_dist = ((oo_dist-xc)/hwidth)**2 + ((angle-yc)/hheight)**2
        if debug:
            print("{:<4d}{:<4d} OO Dist.: {}\tOH Dist.: {}\tAngle: {}\tHBond: {}".format(
                  water_ref.res_uid,water_obs.res_uid,oo_dist,oh1_dist,angle,ellipse_dist<=1))
        return ellipse_dist <= 1
    elif oh2_dist <= oh1_dist:
        angle = math_tools.vec_angle(water_ref.oh2_vec,oo_vec)
        ellipse_dist = ((oo_dist-xc)/hwidth)**2 + ((angle-yc)/hheight)**2
        if debug:
            print("{:<4d}{:<4d} OO Dist.: {}\tOH Dist.: {}\tAngle: {}\tHBond: {}".format(
                  water_ref.res_uid,water_obs.res_uid,oo_dist,oh2_dist,angle,ellipse_dist<=1))
        return ellipse_dist <= 1
    else:
        return False

def is_hydrogenbond_rectangular(water_ref,water_obs,boxlen,dcut,acut,debug=False):

    #Input: A reference Water object and an observed Water object (in that order), as well as the
    #       boxlength. dcut = OO-distance cutoff, acut = angle cutoff
    #
    # Please note that this analysis assumes the reference water is a potential hbond donor, while 
    # the observed water is a potential hbond acceptor

    # O-O distance
    oo_dist = math_tools.pbc_distance(water_ref.o_coords,water_obs.o_coords,boxlen)

    # Unit vector between the two oxygen atoms
    oo_vec = math_tools.unit_vec_between(water_ref.o_coords,water_obs.o_coords)

    # Angle between the two vectors. Must be careful because arccos is not single_valued.
    # Currently avoid this issue by ignoring the O-H pairing if the O_ref-O_obs distance is shorter
    # than the H_ref-O_obs distance (i.e. the hydrogen is pointing away) 
    oh1_dist = math_tools.pbc_distance(water_ref.h1_coords,water_obs.o_coords,boxlen)
    oh2_dist = math_tools.pbc_distance(water_ref.h2_coords,water_obs.o_coords,boxlen)

    # Rectangular H-bond conditions
    if oh1_dist <= oh2_dist:
        angle = math_tools.vec_angle(water_ref.oh1_vec,oo_vec)
        if debug:
            print("{:<4d}{:<4d} OO Dist.: {}\tOH Dist.: {}\tAngle: {}\tHBond: {}".format(
                  water_ref.res_uid,water_obs.res_uid,oo_dist,oh1_dist,angle,angle <= 30 and
                  oo_dist <= 3.5))
        return (oo_dist <= 3.5 and angle <= 40)
    elif oh2_dist <= oh1_dist:
        angle = math_tools.vec_angle(water_ref.oh2_vec,oo_vec)
        if debug:
            print("{:<4d}{:<4d} OO Dist.: {}\tOH Dist.: {}\tAngle: {}\tHBond: {}".format(
                  water_ref.res_uid,water_obs.res_uid,oo_dist,oh2_dist,angle,angle <= 30 and
                  oo_dist <= 3.5))
        return (oo_dist <= 3.5 and angle <= 40)
    else:
        return False

def is_hydrogenbond_distanceonly(water_ref,water_obs,boxlen,cutoff,debug=False):

    #Input: A reference Water object and an observed Water object (in that order), as well as the
    #       boxlength
    #
    # Please note that this analysis assumes the reference water is a potential hbond donor, while 
    # the observed water is a potential hbond acceptor

    # O-O distance
    oo_dist = math_tools.pbc_distance(water_ref.o_coords,water_obs.o_coords,boxlen)

    # OO-distance only criterion
    return oo_dist <= cutoff

def is_hydrogenbond_energetic(water_ref,water_obs,boxlen,ecut):
    # Input: Reference water, observed water, box length (assumes cubic), and energy cutoff in kJ/mol
   
    pair_energy = wpie.TIP4P2005_pairenergy(water_ref,water_obs,boxlen,separate_components=False)    
    
    return pair_energy < ecut

###########################
# Analyzes a single movie #
###########################
    
def analyze_movie(movie,hbond_type,model,frame_name="box",frame_id="Prod: 1, Movie: 1",debug=False):

    # Input: Movie file in pdb format, hbond_type determines which definition is used
    # Output: 1) A list of all water objects
    #         2) A list of all cluster objects
    #         3) Number of hydrogen bonds per water molecule

    # Make all beads from pdb file
    waters = pdb_reader.pdb_to_waters(movie)
    box_size = pdb_reader.get_boxsize(movie)

    # Clusters
    clusters = []
    if debug:
        pair_uids = [] 

    if hbond_type == "ellipse" or hbond_type == "rectangular":

        # Check every water molecule with every other water molecules since this
        # H-bond definition type assumes that the first water molecule passed is
        # the donor
        for i in range(len(waters)):
            for j in range(len(waters)):
                # Check H-bonding criteria
                # If aggregated (as defined by definition in function is_hydrogenbond())
                # either make a new cluster or add it to an existing cluster
                if hbond_type == "ellipse":
                   hbond =  is_hydrogenbond_ellipse(waters[i],waters[j],boxlen=box_size,hwidth=0.70,hheight=52.26,xc=2.90,yc=0)
                elif hbond_type == "rectangular":
                   hbond = is_hydrogenbond_rectangular(waters[i],waters[j],boxlen=box_size,dcut=3.2,acut=40)
                if j != i and hbond:

                    if debug:
                        pair_uids.append((waters[i].res_uid,waters[j].res_uid))

                    # Keep track of donor statistics
                    waters[i].ndonor = waters[i].ndonor + 1
                    waters[j].nacceptor = waters[j].nacceptor + 1

                    #Decide whether to create a new cluster, add a molecule to an existing cluster,
                    #and/or merge two existing clusters

                    # Case 1: The reference water is in a cluster but the observed water is not.
                    # Add the observed water to the cluster
                    if waters[i].incluster and not waters[j].incluster:
                        # Find the cluster that we would like to add the observed water to
                        cluster_index = map_bead_to_cluster(clusters,waters[i]) 
                        # Add observed water to said cluster
                        clusters[cluster_index].append_water(waters[j])
                        # Increment the cluster size by 1
                        clusters[cluster_index].set_size(clusters[cluster_index].size+1)
                        # Change the incluster bool to True for the observed water
                        waters[j].set_incluster(True)
                    # Case 2: Same as Case 1, but now the reference molecule is not in a cluster while
                    #         the observed water molecule is
                    elif not waters[i].incluster and waters[j].incluster:
                        cluster_index = map_bead_to_cluster(clusters,waters[j]) 
                        clusters[cluster_index].append_water(waters[i])
                        clusters[cluster_index].set_size(clusters[cluster_index].size+1)
                        waters[i].set_incluster(True)
                    # Case 3: Both water molecules are already in a cluster. Merge one cluster into another
                    #         and annihilate the extra
                    elif waters[i].incluster and waters[j].incluster:
                        cluster_index1 = map_bead_to_cluster(clusters,waters[i])
                        cluster_index2 = map_bead_to_cluster(clusters,waters[j])
                        # Covers the case that the two water molecules are part of the same cluster
                        # (may indicate something like ring-closing)
                        if cluster_index1 == cluster_index2:
                            continue
                        # Two clusters meet for the first time
                        else:
                            for w in clusters[cluster_index2].waters:
                                clusters[cluster_index1].append_water(w)
                                clusters[cluster_index1].set_size(clusters[cluster_index1].size+1)
                            # Remove extra cluster
                            clusters.remove(clusters[cluster_index2])
                    # Case 4: Neither molecule is yet identified as part of a cluster. Start a new one.
                    else:
                        clusters.append(Cluster(cluster_size=2,waters=[waters[i],waters[j]]))
                        waters[i].set_incluster(True)
                        waters[j].set_incluster(True)
                else:
                    continue

    elif hbond_type == "distance_only" or hbond_type == "energetic":

        # Check unique pairs of water molecules only once because otherwise you double count 
        for i in range(len(waters)-1):
            for j in range(i+1,len(waters)):
                # Check H-bonding criteria
                # If aggregated (as defined by definition in function is_hydrogenbond())
                # either make a new cluster or add it to an existing cluster
                if hbond_type == "distance_only":
                   hbond =  is_hydrogenbond_distance_only(waters[i],waters[j],boxlen=box_size,dcut=3.2)
                elif hbond_type == "energetic":
                   hbond = is_hydrogenbond_rectangular(waters[i],waters[j],boxlen=box_size,ecut=-12.0)
                if hbond:

                    if debug:
                        pair_uids.append((waters[i].res_uid,waters[j].res_uid))

                    # Keep track of donor statistics
                    waters[i].ndonor = waters[i].ndonor + 1
                    waters[j].nacceptor = waters[j].nacceptor + 1

                    #Decide whether to create a new cluster, add a molecule to an existing cluster,
                    #and/or merge two existing clusters

                    # Case 1: The reference water is in a cluster but the observed water is not.
                    # Add the observed water to the cluster
                    if waters[i].incluster and not waters[j].incluster:
                        # Find the cluster that we would like to add the observed water to
                        cluster_index = map_bead_to_cluster(clusters,waters[i]) 
                        # Add observed water to said cluster
                        clusters[cluster_index].append_water(waters[j])
                        # Increment the cluster size by 1
                        clusters[cluster_index].set_size(clusters[cluster_index].size+1)
                        # Change the incluster bool to True for the observed water
                        waters[j].set_incluster(True)
                    # Case 2: Same as Case 1, but now the reference molecule is not in a cluster while
                    #         the observed water molecule is
                    elif not waters[i].incluster and waters[j].incluster:
                        cluster_index = map_bead_to_cluster(clusters,waters[j]) 
                        clusters[cluster_index].append_water(waters[i])
                        clusters[cluster_index].set_size(clusters[cluster_index].size+1)
                        waters[i].set_incluster(True)
                    # Case 3: Both water molecules are already in a cluster. Merge one cluster into another
                    #         and annihilate the extra
                    elif waters[i].incluster and waters[j].incluster:
                        cluster_index1 = map_bead_to_cluster(clusters,waters[i])
                        cluster_index2 = map_bead_to_cluster(clusters,waters[j])
                        # Covers the case that the two water molecules are part of the same cluster
                        # (may indicate something like ring-closing)
                        if cluster_index1 == cluster_index2:
                            continue
                        # Two clusters meet for the first time
                        else:
                            for w in clusters[cluster_index2].waters:
                                clusters[cluster_index1].append_water(w)
                                clusters[cluster_index1].set_size(clusters[cluster_index1].size+1)
                            # Remove extra cluster
                            clusters.remove(clusters[cluster_index2])
                    # Case 4: Neither molecule is yet identified as part of a cluster. Start a new one.
                    else:
                        clusters.append(Cluster(cluster_size=2,waters=[waters[i],waters[j]]))
                        waters[i].set_incluster(True)
                        waters[j].set_incluster(True)

                else:
                    continue

# Output files with detailed statistics
    
    if debug:
        with open("pair_uids.txt","a+") as out:
            out.write("{}\n".format(frame_id))
            for uid in pair_uids:
                out.write("{:>5s}{:>5s}\n".format(str(uid[0]),str(uid[1])))

    with open(somefile,"a+") as out:
        out.write("{}\n".format(frame_id))
        for water in waters:
            out.write("{:>5d}\tndonor: {:<3d}\tnacceptor: {:<3d}\n".format(water.res_uid,water.ndonor,
              water.nacceptor))

    # Add clusters of size 1 (isolated waters)
    for water in waters:
        if not water.incluster:
            clusters.append(Cluster(cluster_size=1,waters=[water]))


    with open("/home/cbunner/BetterH2OH2Mix/DATA/ca{}/{}_{}_clusters.txt".format(model,model,frame_name),"a+") as out:
        out.write("{}\n".format(frame_id))
        for cluster in clusters:
            out.write("{}\n".format(cluster.size))

# Will print coordinates of all clusters in xyz format
#    with open("ca/{}_segmented_clusters.txt".format(frame_name),"a+") as out:
#        out.write("{}\n".format(frame_id))
#        for cluster in clusters:
#            for water in cluster.waters:
#                out.write("{:<4s}{:<12.5f}{:<12.5f}\t{:<12.5f}\n".format("O",*water.o_coords))
#                out.write("{:<4s}{:<12.5f}{:<12.5f}\t{:<12.5f}\n".format("H",*water.h1_coords))
#                out.write("{:<4s}{:<12.5f}{:<12.5f}\t{:<12.5f}\n".format("H",*water.h2_coords))
#            out.write("\n")

    nper = 0

    for water in waters:
        nper += water.ndonor + water.nacceptor

    nper /= len(waters)

    return waters, clusters, nper 

def cluster_statistics(cluster_list):

    # Input:  A list of cluster objects
    # Output: Average cluster size, maximum cluster size, cluster size frequency

    # Make a dictionary to keep track of the number of clusters of each type
    cluster_frequency = {}

    # Track the average cluster size and maximum size
    avg_size = 0.
    max_size = 2

    # Loop over all clusters to find the average cluster size and maximum
    # cluster size. Also, keep track of the number of clusters of each size
    for cluster in cluster_list:
        
        # Add cluster size to total cluster size
        avg_size += cluster.size
        
        # Change the maximum cluster if a larger one is found
        if cluster.size > max_size:
            max_size = cluster.size
            max_cluster = cluster

        # Cluster size distribution
        try:
            cluster_frequency[cluster.size] += 1
        except KeyError:
            cluster_frequency[cluster.size] = 1
       
    # Normalize the cluster list
    for key in cluster_frequency.keys():
        cluster_frequency[key] = cluster_frequency[key]/len(cluster_list)        

    return (avg_size/len(cluster_list)), max_size, cluster_frequency

# Find the cluster object a bead belongs to. Returns the index of the cluster containing the bead
# in the cluster list passed as an argument. Since this analysis program works frame by frame,
# the res_uid should be unique.
def map_bead_to_cluster(cluster_list,test_water):

    unique_id = test_water.res_uid

    for i, cluster in enumerate(cluster_list):
        for water in cluster.waters:
            if water.res_uid == unique_id:
                return i

    print("map_bead_to_cluster_failed")
    return -1

