import numpy as np
import math_tools

# Written by Colin Bunner (bunne043@umn.edu) some time in 2017
 
def lennardjones_potential(sigma,eps,r):
    return 4*eps*((sigma/r)**12-(sigma/r)**6)

def TIP4P2005_pairenergy(h11,h21,o1,m1,h12,h22,o2,m2,box_len):
  
    # Instantiate the interaction energy
    interaction_energy = 0.0

    ### Add Coulombic potential energy ###

    # Partial charges in C
    qh = 0.5564*(1.60217662E-19)
    qo = -2*qh
    # 4 * pi * eps_0
    fourpieps = 4*np.pi*8.854187816E-12
    # Hydrogen and oxygen LJ parameters
    oxygen_eps = 93.20*1.380658E-23
    oxygen_sig = 3.1589E-10

    # Distances between charged beads in Angstroem
    h1_ref_to_h1_obs = math_tools.pbc_distance(h11,h12,box_len)*(1E-10)
    h1_ref_to_h2_obs = math_tools.pbc_distance(h11,h22,box_len)*(1E-10)
    h1_ref_to_m_obs = math_tools.pbc_distance(h11,m2,box_len)*(1E-10)

    h2_ref_to_h1_obs = math_tools.pbc_distance(h21,h12,box_len)*(1E-10)
    h2_ref_to_h2_obs = math_tools.pbc_distance(h21,h22,box_len)*(1E-10)
    h2_ref_to_m_obs = math_tools.pbc_distance(h21,m2,box_len)*(1E-10)

    m_ref_to_h1_obs = math_tools.pbc_distance(m1,h12,box_len)*(1E-10)
    m_ref_to_h2_obs = math_tools.pbc_distance(m1,h22,box_len)*(1E-10)
    m_ref_to_m_obs = math_tools.pbc_distance(m1,m2,box_len)*(1E-10)

    # Coulombic energies
    h1ref_h1obs = (qh*qh)/(fourpieps*h1_ref_to_h1_obs)
    h1ref_h2obs = (qh*qh)/(fourpieps*h1_ref_to_h2_obs)
    h1ref_mobs = (qh*qo)/(fourpieps*h1_ref_to_m_obs)

    h2ref_h1obs = (qh*qh)/(fourpieps*h2_ref_to_h1_obs)
    h2ref_h2obs = (qh*qh)/(fourpieps*h2_ref_to_h2_obs)
    h2ref_mobs = (qh*qo)/(fourpieps*h2_ref_to_m_obs)
 
    mref_h1obs = (qh*qo)/(fourpieps*m_ref_to_h1_obs)
    mref_h2obs = (qh*qo)/(fourpieps*m_ref_to_h2_obs)
    mref_mobs = (qo*qo)/(fourpieps*m_ref_to_m_obs)
   
    interaction_energy += (h1ref_h1obs + h1ref_h2obs + h1ref_mobs + h2ref_h1obs +  
      h2ref_h2obs + h2ref_mobs + mref_h1obs + mref_h2obs + mref_mobs)

    coulombic_energy = interaction_energy

    # Now for the van der Waals energies between oxygens
    oo_dist = math_tools.pbc_distance(o1,o2,box_len)*(1E-10)
    interaction_energy += lennardjones_potential(oxygen_sig,oxygen_eps,oo_dist)
    
    # Go from J to kJ/mol
    interaction_energy *= 6.022E20

    return interaction_energy
 
