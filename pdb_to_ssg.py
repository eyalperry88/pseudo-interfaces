# -*- coding: utf-8 -*-
"""
Transform a PDB file to an SS graph
"""

import copy
import heapq
import numpy as np

from Bio.PDB import PDBParser, DSSP

from data_structures import SSGraph

HELIX_CODES = ['H', 'G', 'I']
BETA_CODES = ['B', 'E']

# minimum lengths for secondary structures
MIN_STRUCT_LENGTH = {'HELIX' : 6, 'BETA' : 4}

# number of residues distances to use in averaged distance
NUM_RESIDUES_FOR_AVG_DIST = 4

# min number of residues threshold between neighboring secondary structures
MIN_NUM_RESIDUES = 12

# below this distance, two residues are considered interacting
INTERACTION_DIST_THRESHOLD = 6.0

# max distance threshold between neighboring secondary structures
MAX_DIST_THRESHOLD = 24.0

  
def check_edge(node1, node2):
    """ 
    Returns information about the relationship of node1 & node2
    
    node1 & node2 are of type SSNode
    
    Returns a tuple of (seq, dist, inter_count). 
    
    seq is True iff node1 and node2 are sequential - the num. of AA between them
    is less than MIN_NUM_RESIDUES
    
    dist is the number of interacting residues between node1 & node2.
    Based on Jiang et al. (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323951/)
    we define two residues to be interacting if the distance between two C_beta
    (C_alpha for Gly) is less than 6.0 Angstrom.
    
    dist is defined as the average of the minimal NUM_RESIDUES_FOR_AVG_DIST distances
    """
    # check residue distance
    seq = False
    if node1.chain == node2.chain:
        if node1.start < node2.start:
            residue_dist = node2.start - (node1.start + len(node1))
        else:
            residue_dist = node1.start - (node2.start + len(node2))
        if residue_dist < MIN_NUM_RESIDUES:
            seq = True
        
    # calculate averaged distance over min(n,m) - 2 C-alphas
    # and also number of interacting atoms
    ca_dists = []    
    interacting_count = 0
    exclude_atoms = ['O', 'N', 'C']
    for i in range(len(node1)):
        ca1 = node1.residues[i]['CA']
        for j in range(len(node2)):
            ca2 = node2.residues[j]['CA']
            ca_dists.append(ca1 - ca2)
            
            found_inter = False
            for atom1 in node1.residues[i]:
                if found_inter: break
                if atom1.get_name() in exclude_atoms:
                    continue
                for atom2 in node2.residues[j]:
                    if found_inter: break
                    if atom2.get_name() in exclude_atoms:
                        continue
                    d = atom1 - atom2
                    if d < INTERACTION_DIST_THRESHOLD:
                        interacting_count += 1
                        found_inter = True
    #n = min(len(node1), len(node2)) - 2
    n = min(NUM_RESIDUES_FOR_AVG_DIST, len(node1), len(node2))
    avg = float(np.average(heapq.nsmallest(n, ca_dists)))     
    
    return (seq, avg, interacting_count)
    

def PDB2SSG(pdb_filename):
    """ 
    Converts a PDB file to a Secondary Structure Graph
    
    Gets the path to a PDB filename, returns an SSGraph    
    """
    # load PDB
    p = PDBParser()
    protein_name = pdb_filename.split('.')[0].split('/')[-1]
    structure = p.get_structure(protein_name, pdb_filename)
    model = structure[0]
    dssp = DSSP(model, pdb_filename)
    print 'Loaded structure for %s (Number of residues: %d)' % (protein_name, len(dssp))
    
    # init graph
    ssg = SSGraph(protein_name, structure)
    
    # first, add nodes
    first = list(dssp)[0][0]
    current_chain = first.get_full_id()[2]
    start_idx = first.get_id()[1]
    current_ss = None  # the current secondary structure
    res_list = []  # residue list for the current structure
    previous_idx = start_idx
    for residue in dssp:
        # struct will be HELIX/BETA/None
        struct = None
        if residue[1] in HELIX_CODES:
            struct = 'HELIX'
        elif residue[1] in BETA_CODES:
            struct = 'BETA'
            
        full_id = residue[0].get_full_id()
        chain = full_id[2]
        current_idx = full_id[3][1]
        
        if chain != current_chain or struct != current_ss:  # if we just changed from 2 different structures
            if current_ss != None:  # if the current structure is meaningful
                res_length = len(res_list)
                if res_length >= MIN_STRUCT_LENGTH[current_ss]: # if the structure length is sufficient
                    # add to graph    
                    node_id = ssg.addNode(current_ss, res_length, start_idx, 
                                          current_chain, copy.copy(res_list))
                    # print('Added node %d (type: %s, length: %d, start: %s, chain: %s)' %
                    # (node_id, current_ss, res_length, start_idx, current_chain))
            
            # set current struct to new one, reset index and residue list
            current_ss = struct
            current_chain = chain
            start_idx = current_idx
            res_list = []
        
        
        if struct:  # if struct is meaningful, append to the building residue list
            res_list.append(residue[0])
       
        previous_idx = current_idx
        
    # second, add edges
    edge_counter = 0
    extended_edge_counter = 0
    for node_id1 in range(len(ssg.nodes)-1):
        for node_id2 in range(node_id1+1, len(ssg.nodes)):
            seq, dist, inter_count = check_edge(ssg.nodes[node_id1], ssg.nodes[node_id2])
            if not seq and dist <= MAX_DIST_THRESHOLD:
                ssg.addEdge(node_id1, node_id2, dist, inter_count)
                edge_counter += 1
            if seq or dist <= MAX_DIST_THRESHOLD:
                ssg.addExtendedEdge(node_id1, node_id2, dist, inter_count)
                extended_edge_counter += 1
    print 'The graph contains %d nodes, %d edges and %d extended edges.' % (len(ssg.nodes), edge_counter, extended_edge_counter)
                
    # return the graph
    return ssg
