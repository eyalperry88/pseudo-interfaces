# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 17:36:34 2015

@author: yot
"""

import os

import Bio
from Bio.PDB import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder

NUM_OF_RESULTS_TO_SAVE = 100000

def pdb_writer(graph, results, dir_name="templates"):
    PDB_CODE = graph.protein_name 
    structure = graph.structure
    
    path = "%s/%s/" % (dir_name, PDB_CODE)
    if len(results) > 0:
        if not os.path.exists(path):
            os.makedirs(path) 
    else:
        print 'No results'
        return
    
    top = min(NUM_OF_RESULTS_TO_SAVE, len(results))
    i=0
    for key, value in results[:top]:
        # TODO - save as different chains
        sb = StructureBuilder()
        sb.init_structure(PDB_CODE)
        sb.init_model(0)
        sb.init_chain('X')
        for node_id in key[0]:
            for res in graph.nodes[node_id].residues:
                sb.init_seg(res.get_segid())
                sb.init_residue(res.get_resname(), res.get_id()[0],
                                res.get_id()[1], res.get_id()[2])
                for atom in res:
                    sb.init_atom(atom.get_name(), atom.get_coord(), atom.get_bfactor(),
                                 atom.get_occupancy(), atom.get_altloc(), atom.get_fullname())
        
        sb.init_chain('Y')
        for node_id in key[1]:
            for res in graph.nodes[node_id].residues:
                sb.init_residue(res.get_resname(), res.get_id()[0],
                                res.get_id()[1], res.get_id()[2])
                for atom in res:
                    sb.init_atom(atom.get_name(), atom.get_coord(), atom.get_bfactor(),
                                 atom.get_occupancy(), atom.get_altloc(), atom.get_fullname())
        filename = path + "interface%d.pdb"%(i)
        io = PDBIO()
        io.set_structure(sb.get_structure())
        io.save(filename)
        i+=1
        

    return 1
def pymol_scripts_creator(graph, results, dir_name="results"):
    PDB_CODE = graph.protein_name   
    
    path = "%s/%s/" % (dir_name, PDB_CODE)
    if len(results) > 0:
        if not os.path.exists(path):
            os.makedirs(path) 
    else:
        print 'No results'
        return
    
    top = min(NUM_OF_RESULTS_TO_SAVE, len(results))
    i=0
    for key, value in results[:top]:
        filename = path + "interface%d.pymol"%(i)
        pymol_file_writer(graph,key,filename)
        i+=1
        

    return 1


def pymol_file_writer(graph, key, filename):
    fd = open(filename,'w')
    text = ["#automated script for pymol selection of interface\n"]
    text+= ["delete Side_A or Side_B\nutil.cbc\n"]
    text+= [pymol_selection(graph,key[0],"Side_A")+"\n"]
    text+= ["print \"Side_! is constructed from nodes %s\"\n"%(str(key[0]))]
    text+= ["color red, Side_A\n"]
    
    text+= [pymol_selection(graph,key[1],"Side_B")+"\n"]
    text+= ["print \"Side_B is constructed from nodes %s\"\n"%(str(key[1]))]
    text+= ["color blue, Side_B\n"]
    
    
    text+= ["show cartoon, Side_A or Side_B\n"]
    text+= ["select None\n"]
    text+= ["deselect"]
    
    fd.writelines(text)
    fd.close()
    




def pymol_selection(graph, triplet_indices, name):
    lst = [graph.nodes[i] for i in triplet_indices]
    inp = []
    result = ['select %s, '%(name)]
    for i in range(3):
        inp.append([lst[i].start,lst[i].start+len(lst[i])-1,lst[i].chain])
        result.append('(resi %d-%d in chain %c)'%(inp[i][0],inp[i][1],inp[i][2]))
    output =  result[0] + result[1] + ' or ' + result[2] + ' or ' + result[3]
    return output
    