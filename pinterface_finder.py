# -*- coding: utf-8 -*-
"""
Pseudo-interface finder

The program receives a protein structure and outputs potential candidates
of pseudo-interfaces in two formats:
1) PDB files containing the structural information of the pseudo-interface, 
   one side of the interface is marked as chain X and the other marked as chain Y.
2) pymol script files which can be executed with the whole protein structure
   in pymol, where the sides of the interface are marked with red and blue.
   
Usage information:

pinterface_finder.py PATH_TO_PDB NUM_OF_RESULTS
 
e.g.
pinterface_finder.py proteins/2ZKM.pdb 100

"""

import argparse

parser = argparse.ArgumentParser(description='Pseudo-interface finder')
parser.add_argument("PDB", help="Path to PDB file")
parser.add_argument("N", default=20, type=int, help="Number of results to output. Default is 20")
parser.add_argument("--core", action="store_true",
                    help="Output only PIs from the protein core (no inter-chain). Default is false")
parser.add_argument("--sd", type=float, default=94.0,
                    help="SD threshold parameter, see \"Revealing pseudo-interfaces\" in paper. D efault is 94")
args = parser.parse_args()

from pdb_to_ssg import PDB2SSG
from MMKKM import MMKKM
from sort_results import sort_results_by_distance_score
from output_pipes import pymol_scripts_creator, pdb_writer

ssg = PDB2SSG(args.PDB)
results = MMKKM(ssg, th=args.sd)
sorted_results = sort_results_by_distance_score(results)

i = 0
final_results = []
print 'Found %d candidates for pseudo-interfaces in %s' % (len(sorted_results), ssg.protein_name)

if args.core:
    core_results = []
    for key, value in sorted_results:
        i += 1
        core_interface = False
        break_flag = False
        chain1 = ssg.nodes[key[0][0]].chain
        for node_id in key[0]:
            node = ssg.nodes[node_id]
            if node.chain != chain1:
                break_flag = True
                break
        chain2 = ssg.nodes[key[1][0]].chain
        for node_id in key[1]:
            node = ssg.nodes[node_id]
            if node.chain != chain2:
                break_flag = True
                break
        if chain1 == chain2 and not break_flag:
            core_results.append((key, value))
    
    print '%d of the pseudo-interface lie in the protein core' % len(core_results)
    final_results = core_results[:args.N]
else:
    final_results = sorted_results[:args.N]
    
print "Interface #1 has SD score %f (IC: %d) " % (final_results[0][1][0], final_results[0][1][1])       
pymol_scripts_creator(ssg, final_results, "scripts")
pdb_writer(ssg, final_results, "templates")
    

