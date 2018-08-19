# Identifying novel pseudo-interfaces in protein cores.
By Yotam Constantini and Eyal Perry.

Identifying novel interface templates by examining the interior structure of resolved proteins.


## Project details
https://github.com/eyalperry88/pseudo-interfaces/blob/master/Pseudo-interfaces_project_paper.pdf

## Dependencies:
- Python 2.7
- BioPython 1.65
- DSSP

## How to run
```python pinterface_finder.py [--sd 94.0] [--core] path_to_pdb num_of_results```

- path_to_pdb - relative path of the PDB file 
- num_of_results - amount of templates to output 
- sd - SD threshold, see paper  
- core - for PDBs containing more than one chain, thi will output only pseudo-interfaces appearing on one single chain.
