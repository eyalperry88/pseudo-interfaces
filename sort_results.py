# -*- coding: utf-8 -*-
"""
Sorting results of MMKKM

@author: eyalp_home
"""

def sort_results_by_distance_score(results):
    return sorted(results.items(), key=lambda e : e[1][0])

def sort_results_by_interaction_score(results):
    return sorted(results.items(), key=lambda e : e[1][1], reverse=True)

def sort_results_by_combined_rank(results):
    by_dist = sort_results_by_distance_score(results)
    d = {}
    for rank in range(len(by_dist)):
        d[by_dist[rank][0]] = rank
    by_inter = sort_results_by_interaction_score(results)
    for rank in range(len(by_inter)):
        d[by_inter[rank][0]] += rank
    return sorted(results.items(), key=lambda e : d[e[0]])
