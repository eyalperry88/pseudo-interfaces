# -*- coding: utf-8 -*-
"""
Many to Many k-k match
Ver 1 - three to three

Created on Tue Sep 29 21:54:10 2015

@author: yot

"""

from scipy.misc import comb



# IN MMKKM
TOTAL_INTERACTIONS_THRESHOLD = 4
calc_weight = lambda x,y,z: x + y + z
# calc_weight = lambda x,y,z: 3*(min(x,y,z)+x+y+z)/4
combine_weights = lambda x, y, z: x + y + z
# combine_weights = lambda x, y, z: 3 * (sum([x, y, z]) * 2 - max(x, y, z)) / 5


def MMKKM(graph, th):
    """given a graph of distances between secondary structures,
    returns a dictionary with pairs of node ID triplets as keys
    and the sum of distances as value, for matches that are 
    less than THRESHOLD apart."""
    dic = dict()
    
    
#    there must be at least 6 nodes to get results
    if len(graph.nodes)<6:
        print 'you must have at least 6 nodes to get results.'
        print '%s has only %d nodes' % (graph.protein_name,len(graph.nodes))
        return dic
    
    #for every set of three points..  
    total = comb(len(graph.nodes), 3)
    global_count = 0
    print 'Running MMKKM...'
    for lst in k_sublists(graph.nodes, 3):
        global_count+=1       
        if (global_count%50000 == 0):
            print '%.1f%% Done' % ((100 * global_count) / total)
#        1. merge list of neighbors.
        a=lst[0]
        b=lst[1]
        c=lst[2]
        lst_tupel = tuple(sorted([a.node_id,b.node_id,c.node_id]))
#        print "a="+str(a.node_id)+"b="+str(b.node_id)+"c="+str(c.node_id)
        
        merged = merge_nodes_neighbors(a.neighbors,b.neighbors,c.neighbors)
        
        
#        2. rank the distances
        ranking =  sorted(merged,lambda x,y:int((x[1]<y[1])-(x[1]>y[1])))
#        print "ranking is: " + str(ranking)
                
#        3. choose only those smaller than the threshold
        cut_index = threshold_index(ranking, th)
#        if len(ranking)>0:
#            print "TH index is: "+str(cut_index)+" Ranking lst len is: "+str(len(ranking))
        ranking = ranking[cut_index:]
        if len(ranking) < 3:
            continue
        
        
#        4. find all combinations of triplets the give suumed weight<threshold
        for trio in k_sublists(ranking,3):
            summed_weight =  combine_weights(trio[0][1], trio[1][1], trio[2][1])
            total_interactions = trio[0][2] + trio[1][2] + trio[2][2]
            if (summed_weight<th):
                strio = tuple(sorted([trio[0][0] , trio[1][0] , trio[2][0]]))
                key = tuple(sorted([lst_tupel,strio]))
                if dic.has_key(key):
                    if (dic[key][0]>summed_weight):
                        dic[key] = (summed_weight, total_interactions)
                dic[key] = (summed_weight, total_interactions)
        
    print '100% Done'
    return dic
    


def merge_nodes_lists(lst1, lst2):
    """returns the intersection of two sorted lists, adding the values"""
    alist = []
    i=0
    j=0
    while i < len(lst1) and j < len(lst2):
        if lst1[i][0] == lst2[j][0]:
#            print "i="+str(i)+"\tj="+str(j)
            alist.append((lst1[i][0],lst1[i][1]+lst2[j][1]))
        if lst1[i][0] < lst2[j][0]:
            i=i+1
        else:
            j=j+1
    return alist

def merge_nodes_neighbors(lst1,lst2,lst3):
    """returns the intersection of three sorted lists, calcing the values"""
#   calc_weight=lambda x,y,z: x+y+z    or    lambda x,y,z: 3*(min(x,y,z)+x+y+z)/4
    alist = []
    i=j=k=0
    while i < len(lst1) and j < len(lst2) and k<len(lst3):
        if lst1[i][0] == lst2[j][0] == lst3[k][0]:
            total_interactions = lst1[i][2]+lst2[j][2]+lst3[k][2]
            if total_interactions>=TOTAL_INTERACTIONS_THRESHOLD:
                calculated_weight = calc_weight(lst1[i][1],lst2[j][1],lst3[k][1])
                alist.append((lst1[i][0],calculated_weight,total_interactions))

        a= lst1[i][0]
        b= lst2[j][0]
        c= lst3[k][0]
        if a< (1+min(b,c)):
            i+=1
        if b< (1+min(a,c)):
            j+=1
        if c< (1+min(a,b)):
            k+=1
       
    return alist
    


def threshold_index(lst,th):
    """Return the index where to find value th in lst, assuming lst is sorted.
    The return value i is such that all e in lst[:i] have e <= th, and all e in
    lst[i:] have e > th.  So if th already appears in the list, 
    lst.insert(th) will insert just after the rightmost x already there."""
    lo=0
    hi = len(lst)
    while lo < hi:
        mid = (lo+hi)//2
        if th > lst[mid][1]: hi = mid
        else: lo = mid+1
    return lo
    
    

# from http://code.activestate.com/recipes/500268-all-k-subsets-from-an-n-set/
def k_subsets_i(n, k):
    '''
    Yield each subset of size k from the set of intergers 0 .. n - 1
    n -- an integer > 0
    k -- an integer > 0
    '''
    # Validate args
    if n < 0:
        raise ValueError('n must be > 0, got n=%d' % n)
    if k < 0:
        raise ValueError('k must be > 0, got k=%d' % k)
    # check base cases
    if k == 0 or n < k:
        yield set()
    elif n == k:
        yield set(range(n))

    else:
        # Use recursive formula based on binomial coeffecients:
        # choose(n, k) = choose(n - 1, k - 1) + choose(n - 1, k)
        for s in k_subsets_i(n - 1, k - 1):
            s.add(n - 1)
            yield s
        for s in k_subsets_i(n - 1, k):
            yield s

def k_subsets(s, k):
    '''
    Yield all subsets of size k from set (or list) s
    s -- a set or list (any iterable will suffice)
    k -- an integer > 0
    '''
    s = list(s)
    n = len(s)
    for k_set in k_subsets_i(n, k):
        yield set([s[i] for i in k_set])

# from http://code.activestate.com/recipes/500268-all-k-subsets-from-an-n-set/
def k_sublists_i(n, k):
    '''
    Yield each subset of size k from the set of intergers 0 .. n - 1
    n -- an integer > 0
    k -- an integer > 0
    '''
    # Validate args
    if n < 0:
        raise ValueError('n must be > 0, got n=%d' % n)
    if k < 0:
        raise ValueError('k must be > 0, got k=%d' % k)
    # check base cases
    if k == 0 or n < k:
        yield list()
    elif n == k:
        yield list(range(n))

    else:
        # Use recursive formula based on binomial coeffecients:
        # choose(n, k) = choose(n - 1, k - 1) + choose(n - 1, k)
        for s in k_sublists_i(n - 1, k - 1):
            s.append(n - 1)
            yield s
        for s in k_sublists_i(n - 1, k):
            yield s

def k_sublists(s, k):
    '''
    Yield all subsets of size k from set (or list) s
    s -- a set or list (any iterable will suffice)
    k -- an integer > 0
    '''
    n = len(s)
    for k_list in k_sublists_i(n, k):
        yield list([s[i] for i in k_list])
