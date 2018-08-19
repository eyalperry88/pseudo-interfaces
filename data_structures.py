# -*- coding: utf-8 -*-

""" Classes for data structures """

class SSNode:
    'Node in the SS graph'

    def __init__(self, node_id, ss_type, res_length, start, chain, residues):
        self.node_id = node_id
        self.ss_type = ss_type
        self.res_length = res_length
        self.start = start
        self.chain = chain
        self.residues = residues
        self.neighbors = []
        self.extended_neighbors = []
    
    def __len__(self):
        return self.res_length

    def addNeighbor(self, node_id, dist, inter_count):
        for i in range(len(self.neighbors)):
            if self.neighbors[i][0] > node_id:
                self.neighbors.insert(i, (node_id, dist, inter_count))
                return
        self.neighbors.append((node_id, dist, inter_count))
        
    def addExtendedNeighbor(self, node_id, dist, inter_count):
        for i in range(len(self.extended_neighbors)):
            if self.extended_neighbors[i][0] > node_id:
                self.extended_neighbors.insert(i, (node_id, dist, inter_count))
                return
        self.extended_neighbors.append((node_id, dist, inter_count))

    def removeNeighbor(self, node_id):
        for i in range(len(self.neighbors)):
            if self.neighbors[i][0] == node_id:
                del self.neighbors[i]
                return
                
    def removeExtendedNeighbor(self, node_id):
        for i in range(len(self.extended_neighbors)):
            if self.extended_neighbors[i][0] == node_id:
                del self.extended_neighbors[i]
                return

class SSGraph:
    """ Secondary Structure Graph """

    def __init__(self, protein_name=None, structure=None):
        self.protein_name = protein_name
        self.nodes = []
        self.residue2node = {}
        self.structure = structure

    def addNode(self, ss_type=None, res_length=-1, start=-1, chain='', residues=[]):
        node_id = len(self.nodes)
        node = SSNode(node_id, ss_type, res_length, start, chain, residues)
        self.nodes.append(node)
        for pos in range(start, start+res_length):
            self.residue2node[(chain, pos)] = node_id
        return node_id

    def addEdge(self, node_id1, node_id2, dist, inter_count):
        node1 = self.nodes[node_id1]
        node1.addNeighbor(node_id2, dist, inter_count)
        node2 = self.nodes[node_id2]
        node2.addNeighbor(node_id1, dist, inter_count)
        
    def addExtendedEdge(self, node_id1, node_id2, dist, inter_count):
        node1 = self.nodes[node_id1]
        node1.addExtendedNeighbor(node_id2, dist, inter_count)
        node2 = self.nodes[node_id2]
        node2.addExtendedNeighbor(node_id1, dist, inter_count)

    def removeEdge(self, node_id1, node_id2):
        node1 = self.nodes[node_id1]
        node1.removeNeighbor(node_id2)
        node2 = self.nodes[node_id2]
        node2.removeNeighbor(node_id1)

    def removeExtendedEdge(self, node_id1, node_id2):
        node1 = self.nodes[node_id1]
        node1.removeExtendedNeighbor(node_id2)
        node2 = self.nodes[node_id2]
        node2.removeExtendedNeighbor(node_id1)