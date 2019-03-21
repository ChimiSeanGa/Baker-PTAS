import sys
import networkx as nx
from networkx import Graph
from networkx.algorithms import approximation, union
from itertools import combinations, chain
import math


class Node:
   def __init__(self, vertices, edges, left, right):
      self.vertices = vertices
      self.edges = edges
      self.left = left
      self.right = right

   def is_leaf():
      return self.left == None and self.right == None


def powerset(iterable):
   "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
   s = list(iterable)
   return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def draw_tree_decomposition(T, name):
   mapping = {}
   for v in list(T.nodes):
      mapping[v] = " ".join(list(v))
   T = nx.relabel_nodes(T, mapping)

   A = nx.nx_agraph.to_agraph(T)
   A.layout()
   A.draw(name + ".png")


def hierarchical_tree(G, T):
   tree_list = dict(nx.bfs_successors(T, list(T.nodes)[0]))
   return hierarchical_tree_recurse(G, tree_list, list(tree_list.keys())[0], [])


def hierarchical_tree_recurse(G, tree_list, node, edge_list):
   edges = []

   for pair in list(combinations(list(node), 2)):
      if G.has_edge(pair[0], pair[1]) and set(pair) not in edge_list:
         edges.append(pair)
         edge_list.append(set(pair))

   left = None
   right = None

   if node in tree_list:
      children = tree_list[node]
      if len(children) > 2:
         left = hierarchical_tree_recurse(G, tree_list, children[0], edge_list)
         tree_list[node].remove(children[0])
         right = hierarchical_tree_recurse(G, tree_list, node, edge_list)
      elif len(children) == 2:
         left = hierarchical_tree_recurse(G, tree_list, children[0], edge_list)
         right = hierarchical_tree_recurse(G, tree_list, children[1], edge_list)
      elif len(children) == 1:
         left = hierarchical_tree_recurse(G, tree_list, children[0], edge_list)
   
   return Node(node, edges, left, right)


def print_bin_tree(tree, n=0):
   if tree == None:
      return
   print("-" * n + str(list(tree.vertices)) + ": " + str(tree.edges))
   print_bin_tree(tree.left, n + 1)
   print_bin_tree(tree.right, n + 1)


def dynamic_prog(B, S, k):
   if B == None or S == set() or k == 0 or k < len(S):
      return 0

   results = [0]
   C1 = B.left
   C2 = B.right

   S1_list = list(powerset(C1.vertices)) if C1 != None else [set()]
   S2_list = list(powerset(C2.vertices)) if C2 != None else [set()]

   if len(S1_list) == 0 and len(S2_list) == 0:
      results.append(len([e for e in B.edges if e[0] in S and e[1] in S]))

   for S1 in S1_list:
      for S2 in S2_list:
         S1 = set(S1)
         S2 = set(S2)
         if not (S1 & B.vertices <= S) or not (S2 & B.vertices <= S):
            continue
         if C1 != None and not S & C1.vertices <= S1:
            continue
         if C2 != None and not S & C2.vertices <= S2:
            continue

         bound = k - len(S - (S1 | S2)) + len(S1 & S2)

         for l1 in range(0, bound + 1):
            l2 = bound - l1
            if l1 < len(S1) or l2 < len(S2):
               continue

            dp1 = dynamic_prog(C1, S1, l1)
            dp2 = dynamic_prog(C2, S2, l2)
            w = len([e for e in B.edges if e[0] in S and e[1] in S])
            results.append(dp1 + dp2 + w)

   return max(results)


def baker_levels(G):
   succ = dict(nx.bfs_successors(G, list(G.nodes)[0]))
   levels = [[list(G.nodes)[0]]]
   done = False
   n = 1

   while not done:
      level = []
      for v in levels[n - 1]:
         if v in succ:
            level += succ[v]

      if level == []:
         done = True
      else:
         levels.append(level)
         n += 1

   return levels


def get_subgraphs(G, levels, b):
   subgraphs = []
   for i in range(0, b):
      for j in range(0, len(levels) // b + 1):
         start = j * b + i
         end = (j + 1) * b + i

         if end >= len(levels):
            end = len(levels)
         if start >= end:
            continue

         vertices = []
         for s in range(start, end):
            vertices += levels[s]

         subgraphs.append(G.subgraph(vertices))

   return subgraphs


def get_tree_decomps(subgraphs):
   trees = []
   for G in subgraphs:
      (w, T) = approximation.treewidth_min_degree(G)
      trees.append(T)

   bin_trees = []
   for i in range(len(subgraphs)):
      G = subgraphs[i]
      T = trees[i]
      bin_trees.append(hierarchical_tree(G, T))

   return bin_trees


def densest_k_subgraph_PTAS(G, k, epsilon):
   levels = baker_levels(G)
   subgraphs = get_subgraphs(G, levels, math.ceil(1 / epsilon))
   trees = get_tree_decomps(subgraphs)
   results = []

   # for T in trees:
   #    S_list = list(powerset(T.vertices))
   #    for S in S_list:
   #       results.append(dynamic_prog(T, set(S), k))

   for T in trees:
      traversal = []
      bin_tree_traversal(T, traversal)
      for B in traversal:
         S_list = list(powerset(B.vertices))
         for S in S_list:
            results.append(dynamic_prog(B, set(S), k))

   return max(results)


def bin_tree_traversal(tree, traversal):
   if tree == None:
      return

   traversal.append(tree)
   bin_tree_traversal(tree.left, traversal)
   bin_tree_traversal(tree.right, traversal)


def main():
   if len(sys.argv) != 4:
      print("Usage: dks_ptas.py [graph file] [k value] [epsilon value]")
      exit(1)
   
   G = nx.read_adjlist(sys.argv[1])
   k = int(sys.argv[2])
   epsilon = float(sys.argv[3])

   print("k = " + str(k) + ", ε = " + str(epsilon))
   print("Treewidth: " + str(approximation.treewidth_min_degree(G)[0]))
   print("Density: " + str(densest_k_subgraph_PTAS(G, k, epsilon)))


if __name__ == '__main__':
   main()