from ete3 import Tree
import growtree

"""
t = Tree()
t.add_child(Tree())
print(t.is_leaf())
print(t)
"""
t = growtree.growtree(0.5,0.3,0.5,20,1,1,1)
print(t)
newick_str = growtree.getNewick(t)
print(newick_str)
growtree.outputNewick(t,"NWtree")
print("done")