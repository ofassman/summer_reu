from ete3 import Tree
import growtree

"""
t = Tree()
t.add_child(Tree())
print(t.is_leaf())
print(t)
"""
t = growtree.growtree(0.5,0.3,0.5,20,1,1,1,2)
print(t)
print(growtree.getNewick(t))
growtree.outputNewick(t,"NWtree")
print("done")