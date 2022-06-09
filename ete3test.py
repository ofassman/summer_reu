from ete3 import Tree
import growtree

"""
t = Tree()
t.add_child(Tree())
print(t.is_leaf())
print(t)
"""
t = growtree.growtree(0.5,0.5,4,20,1,1,1)
print(t)
print("done")