from ete3 import Tree
import growtree

"""
t = Tree()
t.add_child(Tree())
print(t.is_leaf())
print(t)
"""

t = growtree.gen_tree(0.5,0.2,0.5,20,1,1,1,0,100)
print(t)
print(growtree.getNewick(t))
growtree.outputNewick(t,"NWtree")
#growtree.print_seq()
print("done")