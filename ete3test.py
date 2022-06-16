from ete3 import Tree
import growtree

"""
t = Tree()
t.add_child(Tree())
print(t.is_leaf())
print(t)
"""

t = growtree.gen_tree(0.5,0.3,0.5,20,1,1,1,0,100)
print(t)
print(growtree.getNewick(t))
growtree.outputNewick(t,"NWtree")
print("\nsum: ")
print(growtree.tree_branch_sum(t))
print(growtree.tree_branch_mean(t))
print(growtree.tree_branch_variance(t))

print("done")


