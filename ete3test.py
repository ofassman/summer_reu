import growtree

t = growtree.gen_tree(1, 0.01, 0.5, 1, 1, 1, 1, 100, goal_leaves = 10)
print(t)
print(growtree.tree_nleaf(t))
#print(growtree.getNewick(t))
growtree.outputNewick(t, "NWtree")
#growtree.print_seq()
print("done")