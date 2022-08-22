import growtree

t = growtree.gen_tree(0.5, 0.2, 0.5, 1, 1, 1, 1, 100)
print(t)
print(growtree.getNewick(t))
growtree.outputNewick(t, "NWtree")
growtree.print_seq()
print("done")