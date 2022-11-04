import time
import growtree
import matplotlib.pyplot as plt


num_leaves = []
times = []

for i in range(1,91):
    start_time = time.perf_counter()
    t = growtree.gen_tree(i*.1, 0.01, 0.5, 1, 1, 1, 1, 100)
    end_time = time.perf_counter()
    num_leaves.append(growtree.tree_nleaf(t))
    times.append(end_time-start_time)
    print(i)

print(num_leaves)
print(times)
plt.scatter(num_leaves,times)
plt.show()



   


