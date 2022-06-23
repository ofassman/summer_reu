import elfi
import growtree
import matplotlib.pyplot as plt
import scipy.stats as ss
import numpy as np
import scipy



def gen_tree_sims(b, d, n_sims=1):
    arr = []
    for i in range(0, n_sims):
        arr.append(growtree.gen_tree(b,d,1,100,1,1,1,0,100))
    return arr
def treestat(tree):
    return growtree.tree_height(tree)

# true parameters
birth_true = 1.3
death_true = 2.5
sub_true = 1

birth = elfi.Prior(scipy.stats.uniform, 0, 5)
death = elfi.Prior(scipy.stats.uniform, 0, 5)

#y_obs = gen_tree_sims(birth_true, death_true, sub_true)
Y = elfi.Simulator(gen_tree_sims, birth_true, death_true, 1)
X = elfi.Simulator(gen_tree_sims, birth, death, 1)

S1 = elfi.Summary(treestat, Y)
S2 = elfi.Summary(treestat, X)
d = elfi.Distance('euclidean', S1, S2)
rej = elfi.Rejection(d, batch_size=10000)


"""
N = 1000

vis = dict(xlim=[-2,2], ylim=[-1,1])

# You can give the sample method a `vis` keyword to see an animation how the prior transforms towards the
# posterior with a decreasing threshold.
result = rej.sample(N, quantile=0.01, vis=vis)
result.summary()

"""