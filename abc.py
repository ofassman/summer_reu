import elfi
import growtree
import matplotlib.pyplot as plt
import scipy.stats as ss
import numpy as np
import scipy

def gen_tree_sims(b, d, batch_size=1,random_state=None):
    arr = []
    random_state = random_state or np.random
    # should there be a loop here thatcalls gen_tree batch_size times?
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

obs = (gen_tree_sims(birth_true, death_true))[0]
Y = elfi.Simulator(elfi.tools.vectorize(gen_tree_sims),birth,death,observed=obs)

S1 = elfi.Summary(treestat, Y)

dist = elfi.Distance('euclidean', S1)

rej = elfi.Rejection(dist, batch_size=10000)

result = rej.sample(10, threshold=0.2)
result.summary()