import elfi
import growtree
import matplotlib.pyplot as plt
import scipy.stats as ss
import numpy as np
import random

"""
# Define the simulator, the summary and the observed data
def simulator(b, d, s, batch_size=1, random_state=None):
    return growtree.gen_tree(b,d,s,100,1,1,1,0,100)
# Implementation comes here. Return ‘batch_size‘
# simulations wrapped to a NumPy array.
def summary(tree):
    return growtree.tree_height(tree)
# true parameters
birth_true = 1.3
death_true = 2.5
sub_true = 1
y = simulator(birth_true, death_true, sub_true)
birth = random.uniform(0, 5)
death = random.uniform(0, 5)
sub = random.uniform(0, 5)


SIM = elfi.Simulator(simulator, birth, death, sub, observed=y)
S1 = elfi.Summary(summary, SIM)
S2 = elfi.Summary(summary, SIM)
d = elfi.Distance('euclidean', S1, S2)
# Run the rejection sampler
rej = elfi.Rejection(d, batch_size=10000)
result = rej.sample(1000, threshold=0.1)

"""


birth = random.uniform(0, 5)
death = random.uniform(0, 5)
sub = random.uniform(0, 5)

def gen_tree_sims(b, d, s):
    t = growtree.gen_tree(b,d,s,100,1,1,1,0,100)
    # generate stats for t?
    return t

# true parameters
birth_true = 1.3
death_true = 2.5
sub_true = 1

y_obs = gen_tree_sims(birth_true, death_true, sub_true)

Y = elfi.Simulator(gen_tree_sims, birth, death, sub, observed=y_obs)

def treestat(tree):
    return growtree.tree_height(tree)

S1 = elfi.Summary(treestat, Y)
S2 = elfi.Summary(treestat, Y) 
# Finish the model with the final node that calculates the squared distance (S1_sim-S1_obs)**2 + (S2_sim-S2_obs)**2
d = elfi.Distance('euclidean', S1, S2)
rej = elfi.Rejection(d, batch_size=10000)
N = 1000

vis = dict(xlim=[-2,2], ylim=[-1,1])

# You can give the sample method a `vis` keyword to see an animation how the prior transforms towards the
# posterior with a decreasing threshold.
result = rej.sample(N, quantile=0.01, vis=vis)
result.summary()
