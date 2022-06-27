import elfi
import growtree
import numpy as np
import scipy
import ete3

def gen_tree_sims(b = 1, d = 1, random_state = None):
    """
    Returns a simulated phylogenetic tree (using growtree.gen_tree()) with the 
    initial birth rate = 'b', initial death rate = 'd', and initial substitution 
    rate = 1. The tree is returned in a one element array in order to be compatible
    with the ELFI package. 'random_state' is not currently used, but is included
    as a parameter since some ELFI functions pass a value for 'random_state'
    into this function.
    """
    arr = []
    random_state = random_state or np.random # this value is not currently used
    arr.append(growtree.gen_tree(b, d, 1, 100, 1, 1, 1, 0, 100)) # simulate tree and place in 1 element array
    return arr

def tree_height_stat(tree_arr):
    """
    Returns the array of tree heights of the trees in 'tree_arr'. 
    """
    res_arr = [] # array that will hold the tree heights of the trees in 'tree_arr'
    for i in tree_arr: # for each tree in 'tree_arr'
        if(type(i) != ete3.coretype.tree.TreeNode): # if 'tree_arr' is an array of simulated trees
            res_arr.append(growtree.tree_height(i[0])) # calculate height of current tree, 'i'
        else: # 'tree_arr' is a one element array containing only the observed tree ('obs')
            res_arr.append(growtree.tree_height(i.get_tree_root())) # calculate the height of 'obs'
            break
    return res_arr # return array of heights

"""
True parameters for birth and death rates (to be estimated by ABC).
"""
birth_true = 1.3
death_true = 1.5

"""
Prior distributions of rate parameters (uniform distributions).
"""
birth = elfi.Prior(scipy.stats.uniform, 0, 5)
death = elfi.Prior(scipy.stats.uniform, 0, 5)

obs = (gen_tree_sims(birth_true, death_true))[0] # observed tree (tree simulated with true rate parameters)

"""
'sim' is a simulator node with 'gen_tree_sims()' function, the prior distributions of rate 
parameters ('birth' and 'death'), and the obeserved tree ('obs') passed to it as arguments.
"""
sim = elfi.Simulator(elfi.tools.vectorize(gen_tree_sims), birth, death, observed = obs) 

"""
'sum_stat_height' is a summary node with the 'tree_height_stat()' function and the 
simulated trees (includes the observed tree).
"""
sum_stat_height = elfi.Summary(tree_height_stat, sim) 

"""
'dist' is a distance node that calculates the euclidian (squared distance): 
('sum_stat_height_sim' - 'sum_stat_height_obs') * 2
"""
dist = elfi.Distance('euclidean', sum_stat_height)

"""
'rej' is a rejection node used in inference with rejection sampling
using 'dist' values in order to reject. 'batch_size' defines how many 
simulations are performed in computation of 'dist'.
"""
rej = elfi.Rejection(dist, batch_size = 1000)

N = 50 # number of accepted samples needed in 'result' in the inference with rejection sampling below

"""
Below is rejection using a threshold 'thresh'. All simulated trees generated
from rates drawn from the priors that have a generated distance below 'thresh'
will be accepted as samples. The simulator will generate as many trees as it 
takes to accept the specified number of trees ('N' trees).
"""
thresh = 0.5 # distance threshold
result_thresh = rej.sample(N, threshold = thresh)

"""
Below is rejection using quantiles. The quantile of trees size 'quantile' 
with the smallest generated distances are accepted. The simulator will 
generate ('N' / 'quantile') trees and accept 'N' of them.
"""
result_quant = rej.sample(N, quantile = 0.01)

result_thresh.summary() # summary statistics from the inference with rejection sampling