import elfi
import growtree
import numpy as np
import scipy
import ete3
import matplotlib.pyplot as plt

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

def tree_stat(tree_arr, summ_fn):
    """
    Applies function 'summ_fn()' to every element of 'tree_arr' and returns
    the array of results. 'summ_fn()' is a summary function that takes in a single
    tree and returns a numerical value (e.g. tree height, mean branch length,
    colless index).
    """
    res_arr = [] # array that will hold the summary statistic of the trees in 'tree_arr'
    for i in tree_arr: # for each tree in 'tree_arr'
        if(type(i) != ete3.coretype.tree.TreeNode): # if 'tree_arr' is an array of simulated trees
            res_arr.append(summ_fn(i[0])) # calculate the summary statistic of current tree, 'i'
        else: # 'tree_arr' is a one element array containing only the observed tree ('obs')
            res_arr.append(summ_fn(i.get_tree_root())) # calculate the summary statistic for 'obs' tree
            break
    return res_arr # return array of summary statistics

"""
Below are the set of summary statistic functions that can be passed 
into 'tree_stat()'. Each function requires an array of trees ('tree_arr') 
and returns an array where the elements are the summary statistic calculated 
on each tree in 'tree_arr'.
"""
def branch_sum_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_branch_sum)

def branch_mean_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_branch_mean)

def branch_median_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_branch_median)

def branch_variance_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_branch_variance)

def height_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_height)

def depth_mean_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_depth_mean)

def depth_median_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_depth_median)

def depth_variance_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_depth_variance)

def balance_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_balance)

def nleaves_stat(tree_arr):
   return tree_stat(tree_arr, growtree.tree_nleaf)

def root_colless_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_root_colless)

def sum_colless_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_sum_colless)

def mean_colless_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_mean_colless)

def median_colless_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_median_colless)

def variance_colless_stat(tree_arr):
    return tree_stat(tree_arr, growtree.tree_variance_colless)

"""
True parameters for birth and death rates (to be estimated by ABC).
"""
birth_true = .8
death_true = .5

"""
Prior distributions of rate parameters (uniform distributions).
"""
birth = elfi.Prior(scipy.stats.uniform, 0, 1)
death = elfi.Prior(scipy.stats.uniform, 0, 1)

obs = (gen_tree_sims(birth_true, death_true))[0] # observed tree (tree simulated with true rate parameters)

"""
'sim' is a simulator node with 'gen_tree_sims()' function, the prior distributions of rate 
parameters ('birth' and 'death'), and the obeserved tree ('obs') passed to it as arguments.
"""
sim = elfi.Simulator(elfi.tools.vectorize(gen_tree_sims), birth, death, observed = obs) 

"""
Below are summary nodes with each node having a unique tree summary statistic function
and all the simulated trees (includes the observed tree).
"""
summ_branch_sum = elfi.Summary(branch_sum_stat, sim) 
summ_branch_mean = elfi.Summary(branch_mean_stat, sim) 
summ_branch_median = elfi.Summary(branch_median_stat, sim) 
summ_branch_variance = elfi.Summary(branch_variance_stat, sim) 
summ_height = elfi.Summary(height_stat, sim) 
summ_depth_mean = elfi.Summary(depth_mean_stat, sim) 
summ_depth_median = elfi.Summary(depth_median_stat, sim) 
summ_depth_variance = elfi.Summary(depth_variance_stat, sim) 
summ_balance = elfi.Summary(balance_stat, sim) 
summ_nleaves = elfi.Summary(nleaves_stat, sim) 
summ_root_colless = elfi.Summary(root_colless_stat, sim) 
summ_colless_sum = elfi.Summary(sum_colless_stat, sim) 
summ_colless_mean = elfi.Summary(mean_colless_stat, sim) 
summ_colless_median = elfi.Summary(median_colless_stat, sim) 
summ_colless_variance = elfi.Summary(variance_colless_stat, sim) 

"""
Below are distance nodes that calculate the euclidian (squared distance): 
('summ_stat_1_sim' - 'summ_stat_1_obs') * 2 + ... + ('summ_stat_n_sim' - 'summ_stat_n_obs') * 2
for 'n' summary statistics (where 'n' is the number of statistics provided in the creation
of the distance node). 
"""

# 'dist_all_summ' is a distance node containing all the available tree summary statistics
dist_all_summ = elfi.Distance('euclidean', summ_branch_sum, summ_branch_mean, summ_branch_median, 
    summ_branch_variance, summ_height, summ_depth_mean, summ_depth_median, summ_depth_variance, 
    summ_balance, summ_nleaves, summ_root_colless, summ_colless_sum, summ_colless_mean, 
    summ_colless_median, summ_colless_variance)

# 'dist_mmv' is a distance node containing all summary statistics, but exluding tree summary 
# statistics which calculate sums (i.e. mean, median, and variance statistics, but not sum statistics)
dist_mmv = elfi.Distance('euclidean', summ_branch_mean, summ_branch_median, summ_branch_variance, 
    summ_height, summ_depth_mean, summ_depth_median, summ_depth_variance, summ_balance, 
    summ_nleaves,summ_colless_mean, summ_colless_median, summ_colless_variance)

# 'dist_mm' is a distance node containing all summary statistics, but exluding tree summary 
# statistics which calculate sums and variances (i.e. mean and median statistics, but not 
# sum and variance statistics)
dist_mm = elfi.Distance('euclidean', summ_branch_mean, summ_branch_median, summ_height,
    summ_depth_mean, summ_depth_median, summ_balance, summ_nleaves, summ_colless_mean, 
    summ_colless_median)

"""
'rej' is a rejection node used in inference with rejection sampling
using 'dist' values in order to reject. 'batch_size' defines how many 
simulations are performed in computation of 'dist'.
"""
rej = elfi.Rejection(dist_mm, batch_size = 1000)

N = 200 # number of accepted samples needed in 'result' in the inference with rejection sampling below

"""
Below is rejection using a threshold 'thresh'. All simulated trees generated
from rates drawn from the priors that have a generated distance below 'thresh'
will be accepted as samples. The simulator will generate as many trees as it 
takes to accept the specified number of trees ('N' trees).
"""
thresh = 0.5 # distance threshold
result_thresh = rej.sample(N, threshold = thresh) # COMMENT OUT THIS LINE IF QUANTILE METHOD IS USED

result_thresh.summary() # summary statistics from the inference with rejection sampling
result_thresh.plot_marginals() # plotting the marginal distributions of the birth and death rates for the accepted samples
plt.show()

"""
Below is rejection using quantiles. The quantile of trees size 'quant' 
with the smallest generated distances are accepted. The simulator will 
generate ('N' / 'quant') trees and accept 'N' of them.
"""
quant = 0.1 # quantile of accepted trees
result_quant = rej.sample(N, quantile = quant) # COMMENT OUT THIS LINE IF THRESHOLD METHOD IS USED

result_quant.summary() # summary statistics from the inference with rejection sampling
result_quant.plot_marginals() # plotting the marginal distributions of the birth and death rates for the accepted samples
plt.show()

"""
Note that it is not necessary to sample using both types of rejection described 
above (threshold and quantiles). One or the other is sufficient and using both 
will needlessly increase the runtime of the program. (If both types are used, 
the entire inference with rejection sampling process will occur twice.) 
For efficiency, choose one type of sampling per excecution of the file.
"""