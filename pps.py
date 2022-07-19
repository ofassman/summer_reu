import abc_tree
import growtree as gt
import matplotlib.pyplot as plt

from summer_reu.growtree import tree_branch_sum, tree_depth_variance, tree_height, tree_mean_colless, tree_sum_colless

rates_arr = abc_tree.run_main() # get array of inferred rates from ABC

# extract specific arrays of inferred rates
d_rate_arr = rates_arr[0]
r_rate_arr = rates_arr[1]
birth_s_arr = rates_arr[2]
death_s_arr = rates_arr[3]
sub_s_arr = rates_arr[4]
obs_tree = rates_arr[5]

n_sims = len(d_rate_arr) # number of trees that will be simulated from the rates that make up the posterior distributions

i = 0 
branch_mean_arr = [] # will hold the mean branch lengths of the posterior simulated trees 
branch_median_arr = [] # will hold the median branch lengths of the posterior simulated trees 
branch_var_arr = [] # will hold the variance of the branch lengths of the posterior simulated trees 
branch_sum_arr = [] # will hold the sum of the branch lengths of the posterior simulated trees 
height_arr = [] # will hold the heights of the posterior simulated trees
depth_var_arr = [] # will hold the variance of depths of the posterior simulated trees 
sum_colless_arr = [] # will hold the sum of the colless indices of the posterior simulated trees 
mean_colless_arr = [] # will hold the mean of the colless indices of the posterior simulated trees 


branch_mean_obs = gt.tree_branch_mean(obs_tree)
branch_median_obs =  gt.tree_branch_median(obs_tree)
branch_var_obs = gt.tree_branch_variance(obs_tree)
branch_sum_obs = tree_branch_sum(obs_tree)
height_obs = tree_height(obs_tree)
depth_var_obs = tree_depth_variance(obs_tree)
sum_colless_obs = tree_sum_colless(obs_tree)
mean_colless_obs = tree_mean_colless(obs_tree)

while i < n_sims: # generating posterior simulated trees and calculating statistics 
        # generating the posterior simulated tree (using rates from the posterior distributions)
    sim_tree = abc_tree.gen_tree_sims(d = d_rate_arr[i], r = r_rate_arr[i], birth_shape = birth_s_arr[i], death_shape = death_s_arr[i], sub_shape = sub_s_arr[i])[0]
    i += 1

    # calculating the summary statistics on the posterior simulated tree 
    branch_mean_arr.append(gt.tree_branch_mean(sim_tree))
    branch_median_arr.append(gt.tree_branch_median(sim_tree))
    branch_var_arr.append(gt.tree_branch_variance(sim_tree))
    branch_sum_arr.append(gt.tree_branch_sum(sim_tree))
    height_arr.append(gt.tree_height(sim_tree))
    depth_var_arr.append(gt.tree_depth_variance(sim_tree))
    sum_colless_arr.append(gt.tree_sum_colless(sim_tree))
    mean_colless_arr.append(gt.tree_mean_colless(sim_tree))

# plotting the distributions for the 2 statistics calculated on the posterior simulated trees
# with the value of the statistic for the observed tree plotted as a point

fig, axs = plt.subplots(2, 4)
axs[0, 0].hist(branch_mean_arr, bins = 50)
axs[0, 0].plot(branch_mean_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[0, 0].ylabel('Rate frequency')
axs[0, 0].xlabel('Mean branch lengths of posterior simulated trees')
#axs[0, 0].title('Distribution of mean branch lengths for posterior simulated trees compared to the observed tree')
axs[0, 1].hist(branch_median_arr, bins = 50)
axs[0, 1].plot(branch_median_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[0, 1].xlabel('Median branch lengths of posterior simulated trees')
axs[0, 2].hist(branch_var_arr, bins = 50)
axs[0, 2].plot(branch_var_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[0, 2].xlabel('Variance of branch lengths of posterior simulated trees')
axs[0, 3].hist(branch_sum_arr, bins = 50)
axs[0, 3].plot(branch_sum_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[0, 3].xlabel('Sum of branch lengths of posterior simulated trees')
axs[1, 0].hist(height_arr, bins = 50)
axs[1, 0].plot(height_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[1, 0].xlabel('Height of posterior simulated trees')
axs[1, 1].hist(depth_var_arr, bins = 50)
axs[1, 1].plot(depth_var_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[1, 1].xlabel('Variance of depth of posterior simulated trees')
axs[1, 2].hist(sum_colless_arr, bins = 50)
axs[1, 2].plot(sum_colless_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[1, 2].xlabel('Sum of colless indices of posterior simulated trees')
axs[1, 3].hist(mean_colless_arr, bins = 50)
axs[1, 3].plot(mean_colless_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
axs[1, 3].xlabel('Mean colless indices of posterior simulated trees')

plt.show()