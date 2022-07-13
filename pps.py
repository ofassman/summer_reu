import abc_tree
import growtree as gt
import matplotlib.pyplot as plt

rates_arr = abc_tree.run_main(num_accept=3) # get array of inferred rates from ABC

# extract specific arrays of inferred rates
d_rate_arr = rates_arr[0]
r_rate_arr = rates_arr[1]
birth_s_arr = rates_arr[2]
death_s_arr = rates_arr[3]
sub_s_arr = rates_arr[4]
obs_tree = rates_arr[5]

n_sims = len(d_rate_arr) # number of trees that will be simulated from the rates that make up the posterior distributions

i = 0 
depth_mean_arr = [] # will hold the mean depths of the posterior simulated trees 
depth_var_arr = [] # will hold the variance of depths of the posterior simulated trees 

depth_mean_obs = gt.tree_depth_mean(obs_tree) # mean depth of the observed tree
depth_var_obs = gt.tree_depth_variance(obs_tree) # variance of depth of the observed tree

while i < n_sims: # generating posterior simulated trees and calculating statistics 
        # generating the posterior simulated tree (using rates from the posterior distributions)
    sim_tree = abc_tree.gen_tree_sims(d = d_rate_arr[i], r = r_rate_arr[i], birth_shape = birth_s_arr[i], death_shape = death_s_arr[i], sub_shape = sub_s_arr[i])[0]
    i += 1

    print(gt.getNewick(sim_tree))
    print()

    # calculating the summary statistics on the posterior simulated tree 
    depth_mean_arr.append(gt.tree_depth_mean(sim_tree))
    depth_var_arr.append(gt.tree_depth_variance(sim_tree))

# plotting the distributions for the 2 statistics calculated on the posterior simulated trees
# with the value of the statistic for the observed tree plotted as a point
plt.subplot(121)
plt.hist(depth_mean_arr, bins = 30)
plt.plot(depth_mean_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
plt.subplot(122)
plt.hist(depth_var_arr, bins = 30)
plt.plot(depth_var_obs, 1, marker = "o", markersize = 5) # plot observed statistic point
plt.show()





#print(depth_mean_arr)
#print(depth_mean_obs) 

print("----")
#print(depth_var_arr)
#print(depth_var_obs) 

print(gt.getNewick(obs_tree))
