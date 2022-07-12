import abc_tree
import growtree
import matplotlib.pyplot as plt

rates_arr = abc_tree.run_main()
obs_tree = rates_arr[0]
d_rate_arr = rates_arr[1]
r_rate_arr = rates_arr[2]
birth_s_arr = rates_arr[3]
death_s_arr = rates_arr[4]
sub_s_arr = rates_arr[5]

n_sims = len(d_rate_arr)

i = 0 
depth_mean_arr = []
depth_var_arr = []

depth_mean_obs = growtree.tree_depth_mean(obs_tree)
depth_var_obs = growtree.tree_depth_variance(obs_tree)

while i < n_sims:
    sim_tree = abc_tree.gen_tree_sims(d_rate_arr[i], r_rate_arr[i], birth_s_arr[i], death_s_arr[i], sub_s_arr[i])[0]
    i += 1 
    depth_mean_arr.append(growtree.tree_depth_mean(sim_tree))
    depth_var_arr.append(growtree.tree_depth_variance(sim_tree))
    
plt.subplot(121)
plt.hist(depth_mean_arr, bins = 30)
plt.plot(depth_mean_obs, 1, marker = "o", markersize = 5)
plt.subplot(122)
plt.hist(depth_var_arr, bins = 30)
plt.plot(depth_var_obs, 1, marker = "o", markersize = 5)
plt.show()
print(depth_mean_arr)
print(depth_mean_obs) 