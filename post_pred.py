import abc_tree
import growtree as gt
import matplotlib.pyplot as plt

rates_arr = abc_tree.run_main()
d_rate_arr = rates_arr[0]
r_rate_arr = rates_arr[1]
birth_s_arr = rates_arr[2]
death_s_arr = rates_arr[3]
sub_s_arr = rates_arr[4]
obs_tree = rates_arr[5]

n_sims = len(d_rate_arr)

i = 0 
depth_mean_arr = []
depth_var_arr = []

depth_mean_obs = gt.tree_depth_mean(obs_tree)
depth_var_obs = gt.tree_depth_variance(obs_tree)

while i < n_sims:
    sim_tree = abc_tree.gen_tree_sims(d = d_rate_arr[i], r = r_rate_arr[i], birth_shape = birth_s_arr[i], death_shape = death_s_arr[i], sub_shape = sub_s_arr[i])[0]
    print(gt.getNewick(sim_tree))
    print()
    i += 1 
    depth_mean_arr.append(gt.tree_depth_mean(sim_tree))
    depth_var_arr.append(gt.tree_depth_variance(sim_tree))
    
plt.subplot(121)
plt.hist(depth_mean_arr, bins = 30)
#plt.plot(depth_mean_obs, 1, marker = "o", markersize = 5)
plt.subplot(122)
plt.hist(depth_var_arr, bins = 30)
#plt.plot(depth_var_obs, 1, marker = "o", markersize = 5)
plt.show()
#print(depth_mean_arr)
#print(depth_mean_obs) 

print("----")
#print(depth_var_arr)
#print(depth_var_obs) 

print(gt.getNewick(obs_tree))
