import matplotlib.pyplot as plt
import growtree as gt
import abc_tree as abct




d_true = abct.gen_param(abct.d)
r_true = abct.gen_param(abct.r)

rate_arr = abct.calc_rates_bd(d_true, r_true) # calculating the true birth and death parameters 
birth_true = rate_arr[0] # extracting birth rate
death_true = rate_arr[1] # extracting death rate

birth_s_true = abct.gen_param(abct.birth_s)
death_s_true = abct.gen_param(abct.death_s)
sub_s_true = abct.gen_param(abct.sub_s)








