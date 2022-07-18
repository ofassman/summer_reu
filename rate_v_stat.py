import growtree as gt
import abc_tree as abc
import matplotlib.pyplot as plt

d_rate = 0.000005
r_rate = 0.0005
birth_shape = 1
death_shape = 1
sub_shape = 1
i = 0

d_arr = []
r_arr = []
birth_arr = []
death_arr = []
sub_arr = []
stat_arr = []

while(i < 20):
    d_rate += i * 0.0000005
    d_arr.append(d_rate)
    r_rate += i * 0.005
    r_arr.append(r_rate)
    birth_shape += i * 2
    birth_arr.append(birth_shape)
    death_shape = death_shape + i * 2
    death_arr.append(death_shape)
    sub_shape += i * 2
    sub_arr.append(sub_shape)
    bd_rates = abc.calc_rates_bd(d_rate, r_rate)
    t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
    stat_arr.append(gt.tree_mean_colless(t))
    i += 1
    

fig, axs = plt.subplots(2, 3)
axs[0, 0].plot(d_arr, stat_arr, 'ro')
axs[0, 0].set_title('Diversification rate v Branch sum')
axs[0, 1].plot(r_arr, stat_arr, 'ro')
axs[0, 1].set_title('Turnover rate')
axs[1, 0].plot(birth_arr, stat_arr, 'ro')
axs[1, 0].set_title('Birth shape')
axs[1, 1].plot(death_arr, stat_arr, 'ro')
axs[1, 1].set_title('Death shape')
axs[1, 2].plot(sub_arr, stat_arr, 'ro')
axs[1, 2].set_title('Sub shape')
plt.show()

plt.plot(birth_arr, stat_arr, 'ro')
plt.show()




