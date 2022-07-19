import growtree as gt
import abc_tree as abc
import matplotlib.pyplot as plt


def div_rate_v_stats():
    d_rate = 0.000005
    d_arr = []
    b_sum_arr = []
    b_mean_arr = []
    b_median_arr = []
    b_variance_arr = []
    height_arr = []
    d_mean_arr = []
    d_median_arr = []
    d_variance_arr = []
    balance_arr = []
    nleaf_arr = []
    root_colless_arr = []
    sum_colless_arr = []
    mean_colless_arr = []
    median_colless_arr = []
    variance_colless_arr = []

    r_rate = 0.8
    birth_shape = 100
    death_shape = 100
    sub_shape = 100
    i = 0

    while(i < 20):
        d_rate += i * 0.0000005
        d_arr.append(d_rate)
        bd_rates = abc.calc_rates_bd(d_rate, r_rate)
        t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        while(t == None):
            t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        print(t)
        b_sum_arr.append(gt.tree_branch_sum(t))
        b_mean_arr.append(gt.tree_branch_mean(t))
        b_median_arr.append(gt.tree_branch_median(t))
        b_variance_arr.append(gt.tree_branch_variance(t))
        height_arr.append(gt.tree_height(t))
        d_mean_arr.append(gt.tree_depth_mean(t))
        d_median_arr.append(gt.tree_depth_median(t))
        d_variance_arr.append(gt.tree_depth_variance(t))
        balance_arr.append(gt.tree_balance(t))
        nleaf_arr.append(gt.tree_nleaf(t))
        root_colless_arr.append(gt.tree_root_colless(t))
        sum_colless_arr.append(gt.tree_sum_colless(t))
        mean_colless_arr.append(gt.tree_mean_colless(t))
        median_colless_arr.append(gt.tree_median_colless(t))
        variance_colless_arr.append(gt.tree_variance_colless(t))
        i += 1
    

    fig, axs = plt.subplots(3, 5)
    axs[0, 0].plot(d_arr, b_sum_arr, 'ro')
    axs[0, 0].set_title('Div v branch sum')
    axs[0, 1].plot(d_arr, b_mean_arr, 'ro')
    axs[0, 1].set_title('Branch mean')
    axs[0, 2].plot(d_arr, b_median_arr, 'ro')
    axs[0, 2].set_title('Branch median')
    axs[0, 3].plot(d_arr, b_variance_arr, 'ro')
    axs[0, 3].set_title('Branch variance')
    axs[0, 4].plot(d_arr, height_arr, 'ro')
    axs[0, 4].set_title('Height')
    axs[1, 0].plot(d_arr, d_mean_arr, 'ro')
    axs[1, 0].set_title('Depth mean')
    axs[1, 1].plot(d_arr, d_median_arr, 'ro')
    axs[1, 1].set_title('Depth median')
    axs[1, 2].plot(d_arr, d_variance_arr, 'ro')
    axs[1, 2].set_title('Depth variance')
    axs[1, 3].plot(d_arr, balance_arr, 'ro')
    axs[1, 3].set_title('Balance')
    axs[1, 4].plot(d_arr, nleaf_arr, 'ro')
    axs[1, 4].set_title('Nleaf')
    axs[2, 0].plot(d_arr, root_colless_arr, 'ro')
    axs[2, 0].set_title('Root colless')
    axs[2, 1].plot(d_arr, sum_colless_arr, 'ro')
    axs[2, 1].set_title('Sum colless')
    axs[2, 2].plot(d_arr, mean_colless_arr, 'ro')
    axs[2, 2].set_title('Mean colless')
    axs[2, 3].plot(d_arr, median_colless_arr, 'ro')
    axs[2, 3].set_title('Median colless')
    axs[2, 4].plot(d_arr, variance_colless_arr, 'ro')
    axs[2, 4].set_title('Variance colless')
    plt.show()

def turn_rate_v_stats():
    r_rate = 0.000005
    r_arr = []
    b_sum_arr = []
    b_mean_arr = []
    b_median_arr = []
    b_variance_arr = []
    height_arr = []
    d_mean_arr = []
    d_median_arr = []
    d_variance_arr = []
    balance_arr = []
    nleaf_arr = []
    root_colless_arr = []
    sum_colless_arr = []
    mean_colless_arr = []
    median_colless_arr = []
    variance_colless_arr = []

    d_rate = 0.00005
    birth_shape = 100
    death_shape = 100
    sub_shape = 100
    i = 0

    while(i < 20):
        r_rate += i * 0.005
        r_arr.append(r_rate)
        bd_rates = abc.calc_rates_bd(d_rate, r_rate)
        t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        while(t == None):
            t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        print(t)
        b_sum_arr.append(gt.tree_branch_sum(t))
        b_mean_arr.append(gt.tree_branch_mean(t))
        b_median_arr.append(gt.tree_branch_median(t))
        b_variance_arr.append(gt.tree_branch_variance(t))
        height_arr.append(gt.tree_height(t))
        d_mean_arr.append(gt.tree_depth_mean(t))
        d_median_arr.append(gt.tree_depth_median(t))
        d_variance_arr.append(gt.tree_depth_variance(t))
        balance_arr.append(gt.tree_balance(t))
        nleaf_arr.append(gt.tree_nleaf(t))
        root_colless_arr.append(gt.tree_root_colless(t))
        sum_colless_arr.append(gt.tree_sum_colless(t))
        mean_colless_arr.append(gt.tree_mean_colless(t))
        median_colless_arr.append(gt.tree_median_colless(t))
        variance_colless_arr.append(gt.tree_variance_colless(t))
        i += 1
    

    fig, axs = plt.subplots(3, 5)
    axs[0, 0].plot(r_arr, b_sum_arr, 'ro')
    axs[0, 0].set_title('Turn v branch sum')
    axs[0, 1].plot(r_arr, b_mean_arr, 'ro')
    axs[0, 1].set_title('Branch mean')
    axs[0, 2].plot(r_arr, b_median_arr, 'ro')
    axs[0, 2].set_title('Branch median')
    axs[0, 3].plot(r_arr, b_variance_arr, 'ro')
    axs[0, 3].set_title('Branch variance')
    axs[0, 4].plot(r_arr, height_arr, 'ro')
    axs[0, 4].set_title('Height')
    axs[1, 0].plot(r_arr, d_mean_arr, 'ro')
    axs[1, 0].set_title('Depth mean')
    axs[1, 1].plot(r_arr, d_median_arr, 'ro')
    axs[1, 1].set_title('Depth median')
    axs[1, 2].plot(r_arr, d_variance_arr, 'ro')
    axs[1, 2].set_title('Depth variance')
    axs[1, 3].plot(r_arr, balance_arr, 'ro')
    axs[1, 3].set_title('Balance')
    axs[1, 4].plot(r_arr, nleaf_arr, 'ro')
    axs[1, 4].set_title('Nleaf')
    axs[2, 0].plot(r_arr, root_colless_arr, 'ro')
    axs[2, 0].set_title('Root colless')
    axs[2, 1].plot(r_arr, sum_colless_arr, 'ro')
    axs[2, 1].set_title('Sum colless')
    axs[2, 2].plot(r_arr, mean_colless_arr, 'ro')
    axs[2, 2].set_title('Mean colless')
    axs[2, 3].plot(r_arr, median_colless_arr, 'ro')
    axs[2, 3].set_title('Median colless')
    axs[2, 4].plot(r_arr, variance_colless_arr, 'ro')
    axs[2, 4].set_title('Variance colless')
    plt.show()


turn_rate_v_stats()