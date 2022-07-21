import growtree as gt
import abc_tree as abc
import matplotlib.pyplot as plt
import elfi
import scipy
import math


d = elfi.Prior(scipy.stats.expon, 0, .000047) # prior distribution for diversification
r = elfi.Prior(scipy.stats.uniform, 0, 0.999999999999999999) # prior distribution for turnover
birth_s = elfi.Prior(scipy.stats.expon, 0, 100) # prior distribution for birth distribution shape
death_s = elfi.Prior(scipy.stats.expon, 0, 100) # prior distribution for death distribution shape
sub_s = elfi.Prior(scipy.stats.expon, 0, 100) # prior distribution for substitution distribution shape

def zero_log(i):
    if(i == 0):
        return 0
    return math.log10(i)

def div_rate_v_stats(use_prior = True, N = 100):
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

    r_rate = abc.gen_param(r)
    birth_shape = abc.gen_param(birth_s)
    death_shape = abc.gen_param(death_s)
    sub_shape = abc.gen_param(sub_s)
    i = 0

    while(i < N):
        if(use_prior):
            d_rate = abc.gen_param(d)
        else:
            d_rate += i * 0.0000005
        d_arr.append(d_rate)
        bd_rates = abc.calc_rates_bd(d_rate, r_rate)
        t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        n_leaf = gt.tree_nleaf(t)
        while(n_leaf < 1):
            t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
            n_leaf = gt.tree_nleaf(t)
        #print(t)
        b_sum_arr.append(gt.tree_branch_sum(t))
        b_mean_arr.append(gt.tree_branch_mean(t))
        b_median_arr.append(gt.tree_branch_median(t))
        b_variance_arr.append(gt.tree_branch_variance(t))
        height_arr.append(gt.tree_height(t))
        d_mean_arr.append(gt.tree_depth_mean(t))
        d_median_arr.append(gt.tree_depth_median(t))
        d_variance_arr.append(gt.tree_depth_variance(t))
        balance_arr.append(gt.tree_balance(t))
        nleaf_arr.append(n_leaf)
        root_colless_arr.append(gt.tree_root_colless(t))
        sum_colless_arr.append(gt.tree_sum_colless(t))
        mean_colless_arr.append(gt.tree_mean_colless(t))
        median_colless_arr.append(gt.tree_median_colless(t))
        variance_colless_arr.append(gt.tree_variance_colless(t))
        i += 1
    
    """
    b_sum_arr = list(map(zero_log, b_sum_arr))
    b_mean_arr = list(map(zero_log, b_mean_arr))
    b_median_arr = list(map(zero_log, b_mean_arr))
    b_variance_arr = list(map(zero_log, b_variance_arr))
    height_arr = list(map(zero_log, height_arr))
    d_mean_arr = list(map(zero_log, d_mean_arr))
    d_median_arr = list(map(zero_log, d_mean_arr))
    d_variance_arr = list(map(zero_log, d_variance_arr))
    """

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

    axs[0, 0].set_ylim(bottom = 0)
    axs[0, 1].set_ylim(bottom = 0)
    axs[0, 2].set_ylim(bottom = 0)
    axs[0, 3].set_ylim(bottom = 0)
    axs[0, 4].set_ylim(bottom = 0)
    axs[1, 0].set_ylim(bottom = 0)
    axs[1, 1].set_ylim(bottom = 0)
    axs[1, 2].set_ylim(bottom = 0)
    axs[1, 3].set_ylim(bottom = 0)
    axs[1, 4].set_ylim(bottom = 0)
    axs[2, 0].set_ylim(bottom = 0)
    axs[2, 1].set_ylim(bottom = 0)
    axs[2, 2].set_ylim(bottom = 0)
    axs[2, 3].set_ylim(bottom = 0)
    axs[2, 4].set_ylim(bottom = 0)

    plt.show()

def turn_rate_v_stats(use_prior = True, N = 100):
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

    d_rate = abc.gen_param(d)
    birth_shape = abc.gen_param(birth_s)
    death_shape = abc.gen_param(death_s)
    sub_shape = abc.gen_param(sub_s)
    i = 0

    while(i < N):
        if(use_prior):
            r_rate = abc.gen_param(r)
        else: # if use_prior is false, N must be <= 20
            r_rate += i * 0.005
        r_arr.append(r_rate)
        bd_rates = abc.calc_rates_bd(d_rate, r_rate)
        t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        n_leaf = gt.tree_nleaf(t)
        while(n_leaf == 0):
            t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
            n_leaf = gt.tree_nleaf(t)
        #print(t)
        b_sum_arr.append(gt.tree_branch_sum(t))
        b_mean_arr.append(gt.tree_branch_mean(t))
        b_median_arr.append(gt.tree_branch_median(t))
        b_variance_arr.append(gt.tree_branch_variance(t))
        height_arr.append(gt.tree_height(t))
        d_mean_arr.append(gt.tree_depth_mean(t))
        d_median_arr.append(gt.tree_depth_median(t))
        d_variance_arr.append(gt.tree_depth_variance(t))
        balance_arr.append(gt.tree_balance(t))
        nleaf_arr.append(n_leaf)
        root_colless_arr.append(gt.tree_root_colless(t))
        sum_colless_arr.append(gt.tree_sum_colless(t))
        mean_colless_arr.append(gt.tree_mean_colless(t))
        median_colless_arr.append(gt.tree_median_colless(t))
        variance_colless_arr.append(gt.tree_variance_colless(t))
        i += 1

    """
    b_sum_arr = list(map(zero_log, b_sum_arr))
    root_colless_arr = list(map(zero_log, root_colless_arr))
    sum_colless_arr = list(map(zero_log, sum_colless_arr))
    variance_colless_arr = list(map(zero_log, variance_colless_arr))
    """

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
    
    axs[0, 0].set_ylim(bottom = 0)
    axs[0, 1].set_ylim(bottom = 0)
    axs[0, 2].set_ylim(bottom = 0)
    axs[0, 3].set_ylim(bottom = 0)
    axs[0, 4].set_ylim(bottom = 0)
    axs[1, 0].set_ylim(bottom = 0)
    axs[1, 1].set_ylim(bottom = 0)
    axs[1, 2].set_ylim(bottom = 0)
    axs[1, 3].set_ylim(bottom = 0)
    axs[1, 4].set_ylim(bottom = 0)
    axs[2, 0].set_ylim(bottom = 0)
    axs[2, 1].set_ylim(bottom = 0)
    axs[2, 2].set_ylim(bottom = 0)
    axs[2, 3].set_ylim(bottom = 0)
    axs[2, 4].set_ylim(bottom = 0)
    plt.show()

def birth_shape_v_stats(use_prior = True, N = 100):
    birth_shape = 1
    birth_s_arr = []
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

    d_rate = abc.gen_param(d)
    r_rate = abc.gen_param(r)
    death_shape = abc.gen_param(death_s)
    sub_shape = abc.gen_param(sub_s)
    i = 0

    while(i < N):
        if(use_prior):
            birth_shape = abc.gen_param(birth_s)
        else:
            birth_shape += i * 2
        birth_s_arr.append(birth_shape)
        bd_rates = abc.calc_rates_bd(d_rate, r_rate)
        t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        n_leaf = gt.tree_nleaf(t)
        while(n_leaf == 0):
            t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
            n_leaf = gt.tree_nleaf(t)
        #print(t)
        b_sum_arr.append(gt.tree_branch_sum(t))
        b_mean_arr.append(gt.tree_branch_mean(t))
        b_median_arr.append(gt.tree_branch_median(t))
        b_variance_arr.append(gt.tree_branch_variance(t))
        height_arr.append(gt.tree_height(t))
        d_mean_arr.append(gt.tree_depth_mean(t))
        d_median_arr.append(gt.tree_depth_median(t))
        d_variance_arr.append(gt.tree_depth_variance(t))
        balance_arr.append(gt.tree_balance(t))
        nleaf_arr.append(n_leaf)
        root_colless_arr.append(gt.tree_root_colless(t))
        sum_colless_arr.append(gt.tree_sum_colless(t))
        mean_colless_arr.append(gt.tree_mean_colless(t))
        median_colless_arr.append(gt.tree_median_colless(t))
        variance_colless_arr.append(gt.tree_variance_colless(t))
        i += 1
    fig, axs = plt.subplots(3, 5)
    axs[0, 0].plot(birth_s_arr, b_sum_arr, 'ro')
    axs[0, 0].set_title('Birth shape v branch sum')
    axs[0, 1].plot(birth_s_arr, b_mean_arr, 'ro')
    axs[0, 1].set_title('Branch mean')
    axs[0, 2].plot(birth_s_arr, b_median_arr, 'ro')
    axs[0, 2].set_title('Branch median')
    axs[0, 3].plot(birth_s_arr, b_variance_arr, 'ro')
    axs[0, 3].set_title('Branch variance')
    axs[0, 4].plot(birth_s_arr, height_arr, 'ro')
    axs[0, 4].set_title('Height')
    axs[1, 0].plot(birth_s_arr, d_mean_arr, 'ro')
    axs[1, 0].set_title('Depth mean')
    axs[1, 1].plot(birth_s_arr, d_median_arr, 'ro')
    axs[1, 1].set_title('Depth median')
    axs[1, 2].plot(birth_s_arr, d_variance_arr, 'ro')
    axs[1, 2].set_title('Depth variance')
    axs[1, 3].plot(birth_s_arr, balance_arr, 'ro')
    axs[1, 3].set_title('Balance')
    axs[1, 4].plot(birth_s_arr, nleaf_arr, 'ro')
    axs[1, 4].set_title('Nleaf')
    axs[2, 0].plot(birth_s_arr, root_colless_arr, 'ro')
    axs[2, 0].set_title('Root colless')
    axs[2, 1].plot(birth_s_arr, sum_colless_arr, 'ro')
    axs[2, 1].set_title('Sum colless')
    axs[2, 2].plot(birth_s_arr, mean_colless_arr, 'ro')
    axs[2, 2].set_title('Mean colless')
    axs[2, 3].plot(birth_s_arr, median_colless_arr, 'ro')
    axs[2, 3].set_title('Median colless')
    axs[2, 4].plot(birth_s_arr, variance_colless_arr, 'ro')
    axs[2, 4].set_title('Variance colless')

    
    axs[0, 0].set_ylim(bottom = 0)
    axs[0, 1].set_ylim(bottom = 0)
    axs[0, 2].set_ylim(bottom = 0)
    axs[0, 3].set_ylim(bottom = 0)
    axs[0, 4].set_ylim(bottom = 0)
    axs[1, 0].set_ylim(bottom = 0)
    axs[1, 1].set_ylim(bottom = 0)
    axs[1, 2].set_ylim(bottom = 0)
    axs[1, 3].set_ylim(bottom = 0)
    axs[1, 4].set_ylim(bottom = 0)
    axs[2, 0].set_ylim(bottom = 0)
    axs[2, 1].set_ylim(bottom = 0)
    axs[2, 2].set_ylim(bottom = 0)
    axs[2, 3].set_ylim(bottom = 0)
    axs[2, 4].set_ylim(bottom = 0)

    plt.show()

def death_shape_v_stats(use_prior = True, N = 100):
    death_shape = 1
    death_s_arr = []
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

    d_rate = abc.gen_param(d)
    r_rate = abc.gen_param(r)
    birth_shape = abc.gen_param(birth_s)
    sub_shape = abc.gen_param(sub_s)
    i = 0

    while(i < N):
        if(use_prior):
            death_shape = abc.gen_param(death_s)
        else:
            death_shape += i * 2
        death_s_arr.append(death_shape)
        bd_rates = abc.calc_rates_bd(d_rate, r_rate)
        t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        n_leaf = gt.tree_nleaf(t)
        while(n_leaf == 0):
            t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
            n_leaf = gt.tree_nleaf(t)
        #print(t)
        b_sum_arr.append(gt.tree_branch_sum(t))
        b_mean_arr.append(gt.tree_branch_mean(t))
        b_median_arr.append(gt.tree_branch_median(t))
        b_variance_arr.append(gt.tree_branch_variance(t))
        height_arr.append(gt.tree_height(t))
        d_mean_arr.append(gt.tree_depth_mean(t))
        d_median_arr.append(gt.tree_depth_median(t))
        d_variance_arr.append(gt.tree_depth_variance(t))
        balance_arr.append(gt.tree_balance(t))
        nleaf_arr.append(n_leaf)
        root_colless_arr.append(gt.tree_root_colless(t))
        sum_colless_arr.append(gt.tree_sum_colless(t))
        mean_colless_arr.append(gt.tree_mean_colless(t))
        median_colless_arr.append(gt.tree_median_colless(t))
        variance_colless_arr.append(gt.tree_variance_colless(t))
        i += 1
    
    fig, axs = plt.subplots(3, 5)
    axs[0, 0].plot(death_s_arr, b_sum_arr, 'ro')
    axs[0, 0].set_title('Death shape v branch sum')
    axs[0, 1].plot(death_s_arr, b_mean_arr, 'ro')
    axs[0, 1].set_title('Branch mean')
    axs[0, 2].plot(death_s_arr, b_median_arr, 'ro')
    axs[0, 2].set_title('Branch median')
    axs[0, 3].plot(death_s_arr, b_variance_arr, 'ro')
    axs[0, 3].set_title('Branch variance')
    axs[0, 4].plot(death_s_arr, height_arr, 'ro')
    axs[0, 4].set_title('Height')
    axs[1, 0].plot(death_s_arr, d_mean_arr, 'ro')
    axs[1, 0].set_title('Depth mean')
    axs[1, 1].plot(death_s_arr, d_median_arr, 'ro')
    axs[1, 1].set_title('Depth median')
    axs[1, 2].plot(death_s_arr, d_variance_arr, 'ro')
    axs[1, 2].set_title('Depth variance')
    axs[1, 3].plot(death_s_arr, balance_arr, 'ro')
    axs[1, 3].set_title('Balance')
    axs[1, 4].plot(death_s_arr, nleaf_arr, 'ro')
    axs[1, 4].set_title('Nleaf')
    axs[2, 0].plot(death_s_arr, root_colless_arr, 'ro')
    axs[2, 0].set_title('Root colless')
    axs[2, 1].plot(death_s_arr, sum_colless_arr, 'ro')
    axs[2, 1].set_title('Sum colless')
    axs[2, 2].plot(death_s_arr, mean_colless_arr, 'ro')
    axs[2, 2].set_title('Mean colless')
    axs[2, 3].plot(death_s_arr, median_colless_arr, 'ro')
    axs[2, 3].set_title('Median colless')
    axs[2, 4].plot(death_s_arr, variance_colless_arr, 'ro')
    axs[2, 4].set_title('Variance colless')

    
    axs[0, 0].set_ylim(bottom = 0)
    axs[0, 1].set_ylim(bottom = 0)
    axs[0, 2].set_ylim(bottom = 0)
    axs[0, 3].set_ylim(bottom = 0)
    axs[0, 4].set_ylim(bottom = 0)
    axs[1, 0].set_ylim(bottom = 0)
    axs[1, 1].set_ylim(bottom = 0)
    axs[1, 2].set_ylim(bottom = 0)
    axs[1, 3].set_ylim(bottom = 0)
    axs[1, 4].set_ylim(bottom = 0)
    axs[2, 0].set_ylim(bottom = 0)
    axs[2, 1].set_ylim(bottom = 0)
    axs[2, 2].set_ylim(bottom = 0)
    axs[2, 3].set_ylim(bottom = 0)
    axs[2, 4].set_ylim(bottom = 0)

    plt.show()

def sub_shape_v_stats(use_prior = True, N = 100):
    sub_shape = 1
    sub_s_arr = []
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

    d_rate = abc.gen_param(d)
    r_rate = abc.gen_param(r)
    birth_shape = abc.gen_param(birth_s)
    death_shape = abc.gen_param(death_s)
    i = 0

    while(i < N):
        if(use_prior):
            sub_shape = abc.gen_param(sub_s)
        else:
            sub_shape += i * 2
        sub_s_arr.append(sub_shape)
        bd_rates = abc.calc_rates_bd(d_rate, r_rate)
        t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
        n_leaf = gt.tree_nleaf(t)
        while(n_leaf == 0):
            t = gt.gen_tree(bd_rates[0], bd_rates[1], 1, 50000, birth_shape, death_shape, sub_shape, 1, 100)
            n_leaf = gt.tree_nleaf(t)
        #print(t)
        b_sum_arr.append(gt.tree_branch_sum(t))
        b_mean_arr.append(gt.tree_branch_mean(t))
        b_median_arr.append(gt.tree_branch_median(t))
        b_variance_arr.append(gt.tree_branch_variance(t))
        height_arr.append(gt.tree_height(t))
        d_mean_arr.append(gt.tree_depth_mean(t))
        d_median_arr.append(gt.tree_depth_median(t))
        d_variance_arr.append(gt.tree_depth_variance(t))
        balance_arr.append(gt.tree_balance(t))
        nleaf_arr.append(n_leaf)
        root_colless_arr.append(gt.tree_root_colless(t))
        sum_colless_arr.append(gt.tree_sum_colless(t))
        mean_colless_arr.append(gt.tree_mean_colless(t))
        median_colless_arr.append(gt.tree_median_colless(t))
        variance_colless_arr.append(gt.tree_variance_colless(t))
        i += 1
    
    fig, axs = plt.subplots(3, 5)
    axs[0, 0].plot(sub_s_arr, b_sum_arr, 'ro')
    axs[0, 0].set_title('Sub shape v branch sum')
    axs[0, 1].plot(sub_s_arr, b_mean_arr, 'ro')
    axs[0, 1].set_title('Branch mean')
    axs[0, 2].plot(sub_s_arr, b_median_arr, 'ro')
    axs[0, 2].set_title('Branch median')
    axs[0, 3].plot(sub_s_arr, b_variance_arr, 'ro')
    axs[0, 3].set_title('Branch variance')
    axs[0, 4].plot(sub_s_arr, height_arr, 'ro')
    axs[0, 4].set_title('Height')
    axs[1, 0].plot(sub_s_arr, d_mean_arr, 'ro')
    axs[1, 0].set_title('Depth mean')
    axs[1, 1].plot(sub_s_arr, d_median_arr, 'ro')
    axs[1, 1].set_title('Depth median')
    axs[1, 2].plot(sub_s_arr, d_variance_arr, 'ro')
    axs[1, 2].set_title('Depth variance')
    axs[1, 3].plot(sub_s_arr, balance_arr, 'ro')
    axs[1, 3].set_title('Balance')
    axs[1, 4].plot(sub_s_arr, nleaf_arr, 'ro')
    axs[1, 4].set_title('Nleaf')
    axs[2, 0].plot(sub_s_arr, root_colless_arr, 'ro')
    axs[2, 0].set_title('Root colless')
    axs[2, 1].plot(sub_s_arr, sum_colless_arr, 'ro')
    axs[2, 1].set_title('Sum colless')
    axs[2, 2].plot(sub_s_arr, mean_colless_arr, 'ro')
    axs[2, 2].set_title('Mean colless')
    axs[2, 3].plot(sub_s_arr, median_colless_arr, 'ro')
    axs[2, 3].set_title('Median colless')
    axs[2, 4].plot(sub_s_arr, variance_colless_arr, 'ro')
    axs[2, 4].set_title('Variance colless')

    axs[0, 0].set_ylim(bottom = 0)
    axs[0, 1].set_ylim(bottom = 0)
    axs[0, 2].set_ylim(bottom = 0)
    axs[0, 3].set_ylim(bottom = 0)
    axs[0, 4].set_ylim(bottom = 0)
    axs[1, 0].set_ylim(bottom = 0)
    axs[1, 1].set_ylim(bottom = 0)
    axs[1, 2].set_ylim(bottom = 0)
    axs[1, 3].set_ylim(bottom = 0)
    axs[1, 4].set_ylim(bottom = 0)
    axs[2, 0].set_ylim(bottom = 0)
    axs[2, 1].set_ylim(bottom = 0)
    axs[2, 2].set_ylim(bottom = 0)
    axs[2, 3].set_ylim(bottom = 0)
    axs[2, 4].set_ylim(bottom = 0)

    plt.show()

div_rate_v_stats()
turn_rate_v_stats()
birth_shape_v_stats()
death_shape_v_stats()
sub_shape_v_stats()