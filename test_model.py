import matplotlib.pyplot as plt
from sqlalchemy import true
import growtree as gt
import abc_tree as abct
import statistics

i = 0 
N = 30
div_pointx = []
div_pointy = []
turn_pointx = []
turn_pointy = []
births_pointx = []
births_pointy = []
deaths_pointx = []
deaths_pointy = []
subs_pointx = []
subs_pointy = []

while i < N:
    i += 1
    res_arr = abct.run_main(isreal_obs = False, num_accept = 10)
    div_true = res_arr[6]
    turn_true = res_arr[7]
    births_true = res_arr[8]
    deaths_true = res_arr[9]
    subs_true = res_arr[10]

    div_infer_arr = res_arr[0]
    turn_infer_arr = res_arr[1]
    births_infer_arr = res_arr[2]
    deaths_infer_arr = res_arr[3]
    subs_infer_arr = res_arr[4]

    div_pointx.append(div_true)
    div_pointy.append(statistics.mean(div_infer_arr))

    turn_pointx.append(turn_true)
    turn_pointy.append(statistics.mean(turn_infer_arr))

    births_pointx.append(births_true)
    births_pointy.append(statistics.mean(births_infer_arr))

    deaths_pointx.append(deaths_true)
    deaths_pointy.append(statistics.mean(deaths_infer_arr))

    subs_pointx.append(subs_true)
    subs_pointy.append(statistics.mean(subs_infer_arr))

def plot_div_exp_v_true():
    plt.plot(div_pointx, div_pointy, 'ro')
    plt.show()

def plot_turn_exp_v_true():
    plt.plot(turn_pointx, turn_pointy, 'ro')
    plt.show()

def plot_births_exp_v_true():
    plt.plot(births_pointx, births_pointy, 'ro')
    plt.show()

def plot_deaths_exp_v_true():
    plt.plot(deaths_pointx, deaths_pointy, 'ro')
    plt.show()

def plot_subs_exp_v_true():
    plt.plot(subs_pointx, subs_pointy, 'ro')
    plt.show()

def calc_percent(true_arr, interval_arr):
    num_tests = 0
    num_fails = 0
    for i in true_arr:
        num_tests += 1
        if (i < interval_arr[0] or i > interval_arr[1]):
            num_fails += 1
    return (num_tests - num_fails) / num_tests

def plot_coverage(infer_arr, true_arr):
    total_interval = statistics.quantiles(infer_arr, n = 100)
    interval_50 = [] 
    interval_55 = [] 
    interval_60 = [] 
    interval_65 = [] 
    interval_70 = [] 
    interval_75 = [] 
    interval_80 = [] 
    interval_85 = [] 
    interval_90 = [] 
    interval_95 = [] 
    interval_50.append(total_interval[25])
    interval_50.append(total_interval[74])
    interval_55.append(total_interval[23])
    interval_55.append(total_interval[77])
    interval_60.append(total_interval[20])
    interval_60.append(total_interval[79])
    interval_65.append(total_interval[18])
    interval_65.append(total_interval[82])
    interval_70.append(total_interval[15])
    interval_70.append(total_interval[84])
    interval_75.append(total_interval[13])
    interval_75.append(total_interval[87])
    interval_80.append(total_interval[10])
    interval_80.append(total_interval[89])
    interval_85.append(total_interval[8])
    interval_85.append(total_interval[92])
    interval_90.append(total_interval[5])
    interval_90.append(total_interval[94])
    interval_95.append(total_interval[3])
    interval_95.append(total_interval[97])
   
    percent_arr = []
    percent_arr.append(calc_percent(true_arr, interval_50))
    percent_arr.append(calc_percent(true_arr, interval_55))
    percent_arr.append(calc_percent(true_arr, interval_60))
    percent_arr.append(calc_percent(true_arr, interval_65))
    percent_arr.append(calc_percent(true_arr, interval_70))
    percent_arr.append(calc_percent(true_arr, interval_75))
    percent_arr.append(calc_percent(true_arr, interval_80))
    percent_arr.append(calc_percent(true_arr, interval_85))
    percent_arr.append(calc_percent(true_arr, interval_90))
    percent_arr.append(calc_percent(true_arr, interval_95))

    plt.plot([50, 55, 60, 65, 70, 75, 80, 85, 90, 95], percent_arr, 'ro')
    plt.show()

def plot_div_coverage():
    plot_coverage(div_infer_arr, div_pointx)

def plot_turn_coverage():
    plot_coverage(turn_infer_arr, turn_pointx)

def plot_births_coverage():
    plot_coverage(births_infer_arr, births_pointx)

def plot_deaths_coverage():
    plot_coverage(deaths_infer_arr, deaths_pointx)

def plot_subs_coverage():
    plot_coverage(subs_infer_arr, subs_pointx)

plot_div_coverage()
print()