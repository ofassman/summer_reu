import elfi
import growtree
import matplotlib.pyplot as plt
import scipy.stats as ss
import numpy as np
import scipy
import statistics
import ete3

def gen_tree_sims(b=1,d=1,batch_size=1,random_state=None):
    arr = []
    random_state = random_state or np.random
    # should there be a loop here thatcalls gen_tree batch_size times?
    arr.append(growtree.gen_tree(b,d,1,100,1,1,1,0,100))
    return arr

def treestat(tree_arr):
    #print("this should be sim tree array")
    #print(tree_arr)
    #print("done")
    res_arr = []
    for i in tree_arr:
        #print(i)
        if(type(i) != ete3.coretype.tree.TreeNode): # if arr of sim trees
            #print(i[0])
            res_arr.append(growtree.tree_height(i[0]))
        else: # only obs tree
            res_arr.append(growtree.tree_height(i.get_tree_root()))
            break
    #print(res_arr)
    return res_arr

# true parameters
birth_true = 2.5
death_true = 0.3
sub_true = 1

birth = elfi.Prior(scipy.stats.uniform, 0, 5)
death = elfi.Prior(scipy.stats.uniform, 0, 5)

obs = (gen_tree_sims(birth_true, death_true))[0] # observed tree

Y = elfi.Simulator(elfi.tools.vectorize(gen_tree_sims),birth,death,observed=obs)



S1 = elfi.Summary(treestat, Y)

dist = elfi.Distance('euclidean', S1)


rej = elfi.Rejection(dist, batch_size=10)

result = rej.sample(1, threshold=0.2)


result.summary()
