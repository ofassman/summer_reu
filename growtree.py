from ete3 import Tree
import random
import numpy
import statistics

def gen_sequence(length,off_lim = None):
    seq = ""
    if(off_lim == None):
        while(length > 0):
            rnum = random.randint(0, 3)
            if(rnum == 0):
                seq += "A"
            elif(rnum == 1):
                seq += "T"
            elif(rnum == 2):
                seq += "G"
            else:
                seq += "C"
            length -= 1
    else:
        if(off_lim == "A"):
            rnum = random.randint(0, 2)
            if(rnum == 0):
                seq += "T"
            elif(rnum == 1):
                seq += "G"
            elif(rnum == 2):
                seq += "C"
        elif(off_lim == "T"):
            rnum = random.randint(0, 2)
            if(rnum == 0):
                seq += "A"
            elif(rnum == 1):
                seq += "G"
            elif(rnum == 2):
                seq += "C"
        elif(off_lim == "G"):
            rnum = random.randint(0, 2)
            if(rnum == 0):
                seq += "A"
            elif(rnum == 1):
                seq += "T"
            elif(rnum == 2):
                seq += "C"
        else:
            rnum = random.randint(0, 2)
            if(rnum == 0):
                seq += "A"
            elif(rnum == 1):
                seq += "T"
            elif(rnum == 2):
                seq += "G"
    return seq

def gen_event(b,d,s):
    """
    Randomly generates an event based on weighted birth, death, and substitution rates.
    """
    rng = random.Random()
    weights = [b,d,s]
    rnd = rng.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            break
    if(i == 0): 
        return "birth"
    elif(i == 1): 
        return "death"
    return "sub"
    
def gen_rate(mean,shape):
    """
    Samples a new rate based on a gamma distribution given the mean rate and the shape of the distribution. 
    """
    scale_calc = mean/shape
    return numpy.random.gamma(shape, scale=scale_calc, size=None)

def growtree(seq,b,d,s,max_time,shape_b,shape_d,shape_s,branch_info):
    """
    Returns a birth-death tree. All rates (birth, death, and substitution) may change upon a substitution.
    'b', 'd', and 's' are the initial values of the birth, death, and substitution rates (respectively).
    'max_time' is the total amount of time that can be used to generate events to construct the tree
    (if this time is exceeded, the tree returns). Thus the tree stops growing when either all lineages
    go extinct or 'max_time' is exceeded. Only extant lineages are present in the final tree. If there
    are no extant lineages, 'None' will be returned by the function. 'branch_info' is an argument to specify 
    what information is attached to branches. If 'branch_info' is 0, the branch length is a variable of the 
    time for that lineage. If 'branch_info' is 1, the branch length is a variable of the number of substitutions 
    that occurred in that lineage. If 'branch_info' is 2, the branch length is a variable of the expected number 
    of substitutions for that lineage. 'shape_b', 'shape_d', and 'shape_s' are the shapes of the respective 
    gamma distributions from which each rate is sampled from upon a substitution.
    """
    rng = random.Random()
    # initializing the tree and branch length
    t = Tree()
    t.name = seq
    t.dist = 0
    while(True):
        # finding the wait time to any event (b, d, or s) based on rates
        curr_t = 0
        rate_any_event = b + d + s
        wait_t = rng.expovariate(rate_any_event)
        curr_t += wait_t
        if(curr_t <= max_time): # if wait time does not exceed max_time
            # if branch length is a variable of time, add 'wait_time' onto this lineage's branch length
            if(branch_info == 0): 
                t.dist += wait_t
            # if branch length is a variable of expected number of substitutions, add this expected number onto this lineage's branch length
            if(branch_info == 2): 
                t.dist += s * wait_t 
            # calculating weighted rates
            b_weighted = b/rate_any_event
            d_weighted = d/rate_any_event
            s_weighted = s/rate_any_event
            event = gen_event(b_weighted, d_weighted, s_weighted) # generate event based on weighted rates
            if(event == "birth"): # recursively call fn for children, same rates but different max_time
                # max_time for children should be the time remaining (max_time - curr_time) divided by 2 (since 2 children)
                c1 = growtree(seq,b,d,s,(max_time-curr_t)/2,shape_b,shape_d,shape_s,branch_info)
                c2 = growtree(seq,b,d,s,(max_time-curr_t)/2,shape_b,shape_d,shape_s,branch_info)  
                if(c1 == None and c2 == None): # both children are extinct so lineage is extinct (return None)
                    return None
                # children are only concatenated onto the parent tree if they are extant (non-None)
                if(c1 != None): 
                    t.add_child(c1)
                if(c2 != None):
                    t.add_child(c2)
                return t
            elif(event == "sub"): # change current rates based on sampling from a gamma distribution and continue to next event
                # mean of gamma distribution is current rate
                b = gen_rate(b,shape_b)
                d = gen_rate(d,shape_d)
                s = gen_rate(s,shape_s)
                sub_site = random.randint(0, len(seq)-1)
                old_letter = seq[sub_site]
                sub_letter = gen_sequence(1, off_lim=old_letter)
                new_seq = ""
                for i in range(0, sub_site, 1):
                    new_seq += seq[i : i + 1]
                new_seq += sub_letter
                for j in range(sub_site + 1, len(seq), 1):
                    new_seq += seq[j : j + 1]
                seq = new_seq
                # if branch length is a variable of number of substitutions, increase lineage's branch length by 1
                if(branch_info == 1):
                    t.dist += 1
            else: # event is death so return None (lineage goes extinct)
                return None
        else: # wait time exceeded max_time so return tree (lineage ends from timeout)
            return t

def gen_tree(b,d,s,max_time,shape_b,shape_d,shape_s,branch_info,seq_length):
    seq = gen_sequence(seq_length)
    t = growtree(seq,b,d,s,max_time,shape_b,shape_d,shape_s,branch_info)
    return t

def getNewick(t):
    """
    Returns a tree in Newick tree format.
    """
    if (t != None):
        return t.write()
    return ";"

def outputNewick(t,name):
    """
    Writes a tree, 't', in Newick tree format into a file. 'name' specifies the 
    file's name in which the tree is written into. If the tree is empty (i.e. if 
    't' is 'None') no output file is created.
    """
    if (t != None):
        t.write(outfile=name + ".nw")
    else:
        print("Empty tree, no output file created.")

def tree_branch_lst(t,arr):
    if(t == None): # empty tree
        return []
    arr.append(t.dist)
    num_c = len(t.children)  
    if(num_c == 1): # tree with 1 child
        tree_branch_lst(t.children[0], arr) 
    elif(num_c == 2): # tree with 2 children
        tree_branch_lst(t.children[0], arr) 
        tree_branch_lst(t.children[1], arr) 
    return arr

def tree_branch_sum(t):
    """
    Returns the sum of the distances of all the branches in the tree. 
    """
    branch_arr = tree_branch_lst(t,[])
    if(branch_arr == []):
        return 0
    return sum(branch_arr)

def tree_branch_median(t):
    branch_arr = tree_branch_lst(t,[])
    if(branch_arr == []):
        return 0
    return statistics.median(branch_arr)

def tree_branch_mean(t):
    branch_arr = tree_branch_lst(t,[])
    if(branch_arr == []):
        return 0
    return statistics.mean(branch_arr)

def tree_branch_variance(t):
    branch_arr = tree_branch_lst(t,[])
    if(branch_arr == []):
        return 0
    return statistics.variance(branch_arr)

def tree_height(t): # aka max depth
    if t == None:
        return 0 
    left_h = 0
    right_h = 0
    num_c = len(t.children)  
    if(num_c == 1): # tree with 1 child
        left_h = tree_height(t.children[0]) 
    elif(num_c == 2): # tree with 2 children
        left_h = tree_height(t.children[0]) 
        right_h = tree_height(t.children[1]) 
    return max(left_h, right_h) + t.dist

def tree_root_dist(t):
    if t == None or t.up == None:
        return 0 
    return tree_root_dist(t.up) + t.dist

def tree_depth_lst(t,arr):
    if(t == None):
        return []
    if(t.is_leaf()):
        arr.append(tree_root_dist(t))
    else:
        num_c = len(t.children)  
        if(num_c == 1): # tree with 1 child
            tree_depth_lst(t.children[0],arr) 
        elif(num_c == 2): # tree with 2 children
            tree_depth_lst(t.children[0],arr) 
            tree_depth_lst(t.children[1],arr) 
    return arr

def tree_mean_depth(t):
    depth_arr = tree_depth_lst(t,[])
    if(depth_arr == []):
        return 0
    return statistics.mean(depth_arr)

def tree_internal_height_lst(t,arr):
    if(t == None):
        return []
    if(not(t.is_leaf()) and not(t.is_root())):
        arr.append(1/tree_root_dist(t))
    num_c = len(t.children)  
    if(num_c == 1): # tree with 1 child
        tree_internal_height_lst(t.children[0],arr) 
    elif(num_c == 2): # tree with 2 children
        tree_internal_height_lst(t.children[0],arr) 
        tree_internal_height_lst(t.children[1],arr) 
    return arr

def tree_balance(t):
    height_arr = tree_internal_height_lst(t,[])
    if(height_arr == []):
        return 0
    return sum(height_arr)
