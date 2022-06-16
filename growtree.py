from ete3 import Tree
import random
import numpy
import statistics

def gen_sequence(length,off_lim = None):
    """
    Randomly generates a 'length' long genetic sequence of bases. 'off_lim' is by default 'None', but can be used
    to specify a base that is not allowed to appear in the sequence (e.g. upon a substitution, the base that was 
    previously in the substitution site is not a valid newly substituted base).
    """
    seq = ""
    if(off_lim == None): # no restrictions on bases
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
        if(off_lim == "A"): # sequence cannot include 'A'
            rnum = random.randint(0, 2)
            if(rnum == 0):
                seq += "T"
            elif(rnum == 1):
                seq += "G"
            elif(rnum == 2):
                seq += "C"
        elif(off_lim == "T"): # sequence cannot include 'T'
            rnum = random.randint(0, 2)
            if(rnum == 0):
                seq += "A"
            elif(rnum == 1):
                seq += "G"
            elif(rnum == 2):
                seq += "C"
        elif(off_lim == "G"): # sequence cannot include 'G'
            rnum = random.randint(0, 2)
            if(rnum == 0):
                seq += "A"
            elif(rnum == 1):
                seq += "T"
            elif(rnum == 2):
                seq += "C"
        else:
            rnum = random.randint(0, 2) # sequence cannot include 'C'
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
    Returns a birth-death tree. Used as a helper function for 'gen_tree()' that produces the birth-death tree.
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
    gamma distributions from which each rate is sampled from upon a substitution. 'seq_length' specifies the 
    length of the genetic sequence for the cells (the root will have a randomly generated sequence of length 
    'seq_length' and subsequent lineages will carry on this sequence, with modifications upon a substitution).
    """
    seq = gen_sequence(seq_length) # generate random genetic sequence for root cell 
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

#################### TREE STATISTICS ########################################

def __tree_branch_lst(t,arr):
    """
    Returns an array of branch lengths. Private helper function used for calculating summary
    stats regarding branches.
    """
    if(t == None): # empty tree
        return []
    arr.append(t.dist) 
    num_c = len(t.children)  
    if(num_c == 1): # tree with 1 child
        __tree_branch_lst(t.children[0], arr) 
    elif(num_c == 2): # tree with 2 children
        __tree_branch_lst(t.children[0], arr) 
        __tree_branch_lst(t.children[1], arr) 
    return arr

def tree_branch_sum(t):
    """
    Returns the sum of the distances of all the branches in the tree. 
    """
    branch_arr = __tree_branch_lst(t,[]) # get array of branch lengths
    if(branch_arr == []):
        return 0
    return sum(branch_arr)

def tree_branch_mean(t):
    """
    Returns the mean of the distances of all the branches in the tree. 
    """
    branch_arr = __tree_branch_lst(t,[]) # get array of branch lengths
    if(branch_arr == []):
        return 0
    return statistics.mean(branch_arr)

def tree_branch_median(t):
    """
    Returns the median of the distances of all the branches in the tree. 
    """
    branch_arr = __tree_branch_lst(t,[]) # get array of branch lengths
    if(branch_arr == []):
        return 0
    return statistics.median(branch_arr)

def tree_branch_variance(t):
    """
    Returns the variance of the distances of all the branches in the tree. 
    """
    branch_arr = __tree_branch_lst(t,[]) # get array of branch lengths
    if(branch_arr == []):
        return 0
    return statistics.variance(branch_arr)

def tree_height(t): 
    """
    Returns the height (maximum depth) of the tree. 
    """
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

def __tree_root_dist(node):
    """
    Returns the distance from a node to the root. Private helper function used for calculating 
    summary stats regarding tree depth.
    """
    if node == None or node.up == None:
        return 0 
    return __tree_root_dist(node.up) + node.dist

def __tree_depth_lst(node,arr):
    """
    Returns an array of leaf depths. Private helper function used for calculating summary
    stats regarding tree depth.
    """
    if(node == None):
        return []
    if(node.is_leaf()): # if node is leaf, add it's depth to 'arr'
        arr.append(__tree_root_dist(node))
    else: # node is not leaf so recurse to find leaves
        num_c = len(node.children)  
        if(num_c == 1): # tree with 1 child
            __tree_depth_lst(node.children[0],arr) 
        elif(num_c == 2): # tree with 2 children
            __tree_depth_lst(node.children[0],arr) 
            __tree_depth_lst(node.children[1],arr) 
    return arr

def tree_depth_mean(t):
    """
    Returns the mean of leaf depths in the tree. 
    """
    depth_arr = __tree_depth_lst(t,[]) # get array of leaf depths
    if(depth_arr == []):
        return 0
    return statistics.mean(depth_arr)

def tree_depth_median(t):
    """
    Returns the median of leaf depths in the tree. 
    """
    depth_arr = __tree_depth_lst(t,[]) # get array of leaf depths
    if(depth_arr == []):
        return 0
    return statistics.median(depth_arr)

def tree_depth_variance(t):
    """
    Returns the variance of leaf depths in the tree. 
    """
    depth_arr = __tree_depth_lst(t,[]) # get array of leaf depths
    if(depth_arr == []):
        return 0
    return statistics.variance(depth_arr)

def __tree_internal_height_lst(node,arr):
    """
    Returns an array of the reciprocal of the heights (maximum depths) 
    of subtrees of t rooted at internal nodes of t (not including the root). 
    Private helper function used for calculating summary stats regarding 
    tree balance.
    """
    if(node == None):
        return []
    if(not(node.is_leaf()) and not(node.is_root())): # if node is internal (not root or leaf)
        arr.append(1/tree_height(node)) # add reciprocal of height of subtree rooted at 'node' to 'arr'
    # must find all internal nodes to add reciprocal of heights to 'arr'
    num_c = len(node.children)  
    if(num_c == 1): # tree with 1 child
        __tree_internal_height_lst(node.children[0],arr) 
    elif(num_c == 2): # tree with 2 children
        __tree_internal_height_lst(node.children[0],arr) 
        __tree_internal_height_lst(node.children[1],arr) 
    return arr

def tree_balance(t):
    """
    Returns B1 balance index. This is the sum of the reciprocal of the 
    heights (maximum depths) of subtrees of t rooted at internal nodes 
    of t (not including the root).
    """
    height_arr = __tree_internal_height_lst(t,[]) # get array of reciprocal subtree heights
    if(height_arr == []):
        return 0
    return sum(height_arr)
