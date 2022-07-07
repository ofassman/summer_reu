from ete3 import Tree
import random
import numpy
import statistics

__seq_dict = {} # A global dictionary with 'sequence number : sequence' pairs. Populated by 'growtree()'.
__seq_counter = 0 # A global counter for sequence number. Each cell has a unique number. Incremented by 'growtree()'.

def gen_sequence(length, off_lim = None):
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

def gen_event(b, d, s):
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
    
def gen_rate(mean, shape):
    """
    Samples a new rate based on a gamma distribution given the mean rate and the shape of the distribution. 
    """
    scale_calc = mean / shape
    return numpy.random.gamma(shape, scale = scale_calc, size = None)

def growtree(seq, b, d, s, max_time, shape_b, shape_d, shape_s, branch_info):
    """
    Returns a birth-death tree. Used as a recursive helper function for 'gen_tree()' that produces
    the birth-death tree. Populates '__seq_dict' with 'sequence number : sequence' pairs. 
    """
    global __seq_counter # declaring global counter for sequence number in '__seq_dict'
    rng = random.Random()
    # initializing the tree and branch length
    t = Tree()
    key = "SEQUENCE_" + str(__seq_counter) # create sequence number key using the global counter
    t.name = key
    __seq_dict[key] = seq # set the 'sequence number : sequence' pair in the '__seq_dict' dictionary
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
            b_weighted = b / rate_any_event
            d_weighted = d / rate_any_event
            s_weighted = s / rate_any_event
            event = gen_event(b_weighted, d_weighted, s_weighted) # generate event based on weighted rates
            if(event == "birth"): # recursively call fn for children, same rates but different max_time
                # max_time for children should be the time remaining (max_time - curr_time) divided by 2 (since 2 children)
                __seq_counter += 1 # increment global counter so that each child has a unique sequence number
                c1 = growtree(seq, b, d, s, (max_time - curr_t) / 2, shape_b, shape_d, shape_s, branch_info)
                __seq_counter += 1 # increment global counter so that each child has a unique sequence number
                c2 = growtree(seq, b, d, s, (max_time-curr_t)/2, shape_b, shape_d, shape_s, branch_info)  
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
                b = gen_rate(b, shape_b) # generate new birth rate
                d = gen_rate(d, shape_d) # generate new death rate
                s = gen_rate(s, shape_s) # generate new sub rate
                sub_site = random.randint(0, len(seq) - 1) # randomly pick a site to sub a base in the sequence
                old_letter = seq[sub_site] # find old base at this site so that the sub does not change the site to the same base
                sub_letter = gen_sequence(1, off_lim = old_letter) # generate a new base for this site that is not the old base (not 'old_letter')
                # generate the new sequence using the old sequence with one base changed from 'old_letter' to 'new_letter'
                # at index 'sub_site' in the sequence
                new_seq = ""
                for i in range(0, sub_site, 1):
                    new_seq += seq[i : i + 1]
                new_seq += sub_letter
                for j in range(sub_site + 1, len(seq), 1):
                    new_seq += seq[j : j + 1]
                seq = new_seq # update sequence to newly mutated sequence
                __seq_dict[key] = seq # update the 'sequence number : sequence' pair in the '__seq_dict' dictionary
                # if branch length is a variable of number of substitutions, increase lineage's branch length by 1
                if(branch_info == 1):
                    t.dist += 1
            else: # event is death so return None (lineage goes extinct)
                return None
        else: # wait time exceeded max_time so return tree (lineage ends from timeout)
            return t

def gen_tree(b, d, s, max_time, shape_b, shape_d, shape_s, branch_info, seq_length):
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
    Translations from sequence number to sequence can be found in '__seq_dict' (a dictionary populated with
    'sequence number : sequence' pairs).
    """
    seq = gen_sequence(seq_length) # generate random genetic sequence for root cell 
    t = growtree(seq, b, d, s, max_time, shape_b, shape_d, shape_s, branch_info) # generate the tree with a recursive function
    return t

def getNewick(t):
    """
    Returns a tree in Newick tree format.
    """
    if (t != None):
        return t.write()
    return ";"

def outputNewick(t, name):
    """
    Writes a tree, 't', in Newick tree format into a file. 'name' specifies the 
    file's name in which the tree is written into. If the tree is empty (i.e. if 
    't' is 'None') no output file is created.
    """
    if (t != None):
        t.write(outfile = name + ".nw")
    else:
        print("Empty tree, no output file created.")

def get_seq(seq_num = None):
    """
    Returns a genetic sequence for a cell in the simulated phylogeny (simulated using 'gen_tree'). 
    If 'seq_num' is None (which is the default), the whole 'sequence number : sequence' dictionary
    is returned. Otherwise, 'seq_num' can be initialized to an integer argument for the specific
    sequence corresponding to that sequence number to be returned.
    """
    if(seq_num == None): # 'seq_num' is not initialized so return whole dictionary
        return __seq_dict
    else:
        key = "SEQUENCE_" + str(seq_num) # find the specific key corresponding to 'seq_num'
        value = __seq_dict.get(key) # find the sequence in the dictionary corresponding to the key
        if(value == None): # no sequence exists corresponding to 'key'
            print(key, "does not exist.")
            return None
        return value # return this specific sequence

def print_seq(seq_num = None):
    """
    Prints a genetic sequence for a cell in the simulated phylogeny (simulated using 'gen_tree'). 
    If 'seq_num' is None (which is the default), the whole 'sequence number : sequence' dictionary
    is printed in FASTA format. Otherwise, 'seq_num' can be initialized to an integer argument 
    for the specific sequence corresponding to that sequence number to be printed in FASTA format
    (i.e. the single 'sequence number : sequence' pair is printed).
    """
    if(seq_num == None): # 'seq_num' is not initialized so print whole dictionary in FASTA format
        for key, value in __seq_dict.items():
            print('>', key, '\n', value)
    else:
        key = "SEQUENCE_" + str(seq_num) # find the specific key corresponding to 'seq_num'
        value = __seq_dict.get(key) # find the sequence in the dictionary corresponding to the key
        if(value == None): # no sequence exists corresponding to 'key'
            print(key, "does not exist.")
        else:
            print('>', key, '\n', value) # print this 'key' : 'value' pair in FASTA format

#################### TREE STATISTICS ########################################

def __tree_branch_lst(t, arr):
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
    branch_arr = __tree_branch_lst(t, []) # get array of branch lengths
    if(branch_arr == []):
        return 0
    return sum(branch_arr)

def tree_branch_mean(t):
    """
    Returns the mean of the distances of all the branches in the tree. 
    """
    branch_arr = __tree_branch_lst(t,[ ]) # get array of branch lengths
    if(branch_arr == []):
        return 0
    return statistics.mean(branch_arr)

def tree_branch_median(t):
    """
    Returns the median of the distances of all the branches in the tree. 
    """
    branch_arr = __tree_branch_lst(t, []) # get array of branch lengths
    if(branch_arr == []):
        return 0
    return statistics.median(branch_arr)

def tree_branch_variance(t):
    """
    Returns the variance of the distances of all the branches in the tree. 
    """
    branch_arr = __tree_branch_lst(t, []) # get array of branch lengths
    if(branch_arr == [] or len(branch_arr) < 2):
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
    #print(t)
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

def __tree_depth_lst(node, arr):
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
            __tree_depth_lst(node.children[0], arr) 
        elif(num_c == 2): # tree with 2 children
            __tree_depth_lst(node.children[0], arr) 
            __tree_depth_lst(node.children[1], arr) 
    return arr

def tree_depth_mean(t):
    """
    Returns the mean of leaf depths in the tree. 
    """
    depth_arr = __tree_depth_lst(t, []) # get array of leaf depths
    if(depth_arr == []):
        return 0
    return statistics.mean(depth_arr)

def tree_depth_median(t):
    """
    Returns the median of leaf depths in the tree. 
    """
    depth_arr = __tree_depth_lst(t, []) # get array of leaf depths
    if(depth_arr == []):
        return 0
    return statistics.median(depth_arr)

def tree_depth_variance(t):
    """
    Returns the variance of leaf depths in the tree. 
    """
    depth_arr = __tree_depth_lst(t, []) # get array of leaf depths
    if(depth_arr == [] or len(depth_arr) < 2):
        return 0
    return statistics.variance(depth_arr)

def __tree_internal_height_lst(node, arr):
    """
    Returns an array of the reciprocal of the heights (maximum depths) 
    of subtrees of t rooted at internal nodes of t (not including the root). 
    Private helper function used for calculating summary stats regarding 
    tree balance.
    """
    if(node == None):
        return []
    if(not(node.is_leaf()) and not(node.is_root())): # if node is internal (not root or leaf)
        height_subtree = tree_height(node)
        if(height_subtree != 0):
            arr.append(1/tree_height(node)) # add reciprocal of height of subtree rooted at 'node' to 'arr'
        else:
            arr.append(0)
    # must find all internal nodes to add reciprocal of heights to 'arr'
    num_c = len(node.children)  
    if(num_c == 1): # tree with 1 child
        __tree_internal_height_lst(node.children[0], arr) 
    elif(num_c == 2): # tree with 2 children
        __tree_internal_height_lst(node.children[0], arr) 
        __tree_internal_height_lst(node.children[1], arr) 
    return arr

def tree_balance(t):
    """
    Returns B1 balance index. This is the sum of the reciprocal of the 
    heights (maximum depths) of subtrees of t rooted at internal nodes 
    of t (not including the root).
    """
    height_arr = __tree_internal_height_lst(t, []) # get array of reciprocal subtree heights
    if(height_arr == []):
        return 0
    return sum(height_arr)

def tree_nleaf(t):
    """
    Returns number of leaves in a tree.
    """
    if(t == None): # empty tree
        return 0
    if(t.is_leaf()): # node is leaf so increment
        return 1
    num_leaves = 0
    num_c = len(t.children)  
    if(num_c == 1): # tree with 1 child
        num_leaves = tree_nleaf(t.children[0]) 
    elif(num_c == 2): # tree with 2 children
        num_leaves = tree_nleaf(t.children[0]) + tree_nleaf(t.children[1]) # add leaves of both children
    return num_leaves

def __tree_leaf_diff(node):
    """
    Returns the absolute value of the difference between the number of leaves
    on the left side of the tree and the number of leaves on the right side.
    Private helper function used for calculating summary stats regarding the
    colless index.
    """
    if(node == None):
        return 0
    if(not(node.is_leaf())):
        nleaf_l = 0 # will hold number of leaves on left side of 'node'
        nleaf_r = 0 # will hold number of leaves on right side of 'node'
        num_c = len(node.children)  
        if(num_c == 1): # tree with 1 child
            nleaf_l = tree_nleaf(node.children[0]) 
        elif(num_c == 2): # tree with 2 children
            nleaf_l = tree_nleaf(node.children[0]) 
            nleaf_r = tree_nleaf(node.children[1])
        return abs(nleaf_l - nleaf_r)
    return 0

def __tree_leaf_diff_lst(node, arr):
    """
    Returns an array of absolute values of the differences between the number 
    of leaves on the left side of a node and the number of leaves on the right 
    side for all internal nodes including the root. Private helper function 
    used for calculating summary stats regarding the colless index.
    """
    if(node == None):
        return []
    if(not(node.is_leaf())): # node is internal (not a leaf)
        arr.append(__tree_leaf_diff(node)) # add leaf difference to array
    num_c = len(node.children)  
    if(num_c == 1): # tree with 1 child
        __tree_leaf_diff_lst(node.children[0], arr) 
    elif(num_c == 2): # tree with 2 children
        __tree_leaf_diff_lst(node.children[0], arr) 
        __tree_leaf_diff_lst(node.children[1], arr) 
    return arr

def tree_root_colless(t):
    """
    Returns the colless index for the root node of a tree. That is,
    returns the absolute value of the difference between the number of leaves
    on the left side of the tree and the number of leaves on the right side
    for the root node.
    """
    if(t == None):
        return 0
    return __tree_leaf_diff(t)

def tree_sum_colless(t):
    """
    Returns the sum of the colless indices for every internal node of a tree. 
    """
    if(t == None):
        return 0
    arr = __tree_leaf_diff_lst(t, []) # get array of colless indices
    if(arr == []):
        return 0
    return sum(arr)

def tree_mean_colless(t):
    """
    Returns the mean of the colless indices for every internal node of a tree. 
    """
    if(t == None):
        return 0
    arr = __tree_leaf_diff_lst(t, []) # get array of colless indices
    if(arr == []):
        return 0
    return statistics.mean(arr)

def tree_median_colless(t):
    """
    Returns the median of the colless indices for every internal node of a tree. 
    """
    if(t == None):
        return 0
    arr = __tree_leaf_diff_lst(t, []) # get array of colless indices
    if(arr == []):
        return 0
    return statistics.median(arr)

def tree_variance_colless(t):
    """
    Returns the variance of the colless indices for every internal node of a tree. 
    """
    if(t == None):
        return 0
    arr = __tree_leaf_diff_lst(t, []) # get array of colless indices
    if(arr == [] or len(arr) < 2):
        return 0
    return statistics.variance(arr)