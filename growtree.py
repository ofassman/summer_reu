from ete3 import Tree
import random
import numpy

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

def growtree(b,d,s,max_time,shape_b,shape_d,shape_s,branch_info):
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
                c1 = growtree(b,d,s,(max_time-curr_t)/2,shape_b,shape_d,shape_s,branch_info)
                c2 = growtree(b,d,s,(max_time-curr_t)/2,shape_b,shape_d,shape_s,branch_info)  
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
                # if branch length is a variable of number of substitutions, increase lineage's branch length by 1
                if(branch_info == 1):
                    t.dist += 1
            else: # event is death so return None (lineage goes extinct)
                return None
        else: # wait time exceeded max_time so return tree (lineage ends from timeout)
            return t

def getNewick(t):
    """
    Returns a tree in Newick tree format.
    """
    if (t != None):
        return t.write(format=1)
    return ";"

def outputNewick(t,name):
    """
    Writes a tree, 't', in Newick tree format into a file. 'name' specifies the 
    file's name in which the tree is written into. If the tree is empty (i.e. if 
    't' is 'None') no output file is created.
    """
    if (t != None):
        t.write(format=1, outfile=name + ".nw")
    else:
        print("Empty tree, no output file created.")

def tree_sum(t):
    """
    Returns the sum of all the distances of the branches in the tree. 
    """
    if(t == None):
        return 0
    left_h = 0
    right_h = 0
    num_c = len(t.children)  
    if(num_c == 1):
        left_h = tree_sum(t.children[0]) + t.children[0].dist
    elif(num_c == 2):
        left_h = tree_sum(t.children[0]) + t.children[0].dist
        right_h = tree_sum(t.children[1]) + t.children[1].dist
    return left_h + right_h