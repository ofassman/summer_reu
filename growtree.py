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
    Samples a new rate based on a gamma distribution.
    """
    scale_calc = mean/shape
    return numpy.random.gamma(shape, scale=scale_calc, size=None)

def growtree(b,d,s,max_time,shape_b,shape_d,shape_s):
    """
    Returns a birth-death tree. All rates (birth, death, and substitution) may change upon a substitution.
    'b', 'd', and 's' are the initial values of the birth, death, and substitution rates (respectively).
    'max_time' is the total amount of time that can be used to generate events to construct the tree
    (if this time is exceeded, the tree returns). Thus the tree stops growing when either all lineages
    go extinct or 'max_time' is exceeded.
    """
    rng = random.Random()
    t = Tree()
    while(True):
        # finding the wait time to any event (b, d, or s) based on rates
        curr_t = 0
        rate_any_event = b + d + s
        wait_t = rng.expovariate(rate_any_event)
        curr_t += wait_t
        if(curr_t <= max_time): # if wait time did not exceed max_time
            # calculating weighted rates
            b_weighted = b/rate_any_event
            d_weighted = d/rate_any_event
            s_weighted = s/rate_any_event
            event = gen_event(b_weighted, d_weighted, s_weighted) # generate event based on weighted rates
            if(event == "birth"): # recursively call fn for children, same rates but different max_time
                # max_time for children should be the time remaining (max_time - curr_time) divided by 2 (since 2 children)
                c1 = growtree(b,d,s,(max_time-curr_t)/2,shape_b,shape_d,shape_s)
                c2 = growtree(b,d,s,(max_time-curr_t)/2,shape_b,shape_d,shape_s)  
                if(c1 == None and c2 == None):
                    return None
                if(c1 != None):
                    t.add_child(c1)
                if(c2 != None):
                    t.add_child(c2)
                return t
            elif(event == "sub"): # change current rates based on sampling from a gamma dist and continue to next event
                # mean of gamma dist is current rate
                b = gen_rate(b,shape_b)
                d = gen_rate(d,shape_d)
                s = gen_rate(s,shape_s)
            else: # event is death so immediately return t (lineage ends)
                return None
        else: # wait time exceeded max_time (a timeout) so return tree (lineage ends)
            return t

def getNewick(t):
    if (t != None):
        return t.write(format=1)
    return ";"

def outputNewick(t,name):
    if (t != None):
        t.write(format=1, outfile=name + ".nw")
    else:
        print("Empty tree, no output file created.")