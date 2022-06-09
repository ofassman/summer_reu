from ete3 import Tree
import random
import numpy

def gen_event(b,d,s):
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
    scale_calc = mean/shape
    return numpy.random.gamma(shape, scale=scale_calc, size=None)

def growtree(b,d,s,time,shape_b,shape_d,shape_s):
    rng = random.Random()
    t = Tree()
    while(True):
        curr_t = 0
        rate_any_event = b + d + s
        wait_t = rng.expovariate(rate_any_event)
        curr_t += wait_t
        if(curr_t <= time):
            b_weighted = b/rate_any_event
            d_weighted = d/rate_any_event
            s_weighted = s/rate_any_event
            event = gen_event(b_weighted, d_weighted, s_weighted)
            #print(event)
            if(event == "birth"):
                c1 = growtree(b,d,s,(time-curr_t)/2,shape_b,shape_d,shape_s)
                t.add_child(c1)
                c2 = growtree(b,d,s,(time-curr_t)/2,shape_b,shape_d,shape_s)  
                t.add_child(c2)
                #print("returning after b")
                return t
            elif(event == "sub"):
                b = gen_rate(b,shape_b)
                d = gen_rate(d,shape_d)
                s = gen_rate(s,shape_s)
            else:
                #print("returning after d")
                return t
        else: 
            #print("returning after t")
            return t