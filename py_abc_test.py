import pyabc
import scipy

prior = pyabc.Distribution(mu=pyabc.RV("expon", 0, .000047))
print(prior.rvs())
print(prior.rvs()["mu"])