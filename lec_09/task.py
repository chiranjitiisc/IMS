from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

n = 1000
sp = 1
ep = 10
seed = 666
c = 0.90

rng = np.random.default_rng(seed)

def rel_data(n, p_keep=0.9):
    x = np.empty(n)
    x[0] = rng.integers(sp,ep)
    for i in range(1, n):
        if rng.random() < p_keep:
            x[i] = x[i-1] 
        else: 
            x[i] = rng.integers(sp,ep)
    return x

b_size = np.linspace(10,100,10,dtype=int)
# b_size = [10,int(n/10)]
for j in b_size:
    b_data = np.reshape(rel_data(n),[b_size[j],int(n/b_size(j))])
    b_avgs = []
    for i in range (len(b_data)):
        avg = np.mean(b_data[i])
        b_avgs.append(avg)
    b_avgs = np.array(b_avgs)

