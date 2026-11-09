import numpy as np
import matplotlib.pyplot as plt
rn = np.random.default_rng()

ei = np.array([0,1,2,1,0.5])
t = 1
beta = 1/t
wi = []
for i in ei:
    wi.append(np.exp(-beta*i))
z = np.sum(wi)
fi = wi/z
e_avg = np.sum(fi*ei)
# plt.hist(fi,bins=5)
# plt.show()

p0 =np.array([1,0,0,0])
 
