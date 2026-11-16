from scipy import stats
import numpy as np

n = 10
sp = 2
ep = 3
seed = 666
c = 0.90

rn = np.random.default_rng(seed)
# data = rn.uniform(sp,ep,n)
data = np.array([2.21, 2.68, 2.38, 2.52, 2.16, 2.47, 2.30, 2.61, 2.28, 2.39])
s_mean = np.mean(data)
s_std = np.std(data,ddof=1)

# for student t distribution use stats.norm.ppf when to use z distribution  
e = stats.t.ppf((c+1)/2,n-1)

se = s_std/np.sqrt(n)

ci = [s_mean+e*se ,s_mean-e*se]
print("e = ", e)
print("se = ", se)
print("CI : ",s_mean,"±",e*se )