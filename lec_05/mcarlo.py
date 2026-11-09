import numpy as np
rn = np.random.default_rng()

n = 256
rho = 0.776
temp = 0.85
r_c = 3.0 
del_xm = 0.10
n_equil = 1000
n_prod = 2000

vol = n/rho
l = (vol)**(1/3)
posi =[]

def dist(ri,rj):
    rv = ri - rj
    r = np.linalg.norm(rv)
    return r

def u_lj(r):
    coff = (1/r)
    attract =  coff**6
    repel = coff**12
    u = 4*(repel-attract)
    return u

def pbc(r):
    r = r - l*np.floor(r/l)
    return(r)


for i in range(n):
    posi.append(rn.uniform(-l,l,3))

r = []
u = []
rc = 5.5
del_r = 0.01 

for i in range(n):
    for j in range(i+1,n):
        r = dist(posi[i],posi[j])
        if r < rc:
            u.append(u_lj(r))
uto = sum(u)

ridx = np.floor(np.array((rn.uniform(0,n+1,1)))[0])
rmol_posi = posi[int(ridx)]
r_del = rn.uniform(0,del_r,3)
r_new = rmol_posi + r_del
posi[int(ridx)] = r_new

for i in range(n):
    for j in range(i+1,n):
        r = dist(posi[i],posi[j])
        if r < rc:
            u.append(u_lj(r))
uto = sum(u)
