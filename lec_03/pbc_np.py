import numpy as np
import matplotlib.pyplot as plt

box_l = 10
box = [box_l,box_l,box_l]
rn = np.random.default_rng()

ran_posi = rn.uniform(0,box[0]*1.5,3)
posi = np.array(ran_posi)

x_new = posi[0]- box_l*np.floor(posi[0]/box_l)
y_new = posi[1]- box_l*np.floor(posi[1]/box_l)
z_new = posi[2]- box_l*np.floor(posi[2]/box_l)

pbc_posi = np.array([x_new,y_new,z_new])

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(posi[0],posi[1],posi[2],s = 100,label = 'paricle position')
ax.scatter(pbc_posi[0],pbc_posi[1],pbc_posi[2],s = 100, label = 'periodic image')
ax.legend()
plt.show()

