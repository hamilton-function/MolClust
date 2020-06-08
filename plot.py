import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

k = [2,2,2,2,2,4,4,4,4,4,6,6,6,6,6,8,8,8,8,8,10,10,10,10,10]
l = [10,8,6,4,2,10,8,6,4,2,10,8,6,4,2,10,8,6,4,2,10,8,6,4,2]
#rl = [2.256, 3.052, 3.742, 4.042, 8.778, 2.034, 3.247, 4.8, 4.993, 9.434, 1.236, 1.292, 4.154, 6.611, 9.357, 0.439, 1.394, 2.862, 6.934, 11.776, 0.078, 0.648, 2.102, 3.391, 13.109]
#rl = [2.125,2.512,3.651,5.418,10.213,2.019,2.696,4.153,5.121,9.764,1.095,2.076,3.885,6.192,8.831,0.502,1.391,3.511,6.083,10.043,0.18,0.739,2.178,5.453,9.637]
#r total
rl = [8.125,7.512,7.651,8.418,12.213,9.019,8.696,9.153,9.121,12.764,9.095,9.076,9.885,11.192,12.831,9.502,9.391,10.511,12.083,15.043,10.18,9.739,10.178,12.453,15.637]

triang = mtri.Triangulation(k, l)

ax.plot_trisurf(triang, rl, cmap='jet')
ax.scatter(k,l,rl, marker='.', s=10, c="black", alpha=0.5)
ax.view_init(elev=60, azim=-45)


ax.set_xlabel('k')
ax.set_ylabel('l')
ax.set_zlabel('rec. loss')
#ax.scatter(k,l,rl)
plt.show()


