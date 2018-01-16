import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def explode(data):
    size = np.array(data.shape)*2
    data_e = np.zeros(size - 1, dtype=data.dtype)
    data_e[::2, ::2, ::2] = data
    return data_e
    
def viewpiece(piece, dimensions, ax=plt.figure().gca(projection='3d'), voxels=False, color='blue', export=None, axisoff=True):
	# prepare some coordinates
	x, y, z = np.indices(dimensions)
	for i, point in enumerate(piece):
		voxels |= (x >= point[0]) & (x <= point[0]) & (y >= point[1]) & (y <= point[1]) & (z >= point[2]) & (z <= point[2])

	if axisoff:
		plt.axis('off')
	else:
		plt.axis('on')
	ax.voxels(voxels, facecolors=color, edgecolor='#BFAB6E')
	if export is not None:
		plt.savefig(export)
	return voxels
