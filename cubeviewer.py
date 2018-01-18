import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


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
	set_axes_equal(ax)
	if export is not None:
		plt.savefig(export)
	return voxels
