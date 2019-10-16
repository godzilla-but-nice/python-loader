from pyMCDS import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# load cell and microenvironment data
mcds = pyMCDS('output00000275.xml', 'Archive')
# load ECM data
mcds.load_ecm('output00000275_ECM.mat', 'Archive')

cell_df = mcds.get_cell_df()
xx, yy = mcds.get_2D_mesh()
micro = mcds.get_concentrations('ECM anisotropy', 0.0)

# find levels for microenvironment
num_levels = 10
levels = np.linspace(micro.min(), micro.max(), num_levels)

# arrow lengths
dx = mcds.data['ecm']['x_vec'] * micro
dy = mcds.data['ecm']['y_vec'] * micro
print(dx.shape)

# get unique cell types and radii
radii = (cell_df['total_volume'].values * 3 / (4 * np.pi))**(1/3)
types = cell_df['cell_type'].unique()

fig, ax = plt.subplots(figsize=(8, 10))
# add contour layer
cs = plt.contourf(xx, yy, micro, cmap='viridis', levels=levels)
ql = plt.quiver(xx, yy, dx[:, :, 0], dy[:, :, 0])

ax.axis('equal')
ax.set_xlabel('x [micron]')
ax.set_ylabel('y [micron]')
fig.colorbar(cs, ax=ax, shrink=0.75)

plt.savefig('vector.png')
