

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/' \
                 'Fjord/Kangerlussuaq/Figures/Ocean'
vs = [5, 10, 15, 20]

# H = H0 - m*d/v
# (H - H0)*v = - m*d
# (H0-H)*v/d = m

plot_width = 8
plot_height = 8

vmin = 0.03
vmax = 1

fig = plt.figure(figsize=(10,8))
gs = GridSpec(2*plot_height+1, 2*plot_width+3, left=0.08, right=0.92, top=0.93, bottom=0.08)

counter = 1
for i in range(len(vs)):

    v = vs[i]

    h0 = np.arange(30,120,0.1)
    d = np.arange(5,20,0.1)
    H0, D = np.meshgrid(h0, d)

    M = H0*v/D # m * km/yr / km
    M *= 1/365

    if counter==1:
        ax = fig.add_subplot(gs[:plot_height, i*plot_width:(i+1)*plot_width])
        ax.set_xticklabels([])
        ax.set_ylabel('Melange Distance ($d$, km)')
    if counter==2:
        ax = fig.add_subplot(gs[:plot_height, 1+i*plot_width:1+(i+1)*plot_width])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
    if counter == 3:
        ax = fig.add_subplot(gs[plot_height+1:, (i-2) * plot_width:(i - 1) * plot_width])
        ax.set_xlabel('Initial Thickness ($H_0$, m)')
        ax.set_ylabel('Melange Distance ($d$, km)')
    if counter == 4:
        ax = fig.add_subplot(gs[plot_height+1:, 1+(i-2) * plot_width:1+(i - 1) * plot_width])
        ax.set_yticklabels([])
        ax.set_xlabel('Initial Thickness ($H_0$, m)')
    C = ax.pcolormesh(H0, D, M, cmap='turbo', vmin = vmin, vmax = vmax)
    C = ax.contour(H0, D, M, levels=[0.08], colors='w')
    C = ax.contour(H0, D, M, levels=[0.04,0.11], colors='w')
    C = ax.contour(H0, D, M, levels=[0.2], colors='k')
    C = ax.contour(H0, D, M, levels=[0.18, 0.22], colors='k')
    ax.set_title('v = '+str(v)+' km/yr')

    counter +=1

axc = fig.add_subplot(gs[2:-2, -1])
x = np.array([0,1])
y = np.linspace(vmin, vmax, 100)
X, Y = np.meshgrid(x, y)
axc.pcolormesh(X,Y,Y, cmap='turbo', vmin = vmin, vmax = vmax)
axc.set_ylabel('Melt Rate (m/day)')
axc.set_xticklabels([])

plt.suptitle('Melance Melt Rate as a function of distance, velocity, and melt ($m = H_0v/d$)')

plt.savefig(os.path.join(project_folder,'melange_thickness_comparison.png'))
plt.close(fig)








