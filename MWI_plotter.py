"""
Demonstration of MWI imaging system data plotting routines for
K. Verhaegh et al. "Long-legged Divertors and Neutral Baffling as a
Solution to the Tokamak Power Exhaust Challenge"
See for example Fig. 5 therein

Experimental data (MU02):
    46860: SXD
    47079: ED
    46866: CD

Camera numbers:
    cam 2:
    cam 5:
    cam 6:
    cam 7: hydrogen Balmer-alpha
    cam 8:
    cam 9: molecular hydrogen Fulcher band

MWI routines are stored in the mwi_dp package, hosted on UKAEA Gitlab (git.ccfe.ac.uk/twijkamp/mwi_dp) from DOI: 10.1088/1741-4326/acc191.

Routine by T. Wijkamp
"""

"""
import modules
"""
# MWI_inversions_reader is found in /mwi_dp/tools/inversions/MU02
import numpy as np
from MWI_inversions_reader import MWI_inversions_reader
import mwi_dp.grid.grid_tools as grid_tools

"""
user input for loading data
"""

# time range
shot_nb = 46860
t_range = [0.3,0.8]
cam_nbs = [7,9]

# map inversion grid to new grid which is cut to tursted region
re_grid = False

"""
Load the stored inversion data
cio['camX'] stores the results for camera number X
"""
cio = MWI_inversions_reader(shot_nb, cam_nbs,\
      t_range=t_range, quick=False, wait_until_found=True,\
      spectral_cal=False, get_std=False, std_quick=False)

if re_grid:
    print('Interpolate to cut grid.')
    grid_new = '/home/twijkamp/Documents/data/cherab/grids/grid_MAST-U_MU02_5mm_reduced.mat'
    for cam_tag in cio:
        cio[cam_tag].set_new_grid(grid_new)
print('Done loading data')


"""
Plotting tools
"""

# show inversions
CAD_file = '/home/twijkamp/Documents/data/CAD/MAST-U_passive_structure.mat'
cam_nb_plot = 7
plt_time = 0.65
epsilon_max = 4E20
ax, fig = grid_tools.plt_reactor_contour_2D(CAD_file)
cio['cam'+str(cam_nb_plot)].plt_emissivity(time=plt_time,cmap='hot',normalized=False,\
        eps_max=epsilon_max,\
        show_EFIT=True, norm_psi_bound = [0.9,1.1],\
        new_fig = False, ax_plt = ax, fig_plt = fig, ignore_grid_bound = True)

# gif writer
cam_nb_plot = 7
plt_time = 0.65
epsilon_max = 4E18
video_directory = '/home/twijkamp/Documents/script/figures/paper_Kevin_nature_2024/figures/'
video_name = 'emissivity_'+str(shot_nb)+'_cam'+str(cam_nb_plot)+'.gif'
epsilon_max = None
CAD_file = '/home/twijkamp/Documents/data/CAD/MAST-U_passive_structure.mat'
ax, fig = grid_tools.plt_reactor_contour_2D(CAD_file)
cio[cam_tag].save_emissivity_video(video_directory,video_name,video_duration=10,\
    time_bounds=t_range, cmap='hot', eps_max=epsilon_max,\
    R_bound = None, Z_bound = None, show_EFIT=True, log10 = False,\
    new_fig = False, ax_plt = ax, fig_plt = fig, ignore_grid_bound = True)
