import h5py
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap


from matplotlib import rcParams

# Enable LaTeX
rcParams['text.usetex'] = True

# Path to your HDF5 file
project_dir = Path.cwd()
input_folder_name = os.path.join(project_dir,"data")
output_folder_name = os.path.join(project_dir,"plots")

fig_file_name_contours = os.path.join(input_folder_name,"contours_90CI_noise_1237_realization_1.h5")

# Dictionary to store the data
# Structure: data[pno] = [ (x_curve1, y_curve1), (x_curve2, y_curve2), ... ]
data = {}

# Open the HDF5 file in read mode
with h5py.File(fig_file_name_contours, "r") as f:
    # Iterate over PN orders (groups)
    for pno in f.keys():
        pn_group = f[pno]
        curves = []
        
        # Iterate over curves in each PN order
        for curve_name in pn_group.keys():
            curve_group = pn_group[curve_name]
            x = curve_group["x"][:]
            y = curve_group["y"][:]
            curves.append((x, y))
        
        # Save curves for this PN order
        data[pno] = curves

cmap = plt.cm.coolwarm  # You can use 'plasma', 'coolwarm', etc.
# Generate 10 evenly spaced values between 0 and 1
values = np.linspace(0, 1, 10)
# Get colors from the colormap
colors = [cmap(v) for v in values]

pn_props = {
    "-1"       : ("--",colors[0], r"$1000 \times\varphi_{-2}$"),
    "0"        : ("-",colors[1], r"$\varphi_{0}$"),
    "0.5"      : ("-",colors[2], r"$\varphi_{1}$"),
    "1"        : ("-",colors[3], r"$\varphi_{2}$"),  
    "1.5"      : ("-",colors[4], r"$\varphi_{3}$"),
    "2"        : ("-",colors[5], r"$\varphi_{4}$"),
    "3"        : ("-",colors[6], r"$\varphi_{6}$"),
    "3.5"      : ("-",colors[7], r"$\varphi_{7}$"),
    "log(2.5)" : ("-",colors[8], r"$\varphi_{5l}$"),
    "log(3.)"  : ("-",colors[9], r"$\varphi_{6l}$")
}

fig = plt.figure()
ax = fig.add_axes([0,0,0.8,0.5])
labelsize = 12
legendsize = 10

ax_zoom = fig.add_axes([0.75,0.11,0.6,0.45])
ax_zoom.yaxis.tick_right()
ax_zoom.yaxis.set_label_position('right')
#ax_zoom.patch.set_alpha(0.0)  

ax.plot([0,0],[0,0.05], "--k", alpha=0.3)
for key, xy_curves in data.items():
    for c in xy_curves:
        if key == "-1":
            sf = 1000
        else:
            sf = 1
        ax.plot(sf*c[0], sf*c[1], pn_props[key][0], color= pn_props[key][1], label = pn_props[key][2],linewidth = 2)

ax.set_ylim([0., 0.05])
ax.set_xlim([-0.01, 0.02])
ax.set_xlabel(r"$\mu$", fontsize = labelsize)
ax.set_ylabel(r"$\sigma$", fontsize = labelsize)

ax_trans = fig.add_axes([0,0,0.8,0.6])
ax_trans.patch.set_alpha(0.0) 
ax_trans.set_xticklabels([])
ax_trans.set_yticklabels([])
ax_trans.set_yticks([])
ax_trans.set_xticks([])
ax_trans.set_ylim([0., 0.05])
ax_trans.set_xlim([-0.01, 0.02])
ax_trans.set_xlabel("")
ax_trans.set_ylabel("")
ax_trans.spines["top"].set_visible(False)
ax_trans.spines["bottom"].set_visible(False)
ax_trans.spines["right"].set_visible(False)
ax_trans.spines["left"].set_visible(False)
ax_trans.plot([0.001,0.001,-0.001,-0.001],[0,0.002,0.002,0], "-k", linewidth = 0.8)

ax_zoom.plot([0,0],[0,0.05], "--k", alpha=0.3)
for key, xy_curves in data.items():
    for c in xy_curves:
        if key == "-1":
            sf = 1000
        else:
            sf = 1
        ax_zoom.plot(sf*c[0], sf*c[1], pn_props[key][0], color= pn_props[key][1],linewidth = 2)

ax_zoom.set_ylim([0., 0.002])
ax_zoom.set_xlim([-0.001, 0.001])
        
ax_zoom.set_xticks([-0.001, -0.00075, -0.0005, -0.00025, 0, 0.00025, 0.0005, 0.00075, 0.001])
ax_zoom.set_xticklabels(["", "", "-0.005", "","0.00", "", "0.005", "", "0.010"])
ax_zoom.set_yticks([0, 0.0005, 0.001, 0.0015, 0.002])
ax_zoom.set_yticklabels(["0.00","", "0.01", "", "0.02"])
ax_zoom.set_xlabel(r"$10\times \mu$", fontsize = labelsize)
ax_zoom.set_ylabel(r"$10\times \sigma$", fontsize = labelsize)

ax.legend(loc='upper left', bbox_to_anchor=(-0.01, 1.28), ncol=5, frameon=False, fontsize = legendsize )
ax_trans.plot([0.001, 0.0181],[0.002,0.0092], "--k", linewidth = 0.6)
ax_trans.plot([-0.001, 0.0181],[0.002,0.0463], "--k", linewidth = 0.6)

fig.savefig(os.path.join(output_folder_name, "hyperparam_contour"+".pdf"),bbox_inches='tight')
fig.savefig(os.path.join(output_folder_name, "hyperparam_contour"+".png"),bbox_inches='tight')