import os
from pathlib import Path
import h5py
import numpy as np
from scipy.stats import gaussian_kde
import pycbc.conversions as cbc
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from matplotlib import colors

from matplotlib import rcParams
# Enable LaTeX
rcParams['text.usetex'] = True

# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------
project_dir = Path.cwd()
input_folder_name = os.path.join(project_dir,"data")
output_folder_name = os.path.join(project_dir,"plots")

# ---------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------
pno_list = ["pn_three_half"]
data_dir = "results_paper"

waveform_list = ["PhenomD", "PhenomHM"]


# ---------------------------------------------------------------------
# Auxillary functions
# ---------------------------------------------------------------------
def getLevels(density):
    density_sorted = np.sort(density.ravel())[::-1]
    cdf = np.cumsum(density_sorted)
    cdf = cdf/cdf[-1]
    idx_sig1 = np.sum(cdf < 0.30)
    idx_sig2 = np.sum(cdf < 0.90)
    idx_sig3 = np.sum(cdf < 0.98)
    idx_sig4 = np.sum(cdf < 0.99)
    
    levels = [
        density_sorted[idx_sig3],
        density_sorted[idx_sig2],
        density_sorted[idx_sig1],
        1.0
        ]

    return levels

# ---------------------------------------------------------------------
# Load catalog
# ---------------------------------------------------------------------
print("Loading catalog parameters...")
catalog = {}

folder_name = os.path.join(input_folder_name, data_dir)
file_name = os.path.join(folder_name,f"{waveform_list[0]}_{pno_list[0]}.h5")

n_events = np.inf
with h5py.File(file_name, "r") as f:
    grp = f["parameter"]
    z = grp["z"][:]
    idx_redshift = z < 0.5
    idx_redshift = np.logical_and(idx_redshift, np.cumsum(idx_redshift) < n_events) 

    catalog["z"]     = z[idx_redshift]
    catalog["mc"]    = grp["mc"][:][idx_redshift]
    catalog["eta"]   = grp["eta"][:][idx_redshift]
    catalog["chi_1"] = grp["chi_1"][:][idx_redshift]
    catalog["chi_2"] = grp["chi_2"][:][idx_redshift]
    catalog["dl"]    = grp["dL"][:][idx_redshift]
    catalog["theta"] = grp["theta"][:][idx_redshift]
    catalog["phi"]   = grp["phi"][:][idx_redshift]
    catalog["iota"]  = grp["iota"][:][idx_redshift]
    catalog["psi"]   = grp["psi"][:][idx_redshift]

# ---------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------

catalog["m1"] = cbc.mass1_from_mchirp_eta(catalog["mc"], catalog["eta"])
catalog["m2"] = cbc.mass2_from_mchirp_eta(catalog["mc"], catalog["eta"])
catalog["mtot"] = catalog["m1"] + catalog["m2"]
catalog["chi_eff"] = cbc.chi_eff(
    catalog["m1"], catalog["m2"], catalog["chi_1"], catalog["chi_2"]
    )
catalog["invq"] = cbc.invq_from_mass1_mass2(catalog["m1"], catalog["m2"])

# ---------------------------------------------------------------------
# Load Fisher results
# ---------------------------------------------------------------------
print("Loading Fisher results...")
results = {}

for wf_name in waveform_list:
    results[wf_name] = {}
    for pno in pno_list:
        file_name = os.path.join(folder_name, f"{wf_name}_{pno}.h5")
        results[wf_name][pno] = {}
        with h5py.File(file_name, "r") as f:
            grp = f["results"]
            results[wf_name][pno]["delta_k"] = grp["delta_k"][:][idx_redshift]
            results[wf_name][pno]["snr"] = grp["snr"][:][idx_redshift]
            results[wf_name][pno]["isnr"] = grp["isnr"][:][idx_redshift]

# ---------------------------------------------------------------------
# Selection indices
# ---------------------------------------------------------------------
print("Processing indices...")

snr_threshold = 12
index_snr = {wf: results[wf][pno_list[0]]["snr"] > snr_threshold for wf in waveform_list}
isnr_threshold = 6
index_isnr = {wf: results[wf][pno_list[0]]["isnr"] > isnr_threshold for wf in waveform_list}

index_fisher = {
    wf: {
        pno: ~np.isnan(results[wf][pno]["delta_k"])
        for pno in pno_list
    }
    for wf in waveform_list
}

index_select = {
    wf: {
        pno: index_snr[wf] & index_fisher[wf][pno]
        for pno in pno_list
    }
    for wf in waveform_list
}

index_iselect = {
    wf: {
        pno: index_snr[wf] & (index_isnr[wf] & index_fisher[wf][pno])
        for pno in pno_list
    }
    for wf in waveform_list
}

# ---------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------
label_str = {
    "mc": r"$(1+z)\mathcal{M}$",
    "eta": r"$\eta$",
    "chi_1": r"$\chi_1$",
    "chi_2": r"$\chi_2$",
    "dl": r"$\mathrm{d}L$",
    "theta": r"$\theta$",
    "phi": r"$\phi$",
    "iota": r"$\iota$",
    "psi": r"$\psi$",
    "z": r"$z$",
    "m1": r"$m_1$",
    "m2": r"$m_2$",
    "chi_eff": r"$\chi_{eff}$",
    "invq": r"$1/q$",
    "mtot":r"$(1+z)M$"
}


# ---------------------------------------------------------------------
# KDE estimates
# ---------------------------------------------------------------------
print("Do KDE estimates...")

param = [
    "z",
    "invq",
    "iota",
    "chi_eff",
]

kde_estimate = {}

for p in param:
    print(f"   kde estimate: {p}-mc")
    
    # Extract the two parameters from the catalog
    data_x = catalog[p]
    data_y = catalog["mc"]
    
    # Stack them for a 2D KDE
    data = np.vstack([data_x, data_y])
    
    # Compute KDE
    kde_estimate[f"{p}-mc"] = gaussian_kde(data)

    print(f"   kde estimate: {p}-mtot")
    
    # Extract the two parameters from the catalog
    data_x = catalog[p]
    data_y = catalog["mtot"]
    
    # Stack them for a 2D KDE
    data = np.vstack([data_x, data_y])
    
    # Compute KDE
    kde_estimate[f"{p}-mtot"] = gaussian_kde(data)

# ---------------------------------------------------------------------
# Bounds
# ---------------------------------------------------------------------
bounds = {}    
bounds["mtot"]    = (5,200)
bounds["mc"]      = (5,80)
bounds["dl"]      = (0,10)
bounds["iota"]    = (0, np.pi)
bounds["invq"]    = (0,1)
bounds["z"]       = (0,0.5)
bounds["chi_eff"] = (-1,1)

ngrid = 200
valr = {
    "mc"  : np.linspace(bounds["mc"][0],bounds["mc"][1], ngrid),
    "mtot": np.linspace(bounds["mtot"][0],bounds["mtot"][1], ngrid),
    "z"   : np.linspace(bounds["z"][0],bounds["z"][1], ngrid),
    "iota": np.linspace(bounds["iota"][0],bounds["iota"][1], ngrid),
    "invq": np.linspace(bounds["invq"][0],bounds["invq"][1], ngrid),
    "chi_eff": np.linspace(bounds["chi_eff"][0],bounds["chi_eff"][1], ngrid)
}

pno =  pno_list[0]

# ---------------------------------------------------------------------
# selection indices
# ---------------------------------------------------------------------
idx_bwf = np.logical_and(index_select["PhenomD"][pno], index_select["PhenomHM"][pno])
idx_ibwf = np.logical_and(index_iselect["PhenomD"][pno], index_iselect["PhenomHM"][pno])
idx_iprob = np.logical_not(idx_ibwf) & idx_bwf
idx_do  = np.logical_and(index_iselect["PhenomD"][pno], np.logical_not(index_iselect["PhenomHM"][pno]))
idx_hmo = np.logical_and(np.logical_not(index_iselect["PhenomD"][pno]), index_iselect["PhenomHM"][pno])


# Select elements using boolean indexing
delta_k_hm = results["PhenomHM"][pno]["delta_k"]
delta_k_d  = results["PhenomD"][pno]["delta_k"]

# Element-wise log10 of ratio
delta_k_ratio = np.log10(delta_k_hm / delta_k_d)
delta_k_c = delta_k_ratio - np.min(delta_k_ratio)
delta_k_c = delta_k_c/np.max(delta_k_c)

# Compute min and max
clims = (np.min(delta_k_ratio[idx_bwf]), np.max(delta_k_ratio[idx_bwf]))

fig = plt.figure()
plt.rcParams['text.usetex'] = False
labelsize = 12

# set up axes
delta_x = 0.05
delta_y = 0.05
widht   = 0.3
hight   = 0.3
ax_hist = {
    "invq"    : fig.add_axes([0*(widht+delta_x) ,hight+delta_y, widht, hight]) ,
    "iota"    : fig.add_axes([1*(widht+delta_x) ,hight+delta_y, widht, hight]) ,
    "chi_eff" : fig.add_axes([2*(widht+delta_x) ,hight+delta_y, widht, hight]) ,
    "z"       : fig.add_axes([3*(widht+delta_x) ,hight+delta_y, widht, hight]) ,
    "mc"      : fig.add_axes([4*(widht+delta_x) ,0, widht, hight]) 
}

for param, ax in ax_hist.items():
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xlim(bounds[param])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

ax_hist["invq"].set_xticks([0.0, 0.25, 0.5,0.75,1.0])
ax_hist["iota"].set_xticks([0.0, np.pi/4, np.pi/2, 3*np.pi/4,np.pi])
ax_hist["z"].set_xticks([0,0.25,0.5])

ax_hist["mc"].set_xticks([5,20,35,50,65,80])
ax_hist["mc"].set_xticklabels(["5","20","35","50","65","80"])
ax_hist["mc"].set_xlabel(label_str["mc"],fontsize=labelsize)

ax_scatter = {
    "invq"    : fig.add_axes([0*(widht+delta_x) ,0, widht, hight]) ,
    "iota"    : fig.add_axes([1*(widht+delta_x) ,0, widht, hight]) ,
    "chi_eff" : fig.add_axes([2*(widht+delta_x) ,0, widht, hight]) ,
    "z"       : fig.add_axes([3*(widht+delta_x) ,0, widht, hight]) ,
}

for param, ax in ax_scatter.items():
    #ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(bounds[param])
    ax.set_ylim(bounds["mc"])
    ax.set_xlabel(label_str[param], fontsize=labelsize)
    ax.set_ylabel("")
    #ax.spines["top"].set_visible(False)
    #ax.spines["right"].set_visible(False)

ax_scatter["invq"].set_ylabel(label_str["mc"],fontsize=labelsize)
ax_scatter["invq"].set_yticks([5,20,35,50,65,80])
ax_scatter["invq"].set_xticks([0,0.25,0.5,0.75,1]) #[0.0, 0.25, 0.5,0.75,1.0])
ax_scatter["iota"].set_xticks([0.0, np.pi/4, np.pi/2, 3*np.pi/4,np.pi])
ax_scatter["iota"].set_xticklabels([0,"",r"$\pi / 2$","",r"$\pi$"])
ax_scatter["invq"].set_xticklabels(["0.0","","0.5","","1.0"])
ax_scatter["z"].set_xticks([0,0.25,0.5]) #[0.0, 0.25, 0.5,0.75,1.0])
ax_scatter["z"].set_xticklabels(["0.00","0.25","0.50"]);


cmap1c = colors.ListedColormap(['blue', 'blue'])
cmap_nomr = colors.Normalize(vmin=clims[0], vmax=clims[1])
# plot into scatter axes
for param in ["z", "iota", "invq", "chi_eff"]:
    
    param_grid, mc_grid = np.meshgrid(valr[param], valr["mc"])
    xy_vals = np.vstack([param_grid.ravel(), mc_grid.ravel()])
    density = kde_estimate[f"{param}-mc"](xy_vals)
    density = density.reshape(mc_grid.shape)  
    levels = getLevels(density)
    
    ax_scatter[param].contourf(param_grid, mc_grid, density, levels= levels, cmap="Blues")
    ax_scatter[param].contour(param_grid, mc_grid, density, levels= levels, cmap=cmap1c, linewidths = 1, alpha = 0.5)

    scatter = ax_scatter[param].scatter(
        catalog[param][idx_bwf], 
        catalog["mc"][idx_bwf], 
        c=delta_k_ratio[idx_bwf], 
        cmap='viridis',
        norm=cmap_nomr,
        s = 10,
        alpha = 1,
        zorder=10 
        ) 

    ax_scatter[param].scatter(
        catalog[param][idx_do], 
        catalog["mc"][idx_do], 
        s = 10,
        marker = "d", 
        alpha = 1,
        color = "red",
        zorder=10 
        ) 

    ax_scatter[param].scatter(
        catalog[param][idx_hmo], 
        catalog["mc"][idx_hmo], 
        s = 10,
        marker = "^", 
        alpha = 1,
        color = "red",
        zorder=10 
        ) 

    ax_scatter[param].scatter(
        catalog[param][idx_iprob], 
        catalog["mc"][idx_iprob], 
        c=delta_k_ratio[idx_iprob], 
        cmap='viridis',
        norm=cmap_nomr,
        edgecolor='k',      # black edge
        linewidth=0.3,       # very thin line
        s = 20,
        marker = "s", 
        alpha = 1,
        zorder=5 
        ) 

legend_ax = fig.add_axes([4*(widht+delta_x), hight+0.1, widht, 0.2])
legend_ax.set_xticks([])
legend_ax.set_yticks([])
legend_ax.set_xlabel("")
legend_ax.set_ylabel("")
legend_ax.spines["top"].set_visible(False)
legend_ax.spines["right"].set_visible(False)
legend_ax.spines["bottom"].set_visible(False)
legend_ax.spines["left"].set_visible(False)
square1 = patches.Rectangle((0.0, 1.00), 0.1, 0.1, transform=legend_ax.transAxes, color='blue', alpha=0.3)
square2 = patches.Rectangle((0.0, 0.75), 0.1, 0.1, transform=legend_ax.transAxes, color='green', alpha=0.3)
square3 = patches.Rectangle((0.0, 0.50), 0.1, 0.1, transform=legend_ax.transAxes, color='teal', alpha=0.3)
legend_ax.add_patch(square1)
legend_ax.add_patch(square2)
legend_ax.add_patch(square3)

# Add annotations next to the squares
legend_ax.text(0.15, 1.015, 'Events from catalog', transform=legend_ax.transAxes, verticalalignment='center')
legend_ax.text(0.15, 0.765, 'SNR > 12 and FIM invertible' , transform=legend_ax.transAxes, verticalalignment='center')
legend_ax.text(0.15, 0.515, 'Observable for HLV' , transform=legend_ax.transAxes, verticalalignment='center')

cbar_ax = fig.add_axes([4*(widht+delta_x) ,hight+delta_y, widht, 0.05])  # adjust as needed
cbar = fig.colorbar(
    scatter, 
    cax=cbar_ax, 
    orientation='horizontal')
cbar.set_label(
    r"$\mathrm{log}_{10}(\Delta_k^\mathrm{HM}) - \mathrm{log}_{10}(\Delta_k^\mathrm{D})$",
    labelpad = -50,
    fontsize = 10 
    )
cbar.set_ticks([-2.5,-2,-1.5,-1,-0.5])
cbar.set_ticklabels(["-2.5","","-1.5","","-0.5"])

n_bins = 15
for param in ["z","iota", "invq", "mc", "chi_eff"]:

    bins = np.linspace(bounds[param][0],bounds[param][1], n_bins)
    ax_hist[param].hist(
        catalog[param], 
        bins=bins, 
        color='blue', 
        alpha = 0.3,
        edgecolor="blue",
        density=True,
        histtype='stepfilled'
    )
    ax_hist[param].hist(
        catalog[param][idx_bwf], 
        bins=bins, 
        color='green', 
        alpha = 0.3,
        edgecolor="green",
        density=True,
        histtype='stepfilled'
    )
    ax_hist[param].hist(
        catalog[param][idx_ibwf], 
        bins=bins, 
        color='teal', 
        alpha = 0.3,
        edgecolor="teal",
        density=True,
        histtype='stepfilled'
    )

fig.savefig(os.path.join(output_folder_name, pno+".pdf"),bbox_inches='tight')
fig.savefig(os.path.join(output_folder_name, pno+".png"),bbox_inches='tight')

