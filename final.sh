sftp LwKsBCw7QM@io.erda.dk

Password: LwKsBCw7QM


Extract command lines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd /scratch/genomics/hennellyl/feems 

/home/hennellyl/feems

 qrsh -l gpu
 module load bio/feems/1.0.0

python
import feems
%pip install numpy
import numpy as np
import pkg_resources
%pip install scikit-learn
from sklearn.impute import SimpleImputer 
%pip install pandas_plink
from pandas_plink import read_plink
%pip install statsmodels
import statsmodels.api as sm
%pip install matplotlib
import matplotlib.pyplot as plt
%pip install cartopy
import cartopy.crs as ccrs

from feems.utils import prepare_graph_inputs #ModuleNotFoundError: No module named 'feems.utils'

from feems import SpatialGraph, Viz #ImportError: cannot import name 'SpatialGraph' from 'feems' (unknown location)
from feems.cross_validation import run_cv #ModuleNotFoundError: No module named 'feems.cross_validation'
Autosomes_onlyEurasianwolves_filtered_noindels_noastrick_diploid_minQ30_biallelic_maxmiss1_LDprune.bed

%%time
# load data
(bim, fam, G) = read_plink("/scratch/genomics/hennellyl/feems/Autosomes_onlyEurasianwolves_filtered_noindels_noastrick_diploid_minQ30_biallelic_maxmiss1_LDprune".format("./"))
genotypes = (np.array(G)).T #MemoryError: Unable to allocate 3.12 GiB for an array with shape (8822373, 95) and data type float32

# setup graph
coord = np.loadtxt("/scratch/genomics/hennellyl/feems/Autosomes_onlyEurasianwolves_filtered_noindels_noastrick_diploid_minQ30_biallelic_maxmiss0.9_LDprune.coord".format("./"))  # sample coordinates
outer = np.loadtxt("/scratch/genomics/hennellyl/feems/Autosomes_onlyEurasianwolves_filtered_noindels_noastrick_diploid_minQ30_biallelic_maxmiss0.9_LDprune.outer".format("./"))  # outer coordinates
grid_path = "{}/grid_100.shp".format("/scratch/genomics/hennellyl/feems")  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord,
                                             ggrid=grid_path,
                                             translated=False,
                                             buffer=0,
                                             outer=outer)


print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

# plot samples
# need improvements
projection = ccrs.EquidistantConic(central_longitude=82, central_latitude=48)
fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
        edge_alpha=1, edge_zorder=100, sample_pt_size=10,
        obs_node_size=7.5, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_samples()

#!conda install matplotlib=3.5.2 #had to reupdate to a specific version
#import matplotlib.pyplot as plt
#matplotlib to 3.5.2


v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)
fig.savefig('res.eurasia_samples_Sept11_FINAL3.pdf', format="pdf", dpi=500)


# null model
sp_graph.fit_null_model()

%%time
# test lamb
lamb=0.1
sp_graph.fit(lamb = lamb)

%%time
fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=7.5, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()
fig.savefig(f'res.eurasia.lamb_{lamb}.pdf', format="pdf", dpi=500)

%%time
fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=7.5, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
# v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()
fig.savefig(f'res.eurasia.lamb_{lamb}.nosample.pdf', format="pdf", dpi=500)
%%time
# CV
# define grids
# reverse the order of lambdas and alphas for warmstart
lamb_grid = np.geomspace(1e-6, 1e2, 20)[::-1]

# run cross-validation
# mem > 200g
cv_err = run_cv(sp_graph, lamb_grid, n_folds=20, factr=1e10)
#cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)

# average over folds
mean_cv_err = np.mean(cv_err, axis=0)

# argmin of cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])

# plot cv
fig, ax = plt.subplots(dpi=300)
ax.plot(np.log10(lamb_grid), mean_cv_err, ".");
ax.set_xlabel("log10(lambda)");
ax.set_ylabel("L2 CV Error");
ax.axvline(np.log10(lamb_cv), color = "orange")
fig.savefig(f'res.eurasia.cv.pdf', format="pdf", dpi=500)

detached from 381329.pts-191.mjolnirhead01fl

# best lamb
lamb=lamb_cv
sp_graph.fit(lamb = lamb)

%%time
fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=7.5, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()
fig.savefig(f'res.eurasiafinal.lamb_{lamb}.pdf', format="pdf", dpi=500)

#best lambda output
constant-w/variance fit, converged in 154 iterations, train_loss=1179377627.5244088
lambda=100.0000000, alpha=0.8040415, converged in 523 iterations, train_loss=1154953384.3284273


