#cd to directory with data
cd /projects/mjolnir1/people/crq857/Chapter2/03_GeneflowAnalyses/FEEMS/finaldataset

#create a screen session
screen
srun -N 1 -c 1 --mem-per-cpu=64g -t 1-0:0:0 --pty bash -i

# I installed feems following the instructions on the github page
conda create -n=feems_e python=3.8.3 
conda activate feems_e

brew install geos # I didn't need geos since we have a module, I use module load geos/

conda install numpy==1.22.3 scipy==1.5.0 scikit-learn==0.23.1
conda install matplotlib==3.2.2 pyproj==2.6.1.post1 networkx==2.4.0 
conda install shapely==1.7.1 
conda install fiona
conda install pytest==5.4.3 pep8==1.7.1 flake8==3.8.3
conda install click==7.1.2 setuptools pandas-plink
conda install msprime==1.0.0 statsmodels==0.12.2 PyYAML==5.4.1
conda install xlrd==2.0.1 
conda install openpyxl==3.0.7
conda install suitesparse=5.7.2
conda install scikit-sparse=0.4.4 
conda install cartopy=0.18.0


#conda create -n=feems_e python=3.8.3 
conda activate feems_e
ipython

pip install ipython #ipython, it is needed to pip install beforehand so that I can use a conda environment for feems_e instead of for ipython 

ipython #opens ipython

import feems
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer 
from pandas_plink import read_plink
import statsmodels.api as sm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv

#if any did not import I can pip install them without leaving ipython with !conda install scikit-sparse=0.4.4  

# RUN ANALYSIS 

%%time
# load data
(bim, fam, G) = read_plink("{}/Autosomes_onlyEurasianwolves_filtered_noindels_noastrick_diploid_minQ30_biallelic_maxmiss1_LDprune".format("./"))
genotypes = (np.array(G)).T

# setup graph
coord = np.loadtxt("{}/Autosomes_onlyEurasianwolves_filtered_noindels_noastrick_diploid_minQ30_biallelic_maxmiss0.9_LDprune.coord".format("./"))  # sample coordinates
outer = np.loadtxt("{}/Autosomes_onlyEurasianwolves_filtered_noindels_noastrick_diploid_minQ30_biallelic_maxmiss0.9_LDprune.outer".format("./"))  # outer coordinates
grid_path = "{}/grid_100.shp".format("/projects/mjolnir1/people/gnr216/a-software/feems/feems/data")  # path to discrete global grid

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



