#!/usr/local/bin/python

# avoid errors due to no $DISPLAY environment variable available when running sc.pl.paga
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import h5py
import json

import scanpy.api as sc
import anndata

import time
checkpoints = {}

import dynclipy

#   ____________________________________________________________________________
#   Load data                                                               ####
task = dynclipy.main()
# git clone https://github.com/dynverse/dynclipy.git
# cd dynclipy
# pip install git+https://github.com/dynverse/dynclipy.git --upgrade --user
# R -e "devtools::install_github('dynverse/dynutils@devel', dep = F)"
# R -e "devtools::install_github('dynverse/dyncli', dep = F)"
# R -e "devtools::install_github('dynverse/dynwrap@singularity3')"
# task = dynclipy.main(
#   ["--dataset", "/code/example.h5", "--output", "/mnt/output"],
#   "/code/definition.yml"
# )

counts = task["counts"]

params = task["params"]

start_id = task["priors"]["start_id"]
if isinstance(start_id, list):
  start_id = start_id[0]

if "groups_id" in task["priors"]:
  groups_id = task["priors"]['groups_id']
else:
  groups_id = None

# create dataset
if groups_id is not None:
  obs = pd.DataFrame(groups_id)
  obs["louvain"] = obs["group_id"].astype("category")
  adata = anndata.AnnData(counts.values, obs)
else:
  adata = anndata.AnnData(counts.values)

checkpoints["method_afterpreproc"] = time.time()

#   ____________________________________________________________________________
#   Basic preprocessing                                                     ####

n_top_genes = min(2000, counts.shape[1])

# normalisation & filtering
# the recipe_zheng17 only works when > 150 cells because of `np.arange(10, 105, 5)` in filter_genes_dispersion. This should be fixed in the next scanpy release (> 1.2.2) as it is already fixed on github
if counts.shape[1] >= 150:
  sc.pp.recipe_zheng17(adata, n_top_genes=n_top_genes)
else:
  sc.pp.normalize_per_cell(adata)
  sc.pp.scale(adata)

# precalculating some dimensionality reductions
sc.tl.pca(adata, n_comps=params["n_comps"])
sc.pp.neighbors(adata, n_neighbors=params["n_neighbors"])

# denoise the graph by recomputing it in the first few diffusion components
if params["n_dcs"] != 0:
  sc.tl.diffmap(adata, n_comps=params["n_dcs"])

#   ____________________________________________________________________________
#   Cluster, infer trajectory, infer pseudotime, compute dimension reduction ###

# add grouping if not provided
if groups_id is None:
  sc.tl.louvain(adata, resolution=params["resolution"])

# run paga
sc.tl.paga(adata)

# compute a layout for the paga graph
# - this simply uses a Fruchterman-Reingold layout, a tree layout or any other
#   popular graph layout is also possible
# - to obtain a clean visual representation, one can discard low-confidence edges
#   using the parameter threshold
sc.pl.paga(adata, threshold=0.01, layout='fr', show=False)

# run dpt for pseudotime information that is overlayed with paga
adata.uns['iroot'] = np.where(counts.index == start_id)[0][0]
sc.tl.dpt(adata, n_dcs = min(adata.obsm.X_diffmap.shape[1], 10))

# run umap for a dimension-reduced embedding, use the positions of the paga
# graph to initialize this embedding
if params["embedding_type"] != 'fa':
  sc.tl.draw_graph(adata, init_pos='paga')
else:
  sc.tl.umap(adata, init_pos='paga')

checkpoints["method_aftermethod"] = time.time()

#   ____________________________________________________________________________
#   Process & save output                                                   ####
output = {}

# cell ids
output["cell_ids"] = counts.index

# grouping
grouping = pd.DataFrame({"cell_id": counts.index, "group_id": adata.obs.louvain})
output["grouping"] = grouping

# milestone network
milestone_network = pd.DataFrame(
  np.triu(adata.uns["paga"]["connectivities"].todense(), k = 0),
  index=adata.obs.louvain.cat.categories,
  columns=adata.obs.louvain.cat.categories
).stack().reset_index()
milestone_network.columns = ["from", "to", "length"]
milestone_network = milestone_network.query("length >= " + str(params["connectivity_cutoff"])).reset_index(drop=True)
milestone_network["directed"] = False
output["milestone_network"] = milestone_network

# dimred
dimred = pd.DataFrame([x for x in adata.obsm['X_umap'].T]).T
dimred.columns = ["comp_" + str(i) for i in range(dimred.shape[1])]
dimred["cell_id"] = counts.index
output["dimred"] = dimred

# branch progressions: the scaled dpt_pseudotime within every cluster
branch_progressions = adata.obs
branch_progressions["dpt_pseudotime"] = branch_progressions["dpt_pseudotime"].replace([np.inf, -np.inf], 1) # replace unreachable pseudotime with maximal pseudotime
branch_progressions["percentage"] = branch_progressions.groupby("louvain")["dpt_pseudotime"].apply(lambda x: (x-x.min())/(x.max() - x.min())).fillna(0.5)
branch_progressions["cell_id"] = counts.index
branch_progressions["branch_id"] = branch_progressions["louvain"].astype(np.str)
branch_progressions = branch_progressions[["cell_id", "branch_id", "percentage"]]
output["branch_progressions"] = branch_progressions

# branches:
# - length = difference between max and min dpt_pseudotime within every cluster
# - directed = not yet correctly inferred
branches = adata.obs.groupby("louvain").apply(lambda x: x["dpt_pseudotime"].max() - x["dpt_pseudotime"].min()).reset_index()
branches.columns = ["branch_id", "length"]
branches["branch_id"] = branches["branch_id"].astype(np.str)
branches["directed"] = True
output["branches"] = branches

# branch network: determine order of from and to based on difference in average pseudotime
branch_network = milestone_network[["from", "to"]]
average_pseudotime = adata.obs.groupby("louvain")["dpt_pseudotime"].mean()
for i, (branch_from, branch_to) in enumerate(zip(branch_network["from"], branch_network["to"])):
  if average_pseudotime[branch_from] > average_pseudotime[branch_to]:
    branch_network.at[i, "to"] = branch_from
    branch_network.at[i, "from"] = branch_to

output["branch_network"] = branch_network

# timings
timings = pd.Series(checkpoints)
timings.index.name = "name"
timings.name = "timings"
output["timings"] = timings

# save
dynclipy.write_output(
  output,
  task["output"],
  ["branch_trajectory", "timings"]#, "dimred"]
)


dataset = dynclipy.wrap_data(cell_ids = counts.index)
dataset.add_branch_trajectory(
  grouping = output["grouping"], 
  milestone_network = output["milestone_network"],
  branch_progressions = output["branch_progressions"],
  branches = output["branches"],
  branch_network = output["branch_network"]
)
dataset.add_dimred(
  dimred = output["dimred"]
)
dataset.write_output(task["output"])
