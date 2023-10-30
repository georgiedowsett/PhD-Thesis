
import os
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
os.environ["THEANO_FLAGS"] = 'device=0,floatX=float32,force_device=True'
import cell2location
from cell2location.models import RegressionModel
import scvi
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from cell2location.utils import select_slide
from scipy import sparse
from cell2location.utils.filtering import filter_genes
from cell2location.plt import plot_spatial



results_folder = f'./results/230312/C6/'
sp_data_folder = './data/st-test/'
nucseq_folder = './data/'
barcode_folder = './data/barcodes-to-remove/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

def read_and_qc(sample_name, barcode, path=sp_data_folder):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """
    #edited this to include removal of low quality spots
    barcodes = pd.read_csv(barcode_folder + str(barcode))
    adata = sc.read_visium(path + str(sample_name),
                        count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['SYMBOL']
    adata.var.drop(columns='SYMBOL', inplace=True)
    cond = adata.obs_names.isin(barcodes['Barcode'])
    adata = adata[adata.obs_names[~cond].copy()]

    # Calculate QC metrics
    from scipy.sparse import csr_matrix
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['MT'] = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['MT_frac'] = adata[:, adata.var['MT'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                        + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata



#reading in model and andata from here
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
print(inf_aver.iloc[0:5, 0:5])

#to do remove low quality spots  

# Read the list of spatial experiments
sample_data = pd.read_csv(sp_data_folder + 'sample_info.csv')

# Read the data into anndata objects
slides = []
for i, j, k, in zip(sample_data['sample_name'], sample_data['sample_name2'], sample_data['barcodes_to_remove']):
    slides.append(read_and_qc(i, j, k, path=sp_data_folder))

# Combine anndata objects together
t = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['sample_name'],
    index_unique=None
)

# mitochondria-encoded (MT) genes should be removed for spatial mapping
t.obsm['MT'] = t[:, t.var['MT'].values].X.toarray()
t = t[:, ~t.var['MT'].values]

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(t.var_names, inf_aver.index)
t = t[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=t, batch_key="sample")

# create and train the model - ncells and alpha can be modified
mod = cell2location.models.Cell2location(t, cell_state_df=inf_aver, N_cells_per_location=3, detection_alpha=20)
print(mod.view_anndata_setup())

#reduced epochs from 30000 
mod.train(max_epochs=30000,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
        use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
plt.clf()
mod.plot_history(1000)
plt.legend(labels=['full data training']);
plt.savefig(results_folder + "model_training")

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
t = mod.export_posterior(t, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True})

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
t.write(adata_file)
adata_file

plt.clf()
mod.plot_QC()
plt.savefig(results_folder + "model2_qc")

plt.clf()
fig = mod.plot_spatial_QC_across_batches()
fig.savefig(results_folder + "qc-across-batches")




adata_file = f"{run_name}/sp-clusters.h5ad"
t.write(adata_file)
