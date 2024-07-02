import biom
import pandas as pd
import numpy as np
import csv
from sklearn import linear_model
from scipy import stats
from itertools import product
import datetime
import os
import skbio.tree as sktree
import matplotlib.pyplot as plt


def biom2df(_bt):
    """
    Converts biom table (observations x samples) to 
    Pandas dataframe.
    """

    m = _bt.matrix_data
    data = [pd.SparseSeries(m[i].toarray().ravel()) for i in np.arange(m.shape[0])]
    out = pd.SparseDataFrame(data, index=_bt.ids('observation'),
                             columns=_bt.ids('sample'))

    return(pd.DataFrame(out))


def biom_to_df(table):

    return pd.DataFrame(table.matrix_data.todense(), index=table.ids(axis="observation"), columns=table.ids(axis="sample"))


def meta_biom_filter(mymeta, mybiom, anon_name, filter_num):
    mymeta = mymeta.sort_index()
    
    biom_df = biom2df(mybiom)
    meta_df = mymeta.loc[mymeta.ANONYMIZED_NAME == anon_name]
    meta_df = meta_df.loc[set(meta_df.index) & set(biom_df.columns)]
    biom_df = biom_df[list(set(meta_df.index) & set(biom_df.columns))]
    biom_df = biom_df.sort_index()
    biom_df = biom_df.loc[(biom_df.sum(axis=1) > filter_num)]
    
    avg_day = pd.DataFrame(index=biom_df.index)

    for i in set(meta_df.epoch_time):
        samples_day = meta_df.index[meta_df.epoch_time == i]
        avg_day[i] = biom_df.transpose().loc[samples_day].mean()
    
    avg_day = avg_day.reindex_axis(sorted(avg_day.columns), axis=1)
    
    return(meta_df, avg_day)


class claatu:

    def __init__(self, meta, biom, tree, rared):

        self.meta = meta
        self.biom = biom
        self.rared = rared
        self.nodes_thresh = {}
        self.all_regress = {}
        self.all_regress_thresh = {}
        
        self.nodes = list(self.biom.index)

    
    def lag_regress(self, i, j):
        epsilon = 0.000001
        j_relabund = (self.biom.loc[j] + epsilon)/self.rared
        j_term = np.log10(j_relabund).iloc[2:]#.values.reshape(-1,1)
        i_relabund = (self.biom.loc[i] + epsilon)/self.rared
        delta_i_term = (np.log10(i_relabund).iloc[2:] - 
                        np.log10(i_relabund.shift(1).iloc[2:] + epsilon))#.values.reshape(-1,1)

        #reg = LinearRegression()
        #reg.fit(j_term, delta_i_term)
        
        #return([i, j, reg.coef_[0][0], reg.p[0][0]])
        myregression = stats.pearsonr(j_term, delta_i_term)
        return(i, j, myregression[0], myregression[1])
    
        
    def plot_lag(self, nodes):
        i = nodes[0]
        j = nodes[1]
        
        epsilon = 0.000001
        j_relabund = (self.biom.loc[j] + epsilon)/self.rared
        j_term = np.log10(j_relabund).iloc[2:].values.reshape(-1,1)
        i_relabund = (self.biom.loc[i] + epsilon)/self.rared
        delta_i_term = (np.log10(i_relabund).iloc[2:] - 
                        np.log10(i_relabund.shift(1).iloc[2:] + epsilon)).values.reshape(-1,1)
        
        plt.scatter(j_term, delta_i_term)
        

    def plot_abund(self, node):
        
        plt.plot(self.biom.loc[node])
        
    
    def shuffle_biom(self):
        new_cols = np.random.permutation(self.biom.columns)
        self.biom = self.biom[new_cols]
        
    


# In[4]:

if rank == 0:

    myotubiom = biom.load_table("/home/operon/Documents/lozupone_lab/dada2_stoolsal_out/dada2_otu_table_w_tax_no_pynast_failures_rare7500.tsv")

    mymeta = pd.DataFrame.from_csv("/home/operon/Documents/lozupone_lab/dada2_stoolsal_out/full_maps_corrected.txt", sep = "\t")

    mytree = sktree.TreeNode.read("/home/operon/Documents/lozupone_lab/claatu_bangladesh/new_prepped_tree.tre")

    rared = 7500

    meta_df, biom_df = meta_biom_filter(mymeta, mybiom, "F01", 1000)

    f01 = claatu(meta=meta_df, biom=biom_df, tree=mytree, rared=7500)

    blah = [i for i in product(f01.nodes, f01.nodes)]
    blah = np.array_split(blah, size)
    
else:
    f01 = None
    blah = None

#f01 = comm.bcast(f01, root=0)

blah = comm.scatter(blah, root=0)

