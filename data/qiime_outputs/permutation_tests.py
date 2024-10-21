import copy
import biom
import pandas as pd
import numpy as np
import csv
from scipy import stats
from itertools import product
import datetime
import os
import matplotlib.pyplot as plt
from skbio.stats.composition import clr
from scipy.interpolate import PchipInterpolator


def biom2df(_bt):
    """Converts `biom` fromat to pandas dataframe"""
    m = _bt.matrix_data
    data = [pd.SparseSeries(m[i].toarray().ravel()) for i in np.arange(m.shape[0])]
    out = pd.SparseDataFrame(data, index=_bt.ids('observation'),
                             columns=_bt.ids('sample'))

    return(pd.DataFrame(out))


def biom_to_df(table):
    
    return pd.DataFrame(table.matrix_data.todense(), index=table.ids(axis="observation"), columns=table.ids(axis="sample"))

def interp_missing_day(biom_df):
    pchip = PchipInterpolator(biom_df.columns, biom_df.iloc[0])



def meta_biom_filter(mymeta, mybiom, anon_name, filter_num):
    """
    Filters biom table and metadata by subject id based on the prevalence of the
    OTU within the timeseries.
    
    Then return the mean day.
    """
    mymeta = mymeta.sort_index()
    
    biom_df = biom2df(mybiom)
    meta_df = mymeta.loc[mymeta.ANONYMIZED_NAME == anon_name]
    meta_df = meta_df.loc[set(meta_df.index) & set(biom_df.columns)]
    biom_df = biom_df[list(set(meta_df.index) & set(biom_df.columns))]
    biom_df = biom_df.sort_index()
    biom_df = biom_df.loc[(biom_df.sum(axis=1) > filter_num)]
    
    avg_day = pd.DataFrame(index=biom_df.index)
    avg_day = clr(avg_day)

    # for days with > 1 sample, take the mean value for each prevalent OTU
    for i in set(meta_df.epoch_time):
        samples_day = meta_df.index[meta_df.epoch_time == i]
        avg_day[i] = biom_df.transpose().loc[samples_day].mean()
    # sort dataframe columns by chronological order
    avg_day = avg_day.reindex_axis(sorted(avg_day.columns), axis=1)

    avg_day = interp_missing_day(avg_day)
    
    return(meta_df, avg_day)


class MiBiTimeSeries:

    def __init__(self, meta, biom, rared):

        self.meta = meta
        self.biom = biom
        self.rared = rared
        self.nodes_thresh = {}
        self.all_regress = {}
        self.all_regress_thresh = {}
        
        self.biom = biom_df.loc[list(set(biom_df.index) & set([i.name for i in self.tree.traverse()]))]
        self.nodes = list(self.biom.index)
        self.epsilon = 0.000001


    def lag_regress(self, i, j):
        """Get rows i and j (OTUs) from the biom dataframe.
        Convert to relative abundaces and calculate the lagged
        correlation coefficient between j and delta i."""
        # takes 1ms to run.
        # switch to np arrays gives .6ms runtime
        j_relabund = np.array((self.biom.loc[j] + self.epsilon)/self.rared)
        j_term = np.log10(j_relabund)[1:]
        i_relabund = np.array((self.biom.loc[i] + self.epsilon)/self.rared)
        delta_i_term = (np.log10(i_relabund)[1:] - np.log10(np.roll(i_relabund, -1)[1:] + self.epsilon))

        #reg = LinearRegression()
        #reg.fit(j_term, delta_i_term)
        
        #return([i, j, reg.coef_[0][0], reg.p[0][0]])
        myregression = stats.spearmanr(j_term, delta_i_term)
        return(i, j, myregression[0], myregression[1])
        
    
    def shuffle_biom(self):
        new_cols = np.random.permutation(self.biom.columns)
        self.biom = self.biom.reindex_axis(new_cols,1)        
    


# In[4]:


    myotubiom = biom.load_table("./data/qiime_outputs/dada2_stoolsal_out/dada2_otu_table_w_tax_no_pynast_failures_rare7500.tsv")

    mybiom = myMiBiTimeSeriesbiom.merge(myotubiom)

    mymeta = pd.DataFrame.from_csv("/home/operon/Documents/lozupone_lab/dada2_stoolsal_out/full_maps_corrected.txt", sep = "\t")

    mytree = sktree.TreeNode.read("/home/operon/Documents/lozupone_lab/MiBiTimeSeries_bangladesh/new_prepped_tree.tre")

    rared = 7500

    meta_df, biom_df = meta_biom_filter(mymeta, mybiom, "F01", 1000)

    f01 = MiBiTimeSeries(meta_df, biom_df, mytree, 7500)


else:
    f01 = None

f01 = comm.bcast(f01, root=0)
# In[13]:

rank_combs = []
mycombs = f01.non_subclade
# generate all unique permutations of OTU pairs i and j.
for i in range(len(mycombs)):
    if (i % size) == (rank):
        rank_combs.append(mycombs[i])


#x220: 8 hours. 
import time
def permutation_tests(MiBiTimeSeries, i,j, precision):
    #"Calculate"
    MiBiTimeSeries_shuff = copy.copy(MiBiTimeSeries)
    randos = []
    for x in range(precision):
        MiBiTimeSeries_shuff.shuffle_biom()
        randos.append(MiBiTimeSeries_shuff.lag_regress(i,j)[2])
    mypercentile = stats.percentileofscore(randos, MiBiTimeSeries.lag_regress(i,j)[2])/100

    return mypercentile, MiBiTimeSeries.lag_regress(i,j)[2] , np.mean(randos)

perms = []
for i,j in rank_combs:
    percentile, r_coef, rand_r_coef = permutation_tests(f01, i,j, 50000)
    mymontes[(i,j)] = {'pvalue':p_value, 'r_coef':r_coef, 
                       'rand_r_coef':rand_r_coef}    
    
       # print 'mean time: ', np.mean(mytimes), 'max time: ', max(mytimes) ,'\n'
    perms.append(i,j, mymontes[(i,j)]['pvalue'], mymontes[(i,j)]['r_coef'], mymontes[(i,j)]['rand_r_coef'])


