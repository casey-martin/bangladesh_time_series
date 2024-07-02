
import biom
import pandas as pd
import numpy as np
import csv
from sklearn import linear_model
from scipy import stats
from scipy.interpolate import PchipInterpolator
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
        self.tree = tree
        self.rared = rared
        self.nodes_thresh = {}
        self.all_regress = {}
        self.all_regress_thresh = {}
        
        self.biom = biom_df.loc[list(set(biom_df.index) & set([i.name for i in self.tree.traverse()]))]
        self.nodes = list(self.biom.index)
        self.epsilon = 0.000001

    
    def lag_regress(self, i, j):
        # takes 1ms to run.
        # switch to np arrays gives .6ms runtime
        j_relabund = (np.array(self.biom.loc[j]) + self.epsilon)/self.rared
        j_term = np.log10(j_relabund)[1:]
        i_relabund = (np.array(self.biom.loc[i]) + self.epsilon)/self.rared
        delta_i_term = (np.log10(i_relabund)[1:] - np.log10(np.roll(i_relabund, -1)[1:] + self.epsilon))

        #reg = LinearRegression()
        #reg.fit(j_term, delta_i_term)
        
        #return([i, j, reg.coef_[0][0], reg.p[0][0]])
        myregression = stats.pearsonr(j_term, delta_i_term)
        return(i, j, myregression[0], myregression[1])

    def lag_regress_non_subclade(self, nodes):
        i = nodes[0]
        j = nodes[1]
        
        node_i = self.tree.find(i)
        children_i = [child.name for child in node_i.traverse()]
        children_i.remove(i)
        
        node_j = self.tree.find(j)
        children_j = [child.name for child in node_j.traverse()]
        children_j.remove(j)

        if (i not in children_j) & (j not in children_i):
            try:
                return(self.lag_regress(i, j))
            except KeyError:
                return(None)

    def all_non_subclade(self):
        def subclade_test(nodes):
            i = nodes[0]
            j = nodes[1]

            node_i = self.tree.find(i)
            children_i = [child.name for child in node_i.traverse()]
            children_i.remove(i)

            node_j = self.tree.find(j)
            children_j = [child.name for child in node_j.traverse()]
            children_j.remove(j)

            if (i not in children_j) & (j not in children_i):
                try:
                    return(i, j)
                except KeyError:
                    return(None)
        tmp_nodes = [subclade_test([i,j]) for i,j in product(self.nodes,self.nodes)]
        self.non_subclade = [i for i in tmp_nodes if i is not None]    
    
                
    def all_regress_non_subclade2(self):
                
        tmp_regress = [self.lag_regress_non_subclade([i,j]) for i,j in product(self.nodes, self.nodes)]
        
        tmp_regress = [i for i in tmp_regress if i is not None]
        for i in tmp_regress:
            self.all_regress[i[0], i[1]] = {'slope':i[2], 'p_value':i[3]}
    
        
    
    def conserved_interact_rand(self):
        all_correls = pd.DataFrame.from_dict(self.all_regress, orient='index')
        subset_correls = all_correls.sort_values(['p_value'])
        
        x = 0
        drop_set = set()
        while x < len(subset_correls):
            loop_begin = time.time()
            i,j = subset_correls.index[x][0], subset_correls.index[x][1]
            if (i,j) in drop_set:
                x += 1
                continue
            else: 
                i_lineage = [claatu.name for claatu in self.tree.find(i).ancestors()]
                i_lineage += [claatu.name for claatu in self.tree.find(i).traverse()]
                j_lineage = [claatu.name for claatu in self.tree.find(j).ancestors()]
                j_lineage += [claatu.name for claatu in self.tree.find(j).traverse()]
                common_ancest = set(i_lineage) & set(j_lineage)

                for node in common_ancest:
                    i_lineage.remove(node)
                    j_lineage.remove(node)

                to_drop = set(product(i_lineage, j_lineage))
                to_drop.discard((i, j))

                for i in to_drop:
                    drop_set.add(i)
                x += 1
        
        drop_set = drop_set & set(subset_correls.index)
        self.top_regress = subset_correls.drop(drop_set, errors = 'ignore')
        
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
        self.biom = self.biom.reindex_axis(new_cols,1)


myclaatubiom = biom.load_table("/home/operon/Documents/lozupone_lab/claatu_bangladesh/claatu_table.tsv")

myotubiom = biom.load_table("/home/operon/Documents/lozupone_lab/dada2_stoolsal_out/dada2_otu_table_w_tax_no_pynast_failures_rare7500.tsv")

mybiom = myclaatubiom.merge(myotubiom)

mymeta = pd.DataFrame.from_csv("/home/operon/Documents/lozupone_lab/dada2_stoolsal_out/full_maps_corrected.txt", sep = "\t")

mytree = sktree.TreeNode.read("/home/operon/Documents/lozupone_lab/claatu_bangladesh/new_prepped_tree.tre")

rared = 7500

meta_df, biom_df = meta_biom_filter(mymeta, mybiom, "F01", 1000)

f01 = claatu(meta_df, biom_df, mytree, 7500)
f01.all_non_subclade()

rank_combs = []
mycombs = f01.non_subclade
for i in xrange(len(mycombs)):
    if (i % size) == (rank):
        rank_combs.append(mycombs[i])


f01_shuff = f01
#mypercentiles = {}
import time
for i,j in rank_combs:
    t0 = time.time()
    randos = []
    for k in xrange(1000):
        f01_shuff.shuffle_biom()
        randos.append(abs(f01_shuff.lag_regress(i,j)[2]))
    #mypercentiles[(i,j)] = 1 - stats.percentileofscore(randos, abs(f01.lag_regress(i,j)[2]))/100
    randos.append(abs(f01.lag_regress(i,j)[2]))
    p_values = stats.norm.sf(abs(stats.zscore(randos)))
    #normal_p = stats.mstats.normaltest(randos)[1]
    t1 = time.time()
    print t1 - t0
    print i + '\t' + j + '\t' + str(p_values[-1])  + '\t' + str(name) + '_' +  str(rank) + '\n'
