#!/usr/bin/python
#
import numpy
import pandas
from optparse import OptionParser
import os
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
from multiprocessing import Pool
import sys
#
#
global motifposition, source_data, workspace, TFs_file, sample, index, gene_tss
#
motifposition = './Data/motifs'
#source_data = './source_data/'
workspace = sys.argv[3]
TFs_file = sys.argv[4]
inbam1=sys.argv[5]
inbam2=sys.argv[6]
#
#
def gz2dict(gz):
    sites, weighs = {}, {}
    gz = pandas.DataFrame.from_csv(gz, sep='\t')
    chroms, table = gz.index.values, gz.values
    for nindex in index:
        iindex = numpy.where(chroms==nindex)[0]
        position, odd = (table[iindex, 0]+table[iindex,1])//2, table[iindex,3]
        odd[numpy.where(odd>10)] = 10
        prob = (10.0**odd)/(1.0+10.0**odd)
        sites[nindex] = position
        weighs[nindex] = prob
    return sites, weighs
#
#
def gamma(x):
    g = numpy.ones(len(x))
    g[numpy.where(abs(x)>10000)] = 0
    return g
#
#
def score(gz):
    sites, weighs = gz2dict(gz)
    gene_scores = []
    for i in range(0, len(gene_tss)):
        position = sites[gene_tss[i][0]] - gene_tss[i][1]
        weight = gamma(position) * weighs[gene_tss[i][0]]
        gene_scores.append(weight.sum())
    print gz
    return gene_scores
#
#
###### initialize
#
ref=sys.argv[2]
if ref=='hg19':
    index = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY']

elif ref=='mm9':
    index = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
TSS="./Data/Ref/"+ref+'/'+ref+'_refseq_genes_TSS.txt'
gene_tss = [[x.split('\t')[0], int(x.split('\t')[1])] for x in open(TSS).readlines()]


TFs = [x[:-1] for x in open(TFs_file).readlines()]
TFnames = [x.split('.')[0] for x in os.listdir(motifposition) if x[-6:]=='.motif']
TFnames = list(set(TFs).intersection(set(TFnames)))
TFnames.sort()

#sam_names = [x.split('.')[0] for x in os.listdir(source_data) if x[-11:]=='perbase.sam']
#
###### calculate gene score and draw clustermap
#
gz_names1 = [workspace+tf+'_'+inbam1+'.txt.gz' for tf in TFnames]
gz_names2 = [workspace+tf+'_'+inbam2+'.txt.gz' for tf in TFnames]
pool = Pool(int(sys.argv[1]))
matrix1 = numpy.array(pool.map(score, gz_names1))
matrix2 = numpy.array(pool.map(score, gz_names2))
pool.close()
pool.join()
matrix = numpy.vstack((matrix1, matrix2))
numpy.save('tf_matrix.npy', matrix)
#
matrix = numpy.load('tf_matrix.npy')
corr = numpy.corrcoef(matrix)
corr = corr[len(TFnames):, :len(TFnames)]
corr_df = pandas.DataFrame(corr, index=TFs, columns=TFs)
#
seaborn.set_context('notebook', font_scale=1.2)
center = (corr.max()+corr.min())/2
fig1 = seaborn.clustermap(corr_df, cmap='bwr', center=center)
plt.setp(fig1.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.setp(fig1.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
fig1.ax_heatmap.set_xlabel(inbam1)
fig1.ax_heatmap.set_ylabel(inbam2)
fig1.savefig(workspace+'/'+inbam1+'_vs_'+inbam2+'.TF_network.pdf', bbox_inches='tight')
#
#
#
