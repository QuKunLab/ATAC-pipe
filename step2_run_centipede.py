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
import call_binding
import sys
#
#
global motifposition, source_data, workspace, TFs_file, sample_file
#
motifposition = './Data/motifs/'
#source_data = './source_data/'
workspace = sys.argv[4]
if not os.path.exists(workspace): os.makedirs(workspace)
TFs_file = sys.argv[5]
#
#
class options():
    def __init__(self, protocol, model, restarts, mintol, model_file, posterior_file, log_file,
                 window, batch, motif_file, bam_files, bam_file_genomicdna, seed):
        self.protocol = protocol
        self.model = model
        self.restarts = restarts
        self.mintol = mintol
        self.model_file = model_file
        self.posterior_file = posterior_file
        self.log_file = log_file
        self.window = window
        self.batch = batch
        self.motif_file = motif_file
        self.bam_files = bam_files
        self.bam_file_genomicdna = bam_file_genomicdna
        self.seed = seed
#
def runCentipede(tf):
#    sample_file = [source_data+sample+'.prepared.bam']
    sample = sample_file[0].split('/')[-1].split('.')[0]
    TFinfo = motifposition+'/'+tf+'.motif.'+ref+'.bed.txt.gz'
    model = workspace+'/'+tf+'_'+sample+'.pkl'
    log = workspace+'/'+tf+'_'+sample+'.log'
    posterior = workspace+'/'+tf+'_'+sample+'.txt.gz'
    opt = options('ATAC_seq', 'msCentipede', 1, 0.0001, model, posterior, log,
                  128, 5000, TFinfo, sample_file, None, None)
    call_binding.learn_model(opt)
    call_binding.infer_binding(opt)
    return
#
#
###### initialize
#
ref=sys.argv[1]

if ref=='hg19':
    index = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY']

elif ref=='mm9':
    index = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
TSS="./Data/Ref/"+ref+'/'+ref+'_refseq_genes_TSS.txt'
gene_tss = [[x.split('\t')[0], int(x.split('\t')[1])] for x in open(TSS).readlines()]

TFs = sys.argv[5:]
if os.path.exists(TFs_file): TFs = [x[:-1] for x in open(TFs_file).readlines()]
#if len(TFs)<1: TFs = [x[:-1] for x in open(TFs_file).readlines()]
TFnames = [x.split('.')[0] for x in os.listdir(motifposition) if x[-6:]=='.motif']
TFnames = list(set(TFs).intersection(set(TFnames)))
TFnames.sort()
print TFnames

#sam_names = [x.split('.')[0] for x in os.listdir(source_data) if x[-11:]=='perbase.sam']
#
#
###### run Centipede
#
#for sample in sam_names:


sample_file = [sys.argv[3]]
if len(TFnames)==1:
    runCentipede(TFnames[0])
else:
    pool = Pool(int(sys.argv[2]))
    pool.map(runCentipede, TFnames)
    pool.close()
    pool.join()
#
#
#
