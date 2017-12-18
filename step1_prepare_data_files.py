



from optparse import OptionParser
import numpy
import pandas
import os
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
from multiprocessing import Pool
#
#
usage="usage: %prog [options][inputs]"
parser=OptionParser(usage=usage, version="%prog 1.0")

parser.add_option("-r",type='string',help="Set reference genome as hg19 or mm9")
parser.add_option("--tf",type='string',help="Set reference genome as hg19 or mm9")

(options, args) = parser.parse_args()
TFs_file=options.tf
motifposition='./Data/motifs'


def bed2txt():
    if options.r=='hg19':
        index = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY']
    elif options.r=='mm9':
        index = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
    TFs = [x[:-1] for x in open(TFs_file).readlines()]
    TFnames = [x.split('.')[0] for x in os.listdir(motifposition) if x[-6:]=='.motif']
    TFnames = list(set(TFs).intersection(set(TFnames)))
    TFnames.sort()
    bedfiles = [motifposition+'/'+x+'.motif.'+options.r+'.bed' for x in TFnames]
    for x in bedfiles:
        if os.path.exists(x)==False:
            os.system("scanMotifGenomeWide.pl "+motifposition+'/'+x.split('/')[-1].split('.')[0]+'.motif'+' '+options.r+' -bed > '+x)
    for bed in bedfiles:
        outfile ='./Data/motifs/'+ bed.split('/')[-1]+'.txt'
        if os.path.exists(outfile+'.gz')==False:
            with open(outfile,'w') as out, open(bed) as bedfile:
                print >> out, 'Chr'+'\t'+'Start'+'\t'+'Stop'+'\t'+'Strand'+'\t'+'PwmScore'
                for line in bedfile:
                    words = line[:-1].split('\t')
                    if words[0] in index:
                        print >> out, words[0]+'\t'+words[1]+'\t'+words[2]+'\t'+words[5]+'\t'+words[4]
            os.popen('gzip '+outfile)
    return
#
bed2txt()
