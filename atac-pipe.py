import os
import sys
import re
from optparse import OptionParser
from optparse import OptionGroup
import numpy as np
import pandas as pd
import Levenshtein
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import random
from collections import Counter
from pyMappingQC import *
from pyReadCount import *
from pyCorrPCA import *
from pySig import *
from multiprocessing import cpu_count
from multiprocessing import Pool

#######  Options  #######

usage="usage: %prog [options][inputs]"
parser=OptionParser(usage=usage, version="%prog 1.0")

parser.add_option("-r",type='string',help="Set reference genome as hg19 or mm9")
parser.add_option("--group",help="Set group information")
parser.add_option("--project",default='atac',help="Set a name for the project,default='atac'.")
parser.add_option("-o",type="string",default='output',help="Set output directory")
parser.add_option("-t",type='int',default=4,help="Set the max thread,default = 4")

group=OptionGroup(parser,"Mapping Options")
group.add_option("-i",type='string',help="Set input directory of pair-end fastq files.")
group.add_option("-c",type='int',default=0,help="Intercept the first N bases, default=0 means undo interception")
group.add_option("--aq",type="int",default=20,help="Set query_length for adapter trimming,default = 20")
group.add_option("-a",type="string",default="CTGTCTCTTATACACATCTGACGCTGCCGACGA",help="Set adapter sequence,defalut='CTGTCTCTTATACACATCTGACGCTGCCGACGA'")
group.add_option("--MappingQC",action='store_true',default=False,help="ONLY do Mapping & QC, -i, -r, -o are also needed.")

parser.add_option_group(group)

group=OptionGroup(parser,"Merge Options")
group.add_option("--Merge",action='store_true',default=False,help="ONLY do Merge, -r, --group, -o are also needed.")
group.add_option("--bam",help="in bam directory.")

parser.add_option_group(group)

group=OptionGroup(parser,"Peak Calling Options")

group.add_option("--p1",default=2,help="-log10(pvalue),default=2")
group.add_option("--f1",default=1,help="fold enrichment,default=1")
group.add_option("--q1",default=5,help="-log10(qvalue),default=5")
group.add_option("--pipeup",default=10,help="pipeup,default=10")
group.add_option("-u",action="store_true",default=False,help="uncentrized peaks.-u and -w are mutually exclusive.")
group.add_option("-w",default=500,help="centralized peak width, defalut=500,-u and -w are mutually exclusive.")
group.add_option("--PeakCalling",action='store_true',default=False)
group.add_option("--bed",help="in bed directory")

parser.add_option_group(group)

group=OptionGroup(parser,"Significant Analysis")
group.add_option("--p2",default=0.01,help="pvalue,default=0.01")
group.add_option("--f2",default=1,help="log2FoldChange,default=1")
group.add_option("--q2",default=0.05,help="qvalue,default=0.05")
group.add_option("--icv",default=0.8,help="cv within group,default=0.8")
group.add_option("--SigAna",action='store_true',default=False)
group.add_option("--count",help="count file")
group.add_option("--peak",help="peak file")
parser.add_option_group(group)

group=OptionGroup(parser,"Search Motif Options")
group.add_option("--SearchMotif",action='store_true',default=False)
group.add_option("--cluster",help="input cluster list for motifs searching")
group.add_option("--bg",help="input background peak list")
group.add_option("--size",default=200,help="size for homer,default=200")
parser.add_option_group(group)


group=OptionGroup(parser,"Footprint Options")
group.add_option("--Footprint",action='store_true',default=False)
group.add_option("--motif",help="set motif")
group.add_option("--perbase",help="set perbase bed file")
group.add_option("--window",default=200,help="set window width,default=200")
group.add_option("--norm",default=5000000, help="set norm num,default=5000000")
parser.add_option_group(group)



(options, args) = parser.parse_args()

###########################
#  Set global parametres  #
###########################

#######  Set input PE FQ files  #######
if options.i:
    infq_1=[options.i.rstrip('/')+'/'+x for x in os.listdir(options.i) if x[-5:]=='_1.fq']
    infq_2=[options.i.rstrip('/')+'/'+x for x in os.listdir(options.i) if x[-5:]=='_2.fq']
    infqlist=map(list,zip(infq_1,infq_2))
    sn=list(set([x.split('/')[-1].split('_')[0] for x in infq_1+infq_2]))
    cpun=options.t/int(len(sn)) if options.t/len(sn)!=0 else 1



#######  Set Labels  #######
if options.group:
    df=pd.read_table(options.group,header=None)
    kind=[]
    for i in range(df.shape[1]-1):
        group={}
        gn=list(set(df[i+1]))
        gn.sort()
        for x in gn:
            group[x]=list(df[df[i+1]==x][0])
        kind.append(group)
    gn1=kind[0].keys()
    gn1.sort()
    sn=[y for x in list(kind[0].values()) for y in x]

#######  Set Output Dir  #######

outmap=options.o+'/'+'Mapping'
if options.Merge:
    outmap=options.bam
if options.PeakCalling:
    outmap=options.bed
outqc=options.o+'/'+'QC'
outgroup=options.o+'/'+'Group'
outpeak=options.o+'/'+'Peak'
outsig=options.o+'/'+"Sig"
outmotif=options.o+'/'+'Motif'
outfootprint=options.o+'/'+'Footprint'

if options.r:
#######  Set Ref index  #######
    ref_index="./Data/Ref/"+options.r+'/'+options.r
    print ('\tref_index='+ref_index if os.path.exists(ref_index+'.1.bt2') else '\t'+ref_index+' not exists!')

#######  Set Ref size  #######
    ref_size="./Data/Ref/"+options.r+'/'+options.r+'.sizes'
    print ('\tref_size='+ref_size if os.path.exists(ref_size) else '\t'+ref_size+' not exists!')

#######  Set TSS  #######
    TSS="./Data/Ref/"+options.r+'/'+options.r+'_refseq_genes_TSS.txt'
    print ('\tTSS='+TSS if os.path.exists(TSS) else '\t'+TSS+' not exists!')

print "\n"


###########################
#      Mapping & QC       #
###########################

#######   Def Functions for Mapping and QC  #######

def MappingQC_1(sn):
    fq_1=options.i+'/'+sn+'_1.fq'
    fq_2=options.i+'/'+sn+'_2.fq'
    mbase1=outmap+'/'+sn+'_1'
    mbase2=outmap+'/'+sn+'_2'
    mbase=outmap+'/'+sn
    qbase=outqc+'/'+sn
    print "\t\tPretrim for %s"%(sn)
    Pretrim(fq_1,mbase1+".trunc.fastq",options.c)
    Pretrim(fq_2,mbase2+".trunc.fastq",options.c)
    print "\t\ttrim adapters for %s"%(sn)
    trim_adapters(mbase1+".trunc.fastq",mbase2+".trunc.fastq",options.aq,options.a)
    print "\t\tMapping for %s"%(sn)
    Mapping(mbase1+".trim.fastq",mbase2+".trim.fastq",mbase+".sam", cpun, ref_index)
    mdbase=mbase+'.pe.q10.sort'
    print "\t\tDechrM for %s"%(sn)
    DechrM(mbase+".sam")
    mddbase=mdbase+'.rmdup'
    print "\t\tRemove PCR duplicate for %s"%(sn)
    Deduplicate(mdbase+'.bam')
    mddsbase=mddbase+'.shift'
    print "\t\tConvert bam to shifted bed for %s"%(sn)
    Bam2bedshift(mddbase+'.bam')
    print "\t\tConvert bed to bedGraph for %s"%(sn)
    Bed2bedGraph(mddsbase+".bed",ref_size)
    mddsnbase=mddsbase+".norm"
    print "\t\tNormalization by sequence depth for %s"%(sn)
    Normbedgraph(mddsbase+".bedGraph")
    print "\t\tConvert bedGraph to bw for %s"%(sn)
    Sortbedgraph(mddsnbase+'.bedGraph')
    adjust_bedgraph(mddsnbase+".bedGraph",ref_size)
    Bedgraph2bigwig(mddsnbase+".bedGraph",ref_size)
    return

def MappingQC_2(sn):
    mbase=outmap+'/'+sn
    print "\t\tConvert sam to bam for %s"%(sn)
    Sam2bam(mbase+".sam")
    return
def MappingQC_3(sn):
    mbase=outmap+'/'+sn
    print "\t\tConvert chrM.sam to chrM.bam for %s"%(sn)
    ChrM(mbase+'.sam')
    return
def MappingQC_4(sn):
    inbam=outmap+'/'+sn+'.pe.q10.sort.rmdup.bam'
    outfig=outqc+'/'+sn+'_4kb_TSSenrichment.pdf'
    print "\t\tDraw TSS enrichment plot for %s"%(sn)
    TSSEnrichmentplot(TSS, inbam, outfig)
    return
def MappingQC_5(sn):
    inbam=outmap+'/'+sn+'.pe.q10.sort.rmdup.bam'
    outfig=outqc+'/'+sn+'_Fragmentdistribution.pdf'
    print "\t\tDraw fragment distribution plot for %s"%(sn)
    Fragdistribution(inbam,outfig)
    return

def rm(files):
    for x in files:
        os.system('rm %s'%(x))

def MappingQC_6(sn):
    mbase=outmap+'/'+sn
    mdbase=mbase+'.pe.q10.sort'
    mddbase=mdbase+'.rmdup'
    mddsbase=mddbase+'.shift'
    Qc=outqc+'/'+sn+'.qctable'
    print "\t\tQuality control %s"%(sn)
    QC(Qc,sn,mbase+'.sam.map.log',mbase+'.chrM.bam',mddbase+'.bam',mddsbase+'.bed',ref_index+'.blacklist.bed',mbase+'.bam',mdbase+'.bam',mdbase+'.bam.Picard.log')
    rm((mbase+'.sam*',mbase+'.chrM.bam',mbase+'.bam',mdbase+'.bam',mdbase+'*Picard*',mbase+'*.fastq', mbase+'*bai',mddsbase+'*bedGraph'))
    return
def Perbase(sn):
    print "\t\tProcessing per1base files for %s"%(sn)
    inbam=outmap+'/'+sn+'.pe.q10.sort.rmdup.bam'
    Bam2bedshift(inbam,extend=0)
    inbed=outmap+'/'+sn+'.pe.q10.sort.rmdup.shift.per1base.bed'
    Bed2bedGraph(inbed,ref_size)
    outbedGraph=inbed[:-4]+'.bedGraph'
    PerbasebedGraph(outbedGraph)




def MappingQC():
    outmap=options.o+'/'+'Mapping'
    Mkdir(outmap)
    outqc=options.o+'/'+'QC'
    Mkdir(outqc)
    print "######## Step1 Mapping & QC ########\n"

    print '\n\n\t***************************'
    print '\t1.Sequence will be intercepted to: %d base length'%(options.c) if options.c != 0 else '\t1.Use orignal sequence'
    print '\t2.CPU for each sample is set as: %d'%(cpun)
    print '\t3.Reference genome is set as: %s'%(options.r)
    print '\t4.Output directory is set as: %s'%(outmap)
    print '\t5.Query length for adapter trimming is set as: %d'%(options.aq)
    print '\t6.Adapter sequence is set: %s'%(options.a)
    print '\t7.Input fastq file is/are: '
    for x in infqlist:
        print "\t\t%s: %s %s"%(x[0].split('/')[-1].split('.')[0][:-2],x[0],x[1])
    print '\t***************************\n'

    pool = Pool(len(sn))
    pool.map(MappingQC_1, sn)
    pool.close()

    pool = Pool(len(sn))
    pool.map(MappingQC_2, sn)
    pool.close()

    pool = Pool(len(sn))
    pool.map(MappingQC_3, sn)
    pool.close()

    pool = Pool(len(sn))
    pool.map(MappingQC_4, sn)
    pool.close()

    pool = Pool(len(sn))
    pool.map(MappingQC_5, sn)
    pool.close()

    pool = Pool(len(sn))
    pool.map(MappingQC_6, sn)
    pool.close()

    pool =Pool(len(sn))
    pool.map(Perbase,sn)
    pool.close()



###########################
#          Group          #
###########################

#######   Def Functions for Group  #######

def GetMergedBam(gname):
    print "\t\tGet merged bam file for %s"%(gname)
    if kind[0].has_key(gname):
        inbamlist=[outmap+'/'+ sn +'.pe.q10.sort.rmdup.bam' for sn in kind[0][gname]]
        ogname=outgroup+"/"+gname+'.pe.q10.sort.rmdup'
        os.system('samtools merge -1f %s %s'%(ogname+'.bam',' '.join(inbamlist)))
        Bam2bedshift(ogname+'.bam')
        Bed2bedGraph(ogname+".shift.bed",ref_size)
        Normbedgraph(ogname+".shift.bedGraph")
        Sortbedgraph(ogname+'.shift.bedGraph')
        adjust_bedgraph(ogname+".shift.norm.bedGraph",ref_size)
        os.system('sort -k1,1 -k2,2n %s > %s'%(ogname+'.shift.norm.bedGraph',ogname+'.shift.norm.bedGraph.tmp'))
        os.system('mv %s %s'%(ogname+'.shift.norm.bedGraph.tmp',ogname+'.shift.norm.bedGraph'))
        Bedgraph2bigwig(ogname+".shift.norm.bedGraph",ref_size)
    return
def MappingQC_4g(gn):
    inbam=outgroup+'/'+gn+'.pe.q10.sort.rmdup.bam'
    outfig=outgroup+'/'+gn+'_4kb_TSSenrichment.pdf'
    TSSEnrichmentplot(TSS, inbam, outfig)
    return
def MappingQC_5g(gn):
    inbam=outgroup+'/'+gn+'.pe.q10.sort.rmdup.bam'
    outfig=outgroup+'/'+gn+'_Fragmentdistribution.pdf'
    Fragdistribution(inbam,outfig)

def Perbase_g(gn):
    inbam=outgroup+'/'+gn+'.pe.q10.sort.rmdup.bam'
    Bam2bedshift(inbam,extend=0)
    sam2perbasebam(inbam)
    inbed=outgroup+'/'+gn+'.pe.q10.sort.rmdup.shift.per1base.bed'
    Bed2bedGraph(inbed,ref_size)
    outbedGraph=inbed[:-4]+'.bedGraph'
    PerbasebedGraph(outbedGraph)
    Normbedgraph(inbed[:-4]+'.bedGraph')
    Sortbedgraph(inbed[:-4]+'.norm.bedGraph')
    adjust_bedgraph(inbed[:-4]+'.norm.bedGraph',ref_size)
    os.system('sort -k1,1 -k2,2n %s > %s'%(inbed[:-4]+'.norm.bedGraph',inbed[:-4]+'.norm.bedGraph.tmp'))
    os.system('mv %s %s'%(inbed[:-4]+'.norm.bedGraph.tmp',inbed[:-4]+'.norm.bedGraph'))
    Bedgraph2bigwig(inbed[:-4]+'.norm.bedGraph',ref_size)
    return

def TSSReadHeatmap_g(gn):
    inbed=outgroup+'/'+gn+'.pe.q10.sort.rmdup.shift.per1base.bed'
    outfig=outgroup+'/'+gn+'_ReadDenstiyTSS.png'
    ReaddensityTSS(inbed,TSS,outfig)
    return

def Merge():
    Mkdir(outgroup)
    print "######## Step2 Merging samples by groups ########\n"
    print '\n\n\t***************************'
    print "\t1.Samples are groups as:\n"
    print '\t\t',kind[0]
    print '\t2.Reference genome is set as: %s'%(options.r)
    print '\t3.Output directory is set as: %s'%(outgroup)
    print '\n\n\t***************************'
    pool = Pool(options.t)
    pool.map(GetMergedBam, gn1)
    pool.close()

    pool = Pool(options.t)
    pool.map(MappingQC_4g, gn1)
    pool.close()

    pool = Pool(options.t)
    pool.map(MappingQC_5g, gn1)
    pool.close()

    pool=Pool(options.t)
    pool.map(Perbase_g,gn1)
    pool.close()

    pool=Pool(options.t)
    pool.map(TSSReadHeatmap_g,gn1)
    pool.close()


###########################
#Call Peaks and Read Count#
###########################

#######   Def Functions for Peak Calling and Read Count #######
if options.r:
    gsize='hs' if options.r=='hg19' else 'mm' if options.r=='mm9' else 'none'
def callgroup(gn):
    if kind[0].has_key(gn):
        inbedlist=[outmap+'/'+ sn +'.pe.q10.sort.rmdup.shift.bed' for sn in kind[0][gn]]
        ogname=outpeak+"/"+gn
        print "Call Peaks for %s"%(gn)
        os.system("cat %s > %s"%(' '.join(inbedlist),ogname+".combined.bed"))
        os.system("bedSort %s %s"%(ogname+".combined.bed",ogname+".sort.bed"))
        Macs2(ogname+".sort.bed",gsize,outpeak,gn)
        xls_dir=ogname+"_peaks.xls"
        summits_dir=ogname+'_summits.bed'
        bed_dir=ogname+"_peaks.narrowPeak"
        HQpeaks(xls_dir,summits_dir,bed_dir,options.p1,options.f1,options.q1,options.u,options.w,options.pipeup)
        hqbed=ogname+"_peaks.HQ.bed"
        hqrmblbed=ogname+'_peaks.final.bed'
        bl=ref_index+'.blacklist.bed'
        rmBL(hqbed,bl,hqrmblbed)
    return


def GetMergedPeakList():
    mergePeak([outpeak+'/'+x for x in os.listdir(outpeak) if x[-10:]=='.final.bed'],outpeak+'/'+options.project)
    return

def GetCount(sn):
    print "Get Count for %s"%(sn)
    inbed=outmap+'/'+sn+'.pe.q10.sort.rmdup.shift.bed'
    peaklist=outpeak+'/'+options.project+'.merged.peak.list'
    countTable(inbed,sn,peaklist,outpeak)
    return

def GetMergedCount():
    mergecount(outpeak,options.project)
    return

def Read():
    Mkdir(outpeak)
    print "######## Step3 Peak calling & Read count ########\n"
    print '\n\n\t***************************'
    print ('\t1.Use merged files for peak calling.')
    print '\t2.gsize is set as: %s'%(gsize)
    print '\t3.Output directory is set as: %s'%(outpeak)
    print ('\t4.High quality peaks filter: fdc > %f, -log(p) > %f,-log(q) >  %f, pipeup > %f'%(float(options.f1),float(options.p1),float(options.q1),float(options.pipeup)))
    print '\t5.Peak will be centered to summit and extend to %d in both sides'%(int((int(options.w)/2))) if options.u==False else '\t5.Use original peaks'
    print '\n\n\t***************************'
    print "\t Call peaks:"
    pool = Pool(int(options.t))
    pool.map(callgroup, gn1)
    pool.close()

    print "\tGet merged peak list"
    GetMergedPeakList()

    print "\tRead count:"
    pool = Pool(int(options.t))
    pool.map(GetCount, sn)
    pool.close()

    print "\tGet Merged Count table"
    GetMergedCount()


#######################################
#   Significant Analysis with DESeq   #
#######################################

#######   Def Functions for Significant analysis #######

def CallSig():
    Mkdir(outsig)
    print "######## Step4 Get Significant peaks ########\n"
    label=options.group
    if options.SigAna:
        count_file=options.count
        opname=outsig+'/'+count_file.split('/')[-1].split('.')[0]
        norm_count_file=outsig+'/'+count_file.split('/')[-1]+'.log.norm'
        inpeak=options.peak
    else:
        count_file=outpeak+'/'+options.project+'.count'
        opname=outsig+'/'+options.project
        norm_count_file=outsig+'/'+options.project+'.count.log.norm'
        inpeak=outpeak+'/'+options.project+'.merged.peak.list'
    print '\n\n\t***************************'
    print ('\t1.Use DESeq for Significant analysis.')
    print '\t2.Peaks with cv > %f will be set as blacklist'%(float(options.icv))
    print ('\t3.Significant peaks filter: |log2fdc| > %f, p < %f,(q) <  %f'%(float(options.f2),float(options.p2),float(options.q2)))
    print '\t4.Output directory is set as: %s'%(outsig)

    print "\t Caculating Significance with DESeq"
    os.system("Rscript ./deseq.r %s %s %s %s"%(count_file,label,outsig+'/'+count_file.split('/')[-1],outsig))
    count_norm=outsig+'/'+count_file.split('/')[-1]+'.log.norm'
    count_norm_df=pd.read_table(count_norm,header=0,index_col=0)
    CorrFig(count_norm_df,count_norm)
    myPCA(count_norm_df,label,count_norm)

    comlist=[outsig+'/'+x for x in os.listdir(outsig) if x[-4:]=='.com']
    data = pd.read_table(norm_count_file,index_col=0)
    label= pd.Series.from_csv(options.group,index_col=0,sep="\t")[data.columns]
    peaklist=pd.read_table(inpeak,index_col=3,header=None,names=('Chr','Start','End','PeakID','No'))

    print "\t Remove peaks with variance in group"
    BlPeaks=BLpeaks(data,label,outsig,icv=options.icv)
    BlData=PeakID2Data(BlPeaks,data)
    BlData.to_csv(opname+'.BLdata',sep='\t',index=True,header=True)
    BlBed=PeakID2SE(BlPeaks,peaklist)
    BlBed.to_csv(opname+'.BLbed',sep='\t',index=False,header=False)

    SigPeaks=calculateSig(comlist,outsig,f=options.f2,p=options.p2,q=options.q2)
    FinalSigPeaks=rmBLPeaks(SigPeaks, BlPeaks)

    print "\t Get Significant peaks"
    SigData=PeakID2Data(FinalSigPeaks,data)
    SigData.to_csv(opname+'.sigdata',sep='\t',index=True,header=True)
    SigBed=PeakID2SE(FinalSigPeaks,peaklist)
    SigBed.to_csv(opname+'.sigbed',sep='\t',index=False,header=False)

########################
#   Search Motif       #
########################

def SearchMotif():
    Mkdir(outmotif)
    print '\n\n\t***************************'
    print ('\tUse HOMER for Motif Searching.')
    print '\t2.Cluster for motif search is set as :%s'%(options.cluster)
    print ('\t3.Background is set as :%s'%(options.bg))
    print '\t4.Output directory is set as: %s'%(outmotif)
    print '\n\n\t***************************\n'

    cluster=pd.read_table(options.cluster,index_col=0,header=None)
    inpeak=options.peak
    peaklist=pd.read_table(inpeak,index_col=3,header=None,names=('Chr','Start','End','PeakID','No'))
    clsbed=PeakID2SE(cluster.index,peaklist)
    clsbed.to_csv(options.cluster+'.bed',sep='\t',header=False,index=False)
    os.system("findMotifsGenome.pl %s %s %s -size %d -mask"%(options.cluster+'.bed',options.r,outmotif,int(options.size)))
    if options.bg:
        bg=pd.read_table(options.bg,index_col=0,header=None)
        bgbed=PeakID2SE(bg.index,peaklist)
        bgbed.to_csv(options.bg+'.bed',sep='\t',header=False,index=False)
        os.system("findMotifsGenome.pl %s %s %s -size %d -mask -bg %s"%(options.cluster+'.bed',options.r,outmotif,int(options.size),options.bg+'.bed'))
    else:
        os.system("findMotifsGenome.pl %s %s %s -size %d -mask"%(options.cluster+'.bed',options.r,outmotif,int(options.size)))


def Drawfootprint(incsv,motifname,halfwidth):
    lines = open(incsv).readlines()
    length=len(lines[0].split('\t'))
    counts_matrix, genes = [], []
    for iline,line in enumerate(lines):
        if iline>0:
            words = line.split('\t')
            counts_matrix.append(map(float, words[1:]))
            genes.append(words[0])
    counts_matrix = np.asarray(counts_matrix)
    n_ave, n_step, cut_off = 50, 10, 0.7
    n_limit = (len(counts_matrix) - n_ave) // n_step
    ave_matrix = [counts_matrix[i*n_step:i*n_step+n_ave, :].sum(axis=0)/float(n_ave) for i in range(0, n_limit)]
    ave_matrix = np.asarray(ave_matrix)
    log_ave_matrix = np.log(ave_matrix+1)
    max_logAve = log_ave_matrix.max()
    log_ave_matrix[np.where(log_ave_matrix > max_logAve*cut_off)] = max_logAve*cut_off
    nt_label = [' '] * (int(length)-1)
    nt_label[0], nt_label[int(len(nt_label)/2)],nt_label[-1] = '-'+str(halfwidth),motifname,str(halfwidth)
    log_ave_df = pd.DataFrame(log_ave_matrix, index=xrange(0, len(log_ave_matrix)), columns=nt_label)
    sns.set_context('poster', font_scale=1)
    fig=plt.figure(figsize=(6, 14))
    ax11 = plt.subplot2grid((60, 1), (14, 0), rowspan=46)
    sns.heatmap(log_ave_df, xticklabels=True, yticklabels=False, cbar=False, ax=ax11)
    plt.setp(ax11, ylabel=motifname+' motif-sites')
    ax12 = plt.subplot2grid((60, 1), (0, 0), rowspan=9)
    ax12.plot(xrange(-(halfwidth), halfwidth+1), log_ave_df.sum(axis=0), color='firebrick')
    ax12.set_xlim([-(halfwidth), (halfwidth)])
    ax12.set_xticks([-halfwidth, 0, halfwidth])
    ax12.set_xlabel('Distance to motif')
    ax12.set_ylabel('Insertion-site probability')
    ax12.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    fig.savefig(incsv+'.png')
    plt.close(fig)

def Vplot(bed,motifname,peak):
    ref=options.r
    temp1=outfootprint+'/'+peak.split('/')[-1]+'.include.'+motifname
    imname='./Data/motifs/'+motifname+'.motif'
    os.system('bedtools intersect -a %s -b %s -wo > %s'%(peak,imname+'.'+ref+'.bed',temp1))
    inbed=open(bed,'r')
    if  not os.path.exists(bed+'.fragment'):
        outbed=open(bed+'.fragment','w')
        fragid=[]
        for line in inbed:
            items=line.split('\t')
            fchr,fstart,fend,fid,flen=items[0],items[1],str(int(items[1])+int(items[4])),items[3],items[4]
            if fid not in fragid:
                fragid.append(fid)
                outbed.write(fchr+'\t'+fstart+'\t'+fend+'\t'+flen+'\n')
        outbed.close()
    temp2=temp1+'.in.'+bed.split('/')[-1].split('.')[0]+'.bed'

    os.system('bedtools intersect -a %s -b %s -wo > %s'%(bed+'.fragment',temp1,temp2))

    outbed=open(temp2,'r')
    window=500
    max_len=500
    data_matrix=np.zeros((max_len,2*window+1))
    for line in outbed:
        items=line.split('\t')
        mpos=int((int(items[10])+int(items[11]))/2)
        dist=int((int(items[1])+int(items[2]))/2)-mpos
        flen=int(items[3])
        if (abs(dist)<window) & (flen<max_len):
            data_matrix[flen][dist+window]+=1
    norm=5000000

    depth=int(os.popen('wc -l %s'%(bed+'.fragment')).read().split()[0])

    data_matrix=data_matrix*int(norm)/depth
    vplot_data = np.log10(data_matrix + 1)[::-1]
    sns.set_context('poster', font_scale=1)
    fig1 = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot2grid((100, 100), (20, 0), colspan=80, rowspan=80)
    xlabels, ylabels = [' '] * 1000, [' '] * len(vplot_data.sum(axis=1))
    xlabels[100], xlabels[300], xlabels[500], xlabels[700], xlabels[900] = \
    '-400', '-200', motifname, '200', '400'
    ylabels[-1], ylabels[-200], ylabels[-399] = '0', '200', '400'
    sns.heatmap(vplot_data, cmap='Blues', cbar=False,
                    xticklabels=xlabels, yticklabels=ylabels, ax=ax1)
# ax1.imshow(vplot_data, cmap='Blues')
    ax1.set_ylabel('Fragment length')
#
    ax2 = plt.subplot2grid((100, 100), (0, 0), colspan=80, rowspan=20)
    ax2.plot(np.arange(len(data_matrix.sum(axis=0))), data_matrix.sum(axis=0))
    ax2.set_xlim([0, 1000])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax3 = plt.subplot2grid((100, 100), (20, 80), colspan=20, rowspan=80)
    ax3.plot(data_matrix.sum(axis=1), np.arange(len(data_matrix.sum(axis=1))))
    ax3.set_ylim([0, 400])
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    fig1.savefig(temp2[:-4]+'.png')
    plt.close(fig1)



def Footprint():
    peaklist=options.peak
    ref=options.r
    inmotif=options.motif
    Mkdir(outfootprint)
    bg=options.perbase[:-4]+'.bedGraph'
    if not os.path.exists(bg):
        Bed2bedGraph(options.perbase,ref_size)
        PerbasebedGraph(bg)
    imname='./Data/motifs/'+inmotif+'.motif'
    mlen=len(open(imname).readlines())-1
    window=int(options.window)
    omname=outfootprint+'/'+inmotif
    peakname=options.peak.split('/')[-1]
    opeakname=outfootprint+'/'+peakname
    outtxt=opeakname+'.'+inmotif+'_'+bg.split('/')[-1]+'.footprint'
    if (not os.path.exists(imname+'.'+ref+'.bed') or os.path.getsize(imname+'.'+ref+'.bed')==0):
        os.system("scanMotifGenomeWide.pl "+imname+' '+ref+' -bed > '+imname+'.'+ref+'.bed')
    if (not os.path.exists(opeakname+'.'+inmotif+'.in') or os.path.getsize(opeakname+'.'+inmotif+'.in')==0):
        os.system("bedtools intersect -a "+imname+'.'+ref+'.bed'+' -b '+peaklist+' -wo > '+opeakname+'.'+inmotif+'.in')

    inlist=pd.read_table(opeakname+'.'+inmotif+'.in',header=None)
    tf_all=np.asarray((list(inlist.loc[:,0]),map(int,((inlist.loc[:,1]+inlist.loc[:,2])/2)))).T

    chrom={}
    for x in set(tf_all.T[0]):
        chrom[x]=[]
    for i in range(len(tf_all)):
        if chrom.has_key(tf_all[i][0]):
            chrom[tf_all[i][0]].append([tf_all[i][1],[0]*(2*window+1)])

    for site in chrom.keys():
        chrom[site]=dict(chrom[site])

    mysum=0
    inbg=open(bg)
    for line in inbg:
        it=line.rstrip('\n').split('\t')
        mysum+=float(it[3])
        if chrom.has_key(it[0]):
            for site in chrom[it[0]].keys():
                dis=int(it[1])-int(site)
                if abs(dis) <=window:
                    chrom[it[0]][site][dis+window]+=float(it[3])

    df=pd.DataFrame(index=range(-window,window+1))
    for k1 in chrom.keys():
        for k2 in chrom[k1].keys():
            kname=k1+':'+str(int(k2)-int(mlen/2))+'-'+str(int(k2)+int(mlen/2))
            df[kname]=chrom[k1][k2]
    chrom=None
    df=df.T
    df['sum']=df.sum(axis=1)
    df=df.sort_values(by=["sum"],ascending=False).drop(['sum'],axis=1)
    if int(options.norm)!=0:
        df=df*int(options.norm)/mysum

    df.to_csv(outtxt,sep='\t',index=True, header=True)
    df.sum().to_csv(outtxt+'.sum',sep='\t',index=True)
    Drawfootprint(outtxt,inmotif,window)
    Vplot(options.perbase,inmotif,options.peak)
    return




def main():
    if options.MappingQC: MappingQC()
    if options.Merge: Merge()
    if options.PeakCalling: Read()
    if options.SigAna: CallSig()
    if options.SearchMotif: SearchMotif()
    if options.Footprint: Footprint()
    if (options.MappingQC==False)& (options.Merge==False)&(options.PeakCalling==False)&(options.SigAna==False)&(options.SearchMotif==False)&(options.Footprint==False):
        MappingQC()
        Merge()
        Read()
        CallSig()


main()
