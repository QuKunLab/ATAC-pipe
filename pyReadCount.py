import os
from optparse import OptionParser
import re
import pandas as pd
import numpy as np



def Macs2(inbed,gsize,outdir,name):
    os.system("macs2 callpeak -t %s -f BED -g %s --outdir %s -q 0.01 -n %s --nomodel --shift 0"
             %(inbed,gsize,outdir,name))
    return
def HQpeaks(inxls,insummits,inbed,pval,fval,qval,ubool,width,pipeup):
    xls=open(inxls,'r')
    summits=open(inbed,'r')

    outxls=open(inxls[:-4]+".HQ.xls",'w')
    outxls.write("chr\tstart\tend\tlength\tabs_summit\tpileup\t-10*log10(pvalue)\tfold_enrichment\t-10*log10(qvalue)\tname\n")
    outbed=open(inbed[:-11]+".HQ.bed",'w')
    outsummits=open(insummits[:-4]+".HQ.bed",'w')
    for line in xls:
        if re.search(r"^#",line) or re.search(r"^\n",line) or re.search(r"start|end",line):
            continue
        else:
            linesummits=summits.readline()
            items=line.rstrip('\n').split('\t')
            if (float(items[5])>float(pipeup) and float(items[6])>float(pval) and float(items[7])>float(fval) and float(items[8])>float(qval)):
                summit=int(items[4])
                if ubool:
                    start=int(items[1])
                    end=int(items[2])
                else:
                    width=int(width)
                    start=summit-width/2
                    end=summit+width/2
                linebed=items[0]+'\t'+str(start)+'\t'+str(end)+'\t'+items[9]+'\t'+items[8]+'\n'
                outbed.write(linebed)
                outxls.write(line)
                outsummits.write(linesummits)
    xls.close()
    summits.close()
    outxls.close()
    outbed.close()
    outsummits.close()
def rmBL(inbed,blbed,finalbed):
    os.system("bedtools intersect -a %s -b %s -v > %s"%(inbed,blbed,finalbed))
    return


def Add_name(inf,outf):
    named_outf=open(outf,'w')
    i=1
    for line in open(inf):
        data=line.split("\t")
        data.insert(3,"merge_peak_"+str(i))
        i=i+1
        line="\t".join(data)
        named_outf.write(line)
    named_outf.close()
    os.popen('rm ' + inf)

def mergePeak(inlist,pname):
    os.system("cat %s > %s"%(" ".join(inlist),pname+".combined.peak.list"))
    os.system("bedSort %s %s"%(pname+".combined.peak.list ",pname+".sorted.peak.list"))
    os.system("bedtools merge -i %s -c 1 -o count > %s"%(pname+".sorted.peak.list",pname+".tmp.sorted.peak.list"))
    Add_name(pname+".tmp.sorted.peak.list",pname+".merged.peak.list")
    os.system("rm %s"%(pname+".sorted.peak.list"))
    os.system("rm %s"%(pname+".combined.peak.list"))
    return

def Norm(inf, method='DESeq_norm'):
    data0=pd.read_table(inf,header=0,index_col=0)
    if method=="DESeq_norm":
        data=np.log2(data0/(data0.T/data0.apply(gmean,axis=1)).dropna(axis=1,how='any').apply(np.median,axis=1)+1)
    if method=="Quartile_norm":
        rank_mean = data0.stack().groupby(data0.rank(method='first').stack().astype(int)).mean()
        data=np.log2(data0.rank(method='min').stack().astype(int).map(rank_mean).unstack()+1)
    sumPerGroup=data0.apply(lambda x:x.sum())
    sizePerGroup=data0.apply(lambda x:(x>0).sum())
    if method=='Total_norm':
        data=np.log(data0/sumPerGroup*1000000+1)
    if method=="Size_norm":
        data=np.log2(data0*(np.median(sumPerGroup)/sumPerGroup)+1)
    if method=="Mean_norm":
        data=np.log(data0*(np.median(sumPerGroup/sizePerGroup)/(sumPerGroup/sizePerGroup))+1)
    return data


def countTable(inbed,sn,peaklist,outdir):
    os.system("bedtools intersect -a %s -b %s -wo > %s"%(inbed,peaklist,outdir+'/'+sn+'.count'))
    incount=open(outdir+'/'+sn+'.count','r').readlines()
    peakcount={}
    for line in incount:
        items=line.rstrip('\n').split('\t')
        peakcount[items[9]]=0

    for line in incount:
        items=line.rstrip('\n').split('\t')
        peakcount[items[9]]=peakcount[items[9]]+int(items[11])/50
    peakcountSum=open(outdir+'/'+sn+'.readcount','w')
    peakcountSum.write("PeakID\t"+sn+'\n')
    peak=open(peaklist,'r')
    for line in peak:
        items=line.rstrip('\n').split('\t')
        if peakcount.has_key(items[3]):
            count=int(peakcount[items[3]]+0.5)
        else:
            count=0
        peakcountSum.write(items[3]+'\t'+str(count)+'\n')
    peak.close()
    return

def mergecount(inpath,pname):
    peakSumCountList=[inpath+'/'+x for x in os.listdir(inpath) if x.split('.')[-1]=='readcount']
    df0=pd.read_table(peakSumCountList[0],sep='\t',index_col=0)
    df1=pd.DataFrame(index=df0.index)
    for i in range(len(peakSumCountList)):
        name=peakSumCountList[i].split('/')[-1].split(".")[0]
        df=pd.read_table(peakSumCountList[i],sep='\t',index_col=0,names=[name])
        df1[name]=df[name]
    df2=df1.loc[:,sorted(list(df1.columns))]
    df2.to_csv(inpath+'/'+pname+'.count',sep='\t')
    return
