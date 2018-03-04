import os
import sys
import numpy as np
import Levenshtein
import pysam
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
from collections import Counter
import re
import seaborn as sns



#Creat new dir
#************************************************#
def Mkdir(path):
    path=path.strip()
    path=path.rstrip('\\')
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

#Pretrim
#************************************************#
def Pretrim(Infastq,Outfastq,N):
    inf=open(Infastq,'r')
    outf=open(Outfastq,'w')
    lines=inf.readlines()
    if (N==1000):
        os.popen("cp "+Infastq+" "+Outfastq)
    else:
        i=0
        while i<len(lines):
            outf.write(lines[i])
            if len(lines[i+1].rstrip())<=N:
                outf.write(lines[i+1])
            else:           
                outf.write(lines[i+1][0:N])
                outf.write('\n')
            outf.write(lines[i+2])
            if len(lines[i+3].rstrip())<=N:
                outf.write(lines[i+3])
            else:
                outf.write(lines[i+3][0:N])
                outf.write('\n')
            i=i+4
    return

#Adapter Triming
#************************************************#
def mismatch_align(seq1, query_length, read2_rc):
    for s1 in range(len(seq1)-query_length+1, -1, -1):
        temp_read1 = seq1[s1:(s1+query_length)]
        editdist = Levenshtein.distance(temp_read1, read2_rc)
        if editdist<2:
            return s1
    return -1

def rev_comp_dna(read2_rc):
    temp_read2 = ''
    for i in range(len(read2_rc)-1, -1, -1):
        if (read2_rc[i]=='A') | (read2_rc[i]=='a') :
            temp_read2 += 'T'
        elif (read2_rc[i]=='C') | (read2_rc[i]=='c') :
            temp_read2 += 'G'
        elif (read2_rc[i]=='G') | (read2_rc[i]=='g') :
            temp_read2 += 'C'
        elif (read2_rc[i]=='T') | (read2_rc[i]=='t') :
            temp_read2 += 'A'
        elif read2_rc[i]=='N':
            temp_read2 += 'N'
        else:
            return 'error'
    return temp_read2

def trim_adapters(fastq1, fastq2, query_length, adapter_seq):
    trimed1, trimed2 = '.'.join(fastq1.split('.')[:-2])+'.trim.fastq', '.'.join(fastq2.split('.')[:-2])+'.trim.fastq'
    with open(fastq1) as fa1, open(fastq2) as fa2, open(trimed1, 'w') as out1, open(trimed2, 'w') as out2 :
        nReads, mm0_num_read, mm1_num_read = 0, 0, 0
        while 1:
            seq_header1, seq_header2 = fa1.readline()[:-1], fa2.readline()[:-1]
            seq1, seq2 = fa1.readline()[:-1], fa2.readline()[:-1]
            qual_header1, qual_header2 = fa1.readline()[:-1], fa2.readline()[:-1]
            qual1, qual2 = fa1.readline()[:-1], fa2.readline()[:-1]
            nReads += 1
            if ((not seq_header1) | (not seq_header2) | (not seq1) | (not seq2) |
                (not qual_header1) | (not qual_header2) | (not qual1) | (not qual2)): break
            read2_rc = seq2[:query_length]
            read2_rc = rev_comp_dna(read2_rc)
            s1_pos = -1
            s1_pos_find = seq1.rfind(read2_rc)
            if s1_pos_find > 0 :
                s1_pos = s1_pos_find
                mm0_num_read += 1
            else:
                s1_pos = mismatch_align(seq1, query_length, read2_rc)
                if s1_pos>0: mm1_num_read += 1
            if s1_pos >= 0 :
                seq_len = s1_pos + query_length
                trim_seq1 = seq1[seq_len:]
                adapter_trim_seq = adapter_seq[:len(trim_seq1)]
                if adapter_trim_seq==trim_seq1:
                    seq1 = seq1[:seq_len]
                    seq2 = seq2[:seq_len]
                    qual1 = qual1[:seq_len]
                    qual2 = qual2[:seq_len]
            print >> out1, seq_header1
            print >> out1, seq1
            print >> out1, qual_header1
            print >> out1, qual1
            print >> out2, seq_header2
            print >> out2, seq2
            print >> out2, qual_header2
            print >> out2, qual2
    return nReads, mm0_num_read, mm1_num_read

#Mapping
#************************************************#
def Mapping(Infastq1,Infastq2,Outsam,N,ref_index):
    os.system("bowtie2 -p %d --very-sensitive -x %s -1 %s -2 %s -S %s 2>%s.map.log" %(N,ref_index,Infastq1,Infastq2,Outsam,Outsam))
    return

#DechrM
#************************************************#
def DechrM(Insam):
    DechrMbam=Insam[0:-4]+'.pe.q10.sort.bam'
    os.system("awk \'$3!=\"chrM\"\' %s |samtools view -S -b -f 0x2 -q 10 - |samtools sort -o %s" %(Insam,DechrMbam))
    return

#Sam2bam
#************************************************#
def Sam2bam(Insam):
    Outbam=Insam[0:-4]+'.bam'
    os.system("samtools view -Sb %s > %s"%(Insam,Outbam))
    return

#ChrM
#************************************************#
def ChrM(Insam):
    chrMbam=Insam[0:-4]+'.chrM.bam'
    os.system("awk \'$3==\"chrM\"|| NF<10' %s |samtools view -S -b - > %s" %(Insam,chrMbam))
    return

#Deduplicate
#************************************************#
def Deduplicate(Inbam):
    Dedupbam=Inbam[0:-4]+'.rmdup.bam'
    os.system("java -jar ./picard.jar MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s.Picard_Metrics_unfiltered_bam.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true &> %s.Picard.log" %(Inbam,Dedupbam,Inbam,Inbam))
    return
#sam2perbase1bpbam
def sam2perbasebam(inbam,ref_size):
    if os.path.exists(inbam[:-4]+'.per1base.bam')==False:
        os.system('samtools view -h %s >%s'%(inbam, inbam[:-4]+'.sam'))
        os.system('perl shift_sam_bases.pl %s %s %s'%(ref_size,inbam[:-4]+'.sam',inbam[:-4]+'.tmp.sam'))
        os.system('perl sam2perbase1bp.pl %s %s %s'%(ref_size,inbam[:-4]+'.tmp.sam',inbam[:-4]+'.tmp2.sam'))
        header_sam = inbam[:-4]+'.header.sam'
        os.popen('cat '+'./header.sam '+inbam[:-4]+'.tmp2.sam'+' > '+header_sam)
        bam = inbam[:-4]+'.header.bam'
        os.popen('samtools view -bS '+header_sam+' > '+bam)
        sorted_bam = inbam[:-4]+'.per1base.bam'
        os.popen('samtools sort '+bam+' > '+sorted_bam)
        os.popen('samtools index '+sorted_bam)
        os.popen('rm '+header_sam+' '+bam+' '+inbam[:-4]+'.sam'+' '+inbam[:-4]+'.tmp.sam'+' '+inbam[:-4]+'.tmp2.sam')

    return
#Bam2bedshift
#************************************************#
def Bam2bedshift(Inbam,extend=25):
    Outbed=Inbam[:-4]+'.shift.bed' if extend==25 else Inbam[:-4]+'.shift.per'+str(int(extend)*2+1)+'base.bed'
    bed=open(Outbed,'w')
    sam=os.popen("samtools view -h "+Inbam)
    for line in sam:
        if line[0]=='@':
            continue
        else:
            items=line.rstrip('\n').split('\t')
            if (int(items[1])==99 or int(items[1])==163):
                start=int(items[3])
                end=start+abs(int(items[8]))
                start_extl=start+4-int(extend)
                start_extr=start+4+int(extend) if int(extend)!=0 else start+4+int(extend)+1
                end_extl=end-5-int(extend)
                end_extr=end-5+int(extend) if int(extend)!=0 else end-5+int(extend)+1
                if start_extl<0: start_extl=0
                if end_extl<0: end_extl=0
                if int(items[1])==99:
                    bed.write('\t'.join([items[2],str(start_extl),str(start_extr),items[0],str(abs(int(items[8]))-9),'+'])+'\n')
                    bed.write('\t'.join([items[2],str(end_extl),str(end_extr),items[0],str(abs(int(items[8]))-9),'+'])+'\n')

                if int(items[1])==163:
                    bed.write('\t'.join([items[2],str(start_extl),str(start_extr),items[0],str(abs(int(items[8]))-9),'-'])+'\n')
                    bed.write('\t'.join([items[2],str(end_extl),str(end_extr),items[0],str(abs(int(items[8]))-9),'-'])+'\n')
    bed.close()
    return

#Bed2bedGraph
#************************************************#
def Bed2bedGraph(Inbed,ref_size):
    Outbedgraph=Inbed[:-4]+'.bedGraph'
    os.system("genomeCoverageBed -bg -split -i %s -g %s > %s" %(Inbed,ref_size,Outbedgraph))
    return
def PerbasebedGraph(InbedGraph):
    tmp=InbedGraph+'.tmp'
    tmpf=open(tmp,'w')
    lines=open(InbedGraph,'r')
    for line in lines:
        items=line.split('\t')
        for i in range(int(items[2])-int(items[1])):
            tmpf.write(items[0]+'\t'+str(int(items[1])+i)+'\t'+str(int(items[1])+i+1)+'\t'+str(items[3]))
    os.system('mv %s %s'%(tmp,InbedGraph))
    return



#Normbedgraph
#************************************************#
def Normbedgraph(Inbedgraph,Readlength=50,Totalcount=10000000):
    infile=open(Inbedgraph,'r')
    sumOfRead=0
    for line in infile:
        items=line.rstrip('\n').split()
        sumOfRead=sumOfRead+int(items[3])*(int(items[2])-int(items[1]))
    infile.close()
    sumOfRead=abs(sumOfRead)
    infile=open(Inbedgraph,'r')
    Outbedgraph=Inbedgraph[:-9]+'.norm.bedGraph'
    outfile=open(Outbedgraph,'w')
    rawReadLength=int(Readlength)
    for line in infile:
        items=line.rstrip('\n').split()
        value=float(items[3])/(sumOfRead/rawReadLength)*int(Totalcount)
        outfile.write('\t'.join([items[0],items[1],items[2],str(value)])+'\n')
    infile.close()
    outfile.close()
    return

#Sortbedgraph
#************************************************#
def Sortbedgraph(Inbedgraph):
    os.system("bedSort %s %s"%(Inbedgraph,Inbedgraph))
    return

#adjust_bedgraph
#************************************************#
def adjust_bedgraph(Inbedgraph,ref_size):
    size={}
    chrs = open(ref_size,'r')
    for line in chrs:
        items=line.rstrip('\n').split('\t')
        size[items[0]]=items[1]
    chrs.close()

    abgdir=Inbedgraph[:-8]+'tem.bedGraph'

    bg=open(Inbedgraph,'r')
    abg=open(abgdir,'w')
    for line in bg:
        items=line.rstrip('\n').split('\t')
        if (int(items[2])<=(int(size[items[0]])-2) and int(items[1])>=0):
            abg.write(line)
        else:
            print line
    bg.close()
    os.popen("mv "+ abgdir +" "+Inbedgraph)
    return

#Bedgraph2bigwig
#************************************************#
def Bedgraph2bigwig(Inbedgraph,ref_size):
    Outbw=Inbedgraph[:-9]+'.bw'
    os.system("bedGraphToBigWig "+Inbedgraph+' '+ref_size+' '+Outbw)
    return

#TSS Enrichment
#************************************************#
def Indexbam(Inbam):
    os.system("samtools index %s"%Inbam)
    return

def asn_mat(val,mat,s_int,e_int,t,i,weight):
    if float(val)>=s_int and float(val)<e_int-1 and t<1000:
        base = val-s_int
        mat[t][base] += weight
    return mat

def sub_Mat(inbed,inbam):
    mat = np.zeros([1000,4000])
    bedfile=np.loadtxt(inbed,'str')
    bamfile = pysam.Samfile(inbam, "rb")
    end=len(bedfile)
    for i in range(0,end):
        center = int(bedfile[i][1])+(int(bedfile[i][2])-int(bedfile[i][1]))/2
        s_int=center-2000
        e_int=center+2000
        for p2_rds in bamfile.fetch(str(bedfile[i][0]), max(0,s_int-2000), e_int+2000):
            if p2_rds.mapq<30:
                continue
            if p2_rds.is_reverse:
                continue
            else:
                l_pos = p2_rds.pos+4
                ilen = abs(p2_rds.tlen)-9
                r_pos=l_pos+ilen
                c_pos=l_pos+ilen/2
            mat = asn_mat(l_pos,mat,s_int,e_int,ilen,i,1)
            mat = asn_mat(r_pos,mat,s_int,e_int,ilen,i,1)
    return mat

def TSSEnrichmentplot(inbed,inbam,outfig):
    Indexbam(inbam)
    Mat=sub_Mat(inbed,inbam)
    mat0 = np.sum(Mat,0)
    fig=plt.figure(figsize=(8.0, 5.0))
    xran=500
    yran=500
    plt.plot(range(-int(2000),int(2000))[:-1],(mat0/np.mean(mat0[1:200]))[:-1],'k.')
    plt.plot(range(-int(2000),int(2000))[:-1],(np.convolve(mat0,np.ones(20),'same')/20/np.mean(mat0[1:200]))[:-1],'r')
    plt.xlabel('Distance to TSS center')
    plt.ylabel('Enrichment score')
    plt.savefig(outfig)
    return

#Fragment distribution
#************************************************#
def Fragdistribution(Inbam,outfig):
    sam=os.popen("samtools view -f 0x0002 %s"%Inbam)
    fragL=[]
    for line in sam:
        items=line.rstrip("\n").split('\t')
        flag=items[1]
        flag_length=abs(int(items[8]))
        if (int(flag)==99 or int(flag)==163) and int(flag_length)<=600:
            fragL.append(flag_length)
    count=Counter(fragL)
    fig=plt.figure(figsize=(10.0,4.0))
    df=pd.read_table('./Fragment_length_ratio.txt',header=0,index_col=0).iloc[:,1:]
    for i in range(df.shape[1]-1):
        plt.plot(df.index,df.iloc[:,i+1]/1000,color='#DBDBDB',linestyle='-',linewidth=4)
    plt.plot(count.keys(),np.array(count.values())/float(sum(count.values())),'r')
    plt.xlabel('Read length')
    plt.ylabel('Read counts %')

    fig.savefig(outfig)
    plt.close(fig)
    return

#QCTable
#************************************************#
def GetlineNum(inf):
    if inf[-4:]=='.bam':
        return int((os.popen("bedtools bamtobed -i %s |wc -l"%inf)).read().split()[0])
    elif inf[-4:]=='.bed':
        return int((os.popen("wc -l inf")).read().split()[0])
def GetlineNumOfAinB(inf1,inf2):
    return int(os.popen("bedtools intersect -a %s -b %s -u | wc -l"%(inf1,inf2)).read().split()[0])
def GetlineNumOfMapping(inf):
    return int((open(inf,'r').readlines())[0].split()[0]),float((open(inf,'r').readlines())[-1].split()[0][:-1])
def GetlineNumOfPicard(inf):
    return int(re.search(r'Marking (.*) records',open(inf,'r').read()).group(1))

def QC(outqc, name, maplog,chrMbam,dechrMrmdupbam,dechrMrmdupbed,BL,bam,dechrMbam,picardlog):
    qc=open(outqc,'w')
    qc.write("Sample" + "\t" + "TotalRawReads" + "\t" + "OverallAlignmentRate%" + "\t" + "FinalMappedReads" + "\t" + "FinalMapped%" + "\t" +
    "chrM%" + "\t" + "BlackListReads%" + "\t" + "MAPQFiltered%" + "\t" + "Duplicate%" + "\n")
    qc.write(name+'\t')
    totalReads,OverallAlignmentRate=GetlineNumOfMapping(maplog)

    chrMCount=GetlineNum(chrMbam)/2
    chrMPercent=chrMCount/float(totalReads)*100

    mappedReads=GetlineNum(dechrMrmdupbam)/2
    mappedPercent=mappedReads/float(totalReads)*100

    BLCount=GetlineNumOfAinB(dechrMrmdupbed,BL)/2
    BLPercent=BLCount/float(totalReads)*100

#    beforeQC=GetlineNum(bam)
#    afterQC=GetlineNum(dechrMbam)
#    filterQC=(beforeQC-afterQC-chrMCount)/2
#    qcPercent=filterQC/float(totalReads)*100

    DupCount=GetlineNumOfPicard(picardlog)/2
    DupPercent=DupCount/float(totalReads)*100

    qcPercent=OverallAlignmentRate-mappedPercent-chrMPercent-BLPercent-DupPercent

    qc.write("%d\t%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f"
             %(totalReads,OverallAlignmentRate,mappedReads,mappedPercent,chrMPercent,BLPercent,qcPercent,DupPercent))
    return

#####Read density around TSS
def ReaddensityTSS(inbed,TSS,OutFig):
    tss_lines = open(TSS).readlines()
    tss_all = [[x.split('\t')[0], int(x.split('\t')[1])] for x in tss_lines]
    tss_all = np.asarray(tss_all)
    tss_chrs = list(set(tss_all[:, 0]))
    tss_chrs.sort()
    tss_by_chrs, order_by_chrs = {}, {}
    for chrom in tss_chrs:
        order_chrom = np.where(tss_all[:, 0]==chrom)[0]
        tss_chrom = tss_all[order_chrom, 1]
        order_by_chrs[chrom] = order_chrom
        tss_by_chrs[chrom] = np.asarray(map(int, tss_chrom))
    counts_all = np.zeros((len(tss_all), 2001))

    window = 1000
    line=open(inbed,'r').readlines()
    i=0
    while i<len(line):
        chromosome,start,length=line[i].split('\t')[0],int(line[i].split('\t')[1]),int(line[i].split('\t')[4])
        end=start+length
        ranged_diff_order = np.where(abs(tss_by_chrs[chromosome]-start)<window)[0]
        common_order = order_by_chrs[chromosome][ranged_diff_order]
        for order in common_order:
            if abs(length)<=1000:
                counts_all[order,length+1000]+=1
        i=i+1
    i=1
    while i<len(line):
        chromosome,start,length=line[i].split('\t')[0],int(line[i].split('\t')[1]),int(line[i].split('\t')[4])
        end=start-length
        ranged_diff_order = np.where(abs(tss_by_chrs[chromosome]-end)<window)[0]
        common_order = order_by_chrs[chromosome][ranged_diff_order]
        for order in common_order:
            if abs(length)<=1000:
                counts_all[order,-length+1000]+=1
        i=i+1
    matrix = counts_all
    matrix = matrix[np.argsort(matrix.sum(axis=1))[::-1], :]
    matrix = matrix[:, 500:1500]
    n_ave, n_step, cut_off = 50, 50, 0.5
    n_limit = (len(matrix) - n_ave) // n_step
    ave_matrix = [matrix[i*n_step:i*n_step+n_ave, :].sum(axis=0)/float(n_ave) for i in range(0, n_limit)]
    ave_matrix = np.asarray(ave_matrix)
    log_ave_matrix = np.log(ave_matrix+1)
    max_logAve = log_ave_matrix.max()
    log_ave_matrix[np.where(log_ave_matrix > max_logAve*cut_off)] = max_logAve*cut_off
    nt_label = [' '] * 1000
    nt_label[0], nt_label[500], nt_label[999] = '-500', '0', '500'
    log_ave_df = pd.DataFrame(log_ave_matrix, index=xrange(0, len(log_ave_matrix)), columns=nt_label)
    sns.set_context('poster', font_scale=1)
    fig=plt.figure(figsize=(6, 14))
    ax21 = plt.subplot2grid((60, 1), (14, 0), rowspan=46)
    sns.heatmap(log_ave_df, xticklabels=True, yticklabels=False, cmap='Blues', cbar=False, ax=ax21)
    plt.setp(ax21, ylabel='TSSs sorted by expression')
    ax22 = plt.subplot2grid((60, 1), (0, 0), rowspan=9)
    ax22.plot(xrange(0, 500), matrix.sum(axis=0)[500:])
    ax22.set_xlim([0, 490])
    ax22.set_xticks([0, 250, 500])
    #ax22.set_yticks([0, 10000, 20000, 30000])
    ax22.set_xlabel('Fragment length (bp)')
    ax22.set_ylabel('Read density')
    ax22.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    fig.savefig(OutFig)
    plt.close(fig)
