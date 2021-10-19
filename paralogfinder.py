from Bio import SeqIO
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Align
from Bio.Seq import Seq
import swalign
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
def paralogfinder (proj_name,path1,path2):    
#    proj_name='svm'
#    path1='/home/ly/paralog-finder/mutation/'
#    path2='/home/ly/paralog-finder/kpgenome/'
    path3='/home/ly/paralog-finder/mutation/'+proj_name+'/'
    
    #search
    match=2
    mismatch=-1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    
    for tem in SeqIO.parse(path3+proj_name+'.fasta','fasta'):
        template=tem.seq
        templatename=tem.id
    os.chdir(path2)
    filelist=os.listdir()
    data=np.zeros((len(filelist),10000))
    data=pd.DataFrame(data)
    protein_length=len(template)
    
    proteinid=np.zeros((len(filelist),10000))
    proteinid=pd.DataFrame(proteinid)
    proteinid.at[0,:]='a'
    
    for i in range(len(filelist)):
        name=filelist[i]
        m=0
        for seq_record in SeqIO.parse(name,'fasta'):
            seq=str(seq_record.seq)
            seq=seq.replace('*','')
            seq=Seq(seq)
            alig=sw.align(template,seq)
            data.iat[i,m]=alig.score
            proteinid.iat[i,m]=seq_record.id
            m=m+1
        fig,ax=plt.subplots()
        ax.plot(range(10000),data.iloc[i,:],color='g',linewidth=4)
        ax.set_ylim(10,protein_length*2.1)
        plt.savefig(path3+proj_name+str(name)+'.png',dpi=500)
    data.to_csv(path3+proj_name+'_score.csv')
    
    
    #analysis
    grabscore=protein_length*0.6
    sequence=np.zeros((5000,4))
    sequence=pd.DataFrame(sequence)
    sequence.at[0,:]='a'
    n=0
    for i in range(len(filelist)):
        for k in range (10000):
            if data.iat[i,k]>grabscore:
                sequence.iat[n,0]=filelist[i]
                allid= [seq_record.id for seq_record in SeqIO.parse(filelist[i], "fasta")]
                sequence.iat[n,1]=str(allid[k])
                allseq= [seq_record.seq for seq_record in SeqIO.parse(filelist[i], "fasta")]
                sequence.iat[n,2]=str(allseq[k])
                sequence.iat[n,3]=str(filelist[i])+'_'+str(allid[k])
                n=n+1
    sequence.to_csv(path3+proj_name+'_protein.csv')
    
    #to a fasta file
    fastafile=tuple()
    fastafile=(SeqRecord(template, id='template'),)
    for i in range(len(sequence)):
        if sequence.iat[i,0]==0:
            break
        current_gene=(SeqRecord(Seq(sequence.iat[i,2]), id=sequence.iat[i,3]),)
        fastafile=fastafile+current_gene
    SeqIO.write(fastafile, path1+proj_name+'.fasta', "fasta")