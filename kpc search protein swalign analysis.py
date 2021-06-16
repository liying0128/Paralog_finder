import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Align
from Bio import Entrez
from Bio import SeqIO
import os
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


grabscore=100
#grab sequence
sequence=np.zeros((5000,3))
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
            n=n+1
sequence.to_csv('C:\\Users\\Administrator\\Desktop\\kpc_protein\\classA_protein.csv')

#to a fasta file
fastafile=tuple()
fastafile=(SeqRecord(template, id='kpc37 template'),)
for i in range(len(sequence)):
    if sequence.iat[i,0]==0:
        break
    current_gene=(SeqRecord(Seq(sequence.iat[i,2]), id=sequence.iat[i,1]),)
    fastafile=fastafile+current_gene
SeqIO.write(fastafile, "C:\\Users\\Administrator\\Desktop\\kpc_protein\\kpc.fasta", "fasta")


