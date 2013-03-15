#!/usr/bin/env python
'''
Created 2012

Contains various help functions which read or produce an input/ output


@author: Sven Giese
'''
import os
import random
import HTSeq


def readdna(filename):
    """
    Reads in the dna sequence of the given fasta

    @type  filename: string
    @param filename: Fasta-file used as input.
    @rtype:   HTSeq Sequence object
    @return:  Reference Fasta.
    """
    chr = HTSeq.FastaReader(filename)
    for fasta in chr:
        referenz = HTSeq.Sequence(fasta.seq,fasta.name)
    return(referenz)


def writefile(sequenceObject,filename):
    """
    Writes a given sequence object to a fasta file.

    @type  sequenceObject: HTSeq Sequence object
    @param sequenceObject: Reference sequence as fasta.
    """
    
    outfasta = open(filename,"w")
    sequenceObject.write_to_fasta_file(outfasta)
    outfasta.close()


def writeoverview(Ndic_G,aadic_G,Ndic_AR,aadic_AR,filename):
    """
    Creates the "delta" file for the comparison of the two chromosoms. This file contains the differences in nucleotide distribution between reference and artificial.
    input: nucleotid dictionary genom, aa dictionary genome, nucleotid dictionary artificial chromosom, aa dictionary, filename 

    @type  Ndic_G: dictionary
    @param Ndic_G: Nucleotid dictionary genom.
    @type  aadic_G: dictionary
    @param aadic_G: AA dictionary genome.
    @type  Ndic_AR: dictionary
    @param Ndic_AR: Nucleotid dictionary artificial.
    @type  aadic_AR: dictionary
    @param aadic_AR: AA dictionary artificial
    @type  filename: string
    @param filename: Output filename.
    """
    fobj = open(filename,"w")
    fobj.write("NUC /AA \t Genom \t Artificial Reference \t Delta \n")
   
    sum1 =0
    sum2= 0
    for item in Ndic_G.keys():
        fobj.write(item +"\t"+str(Ndic_G[item])+"\t"+str(Ndic_AR[item])+"\t"+str(Ndic_G[item]-Ndic_AR[item])+"\n")
        sum1 +=abs(Ndic_G[item]-Ndic_AR[item])
    fobj.write(str(sum1)+"\n")
    
    for item in aadic_G.keys():
        fobj.write(item +"\t"+str(aadic_G[item])+"\t"+str(aadic_AR[item])+"\t"+str(aadic_G[item]-aadic_AR[item])+"\n")
        sum2 +=abs(aadic_G[item]-aadic_AR[item])
    fobj.write(str(sum2)+"\n")
    
    
    

def nucleotide_dist_seq(seq,txt_file,shallwrite):
    """
    Writes the nucleotide distribution in a file and returns the dictionary. adjust s for % results.
    @type  seq: string
    @param seq: Nucleotide sequence.
    @type  txt_file: string
    @param txt_file: Output compare file.
    @type  shallwrite: Bool
    @param shallwrite: Decides if percentages values are written to the output.
    """
    Nndic={"A":0,"C":0,"G":0,"T":0,"N":0}
    
    for i in range(0,len(seq)):
          Nndic[seq[i]]+=1
    s=len(seq)
    s=1
   
    if (shallwrite==1):
        output_file=open(txt_file,'w')
        for item in Nndic.keys():
            Nndic[item]=Nndic[item]/float(s)
            output_file.write(item + "\t" + str(Nndic[item])+"\n")
            
        output_file.close()
    else:
         for item in Nndic.keys():
            Nndic[item]=Nndic[item]/float(s)
    return (Nndic)    #N can be used for checking: should be the same number in real
                                                                                    # and artificial chromosome


def aa_dist_seq(seq,txt_file,shallwrite):
    """
    Writes the AA distribution in a file and returns the dictionary. adjust s for % results.
    @type  seq: string
    @param seq: Nucleotide sequence.
    @type  txt_file: string
    @param txt_file: Output compare file.
    @type  shallwrite: Bool
    @param shallwrite: Write output in percentages..
    """
    aadic = {"A":0,"R":0,"N":0,"D":0,"C":0,"E":0,"Q":0,"G":0,"H":0,"I":0,"L":0,"K":0,"M":0,"F":0,"P":0,"S":0,"T":0,"W":0,"Y":0,"V":0,"*":0}
    for i in range(0,len(seq)):
        
        '''escape 'n' Sequences '''
        if (seq[i] in aadic):
              aadic[seq[i]]+=1
        else:
            continue
            
    
    n = len(seq)
    n=1
    if (shallwrite==1):
        output_file=open(txt_file,'w')
        for item in aadic.keys():
            aadic[item]=aadic[item]/float(n)
            output_file.write(item + "\t" + str(aadic[item])+"\n")
            
        output_file.close()
    else:
        for item in aadic.keys():
            aadic[item]=aadic[item]/float(n)
            
    return (aadic) 

'''
input: DNA Sequence, outputfilename and 1/0 for writing/not writing outputfile '''

def nucleotide_dist_file(file_fasta,txt_file):
    """
    Writes the DNA distribution in a file and returns the dictionary. adjust n for % results

    @type  file_fasta: string
    @param file_fasta: DNA Sequence
    @type  txt_file: string
    @param txt_file: Filename for output.
    """
    input_file=open(file_fasta,'r')
    output_file=open(txt_file,'a')
    seq=''
    for line in input_file:
        if line[0]!='>':
            line=line.rstrip()
            seq+=line
    output_file.write(str(nucleotide_dist_seq(seq)))
    output_file.write('\n')
    output_file.close()
    input_file.close()


'''gets the number of missmatches between 2 sequences
input: orig sequence, decoy sequence '''
def gethammingdistance(original,artificial):
    """
    Calculates the hamming distances between two sequences.
    @type  original: list
    @param original: Nucleotide sequence from the reference.
    @type  artificial: list
    @param artificial: Nucleotide sequence from the artificial reference.
    """
    hamming = 0
    not_hamming=0
    for i in range(0,len(original)):
        if (original[i]!=artificial[i]):
            hamming +=1
            
        else:
            not_hamming+=1
    print ("#hamming distance REF-ART\t"+ str(hamming))
    print ("avg. distance:\t" + str(len(original)/float(hamming)))
    print("###########################\r\n")