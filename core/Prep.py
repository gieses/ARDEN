# -*- coding: cp1252 -*-
'''
Created 2012

Contains various help functions which initialize / translate /preprocess the data


@author: Sven Giese'''

import cPickle as pickle
import random

''' INIT DICTIONARIES '''
genetic_code={'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
              'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
              'AAT':'N', 'AAC':'N',
              'GAT':'D', 'GAC':'D',
              'TGT':'C', 'TGC':'C',
              'CAA':'Q', 'CAG':'Q',
              'GAA':'E', 'GAG':'E',
              'GGT':'G', 'GGC':'G','GGA':'G', 'GGG':'G',
              'CAT':'H', 'CAC':'H',
              'ATT':'I', 'ATC':'I','ATA':'I',
              'ATG':'M',
              'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
              'AAA':'K', 'AAG':'K',
              'TTT':'F', 'TTC':'F',
              'CCT':'P', 'CCC':'P','CCA':'P', 'CCG':'P',
              'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
              'ACT':'T', 'ACC':'T','ACA':'T', 'ACG':'T',
              'TGG':'W',
              'TAT':'Y', 'TAC':'Y',
              'GTT':'V', 'GTC':'V','GTA':'V', 'GTG':'V',
              'TAA':'*', 'TGA':'*','TAG':'*','NNN':'n'}



def createdic(AAsequence):
    """
    Creates the dictionary for the AA triplets and searches the starting indices 
    of the triplets in the given aminoacid sequence.

    @type  AAsequence: string
    @param AAsequence: aminoacid sequence
    @rtype:   dictionary
    @return:  A dictionary with starting positions of each triplet in the given AA sequence
    
    """
    
    liste = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","*"]
    aa_triplets = {}
    
    # create all possibilities  (triplets)
    for i in range(0,len(liste)):
        for k in range(0,len(liste)):
            for l in range(0,len(liste)):
                aa_triplets[liste[i]+liste[k]+liste[l]]= []
                
    # create lookup dic
    # key = triplet
    # value = list of positions                
    for i in range(1,len(AAsequence),3):
        if i+3 > len(AAsequence):
            break
        if AAsequence[i:i+3] in aa_triplets:
            aa_triplets[AAsequence[i:i+3]].append(i)
    return(aa_triplets)




def isvalidtriplet(codon,dictentry):
    """
    Function which checks if a given triplet has max hamming distance of 1 
    to a other triplet. Used for generation of possible substitutions triplets

    @type  codon: string
    @param codon: nucleotide triplet
    @type  dictentry: string
    @param dictentry: nucleotide triplet
    @rtype:   bool
    @return:  Boolean value. True if max hamming distance 1,else False .
    
    """
    counter = 0
    
    for i in range (0,3):
       
        if codon[i]== dictentry[i]:
            counter+=1
        else:
            continue
        
    if counter == 2:
        return (True)
    else:
        return (False)

def trans_seq(DNA):
    """
    Funtion which translates DNA to AA

    @type  DNA: list
    @param DNA: nucleotide sequence
    @rtype:   prot,rest
    @return:  Translated aminoacid sequence,untranslated nucleotide sequence
    """
    protein=[]
    prot = ""
    rest=""
    
    DNA = "".join(DNA)
    for i in range(0,len(DNA),3):
        # Codon exceeds length 
        if(i+3 > len(DNA)):
            rest +=DNA[i:i+3]
        
            break
        #' found Ns in nucleotid string
        if("N" in DNA[i:i+3]):
            a_a = "n"
            protein.append(a_a)
        else:
            #standard triplet translation
            codon=DNA[i:i+3]
            # look codon up in translation dic
            a_a=genetic_code[codon]
            protein.append(a_a)
            
    # transform to string
    prot = "".join(protein)
    return (prot,rest)

''' DEBUG HELP FUNCTIONS '''


def savepickle(dictionary,outputname):
    """
    basic pickle functions. actually for debugging and to speed up multiple simulations ( possible to load orf lists) 

    @type  dictionary: dictionary
    @param dictionary: Dictionary containg start and end positions of ORFs.
    @type  outputname: string
    @param outputname: Filename for saving.
    
    """
    pickle.dump( dictionary, open(outputname +".p", "wb" ) )
    print("Saved .pickle to: " + outputname +".p")

def loadpickle(inputname):
    """
    basic pickle functions. actually for debugging and to speed up multiple simulations ( possible to load orf lists) 


    @type  inputname: string
    @param inputname: Filename for loading.
    @rtype:   dictionary
    @return:  Dictionary containing start and end positions of ORFs.
    """
    dictionary= pickle.load( open(inputname ))#+".p" ) )
    print("Loaded "+inputname+" pickle!")
    return (dictionary)
