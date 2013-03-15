#!/usr/bin/env python
"""
Created 2012
core Script for the generation of the artificial reference genome

The functions purpose is to to go through a list of positions and find balanced mutations
which fulfill the demands on the artificial reference genome. Once a initial start positions is
randomly selected all possible triplets with hamming distance 1 are generated and looked up 
in a dictionary which contains all triplet positions in the input genome. If a suitable partner
is found for the initial mutation the next start positions is chosen randomly. Else: try all other
triplets with hamming distance 1 until no one is left. This process can be accelerated by allowing
unbalanced mutations, but this will cause differences in the NUC/AA distribution and the AA neighborhood.


@author: Sven Giese
"""

import random as r
import Prep as INI



def getMutation(AA,Codon):
    """
    Returns a random mutation for a given AA and its Codon(DNA). The mutation is done in a way which supports the equilibrium 
    of the nucleotide distribution by only regarding hamming distance=1 Codons as possible mutations 

    @type  AA: string
    @param AA: Single AA.
    @type  Codon: string
    @param Codon: 3-letter dna code.
    @rtype:   list,list
    @return:  A list of all valid mutations (triplet) and the coresponding AA.
    """
    temp_mutationlist = []
    '''create a list of possible triplets within hamming distance 1 '''
    for item in INI.genetic_code.keys():
        isvalid = INI.isvalidtriplet(item,Codon)
        ''' Hamming distance 1, AA is not equal to the given AA,forbid mutation to stopcodon '''
        if (isvalid == True and AA !=INI.genetic_code[item] and INI.genetic_code[item]!="*"):
            temp_mutationlist.append(item)
    
    
    aalist = []
    # generate a list of all possible amino acids resulting from the temp_mutationlist 
    for item in temp_mutationlist:
        if (item in INI.genetic_code):
            aalist.append(INI.genetic_code[item])
        else:
            aalist.append("n")
            
    return(temp_mutationlist,aalist)


def getdifference(triplet_old,triplet_new):
    """
    Given two triplets, returns the differences between them plus the position 

    @type  triplet_old: string
    @param triplet_old: AA triplet.
    @type  triplet_new: string
    @param triplet_new: AA triplet.
    @rtype:   Char,Char,int
    @return:  The new aminoacid, the old aminoacid and the position. 
    """
    for i in range(0,3):
        if (triplet_new[i]!=triplet_old[i]):
            
            return (triplet_new[i],triplet_old[i],i)
        
      
def isvalidposition(pdic,iprime,distance):
    """
    Checks if a position is valid for mutation. It queries all neighboring positions (iprime +-distance) to check whether there already was a mutation in pdic

    @type  pdic: dictionary
    @param pdic: Diciontary containing mutations and start/ stop codons..
    @type  iprime: int
    @param iprime: Position of the prospective mutation (DNA level)
    @type  distance: int
    @param distance: User defined parameter which limits the distance between two mutations.
    @rtype:   Bool
    @return: Boolean which decides if the position is valid (1= yes,0 = no) 
    """
    
    # deal with base shifts 
    distance = distance-2
    
    istforbidden = 0
    for o in range(-distance,distance+2,1):
        if (iprime+o in pdic):
            # E = end of orf
            # S = start of orf
            if((pdic[iprime+o]=="E") or (pdic[iprime+o]=="S")):
                if((o >3) or (o <-3)):
                    pass
                else:
                    istforbidden = 1
                    break
            else:
                istforbidden = 1
                break
    else:
            pass
                
    return(istforbidden)




def mutate_random(DNA,AminoAcid,distance,pdic,rev,header,Random,outputpath):
        """
        Mutates a given DNA(AminoAcid) Genomesequence on several positions (distance based on DISTANCE var. If one mutation is done
        a compareable Triplet is searched to "reverse" the  changes made in AA distribution, N distribution, AA neighborhood
    
        @type  DNA: list
        @param DNA: DNA sequence of the reference genome.
        @type  AminoAcid: list
        @param AminoAcid: AA sequence of the reference genome.
        @type  rev: Bool
        @param rev: Boolean which decides if unbalanced mutations are allowed (only initial mutation is performed)
        @type  pdic: dictionary
        @param pdic: Diciontary containing mutations and start/ stop codons..
        @type  header: string
        @param header: Header for the resulting artificial reference file (fasta format).
        @type  Random: Bool
        @param Random: Boolean for choosing on of the mutation modes (linear = 0,random = 1)
        @type  distance: int
        @param distance: User defined parameter which limits the distance between two mutations.
        @rtype:   list
        @return: Artificial reference genome sequence.
        """
        ##debug vals 
        start = [] # list of start positions of mutations ( start means first mutation in balanced case)
        both = []  # start and end position
        fobj2= open(outputpath+header+"_CompleteLog.txt","a")
        fobj2.write("BalancedMutation"+"\t"+"NewAA" + "\t" + "OldAA"+"\t"+"NewAAPos"+"\t"+"OldAAPos" +"\t"+ "NewDNA"+"\t"+ "OldDNA"+ "\t"+"NewDNAPos"+"\t"+"OldDNAPos"+"\n")
        fobj2.close()
        
        
        # generate start positions for mutation (the samplespace)
        samplespace = []
        for i in range (2,len(AminoAcid),distance/3):
            samplespace.append(i)
        
        
        ##random_modification
        if (Random ==1):
            r.shuffle(samplespace)
        else:
            pass
      
        dna_list = list(DNA)
        AminoAcid_list = list(AminoAcid)
       
        '''the lookup dictionary for the aa triplets '''
        lookup_dic = INI.createdic(AminoAcid)

        #gotit indicator if a possibility was found to revert the initial changes (start of mutation)
        gotit=False
        # stat variables
        succ_counter = 0
        fail_counter = 0 
        skip = 0
       
        ''' Main loop over the AminoAcid'''
        for i in samplespace:
            ''' no triplet left --> break '''
            if(i+2 >len(AminoAcid)):
                print("\t(finished...exceeded length of AA)")
                continue
            
            ''' AA which is going to be mutated'''
            AA = AminoAcid_list[i]
        
            '''index for dna : i*3 --> AminoAcid --> DNA
            #not i*3+3 because i starts at AA 2 since we need a right and left neighbor'''
            iprime = i*3
            
            '''AA and corresponding DNA triplet for the middle AA '''
            AA_triplet= AminoAcid_list[i-1]+AminoAcid_list[i]+AminoAcid_list[i+1]
            DNA_triplet = DNA[iprime:iprime+3]

            # get temporary list of all mutations. Iterate over it to find best possible substitution
            mutationsliste,aaliste = getMutation(AA,DNA_triplet)
            
            
            # isvalidposition returns 1 if the position isforbidden, else 0
            val = isvalidposition(pdic, iprime, distance)
            if (val ==1):
                skip+=1
                fobj2= open(outputpath+header+"_CompleteLog.txt","a")
                fobj2.write(str(0)+"\t"+new_AA_triplet + "\t" + "' '"+"\t"+str(i)+"\t"+"' '" +"\t"+ new_triplet+"\t"+ "' '"+ "\t"+str(iprime+position)+"\t"+"'(skipped)'"+"\n")
                fobj2.close()
                continue
                    
            else:
                pass
            

            for q,item in enumerate(mutationsliste):
                
                if gotit==True:
                    break
                else:
                    pass
                
                ''' old and new variables for before/after the mutation '''
                new_triplet                 = mutationsliste[q]
                new_AA                      = aaliste[q]
                new_N,old_N,position        = getdifference(DNA_triplet,new_triplet)
                new_AA_triplet              = AA_triplet[0]+new_AA+AA_triplet[2]
                tempdic = pdic
                tempdic[iprime+position]="M"
                
                if (new_AA_triplet in lookup_dic):
                    '''templist--> contains all starting positions of the "new_AA_triplet" which we want to substitute back '''
                    templist = lookup_dic[new_AA_triplet]
                    
                    
                    # add potential mutation to dictionary
                    tempposition = [iprime+position,"M"]
                    for l in range(0,len(templist)):
                        posi = templist[l]
                        # i*3 --> protein nach DNA, +3 betrachten IMMER mittlere AA
                        ''' suitable dna position found? '''
                        if (new_triplet == dna_list[posi*3+3]+dna_list[posi*3+3+1]+dna_list[posi*3+3+2]):
                            val = isvalidposition(tempdic, posi*3+3+position, distance)
                            
                            if (val ==1):
                                skip+=1
                                continue
                            else:
                                pass
                            
                            '''back substitution & do subs on 1st position'''
                            pdic[posi*3+3+position]="R"
                            dna_list[posi*3+3+position]= old_N
                            
                            pdic[iprime+position]="M"
                            dna_list[iprime+position]= new_N
                            
                            AminoAcid_list[i]= new_AA
                            AminoAcid_list[posi+1]= AA
                            
                            gotit = True
                            succ_counter+=1
                            #lookup_dic[new_AA_triplet]  = [i for i in lookup_dic[new_AA_triplet] if i!=posi]
                            lookup_dic[new_AA_triplet].remove(posi)
                            
                            '''writing the log file '''
                            fobj= open(outputpath+header+"_CompleteLog.txt","a")
                            fobj.write(str(1)+"\t"+AA_triplet + "\t" + new_AA_triplet+"\t"+str(i)+"\t"+str(posi) +"\t"+ DNA_triplet+"\t"+ str(new_triplet)+ "\t"+str(iprime+position)+"\t"+str(posi*3+3+position)+"\n")
                            fobj.close()
                            
                            ## statistics
                            start.append(iprime+position)
                            both.extend([iprime+position,posi*3+3+position])
                            break
                
                # no possible triplet positions for back substitution in lookup_dic    
                else:
                    continue
                 
            # after loop 
            if (gotit==False):
                fobj2= open(outputpath+header+"_CompleteLog.txt","a")
                fobj2.write(str(0)+"\t"+new_AA_triplet + "\t" + "' '"+"\t"+str(i)+"\t"+"' '" +"\t"+ new_triplet+"\t"+ "' '"+ "\t"+str(iprime+position)+"\t"+"'(tried)'"+"\n")
                fobj2.close()
                fail_counter+=1
                # reverse substitutions on? (=1) off (=0). If one dont change first mutation in the first place. Else: just change it.. 
                if (rev==0):
                    pdic[iprime+position]="M"
                    dna_list[iprime+position]= new_N
                    AminoAcid_list[i]= new_AA
                    start.append(iprime+position)
                    both.extend([iprime+position])        
            elif (gotit==True):
                    gotit = False
            
        # stats (INI.savepickle(pdic,header+"_pdic_e"))
        print("\r\n########Some stats:########")
        print("DNA length:\t" + str(len(DNA)))
        print("max substitutions:\t" + str(len(DNA)/distance))
        print("#Balanced Mutations:\t" + str(succ_counter))
        
        
        return ("".join(dna_list))