#!/usr/bin/env python


"""
FindOrfs, Created 2012

Script to detect start and stop codon of open reading frames in a dna sequence. Note that the algorithm is only rudimentary and does not respect intron-exon structures.

"""


def build_ORF(sequence,file_ORF,pdic):
    """
    Get orf positions (forward/backward) and return them in a dictionary

    @type  sequence: string
    @param sequence: nucleotide sequence
    @type  file_ORF: string
    @param file_ORF: outputfile
    @type  pdic: dictionary
    @param pdic: Used to start start / end positions of ORFs.
    @rtype:   dictionary
    @return:  Stored start start / end positions of ORFs
    """

    
    START,STOP,STARTrev,STOPrev = find_orfs(sequence,pdic)
    t=open(file_ORF,'a')
    print "\t search forward orfs..."
    orf_counter = 0
    orf_name = "forw"
    for i in START:
        for j in STOP:
            #3 conditions
            #first condition: start position needs to be BEFORE the stop position
            #second condition: our proteins need to be longer than 150 aa AND need to be smaller than the biggest protein we know
            #third condition: each ORF needs to give a real aa sequence
            if i<j and (j+3-i)>149 and (j+3-i)<103050 and (j+3-i)%3==0:
                orf_counter+=1
                #if all 3 conditions: write the ORF in a file
                t.write(">"+orf_name+" " + str(orf_counter))
                t.write('\n')
                t.write(sequence[i:j+3])
                t.write('\n')
                pdic[i]="S" # S  = START
                pdic[j]="E" # E  = END
                
    temp = orf_counter
    
    orf_counter = 0
    orf_name = "rev"
    print "\t search backward orfs..."
    for i in STARTrev:
        for j in STOPrev:
            #3 conditions
            #if j<i and (j+2-i+1)>150 and (j+2-i+1)<103050 and (j+2-i+1)%3==0: <--alt
            #first condition: start position needs to be BEFORE the stop position
            #second condition: our proteins need to be longer than 150 aa AND need to be smaller than the biggest protein we know
            #third condition: each ORF needs to give a real aa sequence
            if j<i and (i+3-j)>150 and (j+3-i)<103050 and (j+3-i)%3==0:
                #if all 3 conditions: write the ORF in a file
                orf_counter+=1
                t.write(">"+orf_name+" " + str(orf_counter))
                t.write('\n')
                t.write(sequence[j:i+3][::-1])
                t.write('\n')
                pdic[i]="E"
                pdic[j]="S"
    print(str(temp+orf_counter))
    t.close()
    
    return (pdic)


def findstop_help(posLastStop,sequence,codon):
    """
    return the index of the first position of codon in the dna sequence  

    @type  posLastStop: int
    @param posLastStop: Position of the last found stop codon.
    @type  sequence: string
    @param sequence: Nucleotide sequence.
    @type  codon: string
    @param codon: 3-letter DNA code.
    @rtype:   int
    @return:  The position of the stop codon in the nucleotide sequence.
    """
    try:
        return(sequence.index(codon,posLastStop))
            
    except:
        return(-1)
    
def find_orfs(genomeSequence,pdic):
    """
    function to identify open reading frames in a dna sequence, careful: intron exon structures are not respected!  

    @type  genomeSequence: string
    @param genomeSequence: Nucleotide sequence.
    @type  pdic: dictionary
    @param pdic: Used to store start / end positions of ORFs.
    @rtype:   dictionary
    @return:  Found  start / end positions in the sequence consindering only the ORFs.
    """
    # got a dictionary and a sequence! lets start!
    # Startcodon = ATG
    # Stopcodon = ["TAA","TGA","TAG"]
    
 
    start =[]
    stop =[]
    
    posLastATG = 0
    posLastStop = 3
    orfList = []
    
    # forward
    print("\t..find forward orfs")
    while True:
        foundNew = False
        # catch error if "ATG" not found!
        try:
            
                start.append(genomeSequence.index("ATG",posLastATG))
                # retrieve last element
                posLastATG = start[-1]+1
                foundNew = True
        except:
            pass
        
        stopSub =[]
        stopcodons =["TAA","TGA","TAG"]
        for item in stopcodons:
            stopSub.append(findstop_help(posLastStop, genomeSequence, item))

        stopSub.sort()
        
        if(stopSub[0] > -1):
            stop.append(stopSub[0])
            posLastStop = stop[-1]+1
            foundNew = True
            
        elif(stopSub[1] > -1):
            stop.append(stopSub[1])
            posLastStop = stop[-1]+1 
            foundNew = True;
            
        elif(stopSub[2] > -1):
            stop.append(stopSub[2]); 
            posLastStop = stop[-1]+1 
            foundNew = True;
        
        if(foundNew):
            pass
        else:
            break
        
    # reverse now:
    # start codon: CAT
    # stop codons: TTA;TCA;CTA 
        
    startRev = []
    stopRev = []

    posLastCAT = 3
    posLastStop_rev = 0
       
    print("\t..find reverse orfs")
    while True:
            foundNew_rev  = False
            # catch error if "CAT" not found!
            try:
                
                    startRev.append(genomeSequence.index("CAT",posLastCAT))
                    # retrieve last element
                    posLastCAT = startRev[-1]+1
                    foundNew_rev = True
            except:
                pass
            
            stopSub =[]
            stopcodons =["TTA","TCA","CTA"]
            for item in stopcodons:
                stopSub.append(findstop_help(posLastStop_rev, genomeSequence, item))
    
            stopSub.sort()
            
            if(stopSub[0] > -1):
                stopRev.append(stopSub[0])
                posLastStop_rev = stopRev[-1]+1
                foundNew_rev = True
                
            elif(stopSub[1] > -1):
                stopRev.append(stopSub[1])
                posLastStop_rev = stopRev[-1]+1
                foundNew_rev = True
                
                
            elif(stopSub[2] > -1):
                stopRev.append(stopSub[2])
                posLastStop_rev = stopRev[-1]+1
                foundNew_rev = True
                
            
            if(foundNew_rev):
                pass
            else:
                break
    ###test##################
    print("START codons : "  + str(len(start)))
    print("STOP codons : "  +str(len(stop)))
    print("revSTART codons : "  +str(len(startRev)))
    print("revSTOP codons : "  +str(len(stopRev)))
    
    
    # FORWARD
    #vectors are filled, now test if in the same frame
    # stores start positions without a suitable Stop in frame
    removeList=[]
    print("\t..creating forward orfs")
    for stopPos in stop:
        foundPartner=False
        startPos = 0
        i = 0
        startPos = start[0]
        #= start[i]
        while ((i<len(start)-1) and (start[i] <stopPos)):
                startPos = start[i]
                i+=1
                if((stopPos+3-startPos)% 3 == 0):
                    #Orf found in frame
                    if((foundPartner == False) and (stopPos+3)-startPos >149):
                        #found first matching AT (first ORF+ORF is at least 150bp long
                        foundPartner = True
                        pdic[stopPos+3]="E"
                        pdic[startPos]="S"
                        orfList.append(startPos)
                        orfList.append(stopPos+3)
                        removeList.append(startPos)
                    else:
                        removeList.append(startPos)
        for item in removeList:
            # remove all positions of ATGS leading to "sub-ORFs"
            start.remove(item)
        removeList =[]
        
     
    # BACKWARD
    print("\t..creating reverse orfs")
    removeList_rev=[]
    l = len(stopRev)-1
    
    for r in range(l,-1,-1):#
        stopPos = stopRev[r]
        foundPartner=False
        
        
        
        i = len(startRev)-1
       
        
        while((i >=0) and (startRev[i]>stopPos)):
            startPos=startRev[i]
            i -=1
            if((startPos+3-stopPos)%3 == 0):
                # found ORF in frame
                if((foundPartner != True) and (startPos+3-stopPos > 149)):
                    #found first matching ATG + ORF is at least 150bp long
                    pdic[stopPos]="S"
                    pdic[startPos+3]="E"
                    foundPartner = True
                    orfList.append(startPos)
                    orfList.append(stopPos+3)
                    removeList_rev.append(startPos)
                else:
                    removeList_rev.append(startPos)
                    
        for item in removeList_rev:
        # remove all positions of ATGS leading to "sub-ORFs"
            startRev.remove(item)
        removeList_rev =[]
        
    return(pdic)