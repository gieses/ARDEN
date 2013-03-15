#!/usr/bin/env python
""" This script is used to evalulate the mapping results """

import HTSeq
import cPickle as pickle
import numpy as np
import time
import sys
import re
import subprocess
import random 
import os



alngts = 0
AlignedReadsdic = {}
algnedToRef = 0
algnedToArt = 0

""" Helper debug functions.... """
def savepickle(dictionaary, outputname):
    pickle.dump(dictionaary, open(outputname + ".p", "wb"))
    print("Saved .pickle!")

def loadpickle(inputname):
    dictionary = pickle.load(open(inputname + ".p"))
    print("Loaded " + inputname + ".pickle!")
    return (dictionary)
    

class CustomRead:
    """
    Class for exhaustive reading of the sam alignments.
    
    """
    
    def __init__(self, readname, mr, nm, subs, score, mq, start, end, gaps, mism):
        self.id = readname      # read identifier 
        self.mr = [mr]          # matched reference
        self.subs =[subs]       # substitutions in artificial reference genome
        self.nm = [nm]          # number of mismatches (according to SAM)
        self.rq = score         # read quality score
        self.mq = [mq]          # mapping quality
        self.start = [start]    # alignment start
        self.end = [end]        # alignment end
                                    #mr,nm,subs,score,mq,start,end
        self.gaps = [gaps]
        self.mism = [mism]
        
    def toObjself(self, read): 
        self.mr.append(read.mr[0])          
        self.subs.append(read.subs[0])     
        self.nm.append(read.nm[0])         
        self.mq.append(read.mq[0])          
        self.start.append(read.start[0])     
        self.end.append(read.end[0])       
        self.gaps.append(read.gaps[0])
        self.mism.append(read.mism[0])
        
    def toStr(self, readname): 
        """Converts the object to a string """
        str = ""

        for i, read in enumerate(self.mr):
            #print read
            str += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n" % (readname, self.mr[i], self.subs[i], self.nm[i], self.rq, self.mq[i], self.start[i], self.end[i], self.gaps[i], self.mism[i])
        return(str)
    
    def toStrNoMem(self, identifier): 
        """Converts the object to a string """

        if identifier == "art":
            str = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n" % (self.id, 0, self.subs[0], self.nm[0], self.rq, self.mq[0], self.start[0], self.end[0], self.gaps[0], self.mism[0])
        else:
            str = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n" % (self.id, 1, 0, self.nm[0], self.rq, self.mq[0], self.start[0], self.end[0], self.gaps[0], self.mism[0])
        return(str)
        
        

class TPRead:
    """ Class for exhaustive reading of the sam alignments. Only for TP """
    def __init__(self, nm, score, mq, start, end, gaps, mism):
        self.nm = [nm]       # number of mismatches (according to SAM)
        self.rq = score                  # read quality score
        self.mq = [mq]      # mapping quality
        self.start = [start] # alignment start
        self.end = [end]    # alignment end
        self.gaps = [gaps]
        self.mism = [mism]
        
    def toObjself(self, read):          
        self.nm.append( read.nm[0])         
        self.mq.append(read.mq[0])          
        self.start.append( read.start[0])     
        self.end.append( read.end[0])   
        self.gaps.append( read.gaps[0])
        self.mism.append(read.mism[0])   

    def toStr(self, readname): 
        str = ""
        for i, read in enumerate(self.mr):
            str += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n" % (readname, self.nm[i], self.rq, self.mq[i], self.start[i], self.end[i], self.gaps[i], self.mism[i])
        return(str)
        
    """ returns 1 if an equal alignment is contained in the array """
    def isContained(self, read):
        for i, position in enumerate(self.start):
            # Checks if a given Read (more accurate the alignment) is already in the object. this is used to check against FP
            if position == read.start and self.end[i] == read.end:
                return(1)
        return(0)
    

class ReadID:
    """ Class for efficient addressing for np. arrays! """
    def __init__(self, internalID, quality):
        self.internalID = internalID
        self.quality = quality


############################################# CLASSES END###########################################


def getTPRead(alngt, compareList, readdic):
    """
    Funtion which transforms an HTSeq alignment to a TPRead class.

    @type  alngt: alignment
    @param alngt: alignment from the sam file.
    @type  compareList: list
    @param compareList: list which indicates differences between reference / artificial
    @type  readdic: dictionary
    @param readdic: Contains ranks and qualities for every entry in the fastq file
    
    @rtype:   TPRead
    @return:  Transformed ReadObject.
    """
    cigar = alngt.cigar
    gaps  = sum([i.size for i in cigar if i.type=="I" or i.type=="D"])
    mism = sum([i.size for i in cigar if i.type=="M"]) - sum([int(s) for s in alngt.optional_field("MD") if s.isdigit()])
    nm = alngt.optional_field("NM")
    
    #####->TPRead(nm, score, mq, start, end, gaps, mism, n)
    return (TPRead(nm, readdic[alngt.read.name], alngt.aQual, alngt.iv.start, alngt.iv.end, gaps, mism))
    
def getCustomRead(alngt, compareList, readdic):
    """
    Funtion which transforms an HTSeq alignment to a CustomRead class.

    @type  alngt: alignment
    @param alngt: alignment from the sam file.
    @type  compareList: list
    @param compareList: list which indicates differences between reference / artificial
    @type  readdic: dictionary
    @param readdic: Contains ranks and qualities for every entry in the fastq file
    
    @rtype:   CustomRead
    @return:  Transformed ReadObject.

    """
    cigar = alngt.cigar
    # number of gaps is the sum of "I" (insertions) and "D" deletions.
    gaps  = sum([i.size for i in cigar if i.type=="I" or i.type=="D"])
    
    # mismatches are calculated as cigar XM - MDtag occurrences of numbers
    # which indicate mismatches
    mism = sum([i.size for i in cigar if i.type=="M"]) - getsum(re.findall("(\d+)", alngt.optional_field("MD")))
    
    
    nm = alngt.optional_field("NM")

    ########
    if len(alngt.read.qualstr) <= 1:
       
        return (CustomRead(alngt.read.name,0,nm,getHammingdistance(compareList, alngt.iv.start, alngt.iv.end),readdic[alngt.read.name],alngt.aQual,alngt.iv.start,alngt.iv.end,gaps,mism))
    else:
        return (CustomRead(alngt.read.name,0,nm,getHammingdistance(compareList, alngt.iv.start, alngt.iv.end),(alngt.read.qual.mean() / (alngt.read.qual.max() -  alngt.read.qual.min() + 1)),alngt.aQual,alngt.iv.start,alngt.iv.end,gaps,mism))

def GetOrderDictionary(referenceSAM):
    """
    Function to get a dictionary containing the rank of a read (given a sorted samfile by readname, samtools).

    @type  referenceSAM: string
    @param referenceSAM: Inputfile name for reference SAM file.
    @rtype:   dictionary
    @return:  Internalnaming according to the sorting of samtools. Key = ReadID, Value = rank
    """
    internalDic = {}

    # use ultrafast awk for getting only readnames from aligned fasta file
    readnames = str(subprocess.check_output("""awk < %s -F "\t" '{print  $1 }'""" % referenceSAM, shell=True))

    # make to list
    readnames = readnames.split("\n")

    i = 0
    # skip header files
    while readnames[i][0] == ("@"):
        i += 1
    
    index = 0
    for k in xrange(i, len(readnames)):
        if readnames[k] in internalDic:
            pass
        else:
            internalDic[readnames[k]] = index
            #enumeration of the index
            index += 1
            
            
            
    return(internalDic)


def getNextLine(textfile, compareList, readdic):
    """
    Function which returns the next line from a SAMFile.
    
    @type  textfile: fileobject stream
    @param textfile: SAMfile
    @type  compareList: list
    @param compareList: AA sequence of the reference genome.
    @type  readdic: Dictionary
    @param readdic: Boolean which decides if unbalanced mutations are allowed (only initial mutation is performed)
    @rtype:   Readobj
    @return: Parsed read from text line.
    """
        
    # dummy is needed because nomem and standart reading return "different" values.
    # dummy, equals readname
    line = textfile.readline()
    if not line:
        return(0)
    CRead, dummy = readSAMline(line, "noMem", compareList, readdic)
    return(CRead)



def getMisGap(mdtag, cigar):
    """
    Reads the alignment tag given, the text and a tag to search for
    
    @type  mdtag: string
    @param mdtag: MDTag from alignment
    @type  cigar: string
    @param cigar: Cigar from alignment
    @rtype:   gaps,mismatches
    @return: Parsed gaps and mismatches
    """
    deletions = getsum(re.findall("(\d+)D", cigar))
    insertions = getsum(re.findall("(\d+)I", cigar))
    gaps = int(insertions + deletions)

    # mismatches are:
    # ALL Matches from CIGAR - RealMatches from MDTag
    #

    mismatch = getsum(re.findall("(\d+)M", cigar)) - getsum(re.findall("(\d+)", mdtag))

 
    #print ("%s\t%s\t%s\t%s\t%s\t%s" %(mdtag,cigar,deletions,insertions,mismatch,gaps))
    return(gaps, mismatch)

def getAllID(textfile, read, compareList, readdic):
    """
    Reads all alignments for the current read which are in the SAM file (sorted). If a new read ID is scanned the results are returned.
    
    @type  textfile:  fileobject stream
    @param textfile: SAM file for reading.
    @type  read: read obj
    @param read: the last read obj which defines the current read id
    @type  readdic: dictionary
    @param readdic: Look up for quality values.
    @rtype:   (bool,list,bool)
    @return: a "triple", where 2 bools are defined as indicator variables and a list with all alignments for one read.
    """
    # end of file?
    if read == 0:
        #parse last read and return
        CRead = getNextLine(textfile, compareList, readdic)
        return(getAllID(textfile, CRead, compareList, readdic))
    else:
        # start ne list and append read which was read with the last run of getALLID
        templist = []
        templist.append(read)
    # get all reads with same id as "read" 
    while 1:
        line = textfile.readline()
        # no text left infile
        if not line:
            return(0, templist, 0)
        
        try:
           CRead = getCustomRead(HTSeq.SAM_Alignment.from_SAM_line(line), compareList, readdic)
           #CRead, dummy = readSAMline(line, "noMem", compareList, readdic)
        except:
            CRead, dummy = readSAMline(line, "noMem", compareList, readdic)
        
        if CRead == 0:
            #
            return(1, templist, CRead)
        # same read in current line as was in the last line --> continue reading
        if CRead.id == read.id:
            templist.append(CRead)
        else:
            return(1, templist, CRead)
        
    
def CompareAlignments(reflist, artlist,file):
    """
    Compares alignments for the reference and artificial reference for a specific read id.
    The goal is to identify false positives.

    @type  reflist: read obj list
    @param reflist: list containing alignments witht the same ID (reference)
    @type  artlist: read obj list
    @param artlist: list containing alignments witht the same ID (artificial)
    @rtype:   np.array
    @return:  indices of unique alignments
    """
   
    
    artlen = len(artlist)
    reflen = len(reflist)

    TPcount = reflen
    FPcount = 0
    for RefRead in reflist:
         file.write(RefRead.toStrNoMem("ref"))
    
    # add 1 for every alignment from the artificial that is found on the reference
    for i in xrange(0, artlen):
        # bool which decides if artificial in reference hit
        ArtInRef = 0
        for j in xrange(0,reflen):
            if artlist[i].start == reflist[j].start and artlist[i].end == reflist[j].end:
                ArtInRef = 1
                break
        # hit on artificial not confirmed on reference --> FP
        if (ArtInRef == 0):
            FPcount+=1
            file.write(artlist[i].toStrNoMem("art"))
    
    
    return(TPcount,FPcount)
    
"""def WriteToFile(file, list, identifier):
    
    Writes read obj. to a file.


    @type  file: read obj list
    @param file: list containing alignments witht the same ID (reference)
    @type  list: read obj. list
    @param list: list containing alignments witht the same ID (artificial)
    @type  identifier: string
    @param identifier: indicator if the results are from the reference or artificial
    
    file.write(getBestAlignment(list).toStrNoMem(identifier))"""

def SkipHeader(file, compareList, readdic):
    """
    Skips the header from a SAM file, but reads the first line of the alignment section.

    @type  file: read obj list
    @param file: list containing alignments witht the same ID (reference)
    @type  compareList: list
    @param compareList: list for accumulation of the same read id
    @type  readdic: dictionary
    @param readdic: dictionary containing read ID - read quality mappings.
    @rtype:   read obj.
    @return:  Returns a read obj. from the SAM file.
    """
    while 1:
        temp = file.readline()
        if temp.startswith("@"):
            pass
        else:
            return(readSAMline(temp, "noMem", compareList, readdic))
            break


def getRanks(RefRead,ArtRead,rankdic):
    """
    Function which returns the ranks of 2 given readIDs (read from reference,read from artificial).

    @type  RefRead: read obj 
    @param RefRead: read obj  (reference)
    @type  ArtRead: read obj 
    @param ArtRead: read obj  (artificial)
    @type  rankdic: dictionary
    @param rankdic: dictionary containing ranks of the read IDs (according to the sorted SAM files).
    @rtype:   int,int
    @return:  returns the false positives and true positives for a SAM file pair (reference, artificial)
    """
    # get the ranks resulting from the sorted samfile.
    RefReadRank = rankdic[RefRead[0].id]
    if ArtRead[0].id in rankdic:
        ArtReadRank = rankdic[ArtRead[0].id]
    else:
        # if this ReadID is not in the dictionary then it MUST be a FP!
        ArtReadRank = RefReadRank - 1
        
    return(RefReadRank,ArtReadRank)


def ReadSAMnoMem(ref, art, output, compareList, readdic, rankdic):
    """
    Main function for comparing to mappings. This functions takes the complete alignments for artificial and reference genome and goes through them in parallel.
    Since the mappings are sorted the function alternates the parsing of the samfiles in such a way that no memory is used for comparing these functions.

    @type  ref: string
    @param ref: path to reference alignments (SAM file)
    @type  art: string
    @param art: path to artificial alignments (SAM file)
    @type  output: read obj 
    @param output: read obj  (artificial)
    @type  compareList: read obj 
    @param compareList: read obj  (artificial)
    @type  readdic: dictionary
    @param readdic: dictionary containing read ID - read quality mappings.
    @type  rankdic: dictionary
    @param rankdic: dictionary containing ranks of the read IDs (according to the sorted SAM files).
    @rtype:   ranks
    @return:  Returns the ranks for the 2 read ids.
    """
    
    fobj = open(output, "w")
    fobj.write("#ReadID\tMatchedReference\tSubstitutions\tNMtag\tReadQuality\tMappingQuality\tStart\tEnd\tGaps\tMisM\r\n")
    fileref = open(ref, "r")
    fileart = open(art, "r")
    
    # SKip Header section of SAM file
    CurrentReadRef, dummy = SkipHeader(fileref, compareList, readdic)
    CurrentReadArt, dummy = SkipHeader(fileart, compareList, readdic) 
    
    # get first readIDs for looping    
    hasNextLine, ArtRead, NextReadArt = getAllID(fileart, CurrentReadArt, compareList, readdic)
    hasNextLine, RefRead, NextReadRef = getAllID(fileref, CurrentReadRef, compareList, readdic)
    # variables for counting

    tp = 0
    fp = 0
    MappedReads = 0
    i = 0
    start = time.time()           
    #Indicator variables which file is end first
    ArtBool = 0
    RefBool = 0
    while 1:
        i += 1
        if i == 100000:
            print "\t\t%f" % (time.time() - start),
        
        RefReadRank,ArtReadRank = getRanks(RefRead,ArtRead,rankdic)
        
        #case 1)
        # Keep reading lines from artificial if the id is < than the one from the reference
        # this means unique hits to the artificial!
        #while (RefRead[0].id > ArtRead[0].id and hasNextLine != 0):
        while (RefReadRank > ArtReadRank and hasNextLine != 0):
            # in case of multiple hits, set all mr variables to  0 and get the best alignment
            for l in xrange(0, len(ArtRead)):
                ArtRead[l].mr = [0]
            
            # write FP to file
            for hits in ArtRead:
                fobj.write(hits.toStrNoMem("art"))
                fp += 1

            
            # use akku for read lines
            CurrentReadArt = NextReadArt
            #get all the reads with the same ID into ArtRead
            hasNextLine, ArtRead, NextReadArt = getAllID(fileart, CurrentReadArt, compareList, readdic)
            
            #update rank
            ArtReadRank  = RefReadRank-1 if ArtRead[0].id not in rankdic else rankdic[ArtRead[0].id]
            MappedReads+=1
            
        #case 2)
        # keep reading lines from reference if id is < art.id
        # all these hits are purely TP
        #while (RefRead[0].id < ArtRead[0].id and hasNextLine != 0):
        while (RefReadRank < ArtReadRank and hasNextLine != 0):
            
            for hits in RefRead:
                fobj.write(hits.toStrNoMem("ref"))
                tp += 1

            CurrentReadRef = NextReadRef
            hasNextLine, RefRead, NextReadRef = getAllID(fileref, CurrentReadRef, compareList, readdic)
            RefReadRank = rankdic[RefRead[0].id]
            MappedReads+=1
            
            
        # case 3)
        #if (RefRead[0].id == ArtRead[0].id):
        if (RefReadRank == ArtReadRank):
            MappedReads+=1
            
            #same read.id --> compare all Alignments for one read id (can be multiple alignments)
            # return unique hits to artificial genome
            tempTP,tempFP = CompareAlignments(RefRead, ArtRead,fobj)
            tp +=tempTP
            fp +=tempFP
            # overwrite akku
            CurrentReadRef = NextReadRef
            CurrentReadArt = NextReadArt
            
            
            # get next Read for next ID and check if file has a next line
            hasNextLine, RefRead, NextReadRef = getAllID(fileref, CurrentReadRef, compareList, readdic)
            if hasNextLine == 0:
                #print "REFBREAK"
                RefReadRank = rankdic[RefRead[0].id]
                RefBool = 1
                break
            
            hasNextLine, ArtRead, NextReadArt = getAllID(fileart, CurrentReadArt, compareList, readdic)
            if hasNextLine == 0:
                #print "ARTBREAK"
                ArtReadRank  = RefReadRank-1 if ArtRead[0].id not in rankdic else rankdic[ArtRead[0].id]
                ArtBool = 1
                break

    # if the while loop was broken one read has to be read in the not closed file
    hasNextLine = 1
    # the file which has still some lines is now processed linearly...
    if RefBool == 1:
        #print ("REF FILE FINISHED!!!")
        # reference no lines --> read artificial
        hasNextLine, ArtRead, NextReadArt = getAllID(fileart, CurrentReadArt, compareList, readdic)
        ArtReadRank  = RefReadRank-1 if ArtRead[0].id not in rankdic else rankdic[ArtRead[0].id]
        
        while (RefReadRank > ArtReadRank and hasNextLine != 0):
            MappedReads+=1
            # in case of multiple hits, set all mr variables to  0 and get the best alignment
            for l in xrange(0, len(ArtRead)):
                ArtRead[l].mr = [0]
    
             # write current FP to file
             # write FP to file
            for hits in ArtRead:
                fobj.write(hits.toStrNoMem("art"))
                fp += 1
            # use akku for read lines
            CurrentReadArt = NextReadArt
            #get all the reads with the same ID into ArtRead
            hasNextLine, ArtRead, NextReadArt = getAllID(fileart, CurrentReadArt, compareList, readdic)
    
            #update rank
            ArtReadRank  = RefReadRank-1 if ArtRead[0].id not in rankdic else rankdic[ArtRead[0].id]
            
        
        
        if (RefReadRank == ArtReadRank):
            MappedReads+=1
            #same read.id --> compare all Alignments for one read id (can be multiple alignments)
            # return unique hits to artificial genome
            tempTP,tempFP = CompareAlignments(RefRead, ArtRead,fobj)
            tp +=tempTP
            fp +=tempFP
            
            
        while 1:
            
            if hasNextLine == 0:
                break
            else:
               
                CurrentReadArt = NextReadArt
                # read next artificial read!
                hasNextLine, ArtRead, NextReadArt = getAllID(fileart, CurrentReadArt, compareList, readdic)
                
            
            for hits in ArtRead:
                fobj.write(hits.toStrNoMem("art"))
                fp += 1
            MappedReads+=1
            
    else:
        # read again reference file as long as the final read from he artificial file is reached
        while (RefReadRank < ArtReadRank and hasNextLine != 0):
            MappedReads+=1
        
            for hits in RefRead:
                fobj.write(hits.toStrNoMem("ref"))
                tp += 1
            CurrentReadRef = NextReadRef
            hasNextLine, RefRead, NextReadRef = getAllID(fileref, CurrentReadRef, compareList, readdic)
            RefReadRank = rankdic[RefRead[0].id]
        
        
        if (RefReadRank == ArtReadRank):
            MappedReads+=1
            #same read.id --> compare all Alignments for one read id (can be multiple alignments)
            # return unique hits to artificial genome
            tempTP,tempFP = CompareAlignments(RefRead, ArtRead,fobj)
            tp +=tempTP
            fp +=tempFP

        
        else:
            MappedReads+=1
            for hits in ArtRead:
                fobj.write(hits.toStrNoMem("art"))
                fp += 1
        
        while 1:
                if hasNextLine == 0:
                    MappedReads+=1
                    for hits in RefRead:
                        fobj.write(hits.toStrNoMem("ref"))
                        tp += 1
                    
                    break
                else:
                    CurrentReadRef = NextReadRef
                    MappedReads+=1
                    for hits in RefRead:
                        fobj.write(hits.toStrNoMem("ref"))
                        tp += 1
                hasNextLine, RefRead, NextReadRef = getAllID(fileref, CurrentReadRef, compareList, readdic)

    end = time.time()
    print ("\t\t%f\t\t" %(end - start)),
    fobj.close()
    return(tp, fp,tp+fp,MappedReads)
##################END: FUNCTIONS FOR NO-MEMORY ANALYSIS ####################


def getNumberOf(line, tag):
    """
    Reads the alignment tag, given the text and a tag to search for.

    @type  line: string
    @param line: SAM line
    @type  tag: string
    @param tag: SAM tag. i.e: NM,MD
    @rtype:   int
    @return:  number x behind desired tag tag:i:x
    """

    # gives the number of XX back
    # valid tags: X0,XM (number of missmatches
    m = re.search(tag + ":i:(\d+)", line)
    besthits = m.group(1)
      
    return(besthits)

def getsum(strlist):
    """
    Sum of strings (as ints)

    @type  strlist: list(str)
    @param strlist: SAM line

    @rtype:   int
    @return:  MD tag calculation.
    """

    return(sum([int(i) for i in strlist]))

def getmean(strlist):
    """
    Mean of strings.

    @type  strlist: list(str)
    @param strlist: SAM line

    @rtype:   float
    @return:  mean
    """
    sum = getsum(strlist)
    return(sum/float(len(strlist)))


def ComputeRQScore(quality):
    """
    Computes the ReadQualityScore given a quality string

    @type  quality: string
    @param quality: quality string of a read.

    @rtype:   float
    @return:  ReadQualityScore (RQS)
    """
    # phred --> -33
    qualarray = [ord(i)-33 for i in quality]
    # solexa --> -64
    score =( sum(qualarray)/float(len(qualarray))) / (max(qualarray) -  min(qualarray) + 1)

    return (score)


def getAlignmentLength(cigar):
    """
    Computes the alignment length given a cigar string. Needed for Start + End calculation.

    @type  cigar: string
    @param cigar: Cigar string (SAM)

    @rtype:   int
    @return:  alignmentlength
    """
    # I - insertions do not extend the alignment
    # S - does not extend...
    # P - does not extend...
    # H - does not extend...
    match = re.findall("(\d+)[MDN]", cigar)
    sumOfMatches = 0
    for item in match:
        sumOfMatches += int(item)
    return (sumOfMatches)


def isSaneAlignment(alignment, identifier, compareList, readdic):
    """
    Checks alignment line for unnecessary informations to skip 

    @type  alignment: string
    @param alignment: Line from SAM
    @type  identifier: string
    @param identifier: read id
    @type  compareList: list
    @param compareList: list of alignments with same read id.
    @type  readdic: dictionary
    @param readdic: Dictionary containg a read id; read quality mapping
    @rtype:   readobj, readname
    @return:  Returns the read and it's identifier.
    """
    
    if alignment.startswith("@"):
        return(0, "")
    else:
        read, readname = readSAMline(alignment, identifier, compareList, readdic)
        if read == 0:
            return(0, "")
        else:
            return (read, readname)


def CheckForSameAlignments(readref, readart):
    """
    Function for comparison of artificial alignments and reference alignments. FP are defined such that start and end position must be unique to the artificial reference
    returns 0 if no same read is found (FP found)
    returns 1 if an equal alignment is found

    @type  readref: read obj.
    @param readref: reference
    @type  readart: read obj.
    @param readart: artificial
    @rtype:   bool
    @return:  Indicator if alignment is the same (start & end equal)
    """
    for i, item in enumerate(readref.start):
        if (readref.start[i] == readart.start[0] and readref.end[i] == readart.end[0]):
            return(1)
    return(0)



def getHammingdistance(CompareString, start, end):
    """
    Computes the number of subsitutions in the artificial reference using the CompareString.

    @type  CompareString: string
    @param CompareString: string of 0 and 1s. 1 = hamming 1 between reference and artificial.
    @type  start: int
    @param start: start of alignment
    @type  end: int
    @param end: end ofalignment
    @rtype:   int
    @return:  hamming distance
    """
    return (sum(CompareString[start:end]))
    

def readReadQualities(fastqfile):
    """
    
    Reads a .fastqfile and calculates a defined readscore 
    input: fastq file
    output: fastq dictionary key = readid; value = qualstr 


    @type  fastqfile: string
    @param fastqfile: path to fastq file

    @rtype:   dictionary
    @return:  dictionary containing read ids and read qualities.
    """
    fastq_file = HTSeq.FastqReader(fastqfile , "phred")
    readdictionary = {}
    
    for read in fastq_file:
        readdictionary[read.name] = ComputeRQScore(read.qualstr)
    print("\tReading Fastq file done!")
    return readdictionary




def extendReadDic(readdic):
    """
    Extends a given dictionary with KEY = readid and VALUE = qualstr such that an internal naming is generated which can be used to efficiently create an numpy array


    @type  readdic: dictionary
    @param readdic: dictionary containing read ids and read qualities.

    @rtype:   dictionary
    @return:  extended readdic with KEY = ID, VALUE = READID object with READID.internalid and READID.qulstr = qualstr
    """
    internalnaming = 0
    reverseReadDic = np.zeros(len(readdic), dtype='S50')
    for id in readdic.iterkeys():
        # i
        readdic[id] = ReadID(internalnaming, readdic[id])
        reverseReadDic[internalnaming] = id
        internalnaming += 1
    return(readdic, reverseReadDic)

         

def returnSequence(fasta):
    """
    Returns a sequence string from a fasta file.


    @type  fasta: string
    @param fasta: path to fasta file.

    @rtype:   string
    @return:  sequence
    """
    fastafile = HTSeq.FastaReader(fasta)
    for sequence in fastafile:
        return(sequence.seq)
    
    
   
def CreateCompareList(Reference, ARG):
    """
    Creates a list which is used for comparisons between aligned reads (exact number of mismatches) 



    @type  Reference: string
    @param Reference: reference genome
    @type  ARG: string
    @param ARG: artificial reference genome.
    @rtype:   list
    @return:  list containt 1s, where there is a difference in the genomes and 0s where the nucleotides are equal.
    """
    
    reference = returnSequence(Reference)
    artificialreference = returnSequence(ARG)
    complist = []
    
    
    if (len(reference) != len(artificialreference)):
        print "first 10 letter:\t\t" + reference[0:10] + "..." + reference[-10:]
        print "last 10 letter:\t\t" + artificialreference[0:10] + "..." + artificialreference[-10:]
        print ("Error! Two Sequences have different length! Try to add a line break after the last nucleotide.")
        sys.exit(1)
    else:
        return([ 0 if x == artificialreference[i] else 1 for i, x in enumerate(reference)])
    


def readSAMline(alignment, identifier, compareList, readdic):
    """
    Function for reading SAM alignment text file (one line)



    @type  alignment: string
    @param alignment: SAM alignment
    @type  identifier: string
    @param identifier: read id
    @type  compareList: list
    @param compareList: list containt 1s, where there is a difference in the genomes and 0s where the nucleotides are equal.
    @type  readdic: dictionary
    @param readdic: dictionary containing read ids and read qualities.
    @rtype:   read obj.
    @return:  returns a customRead object
    """
    k = 0

    # ignore header
    columns = alignment.split("\t")
    flag = int(columns[1].strip())
    
    # flag 4 = not mapped to any genome
    if flag == 4:
        return(0, "")
        
    else:
       
        
        readname = columns[0].strip()
        reference = columns[2].strip()
        start = int(columns[3].strip())-1
        mappingquality = columns[4].strip()
        cigar = columns[5].strip()
        qualstr = columns[10].strip()
        tags = " ".join(columns[11:])
        
        mdtag = re.search("MD:Z:([^\s]+)", tags).group(1)
        try:
            mdtag = re.search("MD:Z:([^\s]+)", tags).group(1)
            gaps, mism = getMisGap(mdtag, cigar)
        except:
            gaps = 0
            mism = 0
            print "Error Reading MDtag of %s.Setting gaps = 0,mism = 0" % readname
        try:
            nm = int(getNumberOf(tags, "NM"))
        except:
            nm = 0
        
        leng = getAlignmentLength(cigar)
        
        if len(qualstr) == 1 or (len(qualstr.strip()) == 0):
            score = readdic[readname].quality
        else:
            score = ComputeRQScore(qualstr)
        
        if identifier == "art":
            #tempobj = CustomRead(0,nm, getHammingdistance(compareList, start, start+leng), score, mappingquality,start,start+leng)
            tempobj = CustomRead(readname, 0, nm, getHammingdistance(compareList, start, start + leng), score, mappingquality, start, start + leng, gaps, mism)
        elif identifier == "noMem":
            tempobj = CustomRead(readname, 1, nm, getHammingdistance(compareList, start, start + leng), score, mappingquality, start, start + leng, gaps, mism)
            
        else:
            ##########TPRead(nm, score, mq, start, end, gaps, mism)
            tempobj = TPRead(nm, score, mappingquality, start, start + leng, gaps, mism)



    return(tempobj, readname)





def returnIndex(readdic, readname):
    """
    Returns the index of a read. The index is prescribed by the ordering in the sam file.
    @type  readname: string
    @param readname: read id
    @type  readdic: dictionary
    @param readdic: dictionary containing read ids and read qualities.ArtRead
    @rtype: int
    @return:  index
    """
    return(readdic[readname].internalID)


def ReadArtificialSAMfileHTSeq(art, compareList, RefArray, readdic):
    """
    Function for reading the artificial reference genome using HTSeq.This function is mainly used. Only if no quality string is in the SAM line. The custom
    SAM reading function is used.
    
    @type  art: string
    @param art: artificial file.
    @type  RefArray: array
    @param RefArray: Results from reading the reference SAM file.
    @type  compareList: list
    @param compareList: list containt 1s, where there is a difference in the genomes and 0s where the nucleotides are equal.
    @type  readdic: dictionary
    @param readdic: dictionary containing read ids and read qualities.
    @rtype:   array
    @return:  aligned read objects in an array.
    """
    start = time.time()
    #print ("\tARTIFICIAL:")
    fobj = open(art, "r")
    artdic = {}
    k = 0
    
    read = SkipHeader(fobj,compareList,readdic)
    
    
    for alignment in fobj:
        k += 1
        if k % 1000000 == 0:
            print ("%d.." %(k/1000000)),
        read, readname = isSaneAlignment(alignment, "art", compareList, readdic)
        if read == 0:
            pass
        else:
            # Check by internal naming if spot in array is taken (!= 0) by a read
            index = returnIndex(readdic, readname)
            if RefArray[index] != 0:
                # read already in dic? Check if the alignments are the same
                if (RefArray[index].isContained(read) == 0):
                    if readname in artdic:
                        artdic[readname].toObjself(read)
                    
                    else:
                        artdic[readname] = read
                else:
                    pass
                
            else:
                # new read? Just add it to the dictionary
                artdic[readname] = read
           
 
    fobj.close()
    end = time.time()
    #print("\r\n")
    print ("\t %f " % (end - start)),
    #print ("\tdone in %d seconds" % (end-start))
    return(artdic)



def writeToTabArt(ReadDic, outfile):
    """
    Write function for hits ont he reference genome. Only the best alignment (least mismatches)
    Header :
    ReadID    MatchedReference     Substitutions     NumberOfMismatches     ReadQuality     MappingQuality
    

    @type  outfile: string
    @param outfile: Path to outfile.
    @type  ReadDic: dictionary
    @param ReadDic: dictionary containing read ids and read qualities.
    """

    fobj = open(outfile, "a")
    for read in ReadDic.keys():
       BestReadIndex = np.array(ReadDic[read].nm).argmin()
       fobj.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n" % (read, ReadDic[read].mr[BestReadIndex], ReadDic[read].subs[BestReadIndex], ReadDic[read].nm[BestReadIndex], ReadDic[read].rq, ReadDic[read].mq[BestReadIndex], ReadDic[read].start[BestReadIndex], ReadDic[read].end[BestReadIndex], ReadDic[read].gaps[BestReadIndex], ReadDic[read].mism[BestReadIndex]))
    fobj.close()




def writeToTabRef(RefArray, outfile, reverseReadArray, ReadDic):
    """ 
    Write function for hits ont he reference genome. Only the best alignment (least mismatches)
    
    Header :
    ReadID    MatchedReference     Substitutions     NumberOfMismatches     ReadQuality     MappingQuality
    rows in the infile
    @type  RefArray: string
    @param RefArray: path to inputfile.
    @type  outfile: int
    @param outfile: rows in the infile
    @type  reverseReadArray: dictionary
    @param reverseReadArray: Contains reads = values and ranks = keys
    @type  ReadDic: dictionary
    @param ReadDic: dictionary containing read ids and read qualities.
    """

    fobj = open(outfile, "a")
    for i in xrange(0, len(RefArray)):
        read = RefArray[i]
        # check if entry == 0
        if read == 0:
            continue
        else:
            readname = reverseReadArray[i]
            index = returnIndex(ReadDic, readname)
            BestReadIndex = np.array(RefArray[index].nm,dtype="int").argmin()
            fobj.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n" % (readname, 1, 0, RefArray[index].nm[BestReadIndex], RefArray[index].rq, RefArray[index].mq[BestReadIndex], RefArray[index].start[BestReadIndex], RefArray[index].end[BestReadIndex], RefArray[index].gaps[BestReadIndex], RefArray[index].mism[BestReadIndex]))
    fobj.close()



def ReadFromTab(infile, arraysize):
    """
    Reads the results from writeToTab in for plotting and analysis
    
    Header :
    ReadID    MatchedReference     Substitutions     NumberOfMismatches     ReadQuality     MappingQuality
    
    @type  infile: string
    @param infile: path to inputfile.
    @type  arraysize: int
    @param arraysize: rows in the infile

    @rtype:   np.array(obj)
    @return:  array for unique classified reads.
    """

    
    fobj = open(infile, "r")
    dataArray = np.arange(arraysize, dtype=object)
    for i, line in enumerate(fobj):     
            if line.startswith("#"):
                continue
            else:
                zeile = line.split("\t") 
                id = zeile[0]
                mr = int(zeile[1])
                subs = int(zeile[2])
                nm = int(zeile[3])
                rq = float(zeile[4])
                mq = int(zeile[5])
                gaps = int(zeile[8])
                mism = int(zeile[9])
                #n = int(zeile[10].strip())
                # id is not needed for the analysis 
                dataArray[i - 1] = CustomRead(id, mr, nm, subs, rq, mq, 0, 0, gaps, mism)
                
    return (dataArray)



def readInput(file):
    """
    Function for reading in the controldictionary from the inputfile.

    @type  file: string 
    @param file: Path to the input file.
    
    @rtype:   files
    @return:  Controldictionary which is used to regulate the programs workflow.
    
    """
    fobj = open(file, "r")
    files = {}
    files["mapper"] = {}
    files["fasta"] = {}
    tempref = []
    temparg = []
    reffasta = ""
    artfasta = []
    
    for line in fobj:
        # reference
        if line.startswith("$"):
            files["fasta"]["reffasta"] = GetCompletePath(line.split(":")[1].strip())
        # artificial
        if line.startswith("#"):
            artfasta.append(GetCompletePath(line.split(":")[1].strip()))
        # fastq file
        if line.startswith("&"):
            files["fastqfile"] = GetCompletePath(line.split(":")[1].strip())
        # start sequence for mapper section
        if line.startswith("@"):
            zeile = line.split("\t")
            mapper = zeile[0][1:].strip()
            files["fasta"]["artfasta"] = artfasta
        if line.startswith("ref"):
            tempref.append(GetCompletePath(line.split(":")[1].strip()))
        if line.startswith("art"):
            temparg.append(GetCompletePath(line.split(":")[1].strip()))
        # end of mapper section
        if line.startswith("+"):
            files["mapper"][mapper] = [tempref, temparg]
            tempref = []
            temparg = []
            
            
    fobj.close()
    return(files)

def GetCompletePath(path):
    return(path.replace("~",os.path.expanduser("~")))
    
def initOutFiles(controlDic, mapper, outpath):
    """
    Since we append in the program, we need to make sure no old files remain... 

    @type  controlDic: dictionary 
    @param controlDic: dictionary containing future filenames.
    @type  mapper: string 
    @param mapper: current identifier mapper for which the results are written
    @type  outpath: string
    @param outpath: existing path, where the outfiles will be written.
    """
    for artificial in controlDic["mapper"][mapper][1]:
             filename = artificial.split("/")[-1]
             outfile = mapper + "_" + artificial.split("/")[-1][:-4] + ".esam"
             #easy Sam output 
             fobj = open(outpath + outfile, "w")
             fobj.write("#ReadID\tMatchedReference\tSubstitutions\tNumberOfMismatches\tReadQuality\tMappingQuality\tStart\tEnd\tGaps\tMisM\r\n")
             fobj.close()
     
