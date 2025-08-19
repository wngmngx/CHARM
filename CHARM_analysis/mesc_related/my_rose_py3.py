# re-implementation of my_rose.py from Cai et al. 

import os, sys
import numpy as np
import pysam
import argparse
import re
from collections import defaultdict
import subprocess
import shlex

def parseArgs():
    parser = argparse.ArgumentParser(description="Stitch peaks and call super-enhancer.")

    parser.add_argument("input",
                        help="The input file in BED format.")
    parser.add_argument("ref",
                        help="The reference build. Can be [hg19, hg18, mm9, mm8].")
    parser.add_argument("bam",
                        help="The bam file to calculate the signal.")
    parser.add_argument("-c", "--control", dest="control",
                        help="The control bam file to calculate background signal.")
    parser.add_argument("-o", "--output", dest="output", default="New",
                        help="The output prefix, default is 'New'.")
    parser.add_argument("-f", "--fraction", dest="fraction", default=0, type=float,
                        help="The overlapping fraction (only those with > fraction will be retained.")
    parser.add_argument('-w', '--window', dest='window', default=12500, type=int,
                        help="The window for stitching peaks. Default: 12500")
    parser.add_argument('-e', '--exclude', dest='exclude', default=2000, type=int,
                        help="The size of exclusion zones (on each side) around the TSS. Default: 2000")
    parser.add_argument('--separate', dest='separate', default=False, type=lambda x: (str(x).lower() == 'true'), 
            help="Whether to separate Super-enhancers and typical-enhancers based on their proximity to TSSes. Default: False:")

    args = parser.parse_args()
    return args


class Locus:
    def __init__(self, chrom, start, end, strand, name, commonName=None, signal=None, ctrlSignal=None):
        '''
        strand should be one of ["+", "-", "."]
        '''
        self.__chrom = chrom
        self.__start = int(start)
        self.__end = int(end)
        if strand not in ['+', '-', '.']:
            raise ValueError("strand should be one of [+, -, .]")
        self.__strand = strand
        self.__name = name
        self.__commonName = commonName
        self.__signal = signal  # Average signal
        self.__ctrlSignal = ctrlSignal  # Average ctrl signal

    def getStart(self):
        return self.__start

    def getEnd(self):
        return self.__end

    def getLen(self):
        return self.getEnd() - self.getStart()

    def getChrom(self):
        return self.__chrom

    def getStrand(self):
        return self.__strand

    def getName(self):
        return self.__name

    def setCommonName(self, cn):
        self.__commonName = cn

    def getCommonName(self):
        return self.__commonName

    def overlap(self, otherLocus, forceStrand=False, otherFraction=0):
        if self.getChrom() != otherLocus.getChrom():
            return False
        elif forceStrand and self.getStrand() != otherLocus.getStrand():
            return False
        else:
            # Ensure floating point division for Python 3
            return (min(self.getEnd(), otherLocus.getEnd()) - max(self.getStart(), otherLocus.getStart())) * 1.0 / otherLocus.getLen() > max(otherFraction, 0)

    def setSignal(self, signal):
        '''
        The average signal
        '''
        self.__signal = signal

    def setCtrlSignal(self, ctrlSignal):
        self.__ctrlSignal = ctrlSignal


class StitchedLocus:
    def __init__(self, loci, window, signal=None, ctrlSignal=None, consSignal=None, consCtrlSignal=None):
        if len(loci) <= 0:
            raise ValueError("No Locus")
        self.__loci = loci
        self.__signal = signal
        self.__ctrlSignal = ctrlSignal
        self.__consSignal = consSignal
        self.__consCtrlSignal = consCtrlSignal
        self.__window = window

    def getLoci(self):
        return self.__loci

    def getSignal(self): # This was duplicated, kept one
        return self.__signal

    def getStart(self):
        return self.__loci[0].getStart()

    def getEnd(self):
        return self.__loci[-1].getEnd()

    def setSignal(self, signal):
        self.__signal = signal

    def setCtrlSignal(self, ctrlSignal):
        self.__ctrlSignal = ctrlSignal

    def setConsSignal(self, consSignal):
        self.__consSignal = consSignal

    def setConsCtrlSignal(self, consCtrlSignal):
        self.__consCtrlSignal = consCtrlSignal

    def getConsSignal(self): # Removed unused 'consSignal' argument
        return self.__consSignal

    def getTotalConsSize(self):
        totalSize = 0
        for l in self.getLoci():
            totalSize += l.getLen()
        return totalSize

    def getCtrlSignal(self):
        return self.__ctrlSignal

    def getNumLoci(self):
        return len(self.__loci)

    def getTotalConsSignal(self):
        if self.__consSignal is None: # Changed == None to is None
            return 0
        else:
            return sum(self.__consSignal)

    def getTotalConsCtrlSignal(self):
        if self.__consCtrlSignal is None: # Changed == None to is None
            return 0
        else:
            return sum(self.__consCtrlSignal)

    def getConsCtrlSignal(self):
        return self.__consCtrlSignal

    def getName(self):
        return "%d_%s_%dk_stitched" % (self.getNumLoci(), self.getLoci()[0].getName(), int(self.__window / 1000))


def loadBed(filename, header=False):
    '''
    Returns a dictionary of sorted lists of locis.
    The keys of the dictionary are the chromosome names.
    '''
    loci = {}
    # Use with statement for file handling
    with open(filename, 'r') as f: # Added 'r' for read mode
        if header:
            f.readline()
        for r in f:
            tokens = r.strip().split('\t')
            if tokens[0] not in loci:
                loci[tokens[0]] = []
            loci[tokens[0]].append(Locus(tokens[0], tokens[1], tokens[2], tokens[5], tokens[3]))
    for chrom in loci:
        loci[chrom].sort(key=lambda k: (k.getStart(), k.getEnd()))
    return loci

def loadTSS(filename, header=False, extension=0):
    '''
    Returns a dictionary of sorted lists of locis of TSS.
    The keys of the dictionary are the chromosome names.
    '''
    loci = {}
    # Use with statement for file handling
    with open(filename, 'r') as f: # Added 'r' for read mode
        if header:
            f.readline()
        for r in f:
            tokens = r.strip().split('\t')
            if not re.match(r'^chr(\d+|X|Y)$', tokens[2]): # Added r for raw string
                continue
            if tokens[2] not in loci:
                loci[tokens[2]] = []
            if tokens[3] == '-':
                loci[tokens[2]].append(Locus(tokens[2], int(tokens[5]) - 1 - extension, int(tokens[5]) + extension, tokens[3], tokens[1], commonName=tokens[12]))
            else:
                loci[tokens[2]].append(Locus(tokens[2], int(tokens[4]) - extension, int(tokens[4]) + 1 + extension, tokens[3], tokens[1], commonName=tokens[12]))
    for chrom in loci:
        loci[chrom].sort(key=lambda k: (k.getStart(), k.getEnd()))
    print(list(loci.keys())) # Changed print loci.keys() to print(list(loci.keys())) for Python 3
    return loci

def getStartEndLists(chromLoci):
    '''
    Given a list of loci on the same chromosome sorted by coordinates,
    return one list of the start coordinates and
    one list of end coordinates.
    Note that there should not be a loci in the list being included
    in another loci in the list.
    '''
    starts = []
    ends = []
    for l in chromLoci:
        starts.append(l.getStart())
        ends.append(l.getEnd())
    return starts, ends

def testOverlap(locus, refLoci, starts, ends, fraction):
    '''
    Test whether locus is overlapping with any of refLoci by fraction of the locus.
    The starts and ends are the outputs of getStartEndLists().
    Returns a boolean result
    '''
    import bisect
    leftIdx = bisect.bisect_left(ends, locus.getStart())
    rightIdx = bisect.bisect_right(starts, locus.getEnd())
    if rightIdx <= leftIdx:
        return False
    else:
        for i in range(leftIdx, rightIdx):
            if refLoci[i].overlap(locus, otherFraction=fraction):
                return True
        return False


def removeLociOverlapRef(loci, ref, fraction=0):
    # Use set intersection correctly
    chroms = set(loci.keys()).intersection(set(ref.keys()))
    result = {}
    removed = {}
    for chrom in chroms:
        refStarts, refEnds = getStartEndLists(ref[chrom])
        result[chrom] = []
        removed[chrom] = []
        for locus in loci[chrom]:
            if not testOverlap(locus, ref[chrom], refStarts, refEnds, fraction):
                result[chrom].append(locus)
            else:
                removed[chrom].append(locus)
    return result, removed


def stitchLoci(loci, stitchWindow=0, tss=None):
    '''
    Stitch loci, and if tss is given, will not span any tss.
    That is loci at different sides of tss will not be stitched.
    '''
    chroms = set(loci.keys())
    if tss is not None: # Changed tss != None to tss is not None
        chroms = chroms.intersection(set(tss.keys()))
    stitched = {}
    for chrom in chroms:
        if not loci.get(chrom): # Handle case where chrom might not have loci after intersection
            stitched[chrom] = []
            continue
        chromStitched = []
        stitched[chrom] = []
        cLoci = loci[chrom]
        if not cLoci: # Handle empty cLoci list
            continue

        tempStitched = [cLoci[0], ]
        lastEnd = cLoci[0].getEnd()
        for locus in cLoci[1:]:
            if locus.getStart() < lastEnd + stitchWindow:
                tempStitched.append(locus)
            else:
                chromStitched.append(StitchedLocus(tempStitched, stitchWindow))
                tempStitched = [locus, ]
            lastEnd = locus.getEnd()
        if tempStitched: # Ensure last group is added
            chromStitched.append(StitchedLocus(tempStitched, stitchWindow))

        if tss is not None and chrom in tss: # Check if chrom in tss after intersection
            cTss = tss[chrom]
            i = 0
            j = 0
            while i < len(chromStitched) and j < len(cTss):
                if chromStitched[i].getEnd() <= cTss[j].getStart():
                    stitched[chrom].append(chromStitched[i])
                    i += 1
                elif chromStitched[i].getStart() >= cTss[j].getEnd():
                    j += 1
                else:
                    for locus_item in chromStitched[i].getLoci(): # Renamed locus to locus_item to avoid conflict
                        stitched[chrom].append(StitchedLocus([locus_item, ], stitchWindow))
                    i += 1
            # Add remaining stitched loci if any
            while i < len(chromStitched):
                stitched[chrom].append(chromStitched[i])
                i += 1
        else:
            stitched[chrom] = chromStitched
    return stitched


def calculateSignal(stitched, rankby_path, extension, isCtrl=False): # Renamed rankby to rankby_path
    # Pysam uses AlignmentFile in newer versions
    # Ensure the BAM file has an index (.bai)
    with pysam.AlignmentFile(rankby_path, "rb") as sam: # Use with statement
        # sam.mapped might not be directly available or might be sam.get_index_statistics().mapped
        # For simplicity, assuming total reads might need to be calculated differently if sam.mapped fails
        try:
            total_mapped_reads = sam.mapped
            if total_mapped_reads == 0: # Avoid division by zero if no mapped reads
                print(f"Warning: No mapped reads found in {rankby_path}. Signals will be zero.")
                MMR = 0.0
            else:
                MMR = 1e6 * 1.0 / total_mapped_reads
        except AttributeError:
            print(f"Warning: sam.mapped not available for {rankby_path}. Trying to count reads. This might be slow.")
            # Fallback: count reads (can be slow for large BAMs)
            total_mapped_reads = 0
            for _ in sam.fetch(): # Iterate over all mapped reads
                total_mapped_reads +=1
            sam.reset() # Reset iterator for subsequent fetches
            if total_mapped_reads == 0:
                print(f"Warning: No mapped reads found after counting in {rankby_path}. Signals will be zero.")
                MMR = 0.0
            else:
                 MMR = 1e6 * 1.0 / total_mapped_reads

        print(rankby_path, MMR, total_mapped_reads) # Changed print
        ticker = 0
        processID = 1 # This was for multiprocessing, can be removed if not used

        for chrom in stitched:
            for s in stitched[chrom]:
                ticker += 1
                if ticker % 1000 == 0:
                    print("process %d: %d" % (processID, ticker)) # Changed print

                # Pysam fetch might need `contig` instead of `reference`
                # Ensure start < end for pysam.fetch
                fetch_start = max(0, s.getStart() - extension) # Ensure start is not negative
                fetch_end = s.getEnd()
                if fetch_start >= fetch_end:
                    tags = iter([]) # Return empty iterator if region is invalid
                else:
                    try:
                        tags = sam.fetch(contig=chrom, start=fetch_start, end=fetch_end)
                    except ValueError as e:
                        print(f"Warning: Could not fetch region {chrom}:{fetch_start}-{fetch_end}. Error: {e}")
                        tags = iter([])


                changes = defaultdict(int)
                read_names_processed = set() # To ensure each read is counted once per stitched locus

                for t in tags:
                    if t.is_unmapped or t.query_name in read_names_processed: # Skip unmapped or already processed reads
                        continue
                    read_names_processed.add(t.query_name)

                    # Use reference_start and reference_end for positions
                    # query_alignment_length is the length of the alignment on the query
                    # infer_read_length() might be more robust for read length on reference
                    read_length_on_ref = t.infer_read_length()
                    if read_length_on_ref is None: # Should not happen for mapped reads
                        read_length_on_ref = t.query_length # Fallback

                    if not t.is_reverse:
                        # Read start pos
                        tag_start = t.reference_start
                        # Read end pos (exclusive) for signal calculation; tag end is start + length
                        tag_signal_end = t.reference_start + extension
                        changes[tag_start] += 1
                        changes[tag_signal_end] -= 1
                    else:
                        # For reverse strand, effective start is end of read - extension
                        tag_signal_start = t.reference_end - extension # reference_end is exclusive
                        # Effective end is end of read
                        tag_end = t.reference_end
                        changes[tag_signal_start] += 1
                        changes[tag_end] -= 1

                totalSig = 0
                tracker = 0
                consSignals = []
                tempConsSig = 0
                consI = 0

                # Determine the range for iteration carefully
                # It should cover all coordinates where 'changes' occur and the stitched locus itself
                all_coords = list(changes.keys())
                min_coord = s.getStart()
                max_coord = s.getEnd()
                if all_coords:
                    min_coord = min(min(all_coords), s.getStart())
                    max_coord = max(max(all_coords), s.getEnd())


                if len(changes.keys()) > 0 : # Check if there are any changes to process
                    # Iterate from the minimum coordinate affected by reads up to the end of the stitched locus
                    # The original code iterated from min(changes.keys()) or s.getStart().
                    # Ensure the iteration covers the full range where signal contributes.
                    iter_start = min(min(changes.keys()) if changes else s.getStart(), s.getStart())

                    for i in range(iter_start, max_coord +1): # Iterate up to max_coord where signal might change
                        if i in changes: # Tracker changes only at specific points
                             tracker += changes[i]

                        if i >= s.getStart() and i < s.getEnd(): # Accumulate signal within the stitched locus
                            totalSig += tracker
                            # Constituent signal calculation
                            if consI < len(s.getLoci()):
                                current_locus = s.getLoci()[consI]
                                if current_locus.getStart() <= i < current_locus.getEnd():
                                    tempConsSig += tracker
                                elif i >= current_locus.getEnd():
                                    # Store raw sum of signals for constituent, MMR will be applied later if needed
                                    consSignals.append(tempConsSig) # Original code multiplied by MMR here
                                    tempConsSig = 0
                                    consI += 1
                                    # Check if the new constituent locus starts at current 'i'
                                    if consI < len(s.getLoci()):
                                        current_locus = s.getLoci()[consI]
                                        if current_locus.getStart() <= i < current_locus.getEnd():
                                            tempConsSig += tracker
                    # Add the last constitutional signal
                    if consI < len(s.getLoci()): # and tempConsSig > 0: # only if there's accumulated signal
                        consSignals.append(tempConsSig) # Original code multiplied by MMR here

                # Normalize signals
                # The original code for 'den' was totalSig * MMR / (s.getEnd() - s.getStart())
                # And for consSignals, it was tempConsSig * MMR / (locus_len)
                # Here, we store sums and apply MMR and length normalization at the end if that's the intent.
                # The current script stores total signal (not density) in setSignal
                # and sum of signals for constituents in setConsSignal.
                # The division by length for density was done in the commented out calculateSignalChrom
                # Let's stick to the logic in the active part of the original script: sum * MMR

                final_signal = totalSig * MMR
                final_cons_signals = [cs * MMR for cs in consSignals]


                if not isCtrl:
                    s.setSignal(final_signal)
                    s.setConsSignal(final_cons_signals)
                else:
                    s.setCtrlSignal(final_signal)
                    s.setConsCtrlSignal(final_cons_signals)
    # No need to close sam, 'with' statement handles it


def writeStitched(outfilename, stitched, rankbyName, ctrlName="NONE"):
    '''
    Write the data for ROSE_callSuper.R
    REGION_ID    CHROM    START    STOP    NUM_LOCI    CONSTITUENT_SIZE    SNU16_H3K27AC_merged.bam    SNU16_Total_Input_ATCACG_lane8_read1
    '''
    header = ["REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE", rankbyName, ]
    if ctrlName != "NONE":
        header.append(ctrlName)
    header.append("TOTAL_CONS_SIGNAL")
    if ctrlName != "NONE":
        header.append("TOTAL_CONS_CTRL")

    # Use with statement for file handling
    with open(outfilename, 'w') as out:
        out.write('\t'.join(header))
        out.write('\n')
        for chrom in stitched:
            for s in stitched[chrom]:
                line = [s.getName(), chrom, s.getStart(), s.getEnd(), s.getNumLoci(), s.getTotalConsSize(), s.getSignal()]
                if ctrlName != "NONE":
                    line.append(s.getCtrlSignal())
                line.append(s.getTotalConsSignal()) # This sums up the list from getConsSignal
                if ctrlName != "NONE":
                    line.append(s.getTotalConsCtrlSignal()) # This sums up the list from getConsCtrlSignal
                out.write("\t".join(map(str, line)))
                out.write('\n')

def writeCons(outfilename, cons):
    '''
    Write the peaks to a file.
    '''
    # Use with statement for file handling
    with open(outfilename, 'w') as out:
        for chrom in cons:
            for c in cons[chrom]:
                line = [chrom, c.getStart(), c.getEnd(), c.getName(), "0", c.getStrand()]
                out.write('\t'.join(map(str, line)))
                out.write('\n')


def main():
    args = parseArgs()

    loci = loadBed(args.input)
    script_dir = os.path.dirname(os.path.realpath(__file__)) # Get script directory

    if args.exclude > 0:
        tss_file_path = os.path.join(script_dir, "annotation", f"{args.ref}_refseq.ucsc")
        if not os.path.exists(tss_file_path):
            print(f"Error: TSS annotation file not found: {tss_file_path}")
            sys.exit(1)
        tss2000 = loadTSS(tss_file_path, header=True, extension=args.exclude)
    else:
        tss2000 = None

    # tss50 is loaded but not directly used by processPart if tss2000 is passed to it.
    # If it's meant for a different purpose or for the separatePromoter logic, ensure its path is correct.
    tss50_file_path = os.path.join(script_dir, "annotation", f"{args.ref}_refseq.ucsc")
    if not os.path.exists(tss50_file_path):
        print(f"Error: TSS annotation file for tss50 not found: {tss50_file_path}")
        # Decide if this is critical or can be skipped if separatePromoter is False
        # sys.exit(1) # Potentially exit if this file is always needed
    tss50 = loadTSS(tss50_file_path, header=True, extension=50)


    def processPart(elements, window, tss=None, outputStr="NEW", getCons=False):
        e_stitched = stitchLoci(elements, stitchWindow=window, tss=tss)
        calculateSignal(e_stitched, args.bam, 200, isCtrl=False)
        ctrlName = "NONE"
        if args.control:
            print("has control") # Changed print
            calculateSignal(e_stitched, args.control, 200, isCtrl=True)
            ctrlName = os.path.basename(args.control)
        writeCons(outputStr + "_cons.bed", elements)
        writeStitched(outputStr + ".txt", e_stitched, os.path.basename(args.bam), ctrlName=ctrlName)

        outFolder = "./"
        if "/" in outputStr:
            outFolder = os.path.dirname(outputStr)
            if outFolder and outFolder[-1] != '/': # Ensure outFolder is not empty before checking last char
                outFolder += "/"
            if not outFolder: # If outputStr was like "/file", dirname is "/"
                outFolder = "./" # Or handle as appropriate
        if not os.path.exists(outFolder) and outFolder != "./":
            os.makedirs(outFolder, exist_ok=True) # Create output folder if it doesn't exist


        # R script path needs to be robust
        rose_call_super_r_path = os.path.join(script_dir, "ROSE_callSuper.R")
        if not os.path.exists(rose_call_super_r_path):
            print(f"Error: ROSE_callSuper.R not found at {rose_call_super_r_path}")
            # Optionally, skip this step or exit
        else:
            # Use subprocess.run for better control and error handling
            cmd_str = f'R --no-save --args {outFolder} {outputStr}.txt {os.path.basename(outputStr)} {ctrlName} < {rose_call_super_r_path}'
            print(f"Running R command: {cmd_str}")
            try:
                subprocess.run(cmd_str, shell=True, check=True) 
            except subprocess.CalledProcessError as e:
                print(f"Error running R script: {e}")
            except FileNotFoundError:
                print("Error: R command not found. Is R installed and in PATH?")


        if getCons: 
            print("Performing getCons operations (intersectBed)...")
            super_enhancer_file = os.path.join(outFolder, f"{os.path.basename(outputStr)}_Gateway_SuperEnhancers.bed") 
            cons_bed_file = outputStr + "_cons.bed"

            if not os.path.exists(super_enhancer_file):
                print(f"Warning: Super enhancer file not found for getCons: {super_enhancer_file}")
            else:
                cmd1_str = f'intersectBed -a {cons_bed_file} -b {super_enhancer_file} -u > {outFolder}{os.path.basename(outputStr)}_Super_cons.bed'
                cmd2_str = f'intersectBed -a {cons_bed_file} -b {super_enhancer_file} -v > {outFolder}{os.path.basename(outputStr)}_Typical_cons.bed'
                print(f"Running: {cmd1_str}")
                try:
                    subprocess.run(cmd1_str, shell=True, check=True)
                except (subprocess.CalledProcessError, FileNotFoundError) as e:
                     print(f"Error running intersectBed (cmd1): {e}")
                print(f"Running: {cmd2_str}")
                try:
                    subprocess.run(cmd2_str, shell=True, check=True)
                except (subprocess.CalledProcessError, FileNotFoundError) as e:
                    print(f"Error running intersectBed (cmd2): {e}")

        return outFolder, os.path.basename(outputStr) # Return basename for outputStr as used in R command


    def intersectCons(cons_file, refFile, outFile):
        cmd_str = f'intersectBed -a {cons_file} -b {refFile} -u > {outFile}' 
        print(f"Running intersectBed: {cmd_str}")
        try:
            subprocess.run(cmd_str, shell=True, check=True) 
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Error running intersectBed in intersectCons: {e}")


    # Pass tss2000 directly, which could be None if args.exclude <= 0
    outFolder, baseOutputStr = processPart(loci, args.window, tss=tss2000, outputStr=args.output)
    separatePromoter = args.separate

    if separatePromoter:
        print("Separatzing promoters...") 
        thisRefName = args.output + "_tss.bed" 
        writeCons(thisRefName, tss50) 

        super_enhancer_gateway_file = os.path.join(outFolder, f"{baseOutputStr}_Gateway_SuperEnhancers.bed")
        enhancer_gateway_file = os.path.join(outFolder, f"{baseOutputStr}_Gateway_Enhancers.bed")
        cons_file = args.output + "_cons.bed" # This is outputStr_cons.bed from processPart

        # Check if gateway files exist
        if not os.path.exists(super_enhancer_gateway_file):
            print(f"Warning: File not found, cannot separate promoters (Super Enhancers): {super_enhancer_gateway_file}")
            return 
        if not os.path.exists(enhancer_gateway_file):
            print(f"Warning: File not found, cannot separate promoters (Enhancers): {enhancer_gateway_file}")
            return 

        sp_bed = f"{args.output}_SP.bed"
        se_bed = f"{args.output}_SE.bed"
        tp_bed = f"{args.output}_TP.bed"
        te_bed = f"{args.output}_TE.bed"

        cmd_sp = f'intersectBed -a {super_enhancer_gateway_file} -b {thisRefName} -u > {sp_bed}'
        cmd_se = f'intersectBed -a {super_enhancer_gateway_file} -b {thisRefName} -v > {se_bed}'
        cmd_tp = f'intersectBed -a {enhancer_gateway_file} -b {thisRefName} -u | intersectBed -a stdin -b {sp_bed} -v > {tp_bed}'
        cmd_te = f'intersectBed -a {enhancer_gateway_file} -b {thisRefName} -v | intersectBed -a stdin -b {se_bed} -v > {te_bed}'

        for cmd_str in [cmd_sp, cmd_se, cmd_tp, cmd_te]:
            print(f"Running: {cmd_str}")
            try:
                subprocess.run(cmd_str, shell=True, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"Error during promoter separation: {e}")


        refFiles = [f"{args.output}_{i}.bed" for i in ["SP", "TP", "SE", "TE"]]
        outFiles = [f"{args.output}_{i}_cons.bed" for i in ["SP", "TP", "SE", "TE"]]

        for refFile, outFile in zip(refFiles, outFiles):
            if not os.path.exists(refFile):
                print(f"Warning: refFile for intersectCons not found: {refFile}. Skipping.")
                continue
            intersectCons(cons_file, refFile, outFile)

if __name__ == '__main__':
    main()