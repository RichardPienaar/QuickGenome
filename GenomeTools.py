from SimpleGenome import *
import tests
import time, sys

def get_closest_peak(genomicEntryContainer,entry):
    '''
    returns closest peak in genomicEntryContainer to entry
    '''
    first = 0
    last = len(genomicEntryContainer.entries[entry.chrom])-1
    while last>=first:
        midpoint = (first + last)//2
        closestPeak = genomicEntryContainer.entries[entry.chrom][midpoint]
        
        if distance_between_peaks(entry, closestPeak) == 0:
            return closestPeak 
        elif distance_between_peaks(entry,closestPeak) < 0:
            first = midpoint + 1
        else:
            last = midpoint - 1
        
    if last <0:
        return genomicEntryContainer.entries[entry.chrom][0]

    if first >= len(genomicEntryContainer.entries[entry.chrom])-1:
        return genomicEntryContainer.entries[entry.chrom][-1]
    
    if (abs(distance_between_peaks(entry,genomicEntryContainer.entries[entry.chrom][last])) < abs(distance_between_peaks(entry, genomicEntryContainer.entries[entry.chrom][first]))):
        return  genomicEntryContainer.entries[entry.chrom][last]
    return genomicEntryContainer.entries[entry.chrom][first]

def associate_genes(geneContainer, peaks):
    '''
    Associates each entry in peaks (GenomicEntryContainer) with its closest gene found in geneContainer (GenomicEntryContainer, a gene file .gtf)
    '''
    for peak in peaks: 
        closest_peak = get_closest_peak(geneContainer, peak)
        if abs(distance_between_peaks(peak,closest_peak)) == 0:
            peak.associate_gene(geneContainer.name,closest_peak) # associating genes
        
   
def get_percentage_associated(peaks1, peaks2, threshold):
    '''
    - return[0] the percentage (float) of peak1 peaks (GenomicEntryContainer) that have a peak from peak2 (GenomicEntryContainer) within threshold distance (int)
    - return[1] [peak1 peak, peak1.associatedgene, peak2 peak] for every peak1 peak within threshold
    
    - distance is calculated from entire associated gene region + peak 1 peak region
    '''
    count = 0
    count1 = 0
    Associated_Peaks = []
    for peak in peaks1:
        to_add = []
        if peak.associated_gene:
            count1 +=1
            to_add.append(peak)
            to_add.append(peak.closest_gene)
            if abs(distance_between_peaks(peak, get_closest_peak(peaks2, peak))) < threshold: # peak2 peak within threshold distance of peak1 peak
                to_add.append(get_closest_peak(peaks2, peak))
                count+=1
            elif abs(distance_between_peaks(peak.closest_gene, get_closest_peak(peaks2, peak.closest_gene))) < threshold: # peak2 peak within threshold distance of gene overlapping with peak1 peak
                to_add.append(get_closest_peak(peaks2, peak.closest_gene))
                count+=1

        if len(to_add) > 2:
            Associated_Peaks.append([peak for peak in to_add])

    return count/count1*100, Associated_Peaks
    

    