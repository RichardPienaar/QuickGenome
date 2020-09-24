import bed # bodenlab - binfpy
import GenomeTools,SimpleGenome, getPercentageAssociated, tests
import time, sys

'''
unit testing for GenomeTools - getPecentageAssociated() by comparison to binfpy by boden labs 

TODO: bed.BedFile.getOverlap does not count 1bp shared as overlap 

I think this is affecting the getPecentageAssociated() output in each case
'''


def main(argv):
    start_time = time.time()

    RNA_peaks1 = bed.BedFile(str(argv[0]),format='narrowPeak')
    ATAC_peaks1 = bed.BedFile(str(argv[1]),format='narrowPeak')
    RNA_peaks_no1 = len(RNA_peaks1)

    count=0
    overlappingRNA_ATAC=[]
    for peak in RNA_peaks1:
        if bed.dist(next(iter(ATAC_peaks1.getClosest(peak))),peak) == 0:
            overlappingRNA_ATAC.append((peak,ATAC_peaks1.getOverlap(peak)))
            count+=1

    overlap1a = count/RNA_peaks_no1*100

    print()
    print("--- %s seconds for peak file loading - BINFPY ---" % (time.time() - start_time))
    print()
    start_time = time.time()

    GTF1 = bed.BedFile(str(argv[2]),format='gtf')
    print()
    print("--- %s seconds for geneome file loading - BINFPY  ---" % (time.time() - start_time))
    print()
    start_time = time.time()

    count = 0
    for peak in RNA_peaks1:
        if GenomeTools.distance_between_peaks(next(iter(GTF1.getClosest(peak))),peak) == 0:
            count+=1

    overlap1b = count/RNA_peaks_no1*100
    print()
    print("--- %s seconds for genome file overlap - BINFPY ---" % (time.time() - start_time))
    print()
    start_time = time.time()

    RNA_peaks2 = SimpleGenome.readGeneFile(str(argv[0]),format='narrowPeak',name='RNA_1781')
    ATAC_peaks2 = SimpleGenome.readGeneFile(str(argv[1]),format='narrowPeak',name='ATAC_1781')

    RNA_peaks2_no = len(RNA_peaks2)
    
    overlappingRNA_ATAC2=[]

    count = 0
    for peak in RNA_peaks2:
        if GenomeTools.distance_between_peaks(GenomeTools.get_closest_peak(ATAC_peaks2, peak),peak) == 0:
            overlappingRNA_ATAC2.append((peak, GenomeTools.get_closest_peak(ATAC_peaks2, peak)))
            count+=1

    overlap2a = count/RNA_peaks2_no*100

    print()
    print("--- %s seconds for peak file loading - SimpleGenome ---" % (time.time() - start_time))
    print()
    start_time = time.time()


    GTF2 = SimpleGenome.readGeneFile(str(argv[2]),format='gtf',name='hg19')
    print()
    print("--- %s seconds for geneome file loading - SimpleGenome  ---" % (time.time() - start_time))
    print()
    start_time = time.time()

    count = 0
    for peak in RNA_peaks2:
        if GenomeTools.distance_between_peaks(GenomeTools.get_closest_peak(GTF2, peak),peak) == 0:
            count+=1

    overlap2b = count/RNA_peaks2_no*100

    print()
    print("--- %s seconds for genome file overlap - SimpleGenome ---" % (time.time() - start_time))
    print()
    start_time = time.time()

    count = 0
    count1 = 0
    found_peaks = []
    for RNApeak in RNA_peaks1:
        to_add = []
        found = False
        closest = next(iter(GTF1.getClosest(RNApeak)))
        overlap = GTF1.getOverlap(RNApeak) # only associated genes - empty if kissing bases

        if bed.dist(closest,RNApeak) == 0:  # correcting for kissing bases not counting 
            count1 +=1
        
            for ATAC_peak in ATAC_peaks1.getClosest(RNApeak):
                if bed.dist(ATAC_peak , RNApeak) <2000: 
                    found = True

                for closeGene in GTF1.getClosest(RNApeak):     
                    if bed.dist(ATAC_peak, closeGene) <2000: 
                        found = True
                

            if overlap != []: 
                for overlapping in overlap:
                    for ATAC_peak_overlap in ATAC_peaks1.getClosest(overlapping):
                        if bed.dist(ATAC_peak_overlap , overlapping) <2000: 
                            found = True
                        if bed.dist(ATAC_peak_overlap , RNApeak) <2000: 
                            found = True           
            if found:
                count+=1

                
    overlap3a = count/count1*100
    
    print()
    print("--- %s seconds for getPercentageAssociated() - BINFPY ---" % (time.time() - start_time))
    print()
    start_time = time.time()

    overlap3b = getPercentageAssociated.main(argv)

    print()
    print("--- %s seconds for getPercentageAssociated() - SimpleGenome ---" % (time.time() - start_time))
    print()

    if not tests.test_get_Percentage_Associated(overlap3b[1],2000):
        raise NameError("Test Failed")

    if overlap1a != overlap2a or overlap1b != overlap2b:
        raise NameError("Different outputs")
    else:
        print(overlap1a,overlap1b,overlap2a,overlap2b,overlap3a,overlap3b[0])
        print('All Tests Passed - Yay!')


if __name__=="__main__":
    main(sys.argv[1:])