from GenomeTools import *

def test_distance_between_peaks():
    entry1 = GenomicEntry('chr1',1001,1050, 'entry1')
    entry2 = GenomicEntry('chr1',1050, 1087, 'entry2')

    if distance_between_peaks(entry1, entry2) != 1:
        return False
    
    entry3 = GenomicEntry('chr1',1049, 1087, 'entry3')
    
    if distance_between_peaks(entry1, entry3) != 0:
        return False    
    
    return True
    
def test_get_entries():
    gtf1 = readGeneFile('gtfTest.gtf','gtfTest','gtf')
    gtf2 = readGeneFile('gtfTest2.gtf','gtfTest2','gtf')

    if len(gtf1) != 7:
        return False     

    if len(gtf2) != 1:
        return False
    
    return True
    

def test_get_closest_peak(geneEntryContainer, candidatepeak, requiredpeak):
    '''
    returns True is if candidatepeak is closest peak to requiredpeak
    '''
    first = 0
    last = len(geneEntryContainer.entries[requiredpeak.chrom])-1
    found = False
    while not found:
        midpoint = (first + last)//2
        closestPeak = geneEntryContainer.entries[requiredpeak.chrom][midpoint]
        if candidatepeak == closestPeak:  
            found = True
        elif closestPeak.chromStart < candidatepeak.chromStart:
            first = midpoint + 1
        else:
            last = midpoint - 1
       
    closest = []
    values = (midpoint+0,midpoint+1,midpoint-1,midpoint+2,midpoint-2) # indicies either side to search
    for value in values:
        try:
            closest.append((geneEntryContainer.entries[requiredpeak.chrom][value],abs(distance_between_peaks(requiredpeak,geneEntryContainer.entries[requiredpeak.chrom][value])),value))
        except IndexError:    # fringe cases (no index exists to one side of the found value)    
            continue
        
    closest.sort(key = lambda four_closest : four_closest[1])
    
    # if any peak to either side is actually closer to the peak we're searching for (ie get_cloest_peak returned not the closest peak)
    if closest[0][1] != abs(distance_between_peaks(candidatepeak,requiredpeak)): 
        print([closest[x][1] for x in range(len(closest))])
        print([str(closest[x][0]) for x in range(len(closest))])

        print('candidate peak', str(candidatepeak))
        print('candidate peak distance', abs(distance_between_peaks(candidatepeak,requiredpeak)))
        print('actual peak', closest[0][0])
        print('actual peak distance', closest[0][1])
        print('actual peak index', closest[0][2])
        return False
    else:
        return True
    
def test_get_Percentage_Associated(output, threshold):
    
    # duplication check
    for i in range(len(output[0])):
        try:
            if output[i][0] == output[i+1][0]:
                return False      
        except IndexError:
            continue
    
    # distance checks
    for entry in output:
        if abs(distance_between_peaks(entry[0],entry[1])) >= threshold and abs(distance_between_peaks(entry[0],entry[2])) >= threshold:
            return False
    
    return True

    
