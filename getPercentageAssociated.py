from GenomeTools import *

'''
Arguments

[0] .narrowPeak file
[1] .narrowPeak file
[2] genome file (.gtf)

This fucntion executes:
- For each peak in [0], find the closest gene in [2] and record it in the peak object if they overlap
- For each of these peaks in [0] that HAS an associated gene, record if there exists a peak in [1] within
    2000bps of (this peak in [0] + associated gene boundary). returns this value as a percentage.
       
For example:
[0] is RNA-seq peak data
[1] is ATAC-seq peak data
[2] is the genome assembly used (eg hg19.gtf)

this function would return the perentage of gene-associated RNA-seq peaks that also have an ATAC-seq peak within 2000bp
    of itself or its associated gene region

'''
def main(argv):
    start_time = time.time()
    print()
    print('loading files')
    print()
    
    # get file extensions
    indexes = []
    for i in [0,1,2]:
        for char_index in range(len(str(argv[i]))):
            if str(argv[i])[char_index] == ".":
                indexes.append(char_index)
                break
    
    RNA_peaks = readGeneFile(str(argv[0]),str(argv[0]), str(argv[0])[indexes[0]+1:])
    ATAC_peaks =  readGeneFile(str(argv[1]),str(argv[1]), str(argv[1])[indexes[1]+1:])
    gene_list = readGeneFile(str(argv[2]),str(argv[2])[0:4], str(argv[2])[indexes[2]+1:])

    print("--- %s seconds ---" % (time.time() - start_time))
    
    print()
    print('associating genes')
    print()
    
    print("--- %s seconds ---" % (time.time() - start_time))

    associate_genes(gene_list, RNA_peaks) # only consider peaks that lie within annotated gene regions
        
    print()
    print('calculating percentage of peaks associated')
    print()

    percent_assoc = get_percentage_associated(RNA_peaks,ATAC_peaks, threshold=2000)
    print ("Percentage of " + RNA_peaks.name + " peaks with " + ATAC_peaks.name + " peaks within 2000 bp = ", percent_assoc[0],"%") 

    print()
    print("--- %s seconds ---" % (time.time() - start_time))    
    print()

    return percent_assoc

if __name__ == "__main__":
    main(sys.argv[1:])
