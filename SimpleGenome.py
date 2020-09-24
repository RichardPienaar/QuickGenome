import time,sys

def distance_between_peaks(entry1,entry2):
    '''
    Negative return (float) if entry2 (GeneEntry) to the left of entry1(GeneEntry)
    if any bases of either peak ar eshared, returns 0 
    NOTE: bedFile chromStart indexes at 0 
    '''
    if entry1.chrom == entry2.chrom:
        if (entry1.chromStart >= entry2.chromEnd): return entry2.chromEnd - (entry1.chromStart+1)
        if (entry1.chromEnd <= entry2.chromStart): return (entry2.chromStart+1) - entry1.chromEnd
        
    else:
        raise NameError("different chromosomes")
    return 0 # intersect

def merge_overlapping_peaks(peak1,peak2):
    '''
    merges 2 overlapping peaks (GeneEntry) and combines names
    - adds name of peak2 on to the name of peak1
    '''
    if peak1.chromStart > peak2.chromStart:
        peak1.chromStart = peak2.chromStart
    if peak1.chromEnd < peak2.chromEnd:
        peak1.chromEnd = peak2.chromEnd
    peak1.add_alternate_name(peak2.gene_ids)
    return peak1

class GenomicEntry ():
    '''
    Simple Genomic entry - e.g gene from hg29 assembly
    '''
    def __init__(self, chrom, chromStart, chromEnd, name):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.gene_ids = name

    def add_alternate_name(self, name):
        self.gene_ids += "/ " + name
    
    def __str__(self):
        return self.chrom + " " + str(self.chromStart) + " - " + str(self.chromEnd) + " " + self.gene_ids

    def __eq__(self, other):
        if distance_between_peaks(self,other) == 0: # same chromosome already assumed
            return True
        else:
            return False

class ComplexGenomicEntry (GenomicEntry):
    '''
    Complex Genomic entry - RNA-seq peak, ATAC-seq peak, histone peak etc 
    '''
    def __init__(self, chrom, chromStart, chromEnd,name,dataType, score=None):
        super().__init__(chrom, chromStart, chromEnd,name)
        self.score = score
        self.associated_gene = False # Genomic Entry corresponds to a known gene
        self.type = dataType # type of assay 
    
    def associate_gene(self, genome,gene):
        self.associated_gene = True
        self.associated_genome = genome # Genome assembly name 
        self.closest_gene = gene # associated gene from reference
    
    def __str__(self):
        return self.chrom + " " + str(self.chromStart) + " - " + str(self.chromEnd) + " " + self.gene_ids + " " + str(self.score)

class GenomicEntryContainer ():
    '''
    Container of Genomic Entries

    Represents either GTF file OR narrowPeak file
     - only difference is narrow peak entries will have a score != None
    '''
    def __init__(self, name):
        self.chroms = []
        self.entries = {}
        self.name = name # container name eg "hg19"
    
    def add_chroms(self, chrom):
        if chrom not in self.chroms:
            self.chroms.append(chrom)
    
    def create_containers(self):
        for chromName in self.chroms:
            self.entries[chromName] = []
    
    def add_entry(self, entry):
        self.entries[entry.chrom].append(entry)
    
    def remove_entry(self,entry):
        self.entries[entry.chrom].remove(entry)

    def __len__(self):
        return sum([len(self.entries[x]) for x in self.entries]) 
    
    def __iter__(self):
        for value in self.entries.values():
            for entry in value:
                yield entry
    
def get_chroms(gene_object,filename):
    '''
    parses quicky through a file and records all the chromosomes present to a gene_object (GenomicEntryContainer)
    flename(str) - name of file

    !! Header line not allowed !!
    '''
    f = open(filename)
    row = 0
    
    for line in f:
        row += 1
        words = line.strip().split()
        try:
            chrom = str(words[0])
            gene_object.add_chroms(chrom)
        except RuntimeError:
            raise RuntimeError('Error in BED file at row %d (%s)')

    f.close()
    return gene_object
    

def get_entries(gene_object,filename,format,mergeOverlappingEntries=True):
    """
    parses through a file and records all entries present to a gene_object (GenomicEntryContainer)
    - merges subsequent overlaping entries 

    flename(str) - name of file
    format(str) - file format (.narrowPeaks or .gtf)
    
    !! Header line not allowed !!   
    """
    f = open(filename)
    row = 0
    
    format=format.lower()

    if format == 'gtf':
        if mergeOverlappingEntries:
            previous_entry = None
            for line in f:

                row += 1
                words = line.strip().split()          
                if words[2] == "transcript": # Only full regions considered
                    entry = GenomicEntry(words[0],int(words[3]),int(words[4]),words[9])

                    try:
                        if entry.chrom == previous_entry.chrom:

                            if distance_between_peaks(previous_entry, entry) != 0:
                                gene_object.add_entry(previous_entry)
                                previous_entry = entry
                            else:

                                previous_entry = merge_overlapping_peaks(previous_entry, entry)
                                entry = None

                        else:
                            gene_object.add_entry(previous_entry)
                            previous_entry = entry


                    except AttributeError: # first entry
                        previous_entry = entry
            
            # catch last index
            gene_object.add_entry(previous_entry)

        
        else: # no merging
            for line in f:
                row += 1
                words = line.strip().split()          
                if words[2] == "transcript": # Only full regions considered
                    entry = GenomicEntry(words[0],int(words[3]),int(words[4]),words[9])
                gene_object.add_entry(entry)

    elif format == 'bedraph':
      for line in f:
            row += 1
            words = line.strip().split()
            entry = ComplexGenomicEntry(words[0],int(words[1]),int(words[2]),str(row),float(words[3]))
            gene_object.add_entry(entry)


    elif format == 'narrowpeak':
        if 'rna' in filename.lower():
            dataType = "RNA-seq"
        elif 'atac' in filename.lower():
            dataType = 'ATAC-seq'
        else:
            dataType = input('Please tell me what assay type this is (i.e RNA-seq) - ', filename)
        for line in f:
            row += 1
            words = line.strip().split()
            entry = ComplexGenomicEntry(words[0],int(words[1]),int(words[2]),words[3],dataType,float(words[4]))
            gene_object.add_entry(entry)

    else:
        print('filetype not supported', format)
        raise AttributeError #filetype not supported

    f.close()
    return gene_object    

def readGeneFile(filename, name, format, mergeOverlappingEntries=True):
    '''
    Parses a file of genomic entries and returns a GenomicEntryContainer containing all the entries in the file
    - get_chroms called first to check file completeness quickly
    '''
    new_GTF_object = GenomicEntryContainer(name)
    get_chroms(new_GTF_object,filename)
    new_GTF_object.create_containers()
    get_entries(new_GTF_object, filename, format, mergeOverlappingEntries)
    return new_GTF_object


def writeGenefile(gene_container):
    pass


def main():
    pass

    
if __name__=="__main__":
    main()

