from pybedtools import BedTool

def bedtoolIntersect(gene_file, tad_order_file, output_file):
    '''
        Use bedtools intersect to find which TADs contain at least 90% of a gene's promoter or at least 90% of an enhancer
    '''
    gene_file_bedtoolFile = BedTool(gene_file)
    tad_order_file_bedtoolFile = BedTool(tad_order_file)

    intersection = gene_file_bedtoolFile.intersect(tad_order_file_bedtoolFile, wa=True, wb=True, f=0.9)
    out = open(output_file, 'w')
    for elt in intersection:
        out.write(str(elt))
    out.close()

bedtoolIntersect('TSSgenesbed', 'bTADorder', 'TSSbTAD')
bedtoolIntersect('TSSgenesbed', 'tTADorder', 'TSStTAD')

bedtoolIntersect('CREbedDBenhancers_10092018', 'bTADorder', 'enhancersbTAD')
bedtoolIntersect('CREbedDBenhancers_10092018', 'tTADorder', 'enhancerstTAD')