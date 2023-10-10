from collections import defaultdict

def TADlinks(gene_file, enhancer_tad_file, output_file):
    '''
        Trim down to a file showing just the gene>enhancer>cell>assay
    '''

    tadgenes = defaultdict(dict)
    with open(gene_file, 'r') as tss_file:
        for line in tss_file:
            chr, start, end, geneID, one, two, three, tad, tissue = line.strip().split('\t')
            cell = tissue.split('|')[1]
            tad = tad + '\t' + cell
            tadgenes[tad][geneID] = 1

    tadenh = defaultdict(dict)
    with open(enhancer_tad_file, 'r') as enh_file:
        for line in enh_file:
            chr, start, end, enhID, one, two, three, tad, tissue = line.strip().split('\t')
            cell = tissue.split('|')[1]
            tad = tad + '\t' + cell
            tadenh[enhID][tad] = 1

    n = len(tadenh)
    batch_size = 1000
    start_idx = 0
    count = 1
    while start_idx < n:
        end_idx = min(start_idx + batch_size, n)
        count += 1
        TADlinks_write_to_file(tadenh, tadgenes, output_file, start_idx, end_idx)
        start_idx = end_idx

def TADlinks_write_to_file(tadenh, tadgenes, output_file, start_idx, end_idx):
    with open(output_file, 'a') as out:
        for enh in tadenh.keys():
            for tad in list(tadenh[enh].keys())[start_idx:end_idx]:
                if tad in tadgenes:
                    for gene in tadgenes[tad].keys():
                        panthID = panther_mapping.get(gene, '')
                        cell = tad.split('\t')[1]
                        if enh and panthID and cell:
                            out.write(f"{enh}\t{panthID}\t{cell}\t3\n")

def tissuesReplace(input_file, output_file):
    '''
        Replace the tissues and cell types with their codes
    '''
    with open(output_file, 'w') as out:
        with open(input_file, 'r') as tss_file:
            for line in tss_file:
                line_split = line.strip().split('\t')
                if len(line_split) == 4:
                    enhID, panthID, tiss, assay = line_split
                    cell = tissues.get(tiss, '')
                    out.write(f"{enhID}\t{panthID}\t{cell}\t{assay}\n")

panther_mapping = {}
with open('pantherGeneList.txt', 'r') as pantherGene_file:
    for line in pantherGene_file:
        line_split = line.strip().split('\t')
        panth = line_split[0]
        ensg = line_split[1]
        panther_mapping[ensg] = panth

tissues = {}
with open('tissuetable_10092018.txt', 'r') as tissue_file:
    for line in tissue_file:
        line_split = line.strip().split('\t')
        tissueID = line_split[0]
        tissue = line_split[1]
        tissue = tissue.replace(' ', '_')
        tissues[tissue] = tissueID

TADlinks('TSStTAD', 'enhancerstTAD', 'linkstTAD')
tissuesReplace('linkstTAD', 'linkstTADtissues')