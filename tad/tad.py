from math import floor, ceil
from pybedtools import BedTool
from collections import defaultdict

def bTADoverlap(input_file, output_file, cell_type):
    hash = {}
    with open(input_file, 'r') as tad_file:
        for line in tad_file:
            chr, start1, end1, location, pval = line.strip().split('\t')
            score = pval + '|' + cell_type
            location_split = location.split('___')
            one = location_split[0]
            two = location_split[1]
            three = one.split('|')[2]
            four = two.split('|')[2]
            five = three.split(':')[1]
            six = four.split(':')[1]
            s1, s2 = five.split('-')
            e1, e2 = six.split('-')
            start = floor((int(s1) + int(s2))/2)
            end = ceil((int(e1) + int(e2))/2)
            total = f"{chr}\t{start}\t{end}\t{location}\t{score}"
            hash[total] = 1

    with open(output_file, 'w') as out:
        for key in sorted(hash.keys()):
            out.write(key + '\n')


def tTADOverlap(input_file, output_file, cell_type):
    with open(input_file, 'r') as tad_file, open(output_file, 'w') as out:
        for line in tad_file:
            chr, start, end, location, pval = line.strip().split('\t')
            score = pval + '|' + cell_type
            out.write(f"{chr}\t{start}\t{end}\t{location}\t{score}\n")


def bedtoolIntersect(gene_file, tad_order_file, output_file):
    '''
        Use bedtools intersect to find which TADs contain at least 90% of a gene's promoter or at least 90% of an enhancer
    '''
    gene_file_bedtoolFile = BedTool(gene_file)
    tad_order_file_bedtoolFile = BedTool(tad_order_file)

    intersection = gene_file_bedtoolFile.intersect(tad_order_file_bedtoolFile, wa=True, wb=True, f=0.9)
    with open(output_file, 'w') as out:
        for elt in intersection:
            out.write(str(elt))


def split(file):
    '''
        splits input file by chromosome and creates 25 chromosome specific files
    '''
    chromosomes = ['chr' + str(x) for x in range(1, 23)]
    chromosomes.append('chrX')
    chromosomes.append('chrY')

    output_files = {chromosome: open(chromosome + file, 'w') for chromosome in chromosomes}

    with open(file, 'r') as tad_file:
        for line in tad_file:
            chr = line.strip().split('\t')[0]

            if chr in chromosomes:
                output_files[chr].write(line)

    for output_file in output_files.values():
        output_file.close()


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
    with open(output_file, 'w') as out, open(input_file, 'r') as tss_file:
        for line in tss_file:
            line_split = line.strip().split('\t')
            if len(line_split) == 4:
                enhID, panthID, tiss, assay = line_split
                cell = tissues.get(tiss, '')
                out.write(f"{enhID}\t{panthID}\t{cell}\t{assay}\n")


def concatenate(input_files, output_file):
    '''
        makes a combined file for all tissues
    '''
    with open(output_file, 'wb') as outfile:
        for file_name in input_files:
            with open(file_name, 'rb') as infile:
                for line in infile:
                    outfile.write(line)
    

def cutdowntad(chia_file, eqtl_file, heirarchical_tad_file, tad_file, output_file):
    '''
        only want to include the TAD links that support existing links from ChIA-PET or eQTL data
    '''
    chia = open(chia_file, 'r')
    chia_lines = chia.readlines()
    chia.close()

    hash = {}
    count1 = 0
    count2 = 0

    def line_parse(line, hash):
        line_split = line.strip().split('\t')
        enhancer = line_split[0]
        panthid = line_split[1]
        link = enhacer + '_' + panthid
        hash[link] = 1

    for line in chia_lines:
        line_parse(line, hash)

    with open(eqtl_file, 'r') as eqtl:
        for line in eqtl:
            line_parse(line, hash)

    with open(heirarchical_tad_file, 'r') as heirarchical:
        for line in heirarchical:
            line_parse(line, hash)

    with open(output_file, 'w') as out, open(tad_file, 'r') as tad:
        for line in tad:
            enhancer, panthid, tissue, assay = line.strip().split('\t')
            link = enhacer + '_' + panthid
            rest = f"{enhancer}\t{panthid}\t{tissue}\t{assay}"
            if link in hash:
                count1 += 1
                out.write(f"{rest}\n")
            else:
                count2 += 1

    print(f"Number of TAD links found in TADintxn, eQTL, or ChIA: {count1}")


def orderlinks(input_file, output_file):
    hash = {}
    tissues = {}
    genes = {}
    count = 0

    with open(input_file, 'r') as input:
        for line in input:
            enhID, gene, tissue, assay = line.strip().split('\t')
            hash.setdefault(enhID, {}).setdefault(gene, {})[tissue] = 1

    with open(output_file, 'w') as out:
        for enhancer in sorted(hash.keys()):
            for geen in sorted(hash[enhancer].keys()):
                genes[geen] = 1
                for cell in sorted(hash[enhancer][geen].keys()):
                    tissues[cell] = 1
                    count += 1
                    out.write(f"{enhancer}\t{geen}\t{cell}\t{assay}\n")

    print(f"This file contains {len(hash)} different enhancers linked to {len(genes)} genes in {len(tissues)} tissues.")

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


bTADoverlap('ENCFF274VJU.bed', 'bTADorder', 'Caki2')
tTADOverlap('ENCFF588KUZ.bed', 'tTADorder', 'Caki2')

bedtoolIntersect('TSSgenesbed', 'bTADorder', 'TSSbTAD')
bedtoolIntersect('TSSgenesbed', 'tTADorder', 'TSStTAD')

bedtoolIntersect('CREbedDBenhancers_10092018', 'bTADorder', 'enhancersbTAD')
bedtoolIntersect('CREbedDBenhancers_10092018', 'tTADorder', 'enhancerstTAD')

split('TSSbTAD')
split('enhancersbTAD')

chromosomes = ['chr' + str(x) for x in range(1, 23)]
chromosomes.append('chrX')
chromosomes.append('chrY')
linkbTAD_list = []

for chr in chromosomes:
    TADlinks(chr+'TSSbTAD', chr+'enhancersbTAD', chr+'linksbTAD')

for chr in chromosomes:
    linkbTAD_list.append(chr+'linksbTAD')
concatenate(linkbTAD_list, 'linksbTAD')
    
TADlinks('TSStTAD', 'enhancerstTAD', 'linkstTAD')

tissuesReplace('linksbTAD', 'linksbTADtissues')
tissuesReplace('linkstTAD', 'linkstTADtissues')

concatenate(['linksbTADtissues', 'linkstTADtissues'], 'linksTADtissues')

# cutdowntad('linksDBchia', 'linksDBnumeqtl', 'PSYCHIClinksDB ', 'linkstTADtissues', 'selecttTAD')
# orderlinks('selectTAD', 'linksDBtad')