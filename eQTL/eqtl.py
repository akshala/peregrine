import pandas as pd 
import numpy as np
from pybedtools import BedTool
import glob
import os
import gzip

def split(file):
    '''
        splits input file by chromosome and creates 25 chromosome specific files
    '''
    chromosomes = ['chr' + str(x) for x in range(1, 23)]
    chromosomes.append('chrX')
    chromosomes.append('chrY')

    enh_file = open(file, 'r')
    enh_lines = enh_file.readlines()
    enh_file.close()

    file = os.path.basename(file)
    output_files = {chromosome: open(os.path.join(file.split('.')[0], chromosome + file), 'w') for chromosome in chromosomes}

    for line in enh_lines[1:]:
        variant_id = line.strip().split('\t')[0]
        chr = 'chr' + variant_id.split('_')[0]

        if chr in chromosomes:
            output_files[chr].write(line)

    for output_file in output_files.values():
        output_file.close()

def trim(file):
    '''
        removes undesirable columns from input file 
    '''
    df = pd.read_csv(file, delimiter='\t', header=0)
    new_df = df.filter(['variant_id', 'gene_id', 'pval_nominal'], axis=1)
    new_df['tissue'] = ['Colon_Sigmoid'] * len(new_df)
    new_df.rename(columns={'variant_id': 'variant', 'gene_id': 'gene'}, inplace=True)
    new_df.index = np.arange(1, len(new_df) + 1)
    file = os.path.basename(file)
    new_df.to_csv(os.path.join(file.split('.')[0], 'newdata' + file), sep='\t', index=False)

def eqtl_process(input_file, output_file):
    '''
        rearranges the file into bed format
    '''
    eqtl_file = open(input_file, 'r')
    eqtl_lines = eqtl_file.readlines()
    eqtl_file.close()

    out = []
    for line in eqtl_lines[1:]:
        variant_id, gene_id, pval_nominal, tissue = line.strip().split('\t')
        variant_id_split = variant_id.split('_')
        chr = 'chr' + variant_id[0]
        start = variant_id_split[1]
        key = f"{chr}\t{start}\t{start}\t{variant_id}\t{gene_id}_{pval_nominal}"
        out.append(key + '_' + tissue)

    with open(output_file, 'w') as outfile:
        outfile.write("\n".join(out))

def pantherIDweed(input_file, output_file, leftover_file):
    '''
        removes any genes not found in PANTHER
    '''
    count1 = 0
    count2 = 0
    out = open(output_file, 'w')
    leftover = open(leftover_file, 'w')

    with open(input_file, 'r') as eqtl_file:
        for line in eqtl_file:
            one, two, three, four, five = line.strip().split('\t')
            three = five.split('.')[0]

            if three in panther_mapping:
                count1 += 1
                panther_three = panther_mapping[three]
                out.write(f"{one}\t{two}\t{two}\t{four}\t{five}\n")
            else:
                count2 += 1
                leftover.write(line)
    out.close()
    leftover.close()

def eliminateOverlap(input_file, output_file, overlap_file, unmatched_file):
    '''
        Eliminate any eQTL that are found in exons of their own genes
    '''
    hash = {}
    with open('exons_genes.txt', 'r') as exon_file:
        for line in exon_file:
            line_split = line.strip().split('\t')
            chr = line_split[0]
            start = line_split[1]
            end = line_split[2]
            ensg = line_split[4].split('.')[0]
            line_key = f"{chr}\t{start}\t{end}\t{ensg}"
            if ensg not in hash:
                hash[ensg] = {}
            hash[ensg][line_key] = 1

    count1 = 0
    count2 = 0
    out_dict = {}
    overlap_dict = {}
    unmatched = open(unmatched_file, 'w')

    with open(input_file, 'r') as eqtl_file:
        for line in eqtl_file:
            line = line.strip()
            one, two, three, four, five = line.split('\t')
            key = five.split('.')[0]
            if key in hash:
                for lock in sorted(hash[key].keys()):
                    chr, start, end, ensg = lock.split('\t')
                    two = int(two)
                    start = int(start)
                    end = int(end)

                    if two < start:
                        count1 += 1
                        out_dict[line] = 1
                    elif two > end:
                        count2 += 1
                        out_dict[line] = 1
                    else:
                        overlap_dict[line] = 1
            else:
                unmatched.write(line)

    unmatched.close()
    out = open(output_file, 'w')
    for key in sorted(out_dict.keys()):
        out.write(key + '\n')
    out.close()
    overlap = open(overlap_file, 'w')
    for key in sorted(overlap_dict.keys()):
        overlap.write(key + '\n')
    overlap.close()

def bedtoolIntersect(cred_file, eqtl_file, output_file):
    '''
        see if any of the eQTL intersect with the enhancers
    '''
    cred_file_bedtoolFile = BedTool(cred_file)
    eqtl_file_bedtoolFile = BedTool(eqtl_file)

    intersection = cred_file_bedtoolFile.intersect(eqtl_file_bedtoolFile, wa=True, wb=True)
    out = open(output_file, 'w')
    for elt in intersection:
        out.write(str(elt))
    out.close()

def eqtllinks(input_file, output_file):
    '''
        replace the ENSG IDs with PANTHER long IDs and replace tissue names with tissue IDs
    '''
    match = {}
    with open(input_file, 'r') as eqtl_file:
        for line in eqtl_file:
            chr, start, end, enhID, three, four, five, snp, geee = line.strip().split('\t')
            gene = geee.split('.')[0]
            tissue = geee.split('_')
            tissue.pop(0)
            pval = tissue.pop(0)
            tissue = '_'.join(tissue)
            match.setdefault(enhID, {}).setdefault(gene, {}).setdefault(snp, {}).setdefault(pval, {})[tissue] = 1

    with open(output_file, 'w') as out:
        for enhID in match:
            for gene in match[enhID]:
                for snp in match[enhID][gene]:
                    for pval in match[enhID][gene][snp]:
                        for tissue in match[enhID][gene][snp][pval]:
                            tiss = tissues[tissue]
                            if gene in panther_mapping:
                                out.write(f"{enhID}\t{panther_mapping[gene]}\t{snp}\t{pval}\t{tiss}\t2\n")

def linksDB(input_file, output_file):
    '''
        creates file that summarizes the eQTL information for each link
    '''
    hash = {}
    with open(input_file, 'r') as eqtl_file:
        for line in eqtl_file:
            enhID, panthID, eqtl, pvalue, tissue, assay = line.strip().split('\t')
            hash.setdefault(enhID, {}).setdefault(panthID, {}).setdefault(tissue, {})[eqtl] = 1

    with open(output_file, 'w') as out:
        out.write("enhancer\tgene\ttissue\tnumber_of_eQTL\tassay\n")
        for enh in sorted(hash.keys()):
            for gene in hash[enh]:
                for tis in sorted(hash[enh][gene].keys()):
                    count = len(hash[enh][gene][tis])
                    out.write(f"{enh}\t{gene}\t{tis}\t{count}\t{assay}\n")

def concatenate(input_files, output_file):
    '''
        makes a combined file for all tissues
    '''
    with open(output_file, 'w') as outfile:
        for file_name in input_files:
            infile = open(file_name, 'r')
            outfile.write(infile.read())
            infile.close()

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

gene_pair_gz_files = glob.glob(os.path.join('GTEx_Analysis_v7_eQTL', '*variant_gene_pairs.txt.gz'))
links_files = []
linksnum_files = []

for gene_pair_gz_file in gene_pair_gz_files:

    # get tissue name from gz file and make a folder for it
    folder = os.path.basename(gene_pair_gz_file).split('.')[0]
    os.makedirs(folder)
    gene_pair_file = os.path.basename(gene_pair_gz_file).strip('.gz')

    # extract the txt file from the zip file and place it in the created folder
    with gzip.open(gene_pair_gz_file, 'rb') as gz_file:
        content = gz_file.read()
    with open(os.path.join(folder, gene_pair_file), 'wb') as extracted_file:
        extracted_file.write(content)

    # split(os.path.join(folder, gene_pair_file))

    trim(os.path.join(folder, gene_pair_file))

    tissue = gene_pair_file.split('.')[0]

    eqtl_process(os.path.join(folder, 'newdata' + gene_pair_file), os.path.join(folder, tissue + '_eqtl'))

    eliminateOverlap(os.path.join(folder, tissue + '_eqtl'), os.path.join(folder, tissue + '_eqtl2'), os.path.join(folder, 'overlaps'), os.path.join(folder, 'unmatched'))

    bedtoolIntersect('CREbedDBenhancers_10092018', os.path.join(folder, tissue + '_eqtl2'), os.path.join(folder, tissue + '_intersect'))

    eqtllinks(os.path.join(folder, tissue + '_intersect'), os.path.join(folder, 'links_' + tissue + '_eqtl'))

    links_files.append(os.path.join(folder, 'links_' + tissue + '_eqtl'))

    linksDB(os.path.join(folder, 'links_' + tissue + '_eqtl'), os.path.join(folder, 'linksnum_' + tissue + '_eqtl'))

    linksnum_files.append(os.path.join(folder, 'linksnum_' + tissue + '_eqtl'))

concatenate(links_files, 'linksDBeqtl')
concatenate(linksnum_files, 'linksDBnumeqtl')