from pybedtools import BedTool

def celltypev2(input_file, cell_type, output_file):
    with open(input_file, 'r') as input, open(output_file, 'w') as out:
        for line in input:
            chr1, start1, end1, chr2, start2, end2, number, pval = line.strip().split('\t')
            name = f"{chr1}:{start1}..{end1}-{chr2}:{start2}..{end2}"
            column = f"{tadcell}_{pval}"
            out.write(f"{chr1}\t{start1}\t{end1}\t{name}\t{column}\n")
            out.write(f"{chr2}\t{start2}\t{end2}\t{name}\t{column}\n")


def bedtoolIntersect(gene_file, input_file, output_file):
    '''
        Overlap the enhancers with the regions in this file
    '''
    gene_file_bedtoolFile = BedTool(gene_file)
    input_file_bedtoolFile = BedTool(input_file)

    intersection = gene_file_bedtoolFile.intersect(input_file_bedtoolFile, wa=True, wb=True, F=0.5, f=0.5, e=True)
    with open(output_file, 'w') as out:
        for elt in intersection:
            out.write(str(elt))


def chiageneTSS(input_file, output_file):
    hash = {}
    with open(input_file, 'r') as input:
        for line in input:
            gene, tss, end, start, strand, chr = line.strip().split('\t')
            if strand == 1:
                end = tss
                start = tss - 600
                chr = "chr" + chr
                key = f"{chr}\t{start}\t{end}\t{gene}"
                hash[gene] = key
            elif strand == -1:
                start = tss
                end = tss + 600
                chr = "chr" + chr
                key = f"{chr}\t{start}\t{end}\t{gene}"
                hash[gene] = key

    final = {}
    for ensgene in sorted(panther_mapping.keys()):
        if ensgene in hash:
            line = hash[ensgene]

    with open(output_file, 'w') as out:
        for line in sorted(final.keys()):
            out.write(line + '\n')


def chiagenesenhancers(input_file, gene_file, output_file):
    chia = {}
    enhancer = {}
    with open(input_file, 'r') as input:
        chr, start, end, enhID, one, two, three, chia, cell = line.strip().split()
        thing1, thing2 = chia.split('-')
        cell, pvalue = cell.split('_')
        thing1 = thing1 + '\t' + cell
        thing2 = thing2 + '\t' + cell
        chia.setdefault(thing1, {})[thing2] = pvalue
        chia.setdefault(thing2, {})[thing1] = pvalue
        region = f"{one}:{two}..{three}\t{cell}"
        enhancer.setdefault(region, {})[enhID] = 1

    genes = {}
    with open(gene_file, 'r') as gene:
        for line in gene:
            chr, start, end, geneID, one, two, three, chia_key, cell = line.strip().split('\t')
            thing1, thing2 = chia.split('-')
            cell, pvalue = cell.split('_')
            thing1 = thing1 + '\t' + cell
            thing2 = thing2 + '\t' + cell
            chia.setdefault(thing1, {})[thing2] = pvalue
            chia.setdefault(thing2, {})[thing1] = pvalue
            region = f"{one}:{two}..{three}\t{cell}"
            genes.setdefault(geneID, {})[region] = 1
    
    count1 = 0
    count2 = 0
    final = {}
    for ensg in genes:
        for reg in genes[ensg]:
            if reg in chia:
                for partner in chia[reg]:
                    count1 += 1
                    if partner in enhancer:
                        for enh in enhancer[partner]:
                            count2 += 1
                            tissue = reg.split('\t')[1]
                            cell = tissue.split('-')[0]
                            pvalue = chia[partner][reg]
                            final.setdefault(enh, {}).setdefault(ensg, {})[cell] = 1

    with open(output_file) as out:
        for one in final:
            for two in final[one]:
                for three in final[one][two]:
                    now = panther_mapping[two]
                    tiss = tissues[three]
                    out.write(f"{one}\t{now}\t{tiss}\t1")


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
        tissue.replace('-', '')
        tissues[tissue] = tissueID