import pandas as pd
from pybedtools import BedTool

def idmapping_processing():
    file = open('UP000005640_9606.idmapping', 'r')
    lines = file.readlines()
    file.close()

    mapping = {}
    for line in lines:
        line_split = line.strip().split('\t')
        if line_split[0] not in mapping:
            mapping[line_split[0]] = {}
        if line_split[1] == 'Gene_Name':
            mapping[line_split[0]]['Gene_Name'] = line_split[2]
        elif line_split[1] == 'HGNC':
            mapping[line_split[0]]['HGNC'] = line_split[2]
        elif line_split[1]== 'Ensembl':
            mapping[line_split[0]]['Ensembl'] = line_split[2].split('.')[0]

    df = pd.DataFrame(columns=['HGNC ID', 'symbol', 'UniProt', 'Ensembl gene ID'])
    count = 0
    for key, val in mapping.items():
        if ('Ensembl' in val) and ('HGNC' in val) and ('Gene_Name') in val:
            row = [val['HGNC'], val['Gene_Name'], key, val['Ensembl']]
            df.loc[count] = row
            count += 1
    df.to_csv('resultsHGNC.txt', sep='\t', index=False)

idmapping_processing()

def nonegvalues(input_file, output_file, out_file):
    bed_file = open(input_file, 'r')
    lines = bed_file.readlines()
    bed_file.close()

    tossed_out = open(out_file, 'w')
    out = open(output_file, 'w')
    for line in lines:
        line_split = line.strip().split('\t')
        start = line_split[1]
        end = line_split[2]
        if int(start) < 0 or int(end) < 0:
            tossed_out.write(line)
        else:
            out.write(line)

    tossed_out.close()
    out.close()

nonegvalues('hES.enh_1e-4.bed', 'nn_hES.enh_1e-4.bed', 'out') 

def bedtoolIntersect(cred_file, eqtl_file, output_file):
    cred_file_bedtoolFile = BedTool(cred_file)
    eqtl_file_bedtoolFile = BedTool(eqtl_file)

    intersection = cred_file_bedtoolFile.intersect(eqtl_file_bedtoolFile, wa=True, wb=True, f=0.9)
    out = open(output_file, 'w')
    for elt in intersection:
        out.write(str(elt))
    out.close()

bedtoolIntersect('CREbedDBenhancers_10092018', 'nn_hES.enh_1e-4.bed', 'nn_hES.enh_1e-4_intersect')

def HGNC2PANTH(input_file, output_file, unmatched_file, tissue_name):
    HGNC_mapping = open('resultsHGNC.txt', 'r')
    lines = HGNC_mapping.readlines()
    HGNC_mapping.close()

    hash1 = {}
    hash2 = {}
    hash3 = {}
    for line in lines[1:]:
        hgnc, symbol, uniprot, ensg,  = line.strip().split('\t')
        hgncid = hgnc.split(':')[1]
        hgncid = 'HGNC=' + hgncid
        hash1[symbol] = hgncid
        hash2[symbol] = uniprot
        hash3[symbol] = ensg


    panth_mapping = open('pantherGeneList.txt', 'r')
    lines = panth_mapping.readlines()
    panth_mapping.close()

    panth = {}
    for line in lines:
        line_split = line.strip().split()
        longID = line_split[0]
        ensg = line_split[0]
        longID_split = longID.split('|')
        HGNC = longID_split[1]
        uniprotkb = longID_split[2].split('=')[1]
        panth[HGNC] = longID
        panth[uniprotkb] = longID
        panth[ensg] = longID

    input = open(input_file, 'r')
    lines = input.readlines()
    input.close()

    out = open(output_file, 'w')
    unmatched = open(unmatched_file, 'w')

    for line in lines:
        echr, estart, eend, enhID, pchr, pstart, pend, mix = line.strip().split('\t')
        mix_split = mix.split(':')
        gene = mix_split[0]
        pvalue = mix_split[2]
        if gene in hash1:
            hgncID = hash1[gene]
            if hgncID in panth:
                panthid = panth[hgncID]
                out.write(f'{echr}\t{estart}\t{eend}\t{enhID}\t{pchr}\t{pstart}\t{pend}\t{panthid}\t{pvalue}\t{tissue_name}\n')
            else:
                unmatched.write(f'1\t{line}\n')
        elif gene in hash2:
            uniprot = hash2[gene]
            if uniprot in panth:
                panthid = panth[uniprot]
                out.write(f'{echr}\t{estart}\t{eend}\t{enhID}\t{pchr}\t{pstart}\t{pend}\t{panthid}\t{pvalue}\t{tissue_name}\n')
            else:
                unmatched.write(f'1\t{line}\n')
        elif gene in hash3:
            ensg = hash3[gene]
            if ensg in panth:
                panthid = panth[ensg]
                out.write(f'{echr}\t{estart}\t{eend}\t{enhID}\t{pchr}\t{pstart}\t{pend}\t{panthid}\t{pvalue}\t{tissue_name}\n')
            else:
                unmatched.write(f'1\t{line}\n')
        else:
            unmatched.write(f'2\t{line}\n')

    out.close()
    unmatched.close()

HGNC2PANTH('nn_hES.enh_1e-4_intersect', 'out', 'unmatched', 'hES')

def reformat(input_file, output_file):
    input = open(input_file, 'r')
    lines = input.readlines()
    input.close()

    hash = {}
    for line in lines:
        echr, estart, eend, enhID, pchr, pstart, pend, gene, pval, tissue = line.strip().split('\t')
        break_ = f"{echr}\t{estart}\t{eend}\t{enhID}\t{gene}\t{tissue}\tTADinteractions\t{pval}"
        hash[break_] = 1

    out = open(output_file, 'w')
    for key in sorted(hash.keys()):
        out.write(key + '\n')
    out.close()

reformat('out', 'intTADlinks_hES_1e-4')