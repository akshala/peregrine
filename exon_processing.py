import pandas as pd 
import sys

filename = sys.argv[1]  # exon_file_default.txt
mapping_file = sys.argv[2]   # patch_chromosome_mapping.txt
out = sys.argv[3]   # exons_genes.txt

df = pd.read_csv(filename, delimiter='\t')
new_df = pd.DataFrame()

with open(mapping_file) as f:
    lines = f.readlines()

mapping = {}
for elt in lines:
    elt = elt.split(' ')
    elt = [x for x in elt if x != '']
    mapping[elt[0]] = 'chr' + elt[1]

df['Chromosome/scaffold name'] = df['Chromosome/scaffold name'].apply(lambda x: mapping[x] if x in mapping else 'chr' + x)
new_df['Chromosome number'] = df['Chromosome/scaffold name']
new_df['Exon region start (bp)'] = df['Exon region start (bp)']
new_df['Exon region end (bp)'] = df['Exon region end (bp)']
new_df['blank'] = ['-'] * len(df)
new_df['Gene stable ID'] = df['Gene stable ID']
new_df['Gene type'] = df['Gene type']
new_df['known'] = ['KNOWN'] * len(df)
new_df['lciz'] = ['LCIZ'] * len(df)

new_df.to_csv(out, sep='\t', header=None, index=False)
