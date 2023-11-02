import pandas as pd 
import sys

filename = sys.argv[1]  # TSS_default.txt
mapping_file = sys.argv[2]   # patch_chromosome_mapping.txt
out = sys.argv[3]   # TSSgenesbed

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
new_df['TSS start'] = df['Transcription start site (TSS)']
new_df['TSS end'] = df['Transcription start site (TSS)']
new_df['Gene stable ID'] = df['Gene stable ID']

new_df.to_csv(out, sep='\t', header=None, index=False)
