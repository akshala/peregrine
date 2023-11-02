# Peregrine

To get the final Peregrine data which consists of enhancer-gene links run the following bash script command.

```
./steps.sh
```
The above command creates a python3 virtual environment, installs all the requirements, activates the environment and then runs all the python script for eQTL, heirarchical TAD, TAD and ChIA-PET data. 

This final data is coming from 4 different data sources. 
1. [eQTL](eQTL/)
2. [heirarchical TAD](heirarchicalTAD/)
3. [TAD](tad/)
4. [ChIA-PET](chia_pet/)

The above folders contain the data, code and details for each of the data sources. 


The exon file has been generated using exon_processing.py. Data was downloaded from https://grch37.ensembl.org/biomart/martview/f4c200b649ceea16a7db20ec6e44b8f4 and the file was renamed as exon_file_default.txt. The patch chromosome mapping has been taken from https://www.biostars.org/p/106355/. Then the following command is run to get exons_genes.txt which is used in the eQTL code and placed in that folder. 
```
python3 exon_processing.py exon_file_default.txt patch_chromosome_mapping.txt exons_genes.txt
```

The TSSgenesbed has been generated using TSSgenesbed_generation.py. Data was downloaded from https://grch37.ensembl.org/biomart/martview/ae31c4c3f0b33a8c03a590d93265bf56 and the file was renamed as TSS_default.txt. Then the following command was run to get TSSgenesbed which is used in TAD and the placed in that folder. 
```
python3 TSSgenesbed_generation.py TSS_default.txt patch_chromosome_mapping.txt TSSgenesbed
``` 