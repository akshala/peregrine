## ChIA-PET
The chia_pet folder has the code for generating enhancer-gene links from ChIA-PET data taken from the source - https://www.encodeproject.org/matrix/?type=Experiment&status=released&searchTerm=chia-pet&biosample_ontology.classification=cell+line&assembly=hg19 (fastq files need to be taken from here). This data was processed using [Mango](https://github.com/dphansti/mango) 

The Mango repository would have to be cloned inside the ChIA-PET folder. To run the mango.R script we need an argsfile.txt which will have the bowtieref and bedtoolsgenome path as stated in the Mango README. For generating the refernce files we will first need to download the [hg19 human reference genome in FASTA format](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) (it is the hg19fa.gz file). This file would be placed in the mango folder and unzipped. Then the following command needs to be run to get the bowtie reference. 
```
bowtie-build ./hg19.fasta ./hg19
```
This will create multiple ebwt files in a folder named hg19. After this we need the bedtoolsgenome file which contains the genome name and length. This can be extracted from the FASTA file of the genome by running the following command.
```
python generate_genome_file.py
```
The above command would create the human.hg19.genome file. Now we can make the argsfile.txt as follows.
```
bowtieref         = /peregrine/chia_pet/mango/hg19/hg19_out
bedtoolsgenome    = /peregrine/chia_pet/mango/human.hg19.genome
```
Then the mango.R script can be run to generate the interactions.fdr.mango files which will be used in our ChIA-PET code.


Run the following command to generate enhancer-gene links from the data created above.
```
python ChIA_PET.py
```
