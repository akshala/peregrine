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