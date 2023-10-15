## heirarchical TAD
The heirarchicalTAD folder has the code for generating enhancer-gene links from heirarchical TAD data taken from the source - https://www.cs.huji.ac.il/w~tommy/PSYCHIC/. All the 1e-4 bed files from this source are taken and then filtered on the basis of the tissue file before processing. 

```
python3 heirarchicalTAD.py
```

Running the above command results in PSYCHIClinksDB enhancer-gene links file.

```
$ head PSYCHIClinksDB 
EH37E0105481	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
EH37E0105482	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
EH37E0105483	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
EH37E0105484	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
EH37E0105485	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
EH37E0105486	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
EH37E0105487	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
EH37E0105509	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	4.1e-06	4
66387	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	4.1e-06	4
EH37E0105510	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	4.1e-06	4
```