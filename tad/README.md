## TAD
The tad folder has the code for generating enhancer-gene links from  TAD data taken from the source - https://www.encodeproject.org/files/ENCFF274VJU/, https://www.encodeproject.org/files/ENCFF588KUZ/. 

```
python tad.py
```

If there is any memory issue in running this code then the bash script for TAD can be run which runs different TAD script files (they are just the original ta.py file divided into separate python files to ease with memory release.)

```
./tad.sh
```

Running the above command results in linksDBtad enhancer-gene links file.

#### linksDBtad
```
$ head linksDBtad 
1	HUMAN|HGNC=15846|UniProtKB=Q9NP74	65	3
1	HUMAN|HGNC=321|UniProtKB=P35573	65	3
1000	HUMAN|HGNC=2523|UniProtKB=Q99895	65	3
1000	HUMAN|HGNC=28670|UniProtKB=Q96C19	65	3
10007	HUMAN|HGNC=13433|UniProtKB=Q96MS0	65	3
10007	HUMAN|HGNC=17149|UniProtKB=Q96IQ7	65	3
10007	HUMAN|HGNC=17474|UniProtKB=Q96AP7	65	3
10007	HUMAN|HGNC=20573|UniProtKB=Q96QZ0	65	3
10007	HUMAN|HGNC=26266|UniProtKB=Q6P1R3	65	3
10007	HUMAN|HGNC=29551|UniProtKB=Q3YBR2	65	3

```