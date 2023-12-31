## eQTL
The eQTL folder has the code for generating enhancer-gene links from eQTL data taken from the source - https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz. 

```
python eqtl.py
```

Running the above command generates linksDBeqtl and linksDBnumeqtl files that are the concatention of all the files of respective nature generated from different tissue variant gene pair files. 

#### linksDBeqtl
```
$ head linksDBeqtl 
7	HUMAN|HGNC=23517|UniProtKB=Q8N2H3	10_100169983_A_G_b37	4.57077e-11	23	2
7	HUMAN|HGNC=23517|UniProtKB=Q8N2H3	10_100170037_A_G_b37	7.59471e-11	23	2
7	HUMAN|HGNC=23517|UniProtKB=Q8N2H3	10_100170209_T_C_b37	4.57077e-11	23	2
12	HUMAN|HGNC=28228|UniProtKB=Q96DD0	1_100502868_T_TATAC_b37	4.09143e-08	23	2
16	HUMAN|HGNC=28228|UniProtKB=Q96DD0	1_100555192_A_G_b37	1.5892e-07	23	2
62	HUMAN|HGNC=25613|UniProtKB=Q69YN2	10_101854930_CA_C_b37	5.74856e-06	23	2
63	HUMAN|HGNC=25613|UniProtKB=Q69YN2	10_101863410_C_A_b37	3.37568e-06	23	2
88	HUMAN|HGNC=17452|UniProtKB=Q9HBK9	10_104585405_A_T_b37	2.15079e-08	23	2
92	HUMAN|HGNC=7642|UniProtKB=Q99747	18_10497034_G_A_b37	2.43176e-06	23	2
317	HUMAN|HGNC=28769|UniProtKB=Q6UX65	1_111733724_G_C_b37	1.98817e-07	23	2
```

#### linksDBnumeqtl
```
$ head linksDBnumeqtl 
enhancer	gene	tissue	number_of_eQTL	assay
1026	HUMAN|HGNC=378|UniProtKB=O43823	23	1	2
12	HUMAN|HGNC=28228|UniProtKB=Q96DD0	23	1	2
1212	HUMAN|HGNC=264|UniProtKB=P29275	23	1	2
1213	HUMAN|HGNC=264|UniProtKB=P29275	23	1	2
1295	HUMAN|HGNC=19359|UniProtKB=Q5VWW1	23	1	2
1410	HUMAN|HGNC=28793|UniProtKB=O95568	23	1	2
1411	HUMAN|HGNC=28793|UniProtKB=O95568	23	1	2
16	HUMAN|HGNC=28228|UniProtKB=Q96DD0	23	1	2
1630	HUMAN|HGNC=11984|UniProtKB=Q6ZVM7	23	1	2
```