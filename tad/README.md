## TAD
The tad folder has the code for generating enhancer-gene links from  TAD data taken from the source - https://www.encodeproject.org/files/ENCFF274VJU/, https://www.encodeproject.org/files/ENCFF588KUZ/. 

```
python tad.py
```

If there is any memory issue in running this code then the bash script for TAD can be run which runs different TAD script files (they are just the original ta.py file divided into separate python files to ease with memory release.)

```
./tad.sh
```