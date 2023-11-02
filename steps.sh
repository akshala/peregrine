#!/bin/bash

python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
echo 'venv done'

cd eQTL
python eqtl.py
echo 'eQTL done'

cd ../heirarchicalTAD
python heirarchicalTAD.py
echo 'heirarchical TAD done'

cd ../tad
./tad.sh
echo 'TAD done'

# cd ../chia_pet
# python ChIA_PET.py

deactivate