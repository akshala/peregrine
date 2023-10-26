#!/bin/bash

python3 -m venv venv
pip install -r requirements.txt
source venv/bin/activate

cd eQTL
python eqtl.py

cd ../heirarchicalTAD
python heirarchicalTAD.py

cd ../tad
./tad.sh

cd ../chia_pet
python ChIA_PET.py

deactivate