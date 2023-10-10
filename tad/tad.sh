python3 TADorder.py
echo 'TADorder done'

python3 bedtoolIntersect.py
echo 'bedtoolIntersect done'

python3 split.py
echo 'split done'

python3 bTADlinks.py
echo 'bTADlinks done'

python3 concatenate_linksbTAD.py
echo 'concatenate_linksbTAD done'

python3 tTADlinks.py
echo 'tTADlinks done'

python3 linksbTAD_tissueReplace.py
echo 'linksbTAD_tissueReplace done'

python3 final_concatenate.py
echo 'final_concatenate done'