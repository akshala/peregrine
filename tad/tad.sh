python TADorder.py
echo 'TADorder done'

python bedtoolIntersect.py
echo 'bedtoolIntersect done'

python split.py
echo 'split done'

python bTADlinks.py
echo 'bTADlinks done'

python concatenate_linksbTAD.py
echo 'concatenate_linksbTAD done'

python tTADlinks.py
echo 'tTADlinks done'

python linksbTAD_tissueReplace.py
echo 'linksbTAD_tissueReplace done'

python final_concatenate.py
echo 'final_concatenate done'

python cutdowntad.py
echo 'cutdowntad done'

python orderlinks.py
echo 'orderlinks done'