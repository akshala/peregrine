
with open('hg19.fa', 'r') as bowtieref, open('human.hg19.genome', 'w') as bedtoolsgenome:
    chr = None
    length = 0
    for line in bowtieref:
        if line.startswith('>'):
            if chr != None:
                bedtoolsgenome.write(chr + '\t' + str(length) + '\n')
            chr = line.strip()[1:]
            length = 0
        else:
            length += len(line.strip())
    if chr != None:
        bedtoolsgenome.write(chr + '\t' + str(length) + '\n')