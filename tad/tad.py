from math import floor, ceil

def bTADoverlap(input_file, output_file, cell_type):
    tad_file = open(input_file, 'r')
    tad_lines = tad_file.readlines()
    tad_file.close()

    hash = {}
    for line in tad_lines:
        chr, start1, end1, location, pval = line.strip().split('\t')
        score = pval + '|' + cell_type
        location_split = location.split('___')
        one = location_split[0]
        two = location_split[1]
        three = one.split('|')[2]
        four = two.split('|')[2]
        five = three.split(':')[1]
        six = four.split(':')[1]
        s1, s2 = five.split('-')
        e1, e2 = six.split('-')
        start = floor((int(s1) + int(s2))/2)
        end = ceil((int(e1) + int(e2)))/2
        total = f"{chr}\t{start}\t{end}\t{location}\t{score}"
        hash[total] = 1

    out = open(output_file, 'w')
    for key in sorted(hash.keys()):
        out.write(key + '\n')
    out.close()

bTADoverlap('ENCFF274VJU.bed', 'bTADorder', 'Caki2')

