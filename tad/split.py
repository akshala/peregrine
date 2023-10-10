def split(file):
    '''
        splits input file by chromosome and creates 25 chromosome specific files
    '''
    chromosomes = ['chr' + str(x) for x in range(1, 23)]
    chromosomes.append('chrX')
    chromosomes.append('chrY')

    tad_file = open(file, 'r')
    tad_lines = tad_file.readlines()
    tad_file.close()

    output_files = {chromosome: open(chromosome + file, 'w') for chromosome in chromosomes}

    for line in tad_lines:
        chr = line.strip().split('\t')[0]

        if chr in chromosomes:
            output_files[chr].write(line)

    for output_file in output_files.values():
        output_file.close()

split('TSSbTAD')
split('enhancersbTAD')