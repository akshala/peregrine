def concatenate(input_files, output_file):
    '''
        makes a combined file for all tissues
    '''
    with open(output_file, 'wb') as outfile:
        for file_name in input_files:
            with open(file_name, 'rb') as infile:
                for line in infile:
                    outfile.write(line)

chromosomes = ['chr' + str(x) for x in range(1, 23)]
chromosomes.append('chrX')
chromosomes.append('chrY')
linkbTAD_list = []

for chr in chromosomes:
    linkbTAD_list.append(chr+'linksbTAD')
concatenate(linkbTAD_list, 'linksbTAD')