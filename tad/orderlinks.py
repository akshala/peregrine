def orderlinks(input_file, output_file):
    hash = {}
    tissues = {}
    genes = {}
    count = 0

    with open(input_file, 'r') as input:
        for line in input:
            enhID, gene, tissue, assay = line.strip().split('\t')
            hash.setdefault(enhID, {}).setdefault(gene, {})[tissue] = 1

    with open(output_file, 'w') as out:
        for enhancer in sorted(hash.keys()):
            for geen in sorted(hash[enhancer].keys()):
                genes[geen] = 1
                for cell in sorted(hash[enhancer][geen].keys()):
                    tissues[cell] = 1
                    count += 1
                    out.write(f"{enhancer}\t{geen}\t{cell}\t{assay}\n")

    print(f"This file contains {len(hash)} different enhancers linked to {len(genes)} genes in {len(tissues)} tissues.")

orderlinks('selectTAD', 'linksDBtad')