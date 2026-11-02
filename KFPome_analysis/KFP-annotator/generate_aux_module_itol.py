import sys

# A selection of Auxiliary modules to be visualised
display_colours = {
                    'NO FUSED DOMAIN': '#000000', 
					'Kinase': '#e6194b', 
					'MSP': '#3cb44b', 
					'MSP-WD40': '#bcf60c',
					'NBARC-LRR': '#4363d8', 
					'NLR': '#f58231', 
					'C2': '#911eb4', 
					'WD40': '#46f0f0', 
					'VWA': '#f032e6', 
					'TM': '#ffe119', 
					'PP2': '#aaffc3', 
					'HMA': '#008080', 
					'Nakanori': '#e6beff', 
					'PSK-like': '#9a6324', 
					'MLKL': '#fabebe', 
					'Jacalin': '#808000', 
}

args = sys.argv[1:]

annotations = args[0]
basename = args[1]

output_annotation= f"{basename}_auxiliary_modules.itol"

kfp_dict = dict()

with open(annotations, 'rt') as file:

     for line in file:

        data = line.strip().split('\t')
        kfp = data[0]
        aux = data[1]

        kfp_dict[kfp] = aux

with open(output_annotation, "w+") as out:

    out.write("DATASET_COLORSTRIP\n")
    out.write("SEPARATOR COMMA\n")
    out.write("DATASET_LABEL,KFP_Annotations\n")
    out.write("HEIGHT_FACTOR,10\n")
    out.write("BORDER_WIDTH,0.5\n")
    out.write("BORDER_COLOR,#bbbbbb\n")
    out.write("COMPLETE_BORDER,1\n")
    out.write("STRIP_WIDTH,100\n")
    out.write("STRIP_LABEL_AUTO_COLOR,1\n")
    out.write("DATA\n")

    for kfp, aux in kfp_dict.items():
        try:
            out.write(f"{kfp},{display_colours[aux]},{aux}\n")
        except KeyError:
            out.write(f"{kfp},#ffffff,{aux}\n")

