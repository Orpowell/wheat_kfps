import pandas as pd

primary_modules = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation/T2T_KFP_primary_modules.tsv"
annotations = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_annotation/T2T_KFP_full_length.faa.tsv"
output_annotation = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_auxiliary_modules.itol"
output_table = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_auxiliary_modules.tsv"
output_count = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/T2T_chinese_spring_analysis/analysis/KFP_phylogenetics/T2T_KFP_auxiliary_module_counts.tsv"
#pd.options.mode.copy_on_write = True
class Sorter:

    def __init__(self, name, domains, length, primary_length) -> None:
        self.name: str = name
        self.pa_ratio = primary_length / length
        self.data: pd.DataFrame = domains
        self.missing = set()  # type: ignore
        self.aux = set()  # type: ignore
        self.annotation: str = "***TBD***"
  
    def __repr__(self) -> str:
        return f"{self.name} {self.annotation}"

    def get_aux(self):

        aux_only: pd.DataFrame = self.data[(~self.data.overlap) & (~self.data.analysis.isin(['MobiDBLite', 'Coils', "PANTHER"]))]
        
        if (len(aux_only) == 0):

            primary_domains_df: pd.DataFrame = self.data[self.data.overlap]

            primary_domains = primary_domains_df.interpro_accession.unique()

            if 'G3DSA:3.30.200.20:FF:000465' in primary_domains:
                self.aux.add('Singleton_crRLK-6')
                return
                
            else:

                self.aux.add("Singleton_kinase")

            return

        domain_set = aux_only.interpro_accession.unique()

        for domain in domain_set:

            if domain == 'IPR027806':
                self.aux.add("Endonuclease")
            
            if domain == 'IPR021094':
                self.aux.add("TM")
        
            if domain in ("IPR002182", "IPR04219", "IPR042197", "IPR027417", "IPR003593", "PF23559"):
                self.aux.add("NBARC")
                continue

            if domain in ("IPR032675", "IPR055414", "IPR056789"):
                self.aux.add("LRR")
                continue

            if domain in ("IPR041118", "IPR038005"):
                self.aux.add("CC")
                continue

            if domain in ("IPR008962", "IPR013783", "IPR000535"):
                self.aux.add("MSP")
                continue

            if domain in ("IPR006121", "IPR036163", "G3DSA:3.30.70.100"):
                self.aux.add("HMA")
                continue

            if domain == "IPR049065":
                self.aux.add("Nakanori")
                continue

            if domain in ("IPR011009", "IPR000719", "IPR001245"):
                self.aux.add("Kinase")
                continue

            if domain in ("IPR015943", "IPR036322", "IPR001680", "IPR020472", "IPR011659", "IPR050844", "IPR011047", "IPR011042"):
                self.aux.add("WD40")
                continue

            if domain in ("IPR036465", "IPR002035"):
                self.aux.add("VWA")
                continue

            if domain in ("IPR025886", "IPR001932"):
                self.aux.add("PP2")
                continue

            if domain in ("IPR036404", "IPR001229", "IPR033734", "PTHR46506"):
                self.aux.add("Jacalin")
                continue

            if domain in ("IPR000008", "IPR035892", "IPR013583", "IPR014020"):
                self.aux.add("C2")
                continue

            if domain in ("IPR017927", "IPR017938", "IPR039261", "IPR000778", "IPR017938", "IPR013121", "IPR013112", "IPR036188", "IPR020946"):
                self.aux.add("FAD")
                continue

            if domain == "IPR021864":
                self.aux.add("DUF3475")
                continue

            if domain == "IPR007700":
                self.aux.add("DUF668")
                continue

            if domain == "IPR005202":
                self.aux.add("GRAS")
                continue

            if domain in ("IPR046829", "IPR046830"):
                self.aux.add("Calmodulin")
                continue

            if domain == "IPR036537":
                self.aux.add("MLKL")
                continue

            if domain in ("IPR013083", "IPR000306", "IPR011011", "IPR017455", "IPR057024"):
                self.aux.add("zinc_finger")
                continue

            if domain == "IPR005559":
                self.aux.add("CG1")
                continue

            if domain == "IPR056284":
                self.aux.add("A9")
                continue

            if domain == "IPR055474":
                self.aux.add("DUF7046")
                continue

            if domain in ("IPR018181", "IPR029047", "IPR013126"):
                self.aux.add("HSP70")
                continue

            if domain in ("IPR006189", "IPR042240"):
                self.aux.add("CHASE")
            
            if domain in ("IPR004158"):
                self.aux.add("CHASE")
            
            if domain in ("IPR002902", "IPR038408"):
                self.aux.add("GNK2")
            
            if domain in ("IPR026961"):
                self.aux.add("PGG")
            
            if domain in ("IPR011990"):
                self.aux.add("Tetratricopeptide")
            
            if domain == "TMhelix":
                self.aux.add("TM")
   			
            else:
                self.missing.add(f"Missing - {domain}")

        if self.name == "XBI87574.1":
        	print(self.name, self.aux, self.data.sort_values(by='start'))
        
    

    def get_annotation(self):
  
        if len(self.aux) == 1:
            
            if {"G3DSA:6.10.140.890"} == self.aux:
                self.annotation = "Singleton_kinase"

            if {"G3DSA:6.10.140.890"} == self.aux:
                self.annotation = "Singleton_kinase"
            
            if {"NBARC"} == self.aux:
                self.annotation = "NBARC-LRR"

            else:
                self.annotation: str = list(self.aux)[0]

        if len(self.aux) == 2:

            if ('MSP' in self.aux) and ('WD40' in self.aux):
                self.annotation = "MSP-WD40"

            if ('MSP' in self.aux) and ('Jacalin' in self.aux):
                self.annotation = "MSP-Jacalin"

            if ('GRAS' in self.aux) and ('Calmodulin' in self.aux):
                self.annotation = "GRAS"
            
            if ('A9' in self.aux):
                self.annotation = "A9"
            
            if ('LRR' in self.aux) and ('NBARC' in self.aux):
                self.annotation = "NBARC-LRR"
            
            if ('DUF3475' in self.aux) and ('DUF668' in self.aux):
                self.annotation = "PSK"
            
            if {'CC', 'NBARC'} == self.aux:
                self.annotation = "NLR"
            
            if {'MSP', 'CC'} == self.aux:
                self.annotation = "MSP"
            
            if ({'VWA', 'Jacalin'} == self.aux) or ({'Kinase', 'Jacalin'} == self.aux):
                self.annotation = 'Jacalin'
            
            if ({'C2', 'Kinase'} == self.aux) or ({'C2', 'TM'} == self.aux):
                self.annotation = "C2"
            
            if ( {'MSP', 'Kinase'} == self.aux):
                self.annotation = "MSP"
            
            if ( {'MSP', 'VWA'} == self.aux):
            	self.annotation = "MSP-VWA"
            
            if {'FAD', 'TM'} == self.aux:
                self.annotation = "FAD"
            
            if {'TM', 'Kinase'} == self.aux:
                self.annotation = "Kinase"
            
            if {'PGG', 'TM'} == self.aux:
                self.annotation = "TM"

            if {'TM', 'HMA'} == self.aux:
                self.annotation = "HMA"

        if len(self.aux) == 3:
            if len(self.aux.intersection({"CC", "LRR", "NBARC"})):
                self.annotation = "NLR"

            if {'Kinase', 'WD40', 'MSP'} == self.aux:
                self.annotation = "MSP-WD40"
                
            if {'MSP', 'TM', 'CC'} == self.aux:
            	self.annotation = "MSP"
        
        if len(self.aux) == 4:
            if {'LRR', 'MSP', 'NBARC', 'CC'} or {'CC', 'MSP', 'NBARC', 'LRR'} == self.aux:
                self.annotation = "MSP"

            if {'LRR', 'NBARC', 'Nakanori', 'HMA'} == self.aux:
                self.annotation = "Nakanori"
            
            if {'PGG', 'CC', 'LRR', 'NBARC'} or {'NBARC', 'TM', 'CC', 'LRR'} == self.aux:
                self.annotation = "NLR"

            if {'CC', 'NBARC', 'LRR', 'Jacalin'} == self.aux:
                self.annotation = "Jacalin"
        
        if len(self.aux) == 5:
        	if {'LRR', 'CG1', 'NBARC', 'MSP', 'CC'} == self.aux:
        		self.annotation = "NLR"

        if self.name in {"XBI24482.1", 'XBH81912.1', 'XBI34015.1'}:
            self.annotation = "Singleton_kinase"
        
        if self.name in {"XBI77658.1", "XBJ13330.1", " XBI77634.1", "XBI05726.1", "XBJ13330.1", "XBI33615.1", "XBI33411.1", "XBI33413.1", "XBI77634.1"}:
        	self.annotation = "MSP"
        
        if self.name == "XBI87574.1":
        	self.annotation = "Nakanori"

df_primaries = pd.read_csv(primary_modules, sep="\t", header=None)

df_primaries.columns = ['target', 'begin', 'end']

df_primaries['p_coords'] = pd.IntervalIndex.from_arrays(df_primaries.begin, df_primaries.end)
df_primaries['p_len'] = df_primaries.end - df_primaries.begin


df_annotations = pd.read_csv(annotations, sep="\t", header=None)

df_annotations.columns = ['target', 
                          'md5', 
                          'sequence_len', 
                          'analysis', 
                          'signature_accession',
                          'signature_description',
                          'start',
                          'stop',
                          'score',
                          'status',
                          'date',
                          'interpro_accession',
                          'interpro_description',
                          'go_annotations',
                          'pathway']

df_annotations['h_coords'] = pd.IntervalIndex.from_arrays(df_annotations.start, df_annotations.stop, closed='both')

merged_df = df_annotations.merge(df_primaries)
merged_df['overlap'] = [data.p_coords.overlaps(data.h_coords) for index, data in merged_df.iterrows()]
merged_df.interpro_accession[merged_df.interpro_accession == '-'] = merged_df.signature_accession

aux_info = []

for index, data in merged_df.groupby(by='target'):

    kfp = Sorter(name=index, domains=data, length=data.sequence_len.to_list()[0], primary_length=data.p_len.to_list()[0])
    kfp.get_aux()
    kfp.get_annotation()
    aux_info.append(kfp)

annotated_auxs = dict()

for kfp in aux_info:

    annotation = kfp.annotation

    if annotation in annotated_auxs.keys():
        annotated_auxs[annotation] += 1
    
    else:
        annotated_auxs[annotation] = 1

counter = 0
total_kfps = 0
for k, v in annotated_auxs.items():
    
    if v >= 1:
        counter += 1
    total_kfps += v

df = pd.DataFrame(annotated_auxs, index=[0]).transpose()
df.columns = ['counts']
df.sort_values(by="counts", ascending=False, inplace=True)

colours = [
    "#808080",
    "#e6194b",
    "#000000",
    "#3cb44b",
    "#ffe119",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#46f0f0",
    "#f032e6",
    "#bcf60c",
    "#fabebe",
    "#008080",
    "#e6beff",
    "#9a6324",
    "#fffac8",
    "#800000",
    "#aaffc3",
    "#808000",
    "#ffd8b1",
    "#000075",
]

display = df[(df.counts >= 0)]

display.sort_values(by="counts", ascending=False, inplace=True)

display_colours = {'Singleton_kinase': '#000000', 
					 'Kinase': '#e6194b', 
					 'Singleton_crRLK-6': '#000000', 
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
					 'PSK': '#9a6324', 
					 'A9': '#fffac8', 
					 'MSP-Jacalin': '#800000', 
					 'MLKL': '#fabebe', 
					 'Jacalin': '#808000', 
					 'NBARC': '#ffd8b1', 
					 'zinc_finger': '#000075'}


display['colour'] = display_colours
display['colour'].fillna('#FFFFFF', inplace=True)
colour_dict = dict(zip(display.index, display.colour))
print(display)
# Remove cloned KFPs from counts
display.at["Kinase", 'counts'] = (display.loc["Kinase"]['counts'] - 4)
display.at["HMA", 'counts'] = display.loc["HMA"]['counts'] - 2
display.at["MSP", 'counts'] = display.loc["MSP"]['counts'] - 2
display.at["MSP-WD40", 'counts'] = display.loc["MSP-WD40"]['counts'] - 1
display.at["VWA", 'counts'] = display.loc["VWA"]['counts'] - 2
display.at["NBARC-LRR", 'counts'] = display.loc["NBARC-LRR"]['counts'] - 1
display.at["PSK", 'counts'] = display.loc["PSK"]['counts'] - 1
display.at["TM", 'counts'] = display.loc["TM"]['counts'] - 1
display.at["C2", 'counts'] = display.loc["C2"]['counts'] - 4
display.at["MLKL", 'counts'] = display.loc["MLKL"]['counts'] - 1
print(display)

display.to_csv(output_count, header=None, index=True, columns=['counts'])

with open(output_annotation, "w+") as out:

    legend_shapes = list(map(lambda x: 1, display.index))
    legend_labels = ','.join(display.index.to_list())
    legend_colour = ','.join(display.colour.to_list())

    out.write("DATASET_COLORSTRIP\n")
    out.write("SEPARATOR COMMA\n")
    out.write("DATASET_LABEL,KFP_Annotations\n")
    out.write("HEIGHT_FACTOR,5\n")
#    out.write(f"FIELD_COLORS,{','.join(display.colour.to_list())}\n")
#    out.write(f"FIELD_LABELS,{','.join(display.index.to_list())}\n")
    out.write("BORDER_WIDTH,0.5\n")
    out.write("BORDER_COLOR,#000000\n")
    out.write("COMPLETE_BORDER,1\n")
    out.write("STRIP_WIDTH,100\n")
    out.write("STRIP_LABEL_AUTO_COLOR,1\n")
    # Legend
    #out.write("LEGEND_TITLE,Dataset_legend\n")
    #out.write("LEGEND_HORIZONTAL,1\n")
    #out.write(f"LEGEND_COLORS,{legend_colour}\n")
    #out.write(f"LEGEND_LABELS,{legend_labels}\n")
    #out.write(f"LEGEND_SHAPES,{legend_shapes}\n")
    out.write("DATA\n")

    for kfp in aux_info:
        try:
            out.write(f"{kfp.name},{colour_dict[kfp.annotation]},{kfp.annotation}\n")
        except KeyError:
            pass

with open(output_table, "w+") as out:
	for kfp in aux_info:
		out.write(f"{kfp.name}\t{kfp.annotation}\n")


