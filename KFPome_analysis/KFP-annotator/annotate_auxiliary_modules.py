import pandas as pd
import sys

class Sorter:

    def __init__(self, name, domains) -> None:
        self.name: str = name
        self.data: pd.DataFrame = domains
        self.missing = set()  # type: ignore
        self.aux = set()  # type: ignore
        self.annotation: str = "***TBD***"
  
    def __repr__(self) -> str:
        return f"{self.name} {self.annotation}"

    def get_aux(self):

        aux_only: pd.DataFrame = self.data
        
        if len(aux_only) == 0:

            self.aux.add("NO FUSED DOMAIN")

            return

        domain_set = aux_only.interpro_accession.unique()

        for domain in domain_set:

            if domain in ('IPR013194', 'IPR031693'):
                self.aux.add('Sin3-like')
                continue

            if domain == 'IPR001554':
                self.aux.add('Glycoside Hydrolase')

            if domain in ("IPR006652", "IPR015915"):
                self.aux.add('Kelch')
                continue

            if domain == 'IPR008480':
                self.aux.add('DUF761')
                continue

            if domain in ("IPR017923", "IPR035441"):
                self.aux.add("TFIIS")
                continue

            if domain in ("IPR036388", "IPR058922"):
                self.aux.add("Winged Helix")
                continue

            if domain == 'IPR027806':
                self.aux.add("Endonuclease")
            
            if domain == 'IPR021094':
                self.aux.add("TM")
                continue
        
            if domain in ("IPR002182", "IPR04219", "IPR042197", "IPR027417", "IPR003593", "PF23559"):
                self.aux.add("NBARC")
                continue

            if domain in ("IPR032675", "IPR055414", "IPR056789"):
                self.aux.add("LRR")
                continue

            if domain in ("IPR041118", "IPR038005", 'G3DSA:1.20.5.4130'):
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

            if domain in ("IPR013083", "IPR000306", "IPR011011", "IPR017455", "IPR057024", "IPR001878"):
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

            if domain in ('G3DSA:3.10.20.90', 'IPR000626'):
                self.aux.add('Ubiquitin')
                continue

            if domain == 'IPR033389':
                self.aux.add('AUX/IAA')

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

            if domain in ('IPR023175','IPR039431'):
                self.aux.add("Callose Synthase")
                continue

            if domain in ('IPR032881'):
                self.aux.add('PHD Finger')

            if domain in ('IPR017438', 'IPR001206','IPR000756'):
                self.aux.add('ATP-NAD kinase')

            if domain in ('IPR012340', 'IPR007199'):
                self.aux.add('Replication Factor A')
                continue

            if domain in ('IPR036770', 'IPR002110', 'IPR002110', 'IPR036770', 'IPR002110'):
                self.aux.add('Ankyrin')
                continue

            if domain == 'IPR048628':
                self.aux.add('SEC3-like')
            
            if domain in ('IPR019312', 'IPR007207'):
                self.aux.add('CCR4-NOT Component')
                continue

            if domain == 'IPR022228':
                self.aux.add('DUF3755')
                continue

            if domain in ('IPR004883', 'IPR004883'):
                self.aux.add('LOB')

            if domain in ('IPR046342', 'IPR000644'):
                self.aux.add('CBS')

            if domain == 'IPR001005':
                self.aux.add('MYB')

            if domain == 'IPR036397':
                self.aux.add('Ribonuclease H')

            if domain == 'PF23733':
                self.aux.add('Glutaredoxin')

            if domain in ('IPR012677', 'IPR000504'):
                self.aux.add('RNA recognition motif')

            if domain == 'G3DSA:6.10.140.890':
                self.aux.add('Proton transporter')
            
            else:
                self.missing.add(f"Missing - {domain}")

        if self.name == "XBI87574.1":
        	print(self.name, self.aux, self.data.sort_values(by='start'))

    def get_annotation(self):

        if len(self.aux) == 1:
            self.annotation = list(self.aux)[0]

            if self.annotation == 'DUF668':
                self.annotation = 'PSK-like'

            return
        
        match sorted(list(self.aux)):

            case ['DUF761', 'MSP', 'TM']:
                self.annotation = 'MSP-DUF761'

            case ['MSP', 'NBARC']:
                self.annotation = 'MSP'

            case ['MLKL']:
                self.annotation = 'MLKL'
                return

            case ['Jacalin']:
                self.annotation = 'Jacalin'
                return

            case ['Kinase', 'MLKL']:
                self.annotation = 'MLKL'
                return

            case ['DUF3475', 'DUF668']:
                self.annotation = 'PSK-like'
                return

            case ['Kinase']:
                self.annotation = 'Kinase'
                return
        
            case ['WD40']:
                self.annotation = 'WD40'
                return

            case ['HMA']:
                self.annotation = 'HMA-like'
                return
            
            case ["NO FUSED DOMAIN"]:
                self.annotation = "No Fused Domain"
                return

            case ['Nakanori']:
                self.annotation = 'Nakanori'
                return

            case ['MSP']:
                self.annotation = 'MSP'
                return

            case ['C2']:
                self.annotation = "C2"
                return

            case ['C2', 'Kinase']:
                self.annotation = "C2"
                return
                    
            case ['MSP', 'Kinase']:
                self.annotation = "MSP"
                return

            case ['MSP', 'VWA']:
                self.annotation = "MSP-VWA"
                return

            case ['FAD', 'TM']:
                self.annotation = "FAD"
                return

            case ['MSP', 'WD40']:
                self.annotation = "MSP-WD40"
                return

            case ['LRR', 'Winged Helix']:

                self.annotation = "NBARC-LRR"
                return

            case ['LRR', 'NBARC', 'Winged Helix']:

                self.annotation = "NBARC-LRR"
                return

            case ['CC', 'LRR', 'NBARC', 'Winged Helix']:

                self.annotation = "NLR"
                return

            case ['PP2']:

                self.annotation = "PP2"
                return

            case ['CC', 'LRR', 'MSP', 'NBARC', 'Winged Helix']:

                self.annotation = "MSP"
                return


            case ['C2', 'TM']:

                self.annotation = "C2"
                return

            case ['A9', 'DUF7046']:

                self.annotation = 'A9'
                return

            case ['Calmodulin', 'GRAS']:

               self.annotation = 'Calmodulin-GRAS'
               return

            case ['HMA', 'TM']:

                self.annotation = 'HMA'
                return

            case ['Jacalin', 'MSP']:

                self.annotation = 'Jacalin-MSP'
                return

            case ['CC', 'LRR', 'NBARC', 'TM', 'Winged Helix']:

                self.annotation = 'TM'
                return

            case ['CC', 'MSP', 'NBARC']:

                self.annotation = 'MSP'
                return

            case ['CC', 'NBARC', 'Winged Helix']:

                self.annotation = 'NLR'
                return

            case ['CC', 'CG1', 'LRR', 'MSP', 'NBARC', 'Winged Helix']:

                self.annotation = 'NLR'
                return

            case ['HMA', 'LRR', 'NBARC', 'Nakanori', 'Winged Helix']:

                self.annotation = 'Nakanori'
                return

            case ['Kinase', 'LRR', 'Winged Helix']:

                self.annotation = 'WHD-LRR'
                return

            case ['MSP', 'Winged Helix']:
                self.annotation = 'MSP'
                return

            case ['Kinase', 'MSP', 'WD40']:
                self.annotation = 'MSP-WD40'
                return

            case ['MSP', 'TM', 'WD40']:
                self.annotation = 'MSP-WD40'

            case ['LRR', 'NBARC']:
                self.annotation = 'NBARC-LRR'

            case ['CC', 'MSP', 'TM']:
                self.annotation = 'MSP'

            case ['Kinase', 'VWA']:

                self.annotation = 'VWA'

            case ['MSP', 'Tetratricopeptide']:

                self.annotation = 'MSP-Tetratricopeptide'

            case ['LRR', 'MSP', 'NBARC', 'Winged Helix']:
                self.annotation = 'MSP'

            case ['CC', 'LRR', 'NBARC', 'PGG', 'TM', 'Winged Helix']:
                self.annotation = 'NLR'

            case ['CC', 'Jacalin', 'LRR', 'NBARC', 'Winged Helix']:
                self.annotation = 'Jacalin'

            case ['CC', 'NBARC']:

                self.annotation = 'CC-NBARC'

            case ['LRR', 'NBARC', 'TM', 'Winged Helix']:

                self.annotation = 'NLR'

            case ['PGG', 'TM']:

                self.annotation = 'PGG'

            case ['NBARC', 'Winged Helix']:

                self.annotation = 'NBARC'

            case ['Kinase', 'WD40']:

                self.annotation = 'WD40'

            case ['Jacalin', 'Kinase']:

                self.annotation = 'Jacalin'

            case ['Kinase', 'MSP']:

                self.annotation = 'MSP'

            case ['MSP', 'TM']:

                self.annotation = 'MSP'
            
            case ['Jacalin', 'VWA']:

                self.annotation = 'Jacalin-VWA'

            case ['LRR', 'MSP', 'NBARC']:

                self.annotation = 'MSP'

            case ['CC', 'Kinase', 'LRR', 'NBARC', 'Winged Helix']:

                self.annotation = 'NLR'

            case ['CC', 'MSP']:

                self.annotation = 'MSP'

            case ['Kinase', 'TM']:

                self.annotation = 'Kinase'


            case ['C2', 'MSP']:

                self.annotation = 'C2'

            case ['CC', 'MSP', 'NBARC', 'Winged Helix']:

                self.annotation = 'MSP'

            case ['CC', 'Kinase', 'NBARC', 'Winged Helix']:

                self.annotation = 'NLR'

            case ['Callose Synthase', 'MSP']:

                self.annotation = 'MSP-Callose Synthase'

            case ['MSP', 'WD40', 'zinc_finger']:

                self.annotation = 'MSP-WD40'

            case ['CC', 'LRR', 'Winged Helix']:

                self.annotation = 'NLR'

            case ['Ankyrin', 'TM']:

                self.annotation = 'Ankyrin'

            case ['CHASE', 'TM']:

                self.annotation = 'CHASE'

            case ['MSP', 'Ubiquitin']:

                self.annotation = 'MSP-Ubiquitin'

            case ['MSP', 'PGG', 'TM']:

                self.annotation = 'MSP-PGG'

            case ['AUX/IAA', 'MSP', 'Ubiquitin']:
                self.annotation = 'MSP-AUX/IAA-Ubiquitin'

            case ['CBS', 'MSP']:
                self.annotation = 'CBS-MSP'

            case ['Ankyrin', 'CC', 'CG1', 'LRR', 'MSP', 'NBARC', 'Winged Helix']:
                self.annotation = 'NLR'
            
            case _:
                self.annotation = "Missing!"


args = sys.argv[1:]

annotations = args[0]
primary_modules = args[1]
basename = args[2]

output_annotation = f"{basename}_KFP_auxiliary_modules.itol"
output_table = f"{basename}_KFP_auxiliary_modules.tsv"
output_count = f"{basename}_auxiliary_module_counts.tsv"


with open(primary_modules, "rt") as file:
    primary_module_data = dict()

    for line in file:
        data = line.strip().split("\t")

        primary_module_data[data[0]] = (int(data[1]), int(data[2]))

df_annotations = pd.read_csv(annotations, sep="\t", header=None)

df_annotations.columns = [
    "target",
    "md5",
    "sequence_len",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_annotations",
    "pathway",
]

# Only Pfam Gene3D and THMM are used to annotate auxiliary modules
df_annotations = df_annotations[df_annotations.analysis.isin(["Pfam", "Gene3D", "TMHMM"])]

# Remove kinase annotations from Gene3D for simplicity
df_annotations = df_annotations[~df_annotations.signature_accession.isin(["G3DSA:1.10.510.10", "G3DSA:3.30.200.20"])]

# Replace interpro acession with signature acession if no interpro accession is available
df_annotations.interpro_accession[df_annotations.interpro_accession == "-"] = df_annotations.signature_accession

aux_info = []

for index, df_protein in df_annotations.groupby("target"):
    primary_start = primary_module_data[index][0]
    primary_end = primary_module_data[index][1]

    aux = df_protein[((df_protein.stop <= primary_start) | (df_protein.start >= primary_end))]

    kfp = Sorter(index, aux)
    kfp.get_aux()
    kfp.get_annotation()
    aux_info.append(kfp)

annotated_auxs = dict()

for kfp in aux_info:

    annotation = kfp.annotation

    if annotation in annotated_auxs.keys():
        annotated_auxs[annotation] += 1
        continue
    
    annotated_auxs[annotation] = 1


df = pd.DataFrame(annotated_auxs, index=[0]).transpose()
df.columns = ['counts']
df.sort_values(by="counts", ascending=False, inplace=True)
df.to_csv(output_count, header=None, index=True, columns=['counts'])

with open(output_table, "w+") as out:
	for kfp in aux_info:
		out.write(f"{kfp.name}\t{kfp.annotation}\n")


