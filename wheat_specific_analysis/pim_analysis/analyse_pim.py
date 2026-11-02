import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_pim(pim) -> dict[str, list[float]]:
    with open(pim) as input:

        pim_dict: dict[str, list[float]] = {}

        n_keys = int(next(input))

        for line in input:
            data = line.split()

            key = data[0]
            value = list(map(float, data[1:]))

            pim_dict[key] = value

        if n_keys == len(pim_dict):
            print("All data loaded..")
    
    return pim_dict

def load_kfp_subset(itol, colour):
    
    kfps = []

    switch = False
    
    with open(itol) as input:
        for line in input:
            
            if switch:
            	if colour in line:
                	kfps.append(line.split(',')[0])
                
            
            if line.startswith('DATA\n'):
                switch = True
    
    return kfps
            
def generate_pim_matrix(pim):
    pim_dict = load_pim(pim)
    pim_df = pd.DataFrame(pim_dict)
    pim_df.index = pim_df.columns

    return pim_df

def plot_pim_by_aux(itol, aux_colour, PIM, output) -> None:

	aux_label = load_kfp_subset(itol, aux_colour)
	df = generate_pim_matrix(PIM)
	kinase_df = df[(df.index.isin(aux_label))]
	kinase_df = kinase_df[(aux_label)]
	ax = sns.clustermap(kinase_df, vmin=0, vmax=100, xticklabels=aux_label, yticklabels=aux_label, cmap='crest_r', method='single')
	fontsize = 8
	ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels(), fontsize = fontsize)
	ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize = fontsize)
	plt.savefig(output, dpi=300)
	

if __name__ == "__main__":
    PIM = "wheat_chinese_spring_primary_modules.pim" # from clustal Omega
    itol = "wheat_chinese_spring_auxiliary_modules.itol" #from KFP-annotator
    kinase = '#e6194b'
    msp = '#3cb44b'
    
    outdir="outdir"
    msp_out=f"{outdir}/wheat_chinese_spring_MSP_KFP_PIM.svg"
    kinase_out=f"{outdir}/wheat_chinese_spring_kinase_KFP_PIM.svg" 
    
    plot_pim_by_aux(itol, kinase, PIM, kinase_out)
    plot_pim_by_aux(itol, msp, PIM, msp_out)
	
    	
    	
