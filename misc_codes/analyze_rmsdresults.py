import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
homomers_2mers_rmsd_data =pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/homomers/homomers_rmsd_analysis.csv")
heteromers_2mers_rmsd_data =pd.read_csv("/fs/project/PAS1146/turzo.1/complex/ideal_complex/2mers/heteromers/heteromers_rmsd_analysis_B.csv")

def plt_rmsd_plots(type_,rmsd_data):
	subunit_v222_rmsd = rmsd_data["SubV222"]
	subunit_v231_rmsd = rmsd_data["SubV231"]
	complex_v222_rmsd = rmsd_data["CompV222"]
	complex_v231_rmsd = rmsd_data["CompV231"] 
	protein_name      = rmsd_data["Protein"]
	comp_diff = complex_v231_rmsd - complex_v222_rmsd
	subunit_diff = subunit_v231_rmsd - subunit_v222_rmsd
	plt.figure()
	barWidth = 0.25
	# Set position of bar on X axis
	br1 = np.arange(len(protein_name))
	br2 = [x + barWidth for x in br1]
	plt.bar(br1,subunit_diff, width = 0.25, edgecolor ='black', label ='Subunit')
	plt.bar(br2,comp_diff,    width = 0.25, edgecolor ='black', label ='Complex')
	# Set the xticks to the center of the bars
	plt.xticks([r + barWidth/2 for r in range(len(protein_name))], protein_name,rotation=45)
	plt.xlabel("PDB ID")
	plt.ylabel("RMSD Diff")
	plt.legend()
	plt.tight_layout()
	plt.savefig("../"+str(type_)+"/Complex_Subunit_Diff.png",dpi=300)
	plt.close()
	
	plt.figure()
	plt.scatter(subunit_v222_rmsd,complex_v222_rmsd,color='crimson',alpha=0.50, label='V222')
	plt.scatter(subunit_v231_rmsd,complex_v231_rmsd,color='black',alpha=0.70, label='V231')
	plt.legend()
	plt.xlabel("Subunit RMSD")
	plt.ylabel("Complex RMSD")
	plt.tight_layout()
	plt.savefig("../"+str(type_)+"/Complex_vs_Subunit.png",dpi=300)
	plt.close()


plt_rmsd_plots("heteromers",heteromers_2mers_rmsd_data)
