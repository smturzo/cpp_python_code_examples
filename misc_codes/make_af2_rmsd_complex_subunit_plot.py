import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def violin_plot_r1_r5_complex(path_csv_file):
	complex_rmsd = pd.read_csv(path_csv_file,header=0)
	# Melt the DataFrame so that the RMSD_R0, RMSD_R1, RMSD_R2, RMSD_R3, and RMSD_R4 columns are in a single column
	df_melted = pd.melt(complex_rmsd, id_vars=['PDBID_CHAINID'], value_vars=['RMSD_R0', 'RMSD_R1', 'RMSD_R2', 'RMSD_R3', 'RMSD_R4'], var_name='RMSD', value_name='Value')
	print(df_melted)
	plt.figure()
	sns.violinplot(x='RMSD', y='Value', data=df_melted, split=True,cut=0)
	plt.xlabel("AF Ranks")
	plt.ylabel("RMSD")
	plt.tight_layout()
	output_file_name = path_csv_file.split(".csv")[0]
	plt.savefig(f"{output_file_name}_plot.png",dpi=300)
	plt.close()

def violin_plot_r1_r5_subunit(path_csv_file):
	# Read the CSV file into a pandas dataframe
	df = pd.read_csv(path_csv_file)

	# Split the dataframe into two based on chain id
	df_A = df[df['CHAINID'] == 'A']
	df_B = df[df['CHAINID'] == 'B']
	
	# Set up the plot
	plt.figure()
	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

	# Plot the violin plot for chain id A
	sns.violinplot(data=df_A[['RMSD_R0', 'RMSD_R1', 'RMSD_R2', 'RMSD_R3', 'RMSD_R4']], ax=axs[0], split=True,cut=0)
	axs[0].set_title('Chain A')

	# Plot the violin plot for chain id B
	sns.violinplot(data=df_B[['RMSD_R0', 'RMSD_R1', 'RMSD_R2', 'RMSD_R3', 'RMSD_R4']], ax=axs[1], split=True,cut=0)
	axs[1].set_title('Chain B')

	# Add a common y-axis label
	axs[0].set_ylabel("RMSD")
	axs[1].set_ylabel("RMSD")
	#fig.text(0.04, 0.5, 'RMSD', ha='center', va='center', rotation='vertical')

	
	# Save the plot
	output_file_name = path_csv_file.split(".csv")[0]
	plt.tight_layout()
	plt.savefig(f"{output_file_name}_plot.png",dpi=300)
	plt.close()

#violin_plot_r1_r5_complex("../heteromers/V222_full_complex_rmsd.csv")
#violin_plot_r1_r5_complex("../heteromers/V231_full_complex_rmsd.csv")
#violin_plot_r1_r5_subunit("../heteromers/V222_single_chain_subunit_rmsd.csv")
#violin_plot_r1_r5_subunit("../heteromers/V231_single_chain_subunit_rmsd.csv")

violin_plot_r1_r5_complex("../homomers/V222_full_complex_rmsd.csv")
violin_plot_r1_r5_complex("../homomers/V231_full_complex_rmsd.csv")
violin_plot_r1_r5_complex("../homomers/V222_single_chain_subunit_rmsd.csv")
violin_plot_r1_r5_complex("../homomers/V231_single_chain_subunit_rmsd.csv")






