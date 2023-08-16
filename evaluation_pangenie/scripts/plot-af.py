import sys
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import upsetplot
from scipy.stats import pearsonr
from pylab import *

filename = sys.argv[1]
outname = sys.argv[2]

df = pd.read_csv(filename, sep='\t')

df_intersection = df[df.in_gnomAD & df.in_assemblies]
df_only_gnomad = df[df.in_gnomAD & ~(df.in_assemblies)]
df_only_assembly = df[df.in_assemblies & ~(df.in_gnomAD)]



print('Intersection: ' + str(len(df_intersection)))
print('Only gnomAD: ' + str(len(df_only_gnomad)))
print('Only assembly: ' + str(len(df_only_assembly)))


with PdfPages(outname) as pdf:
	# plot AF for only assemblies
	plt.figure()
	fig, ax = plt.subplots()
#	fig.title('allele frequency')
	weights1 = np.ones_like(df_only_assembly.allele_freq) / len(df_only_assembly)
	weights2 = np.ones_like(df_intersection.allele_freq) / len(df_intersection)
	df_only_assembly.hist(ax=ax, column="allele_freq", bins=64, weights=weights1, alpha=0.45, label="only PanGenie")
	df_intersection.hist(ax=ax, column="allele_freq", bins=64, weights=weights2, alpha=0.45, label="Intersection")
	ax.set_yscale('log')
	ax.legend()
#	ax.set_xlabel('allele frequency')
	pdf.savefig()
	plt.close()
