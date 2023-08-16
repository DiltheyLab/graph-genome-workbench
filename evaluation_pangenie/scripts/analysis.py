#!/usr/bin/env python
import sys
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import upsetplot
from scipy.stats import pearsonr
from pylab import *

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import AdaBoostClassifier
from sklearn.svm import SVR
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer



if __name__ == "__main__":
	filename = sys.argv[1]
	outname = sys.argv[2]
	med_svs = sys.argv[3]

	df = pd.read_csv(filename, sep='\t')
	df = df.assign(pangenie_mendelian_consistency=lambda d: d['pangenie_mendelian_consistent_trios'] / d['pangenie_considered_trios'])

	df = df.assign(ac0_fail = lambda df: df['pangenie_allele_freq'] == 0)
	df = df.assign(mendel_fail = lambda df: (df.pangenie_mendelian_consistency < 0.9) & (df['pangenie_considered_trios']>=5))
	df = df.assign(gq_fail = lambda df: df['pangenie_GQ>=200'] < 5)
	df = df.assign(self_fail = lambda df: (df['pangenie_self-genotyping_correct [%]'] < 90.0) )
	df = df.assign(nonref_fail = lambda df: (df['pangenie_self-genotyping_0/1_typed_0/1'] + df['pangenie_self-genotyping_1/1_typed_1/1'])==0)
	df = df.assign(all_pass = lambda df: ~(df.ac0_fail | df.mendel_fail | df.gq_fail | df.nonref_fail | df.self_fail) )
	df = df.assign(negative = lambda df: ~df.ac0_fail & (((1*df.gq_fail) + (1*df.mendel_fail) + (1*df.nonref_fail) + (1*df.self_fail)) >= 3))

	snps = set(id for id in df.variant_id if 'SNV' in id)
	indels = set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])<50))
	small_indels = set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])<=19))
	midsize_indels = set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])>19) and (int(id.split('-')[-1])<50))
	svs = set(id for id in df.variant_id if (not 'snp' in id) and (int(id.split('-')[-1])>=50))

	small_insertions = set(id for id in small_indels if id.split('-')[2]=="INS")
	small_deletions = set(id for id in small_indels if id.split('-')[2]=="DEL")
	midsize_insertions = set(id for id in midsize_indels if id.split('-')[2]=="INS")
	midsize_deletions = set(id for id in midsize_indels if id.split('-')[2]=="DEL")
	large_insertions = set(id for id in svs if id.split('-')[2]=="INS")
	large_deletions = set(id for id in svs if id.split('-')[2]=="DEL")
	med = set()

	for line in open(med_svs, 'r'):
		if not 'allele-' in line:
			med.add(line.strip())

	metrics = ['panel_allele_freq', 'pangenie_allele_freq', 'pangenie_heterozygosity']

	# features used for regression
	features = [
		'pangenie_self-genotyping_correct [%]',
		'pangenie_self-genotyping_wrong [%]',
		'pangenie_self-genotyping_not_typed [%]',
		'pangenie_self-genotyping_correct',
		'pangenie_self-genotyping_wrong',
		'pangenie_self-genotyping_not_typed',
		'pangenie_self-genotyping_absent_in_truth',
		'pangenie_self-genotyping_0/0_typed_0/0',
		'pangenie_self-genotyping_0/0_typed_0/1',
		'pangenie_self-genotyping_0/0_typed_1/1',
		'pangenie_self-genotyping_0/0_typed_./.',
		'pangenie_self-genotyping_0/1_typed_0/0',
		'pangenie_self-genotyping_0/1_typed_0/1',
		'pangenie_self-genotyping_0/1_typed_1/1',
		'pangenie_self-genotyping_0/1_typed_./.',
		'pangenie_self-genotyping_1/1_typed_0/0',
		'pangenie_self-genotyping_1/1_typed_0/1',
		'pangenie_self-genotyping_1/1_typed_1/1',
		'pangenie_self-genotyping_1/1_typed_./.',
		'panel_allele_freq',
		'panel_alternative_alleles',
		'panel_total_alleles',
		'pangenie_alternative_alleles',
		'pangenie_total_alleles',
		'pangenie_heterozygosity',
		'pangenie_heterozygous_genotypes',
		'pangenie_total_genotypes',
		'pangenie_unique_kmers',
		'pangenie_GQ>=200',
		'pangenie_allele_freq',
		'pangenie_mendelian_consistent_trios',
		'pangenie_alternative_transmitted',
		'pangenie_considered_trios',
	]

	variant_types = {
		'snps': [snps, ['unfiltered', 'strict']],
		'small_insertions': [small_insertions, ['unfiltered', 'strict']],
		'small_deletions': [small_deletions, ['unfiltered', 'strict']],
		'midsize_deletions': [midsize_deletions, ['unfiltered', 'strict']],
		'midsize_insertions': [midsize_insertions, ['unfiltered', 'strict']],
		'large_deletions': [large_deletions, ['unfiltered', 'lenient_-0.5', 'strict']],
		'large_insertions': [large_insertions, ['unfiltered', 'lenient_-0.5', 'strict']],
		'giab_med_svs': [med, ['unfiltered', 'lenient_-0.5', 'strict']]
	}

	df['score_SVR'] = np.nan

	for name, info in variant_types.items():
		if 'large' in name:
			df_sub = df[df.variant_id.isin(info[0]) & ~df.ac0_fail].copy()
			df_sub.sort_values(by=['variant_id'], inplace=True)
			df_sub.set_index('variant_id')

			# impute missing values
			imp = IterativeImputer(max_iter=10, verbose=0, random_state=0)
			imp.fit(df_sub[features])
			imputed_df = pd.DataFrame(imp.transform(df_sub[features]), columns=df_sub[features].columns, index=df_sub.index)
			df_sub.update(imputed_df, overwrite=True)

			# Fit transform using all data points (labeled + unlabeled)
			df_sub = df_sub.assign(target = lambda df: (-1*df.negative) + (1*df.all_pass))
			x = df_sub.loc[:,features].values
			scaler = StandardScaler()
			scaler.fit_transform(x)

			# Train model only on labeled points
			autosomal_ids = [i for i in info[0] if not 'chrX' in i]
			df_labeled = df_sub[df.variant_id.isin(autosomal_ids) & (df_sub.target != 0)]
			x = scaler.transform(df_labeled.loc[:,features].values)
			y = df_labeled.loc[:,['target']].values
			print('Training regression model')
			regressor = SVR(C=50)
#			regressor = GradientBoostingRegressor(n_estimators=20, random_state=0)
			regressor.fit(x,y.ravel())

			# Apply to unlabeled data
			y_pred = regressor.predict(scaler.transform(df_sub.loc[:,features].values))
			column_label = "score_" + name
			df_score = pd.DataFrame({"variant_id":df_sub.variant_id, column_label:y_pred})

			# Add column with variant specific scores to table
			df = df.merge(df_score, on='variant_id', how='left')

		for filter in info[1]:
			variant = info[0]
			with PdfPages(outname + '_' + name + '_' + filter + '.pdf') as pdf:
				ids_autosomes = set([i for i in variant if not 'chrX' in i])
			
				if filter == 'unfiltered':
					df_sub = df[df.variant_id.isin(variant)]
					df_sub_autosomes = df[df.variant_id.isin(ids_autosomes)]
				elif filter.startswith('lenient'):
					score_column = 'score_' + name
					plt.figure()
					fig, ax = plt.subplots()
					df[df.variant_id.isin(variant) & df.negative].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='red')
					df[df.variant_id.isin(variant) & df.all_pass].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='blue')
					df[df.variant_id.isin(variant) & ~df.negative & ~df.all_pass & ~df.ac0_fail].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='grey')
					ax.set_yscale('log')
					pdf.savefig()
					plt.close()
	
					score_cutoff = float(filter.split('_')[1])
					df_sub = df[df.variant_id.isin(variant) & ( df.all_pass | ((~df.ac0_fail) & (df[score_column] > score_cutoff)) ) ]
					df_sub_autosomes = df[df.variant_id.isin(ids_autosomes) & ( df.all_pass | ((~df.ac0_fail) & (df[score_column] > score_cutoff)) ) ]
				else:
					assert(filter == 'strict')
					df_sub = df[df.variant_id.isin(variant) & df.all_pass]
					df_sub_autosomes = df[df.variant_id.isin(ids_autosomes) & df.all_pass]

				print('  variant count ' + filter + ' ' + name + ':', len(df_sub))

				# create upset plot with filters
				filter_counts = df_sub.groupby(by=['ac0_fail','mendel_fail','gq_fail','nonref_fail', 'self_fail']).size()
				plt.figure()
				upsetplot.plot(filter_counts, sort_by='cardinality')
				pdf.savefig()
				plt.close()

				plt.figure()
				fig, ax = plt.subplots()
				metric = 'pangenie_mendelian_consistency'
				df_sub.hist(ax=ax, column=metric, bins=64, bottom = 0.1)
				ax.set_yscale('log')
				fig.suptitle('min={}, max={}, mean={}, median={}'.format(df_sub[metric].min(),df_sub[metric].max(),df_sub[metric].mean(),df_sub[metric].median()) )
				pdf.savefig()
				plt.close()

				# heatmaps
				for i in range(len(metrics)):
					for j in range(i+1, len(metrics)):
							plt.figure()
							fig, ax = plt.subplots()
							x_values = []
							y_values = []
							for l,f in zip(df_sub_autosomes[metrics[i]], df_sub_autosomes[metrics[j]]):
								if not pd.isnull(f) and not pd.isnull(l):
									x_values.append(l)
									y_values.append(f)
							assert len(x_values) == len(y_values)
							if len(x_values) == 0:
								continue
							cmap = cm.get_cmap('Greys', 6) 
							joint_kws=dict(gridsize=35, cmap=cmap)
							ax = sns.jointplot(x=x_values, y=y_values, xlim=(-0.05, 1.05), ylim=(-0.05,1.05), bins='log', kind='hex', joint_kws=joint_kws, marginal_ticks=True, color="grey")
							ax.set_axis_labels(metrics[i], metrics[j])
							if 'allele_freq' in metrics[i] and 'allele_freq' in metrics[j]:
								pearson_corr, p_value = pearsonr(x_values, y_values)
								ax.fig.suptitle(name + " (r=" + str(pearson_corr) + ")")
								print('  pearson correlation ' + name + ':', metrics[i], metrics[j], pearson_corr)
							if 'pangenie_allele_freq' in [metrics[i], metrics[j]] and 'pangenie_heterozygosity' in [metrics[i], metrics[j]]:
								# plot theoretical line
								t = np.arange(0.0, 1.01, 0.01)
								s = [2*j*(1.0-j) for j in t]
								ax.fig.suptitle(name)
								ax.ax_joint.plot(t,s,'r-')
							plt.colorbar()
							plt.tight_layout()
							pdf.savefig()
							plt.close()
			
				# boxplots
				pop_metrics = ['pangenie_allele_freq', 'panel_allele_freq']
				length_intervals = []
				n_intervals = 10
				previous = -0.1
				length_intervals = []
				for i in range(1, n_intervals + 1):
					length_intervals.append([previous, i/n_intervals + 0.0001])
					previous = i/n_intervals
				print(length_intervals)
				length_afs = [[] for i in range(n_intervals)]
				for l,f in zip(df_sub_autosomes[pop_metrics[0]], df_sub_autosomes[pop_metrics[1]]):
					if not pd.isnull(l):
						for i,interval in enumerate(length_intervals):
							if interval[0] < f <= interval[1]:
								length_afs[i].append(l)
								break
	
				print(sum([len(b) for b in length_afs]))
				# create boxplots
				fig = plt.figure()
				ax = plt.axes()		
				bp = ax.boxplot(length_afs)
				length_labels = []
				for i,x in enumerate(length_intervals):
					label = '<=' + str(i+1) + '/' + str(n_intervals)
					length_labels.append(label)
#				ax.set_yticks([0] + [i[1] for i in length_intervals])
#				ax.set_yticklabels(["0"] + [str(i) + '/' + str(n_intervals) for i in range(1,n_intervals+1)])
				ax.set_xticklabels(length_labels, rotation=45)
		#		ax.set_yscale("log")
				plt.tight_layout()
				pdf.savefig()
				plt.close()

	for variant in variant_types.keys():
		if 'large' in variant:
			column_name = 'score_' + variant
			df['score_SVR'].fillna(df[column_name], inplace=True)
	
	df['confidence_level'] = 0
	df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<-0.5)), 1, inplace=True )
	df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<0.0)), 2, inplace=True )
	df.confidence_level.where( (df.score_SVR.isnull() | (df.score_SVR<0.5)), 3, inplace=True )
	df.confidence_level.where(~df.all_pass, 4, inplace=True)

	df = df.assign(is_snp = lambda df: df.variant_id.isin(snps))
	df = df.assign(is_small_insertion = lambda df: df.variant_id.isin(small_insertions))
	df = df.assign(is_small_deletion = lambda df: df.variant_id.isin(small_deletions))
	df = df.assign(is_midsize_insertion = lambda df: df.variant_id.isin(midsize_insertions))
	df = df.assign(is_midsize_deletion = lambda df: df.variant_id.isin(midsize_deletions))
	df = df.assign(is_large_insertion = lambda df: df.variant_id.isin(large_insertions))
	df = df.assign(is_large_deletion = lambda df: df.variant_id.isin(large_deletions))

	header = ["variant_id", "ac0_fail", "mendel_fail", "gq_fail", "nonref_fail", "score_SVR", "all_pass", "confidence_level", "is_snp", "is_small_insertion", "is_small_deletion", "is_midsize_insertion", "is_midsize_deletion", "is_large_insertion", "is_large_deletion"]
	df.to_csv(outname + '_filters.tsv', columns=header, sep='\t', index=False, na_rep='nan')

