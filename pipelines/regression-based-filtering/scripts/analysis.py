#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
from pylab import *

import sklearn
from sklearn.experimental import enable_iterative_imputer

SV_min_size = 50
small_indel_max_size = 20

# sklearn settings
seed = 42
verbose = False
impute_iterations = 10
Epsilon_Support_Vector_Regression_regularization_parameter = 50
Epsilon_Support_Vector_Regression_epsilon = 0.1

# sklearn generic processing objects
imputor = sklearn.impute.IterativeImputer(max_iter=impute_iterations, verbose=verbose, random_state=seed)
scaler = sklearn.preprocessing.StandardScaler()
regressor = sklearn.svm.SVR(C=Epsilon_Support_Vector_Regression_regularization_parameter 
							epsilon=Epsilon_Support_Vector_Regression_epsilon,
							verbose=verbose)
#regressor = GradientBoostingRegressor(n_estimators=20, random_state=0)

if __name__ == "__main__":
	filename = sys.argv[1]
	outname = sys.argv[2]
#	med_svs = sys.stgv[3]

	df = pd.read_csv(filename, sep='\t')
	d['pan']
	df = df.assign(pangenie_mendelian_consistency=lambda d: d['pangenie_mendelian_consistent_trios'] / d['pangenie_considered_trios'])

	# consider allele frequency computed across all genotyped samples
	df = df.assign(ac0_fail = df['pangenie-all_allele_freq'] == 0,
				   mendel_fail = (df.pangenie_mendelian_consistency < 0.8) & (df['pangenie_considered_trios']>=5),
	               gq_fail = df['pangenie-all_GQ>=200'] < 50,
	               self_fail = df['pangenie_self-genotyping_correct [%]'] < 90.0,
	               nonref_fail = ((df[['pangenie_self-genotyping_0/1_typed_0/1',
				                       'pangenie_self-genotyping_1/1_typed_1/1']] == 0).all(axis=1) & 
				                  (df[['pangenie_self-genotyping_0/1_typed_0/0', 
								       'pangenie_self-genotyping_0/1_typed_1/1',
								       'pangenie_self-genotyping_1/1_typed_0/0',
									   'pangenie_self-genotyping_1/1_typed_0/1']] != 0).any(axis=1)),
	               likely_true = ~(df.ac0_fail | df.mendel_fail | df.gq_fail | df.nonref_fail | df.self_fail),
	               likely_false = ~(df.ac0_fail & ((df.gq_fail + df.mendel_fail + df.nonref_fail + df.self_fail) >= 3))
				   test_set = ~df.ac0_fail
				   )
	df = df.assign(training_set_category = 0)
	df.loc[likely_true, 'training_set_category'] = 1
	df.loc[likely_false, 'training_set_category'] = -1

	snps = df.variant_id.str.contains('SNV') #set(id for id in df.variant_id if 'SNV' in id)
	indel_size = df.variant_id.str.split('-').str[-1].astype(int)
	svs = ~snps & (indel_size >= SV_min_size)
	indels = ~snps & ~svs
	small_indels = ~snps & ~svs & (indel_size < small_indel_max_size) #set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])<=19))
	midsize_indels = ~snps & ~svs & ~small_indels #set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])>=small_indel_max_size) and (int(id.split('-')[-1])<SV_min_size))
	
	insertions = df.variant_id.str.split('-')[2] == "INS"
	deletions = df.variant_id.str.split('-')[2] == "DEL"
	complex_indels = df.variant_id.str.split('-')[2] == "COMPLEX"

	small_insertions = small_indels & insertions
	small_deletions = small_indels & deletions
	small_complex = small_indels & complex_indels
	assert sum(small_insertions) + sum(small_deletions) + sum(small_complex) == sum(small_indels)

	midsize_insertions = midsize_indels & insertions
	midsize_deletions = midsize_indels & deletions
	midsize_complex = midsize_indels & complex_indels
	assert sum(midsize_insertions) + sum(midsize_deletions) + sum(midsize_complex) == sum(midsize_indels)

	large_insertions = svs & insertions
	large_deletions = svs & deletions
	large_complex = svs & complex_indels
	assert sum(large_insertions) + sum(large_deletions) + sum(large_complex) == sum(svs)

	df.assign(mutation_category = '', mutation_size = '')
	df.loc[snps, 'mutation_category'] = 'SNP'
	df.loc[insertions, 'mutation_category'] = 'insertion'
	df.loc[deletions, 'mutation_category'] = 'deletion'
	df.loc[complex_indels, 'mutation_category'] = 'complex_indel'
	df.loc[snps, 'mutation_size'] = 'SNP'
	df.loc[svs, 'mutation_size'] = 'SV'
	df.loc[small_indels, 'mutation_size'] = 'small_indel'
	df.loc[midsize_indels, 'mutation_size'] = 'midsize_indel'


#	med = set()
#	for line in open(med_svs, 'r'):
#		if not 'allele-' in line:
#			med.add(line.strip())

	# for pangenie, consider allele frequencies computed across unrelated samples (+ panel samples) only, and skip children (unless they are panel samples)
	metrics = ['panel_allele_freq', 'pangenie-unrelated_allele_freq', 'pangenie-unrelated_heterozygosity']

	# features used for regression
	regression_features = [
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
		'pangenie-all_alternative_alleles',
		'pangenie-all_total_alleles',
		'pangenie-all_heterozygosity',
		'pangenie-all_heterozygous_genotypes',
		'pangenie-all_total_genotypes',
		'pangenie-all_unique_kmers',
		'pangenie-all_GQ>=200',
		'pangenie-all_allele_freq',
		'pangenie_mendelian_consistent_trios',
		'pangenie_alternative_transmitted',
		'pangenie_considered_trios',
	]

	variant_types = {
		'snps': [snps, ['unfiltered', 'strict'], ['all-regions']],
#		'small_insertions': [small_insertions, ['unfiltered', 'strict'], ['all-regions']],
#		'small_deletions': [small_deletions, ['unfiltered', 'strict'], ['all-regions']],
#		'small_complex': [small_complex, ['unfiltered', 'strict'], ['all-regions']],
#		'midsize_deletions': [midsize_deletions, ['unfiltered', 'strict'], ['all-regions']],
#		'midsize_insertions': [midsize_insertions, ['unfiltered', 'strict'], ['all-regions']],
#		'midsize_complex': [midsize_complex, ['unfiltered', 'strict'], ['all-regions']],
		'indels': [indels, ['unfiltered', 'strict'], ['all-regions']],
		'large_deletions': [large_deletions, ['unfiltered', 'lenient_-0.5', 'strict'], ['all-regions']],
		'large_insertions': [large_insertions, ['unfiltered', 'lenient_-0.5', 'strict'], ['all-regions']],
		'large_complex': [large_complex, ['unfiltered', 'lenient_-0.5', 'strict'], ['all-regions']]
#		'giab_med_svs': [med, ['unfiltered', 'lenient_-0.5', 'strict']]
	}

	df['score_SVR'] = np.nan

	df = df.sort_values(by='variant_id').set_index('variant_id')

	for variant_category, (variant_category_mask, categories, regions) in variant_types.items():
		if 'large' in variant_category:
			variants_to_test = variant_category_mask & test_set
			variant_category_subset = df.loc[variants_to_test, regression_features].copy()
			# fill in any nans of features we will be regressing on
			imputed_regression_features = imputor.fit_transform(variant_category_subset)
			imputed_df = pd.DataFrame(imputed_regression_features,
							          columns=regression_features,
									  index=df.index[variants_to_test])
			variant_category_subset.update(imputed_df, overwrite=True) # Why would this be True? Wouldn't you want to regain original values?

			# Fit transform using all data points (labeled + unlabeled)
			scaler.fit_transform(variant_category_subset[regression_features].values)

			# Train model only on labeled points
			autosomal_variants = ~variant_category_subset.index.str.contains('chrX')
			autosomal_training_set = autosomal_variants & (~test_set)
			x = scaler.transform(variant_category_subset.loc[autosomal_training_set, regression_features].values)
			y = variant_category_subset.loc[autosomal_training_set, ['training_set_category']].values
			print('Training regression model')
			regressor.fit(x,y.ravel())

			# Apply to unlabeled data
			y_pred = regressor.predict(scaler.transform(variant_category_subset.loc[:,features].values))
			column_label = "regression_score"
			df_score = pd.DataFrame({"variant_id":variant_category_subset.index, column_label:y_pred})

			# Add column with variant specific scores to table
			df = df.merge(df_score, on='variant_id', how='left')

		for filter in categories:
			variant = df.index[variant_category_mask]
			for region in regions:
				plot_results(variant, filter, region, df)


	for variant in variant_types.keys():
		if 'large' in variant:
			column_name = 'score_' + variant
			df['score_SVR'].fillna(df[column_name], inplace=True)
	
	df['confidence_level'] = 0
	df.confidence_level.where( svs & (df.regression_score <  0.5), 3, inplace=True )
	df.confidence_level.where( svs & (df.regression_score <  0.0), 2, inplace=True )
	df.confidence_level.where( svs & (df.regression_score < -0.5), 1, inplace=True )
	df.confidence_level.where(df.likely_true, 4, inplace=True)
	
	header = ["variant_id", "ac0_fail", "mendel_fail", "self_fail", "gq_fail", "nonref_fail", 
		      "score_SVR", "likely_true", "regression_score", "confidence_level", 
			  "mutation_category", "mutation_size"]
	df.to_csv(outname + '_filters.tsv', columns=header, sep='\t', index=False, na_rep='nan')

