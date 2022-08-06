import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import upsetplot

def plot_results(variant, filter, region, df)
    # collect all variants in the region
    variant = df.index[variant]
    if region == "all-regions":
        ids_region = variant
    elif region == "repeat-regions":
#					ids_region = set([i for i in info[0] if (df.loc[i, 'ucsc_overlaps'] > 0)])
        ids_region = set(id for id in df[df["ucsc-repeats_overlaps"]>=0.5].variant_id if id in variant)
    elif region == "nonrepeat-regions":
#					ids_region = set([i for i in info[0] if (df.loc[i, 'ucsc_overlaps'] <= 0)])
        ids_region = set(id for id in df[df["ucsc-repeats_overlaps"]<0.5].variant_id if id in variant)

    # create plots
    with PdfPages(outname + '_' + variant_category + '_' + filter + '_' + region + '.pdf') as pdf:
        ids_autosomes = set([i for i in ids_region if not 'chrX' in i])
    
        if filter == 'unfiltered':
            df_sub = df[df.variant_id.isin(ids_region)]
            df_sub_autosomes = df[df.variant_id.isin(ids_autosomes)]
        elif filter.startswith('lenient'):
            score_column = 'score_' + variant_category
            plt.figure()
            fig, ax = plt.subplots()
            df[df.variant_id.isin(variant) & df.likely_false].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='red')
            df[df.variant_id.isin(variant) & df.likely_true].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='blue')
            df[df.variant_id.isin(variant) & ~df.likely_false & ~df.likely_true & ~df.ac0_fail].hist(ax=ax, column=score_column, bins=64, bottom = 0.1, alpha=0.5, color='grey')
            ax.set_yscale('log')
            pdf.savefig()
            plt.close()

            score_cutoff = float(filter.split('_')[1])
            df_sub = df[df.variant_id.isin(ids_region) & ( df.likely_true | ((~df.ac0_fail) & (df[score_column] > score_cutoff)) ) ]
            df_sub_autosomes = df[df.variant_id.isin(ids_autosomes) & ( df.likely_true | ((~df.ac0_fail) & (df[score_column] > score_cutoff)) ) ]
        else:
            assert(filter == 'strict')
            df_sub = df[df.variant_id.isin(ids_region) & df.likely_true]
            df_sub_autosomes = df[df.variant_id.isin(ids_autosomes) & df.likely_true]

        print('  variant count ' + filter + ' ' + region + ' ' + variant_category + ':', len(df_sub))

        # create upset plot with filters
        filter_counts = df_sub.groupby(by=['ac0_fail','mendel_fail','gq_fail','nonref_fail', 'self_fail']).size()
        plt.figure()
        upsetplot.plot(filter_counts, sort_by='cardinality')
        pdf.savefig()
        plt.close()


        for metric in ['pangenie_mendelian_consistency', 'pangenie-unrelated_allele_freq']:
            plt.figure()
            fig, ax = plt.subplots()
            df_sub.hist(ax=ax, column=metric, bins=64, bottom = 0.1)
            ax.set_yscale('log')
            fig.suptitle('min={}, max={}, mean={}, median={}'.format(df_sub[metric].min(),df_sub[metric].max(),df_sub[metric].mean(),df_sub[metric].median()) )
            pdf.savefig()
            plt.close()

        # plot panel allele frequencies for variants typed with AF=0
        df_absent = df_sub[df_sub['pangenie-all_allele_freq'] == 0.0]
        print(' typed AF=0 (allowing missing genotypes) ' + filter + ' ' + region + ' ' + variant_category + ':', len(df_absent))
        print(' typed 0/0 in all samples ' + filter + ' ' + region + ' ' + variant_category + ':', len(df_absent[df_absent['pangenie-all_unknown_alleles'] == 0]))
        plt.figure()
        fig, ax = plt.subplots()
        df_absent.hist(ax=ax, column='panel_allele_freq', bins=64, bottom = 0.1)
        ax.set_yscale('log')
        fig.suptitle('panel allele freq. for variants typed with AF=0')
#					fig.suptitle('min={}, max={}, mean={}, median={}'.format(df_absent['panel_allele_freq'].min(),df_absent['panel_allele_freq'].max(),df_absent['panel_allele_freq'].mean(),df_absent['panel_allele_freq'].median()) )
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
#							cmap = cm.get_cmap('Greys', 6) 
#							joint_kws=dict(gridsize=35, cmap=cmap)
                    joint_kws=dict(gridsize=35, cmap="hot_r")
                    ax = sns.jointplot(x=x_values, y=y_values, xlim=(-0.05, 1.05), ylim=(-0.05,1.05), bins='log', kind='hex', joint_kws=joint_kws, marginal_ticks=True, color="red")
                    ax.set_axis_labels(metrics[i], metrics[j])
                    if 'allele_freq' in metrics[i] and 'allele_freq' in metrics[j]:
                        pearson_corr, p_value = pearsonr(x_values, y_values)
                        ax.fig.suptitle(variant_category + " (r=" + str(pearson_corr) + ")")
                        print('  pearson correlation ' + variant_category + ' ' + filter + ' ' + region + ':', metrics[i], metrics[j], pearson_corr)
                    if 'pangenie-unrelated_allele_freq' in [metrics[i], metrics[j]] and 'pangenie-unrelated_heterozygosity' in [metrics[i], metrics[j]]:
                        # plot theoretical line
                        t = np.arange(0.0, 1.01, 0.01)
                        s = [2*j*(1.0-j) for j in t]
                        ax.fig.suptitle(variant_category)
                        ax.ax_joint.plot(t,s,'r-')
                    plt.colorbar()
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()

        # boxplots
        pop_metrics = ['pangenie-unrelated_allele_freq', 'panel_allele_freq']
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


        # plot distribution of AF=0 variants along the chromosomes
        # do this for large variants only
        if (filter == 'unfiltered') and ('large' in variant_category):
            chromosome_to_total = defaultdict(list)
            chromosome_to_absent = defaultdict(list)

            for var_id in df_sub.variant_id:
                fields = var_id.split('-')
                chromosome_to_total[fields[0]].append(int(fields[1]))
            for var_id in df_absent.variant_id:
                fields = var_id.split('-')
                chromosome_to_absent[fields[0]].append(int(fields[1]))

            # create histogram for each chromosome
            for chrom in sorted(chromosome_to_total.keys()):
                fig = plt.figure(figsize=(80,10))
                ax = plt.axes()
                fig.suptitle(chrom + ' ' + variant_category)

                _, bins, _ = plt.hist(chromosome_to_total[chrom], bins=500, alpha=0.5, label="total alleles")
                _ = plt.hist(chromosome_to_absent[chrom], bins=bins, alpha=0.5, label="AF=0 alleles")

                plt.xlabel('genomic position (bp)')
                plt.ylabel('allele count')
#							plt.yscale('log', nonposy='clip')
                plt.legend(loc='upper right')
                plt.tight_layout()
                pdf.savefig()
                plt.close()
