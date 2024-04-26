import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Plot(object):

    def __init__(self, args):
        self.args = args
        self.out = args.out
        self.snpEff = args.snpEff
        self.fig_width = args.fig_width
        self.fig_height = args.fig_height
        self.white_space = args.white_space
        self.plot_with_indel = args.indel
        try:
            import seaborn as sns
            sns.set_style('whitegrid')
        except ModuleNotFoundError:
            pass

    def read_files(self):
        if self.snpEff is None:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                    sep='\t',
                                    names=['CHROM',
                                           'POSI',
                                           'variant',
                                           'depth',
                                           'p99',
                                           'p95',
                                           'SNPindex'])
        else:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                    sep='\t',
                                    names=['CHROM',
                                           'POSI',
                                           'variant',
                                           'impact',
                                           'depth',
                                           'p99',
                                           'p95',
                                           'SNPindex'])

        sliding_window = pd.read_csv('{}/sliding_window.tsv'.format(self.out),
                                     sep='\t',
                                     names=['CHROM',
                                            'POSI',
                                            'mean_p99',
                                            'mean_p95',
                                            'mean_SNPindex'])

        snp_index['CHROM'] = snp_index['CHROM'].astype('str')
        snp_index['POSI'] = snp_index['POSI']/1000000
        sliding_window['CHROM'] = sliding_window['CHROM'].astype('str')
        sliding_window['POSI'] = sliding_window['POSI']/1000000
        sliding_window = sliding_window.groupby('CHROM')

        return snp_index, sliding_window

    def get_significant_contigs(self, N_chr, snp_index, sliding_window):

        print(('!!WARNING!! Your reference genome has too many contigs (>50). '
               'Therefore, only significant contigs will be plotted.'), file=sys.stderr)
        sliding_window=sliding_window.obj
        significant_windows = sliding_window[abs(sliding_window['mean_p95']) <= \
                                             abs(sliding_window['mean_SNPindex'])]

        significant_contigs = list(significant_windows['CHROM'].drop_duplicates())

        N_chr = len(significant_contigs)
        snp_index = snp_index[snp_index['CHROM'].isin(significant_contigs)]
        sliding_window = sliding_window[sliding_window['CHROM'].isin(significant_contigs)]
        sliding_window = sliding_window.groupby('CHROM')
        return N_chr, snp_index, sliding_window

    def set_plot_style(self, N_chr):
        if N_chr == 1:
            style = (1, 1)
        elif N_chr == 2:
            style = (1, 2)
        else:
            if N_chr % 2 == 0:
                style = (2, int(N_chr/2))
            else:
                style = (2, int(N_chr/2) + 1)
        return style

    def run(self):
        snp_index, sliding_window = self.read_files()

        N_chr = len(sliding_window['CHROM'].unique())

        if N_chr > 50:
            N_chr, self.snp_index, sliding_window = self.get_significant_contigs(N_chr, 
                                                                                      snp_index, 
                                                                                      sliding_window)

        N_col, N_raw = self.set_plot_style(N_chr)

        fig = plt.figure(figsize=(self.fig_width*N_col, self.fig_height*N_raw))
        plt.subplots_adjust(hspace=self.white_space)

        xmax = snp_index['POSI'].max()

        for i, (chr_name, chr_sliding_window) in enumerate(sliding_window):
            chr_snp_index = snp_index[snp_index['CHROM']==chr_name]
            if not self.plot_with_indel:
                chr_snp_index = chr_snp_index[chr_snp_index['variant']=='snp']

            ax = fig.add_subplot(N_raw, N_col, i+1)
            ax.plot(chr_sliding_window['POSI'],
                    chr_sliding_window['mean_p99'],
                    color='orange',
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'],
                    chr_sliding_window['mean_p95'],
                    color='lime',
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'],
                    chr_sliding_window['mean_SNPindex'],
                    color='red',
                    linewidth=3)

            if self.snpEff is None:
                ax.scatter(chr_snp_index['POSI'],
                           chr_snp_index['SNPindex'],
                           marker='.',
                           color='navy')
            else:
                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODIFIER']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODIFIER']['SNPindex'],
                           marker='.',
                           color='navy')

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='LOW']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='LOW']['SNPindex'],
                           marker='.',
                           color='navy')

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODERATE']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODERATE']['SNPindex'],
                           marker='.',
                           color='navy')

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='HIGH']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='HIGH']['SNPindex'],
                           marker='x',
                           color='firebrick')

            ax.hlines([0.5], 0, xmax, linestyles='dashed')
            ax.set_xlabel('position (Mb)', fontsize=15)
            ax.set_ylabel('SNP-index', fontsize=15)
            ax.set_xlim(0, xmax)
            ax.set_ylim(0, 1.05)
            ax.set_title(chr_name, fontsize=17)

        plt.savefig('{}/mutmap_plot.{}'.format(self.out, self.args.format))
