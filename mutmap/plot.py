import os
import sys
import re
import gzip
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 13


class Plot(object):

    def __init__(self, args):
        self.args = args
        self.out = args.out
        self.vcf = args.vcf
        self.snpEff = args.snpEff
        self.line_colors = args.line_colors.split(',') #SNP-index, p95, and p99
        self.dot_color = args.dot_color
        self.fig_width = args.fig_width
        self.fig_height = args.fig_height
        self.white_space = args.white_space
        self.plot_with_indel = args.indel
        try:
            import seaborn as sns
            sns.set_style('whitegrid')
        except ModuleNotFoundError:
            pass

        if len(self.line_colors) != 3 or '' in [color.strip() for color in self.line_colors]:
            raise ValueError(('Invalid "--line-colors" value: "{}"\n'
                              'Please specify exactly 3 comma-separated colors in this order:\n'
                              'SNP-index,p95,p99 (e.g. "#D55D00,#009E72,#FDB003").').format(args.line_colors))

        self.line_colors = [color.strip() for color in self.line_colors]
        self.dot_color = self.dot_color.strip()
        if self.dot_color == '':
            raise ValueError('Invalid "--dot-color" value: empty string.')

    def read_files(self):
        snp_index_path = '{}/snp_index.tsv'.format(self.out)
        sliding_window_path = '{}/sliding_window.tsv'.format(self.out)
        if not os.path.isfile(snp_index_path):
            raise FileNotFoundError(('snp_index.tsv was not found.\n'
                                     'Please check "{}".').format(snp_index_path))
        if not os.path.isfile(sliding_window_path):
            raise FileNotFoundError(('sliding_window.tsv was not found.\n'
                                     'Please check "{}".').format(sliding_window_path))

        if self.snpEff is None:
            cols = ['CHROM',
                    'POSI',
                    'variant',
                    'depth',
                    'p99',
                    'p95',
                    'SNPindex']
        else:
            cols = ['CHROM',
                    'POSI',
                    'variant',
                    'impact',
                    'depth',
                    'p99',
                    'p95',
                    'SNPindex']

        try:
            snp_index = pd.read_csv(snp_index_path, sep='\t', names=cols)
            sliding_window = pd.read_csv(sliding_window_path,
                                         sep='\t',
                                         names=['CHROM',
                                                'POSI',
                                                'mean_p99',
                                                'mean_p95',
                                                'mean_SNPindex'])
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as err:
            raise ValueError(('Failed to parse MutPlot intermediate tables.\n'
                              'Please check "{}" and "{}".\n'
                              'Reason: {}').format(snp_index_path,
                                                   sliding_window_path,
                                                   err)) from err

        snp_index['CHROM'] = snp_index['CHROM'].astype('str')
        snp_index['POSI'] = pd.to_numeric(snp_index['POSI'], errors='coerce')/1000000
        sliding_window['CHROM'] = sliding_window['CHROM'].astype('str')
        sliding_window['POSI'] = pd.to_numeric(sliding_window['POSI'], errors='coerce')/1000000
        return snp_index, sliding_window

    def read_contig_length(self):
        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')
        contig_length = {}
        for line in vcf:
            if re.match(r'^##contig', line):
                contig = line.split('=')[2].split(',')[0]
                length = int(line.split('=')[3].replace('>', ''))
                contig_length[contig] = length
        vcf.close()
        return contig_length

    def get_50_contigs(self, N_chr, snp_index, sliding_window):
        contig_length = self.read_contig_length()
        contig_50_names = [k for k, v in sorted(contig_length.items(), key=lambda x:x[1], reverse=True)[:50]]
        N_chr = len(contig_50_names)
        snp_index = snp_index[snp_index['CHROM'].isin(contig_50_names)]
        sliding_window = sliding_window[sliding_window['CHROM'].isin(contig_50_names)]
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
        if snp_index.empty:
            raise RuntimeError(('No records were found in "snp_index.tsv".\n'
                                'This usually means all variants were filtered out.\n'
                                'Please review -d/-D/-m/-S parameters.'))
        if sliding_window.empty:
            raise RuntimeError(('No records were found in "sliding_window.tsv".\n'
                                'Please review upstream SNP-index filtering results.'))

        N_chr = len(sliding_window['CHROM'].unique())
        if N_chr == 0:
            raise RuntimeError(('No contigs are available for plotting.\n'
                                'Please check input VCF content and filtering parameters.'))

        if N_chr > 50:
            print(('!!WARNING!! Your reference genome has too many contigs (>50). '
                   'Therefore, only the 50 longest contigs will be used for plotting.'), file=sys.stderr)
            N_chr, snp_index, sliding_window = self.get_50_contigs(N_chr, 
                                                                   snp_index, 
                                                                   sliding_window)
            if N_chr == 0:
                raise RuntimeError('No contigs remained after selecting the top 50 longest contigs.')

        N_col, N_raw = self.set_plot_style(N_chr)

        fig = plt.figure(figsize=(self.fig_width*N_col, self.fig_height*N_raw))
        plt.subplots_adjust(hspace=self.white_space)

        xmax = snp_index['POSI'].max()
        if pd.isna(xmax):
            raise RuntimeError(('Failed to determine x-axis range from SNP positions.\n'
                                'Please check whether POSI values are numeric in "snp_index.tsv".'))
        sliding_window = sliding_window.groupby('CHROM')
        for i, (chr_name, chr_sliding_window) in enumerate(sliding_window):
            chr_snp_index = snp_index[snp_index['CHROM']==chr_name]
            if not self.plot_with_indel:
                chr_snp_index = chr_snp_index[chr_snp_index['variant']=='snp']

            ax = fig.add_subplot(N_raw, N_col, i+1)
            ax.plot(chr_sliding_window['POSI'],
                    chr_sliding_window['mean_p99'],
                    color=self.line_colors[2],
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'],
                    chr_sliding_window['mean_p95'],
                    color=self.line_colors[1],
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'],
                    chr_sliding_window['mean_SNPindex'],
                    color=self.line_colors[0],
                    linewidth=3)

            if self.snpEff is None:
                ax.scatter(chr_snp_index['POSI'],
                           chr_snp_index['SNPindex'],
                           marker='.',
                           color=self.dot_color)
            else:
                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODIFIER']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODIFIER']['SNPindex'],
                           marker='.',
                           color=self.dot_color)

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='LOW']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='LOW']['SNPindex'],
                           marker='.',
                           color=self.dot_color)

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODERATE']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODERATE']['SNPindex'],
                           marker='+',
                           color=self.dot_color)

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='HIGH']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='HIGH']['SNPindex'],
                           marker='x',
                           color='#772D8B')

            ax.hlines([0.5], 0, xmax, linestyles='dashed', colors='black')
            ax.set_xlabel('Position (Mbp)')
            ax.set_ylabel('SNP-index')
            ax.set_xlim(0, xmax)
            ax.set_ylim(0, 1.05)
            ax.set_title(chr_name, fontsize=17)

        try:
            # Prevent labels/titles from being clipped in the exported figure.
            fig.tight_layout(pad=0.8, h_pad=self.white_space, w_pad=0.8)
            plt.savefig('{}/mutmap_plot.{}'.format(self.out, self.args.format),
                        bbox_inches='tight')
        except Exception as err:
            raise RuntimeError(('Failed to create/save the plot image.\n'
                                'Try increasing "--fig-width"/"--fig-height" or reducing contigs.\n'
                                'Reason: {}').format(err)) from err
        finally:
            plt.close(fig)
