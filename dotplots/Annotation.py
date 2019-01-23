from helper import Nearest_gene, size, cordinate_converter, Coordinate
import matplotlib.patches as patches
import re
import numpy as np
class Annotation:
    def __init__(self, inputs, ax, fs, gene_df, red_box_size, search_pad):
        self.__inputs = inputs
        self.__ax = ax
        self.__red_box_size = red_box_size
        self.__fs = fs
        self.__coordinate = self.__inputs.get_coordinate()
        self.__search_pad = search_pad
        if len(gene_df) > 1:
            self.__gene_df1 = gene_df[0]
            self.__gene_df2 = gene_df[1]
            self.__is_genic = False
        else:
            self.__gene_df = gene_df[0]
            self.__is_genic = True
        self.__head_length = 0.04
        self.__head_width = 0.03

    def plot_annotation(self, gene_color='k'):
        if self.__is_genic:
            self.__plot_genic(gene_color)

        else:
            kwargs = {
            'fc':'k',
            'ec':'k',
            'head_width':self.__head_width,
            'head_length':self.__head_length,
            'length_includes_head':True,
            'zorder' : 4
            }
            self.__plot_intergenic(kwargs)

    def __plot_intergenic(self, kwargs):
        self.__ax[0].set_xlim([-self.__head_length, 1 + self.__head_length])
        self.__ax[0].set_ylim([-self.__head_width/1.5, self.__head_width/1.5])            
        self.__ax[0].arrow(0.0001, 0, -self.__head_length/2, 0, **kwargs)
        self.__ax[0].arrow(self.__head_length/2, 0, -self.__head_length/2, 0, **kwargs)
        self.__ax[0].arrow(1, 0, self.__head_length/2, 0, **kwargs)
        self.__ax[0].arrow(1 - self.__head_length/2, 0, self.__head_length/2, 0, **kwargs)
        self.__ax[0].set_title(
            'there is no gene in region\n %s:%d-%d'% (
                self.__coordinate.get_chrom(), 
                self.__coordinate.get_b() - self.__search_pad - self.__inputs.get_sequence_padding(), 
                self.__coordinate.get_e() + self.__search_pad + self.__inputs.get_sequence_padding()
                ),
            fontsize = size(self.__fs, 30), y = 0.03
            )
        gene_us = Nearest_gene()
        gene_ds = Nearest_gene()
        if self.__gene_df1['strand'] == '+':
            gene_us.utr = '3\''
        else:
            gene_us.utr = '5\''
        gene_us.coordinate = self.__coordinate.get_b() - int(self.__gene_df1['txEnd'])
        gene_us.name = self.__gene_df1['name2']
        if self.__gene_df2['strand'] == '-':
            gene_ds.utr = '3\''
        else:
            gene_ds.utr = '5\''
        gene_ds.coordinate = int(self.__gene_df2['txStart']) - self.__coordinate.get_e()
        gene_ds.name = self.__gene_df2['name2']
        self.__ax[0].text(
            self.__head_length + 0.02, 
            -self.__head_width/2 - 0.005, 
            '%s\n%s\n%d'%(gene_us.name, gene_us.utr, gene_us.coordinate), 
            fontsize = size(self.__fs, 20), 
            horizontalalignment='left', 
            zorder = 2)
        self.__ax[0].text(
            1-self.__head_length - 0.02, 
            -self.__head_width/2 -0.005, 
            '%s\n%s\n%d'%(gene_ds.name, gene_ds.utr, gene_ds.coordinate), 
            fontsize = size(self.__fs, 20), 
            zorder = 2, 
            horizontalalignment='right')

    def __plot_arrows(self, i):
        # Patches of white rectangles to block the background for the arrows
        arrow1 = patches.Rectangle((-self.__head_length, -self.__head_width / 2),     #width of arrow 1
                self.__head_length, self.__head_width, fc='w', zorder=3.5)
        arrow2 = patches.Rectangle((1, -self.__head_width / 2), 1             #width of arrow 2
                + self.__head_length, self.__head_width, fc='w', zorder=3.5)
        self.__ax[i].add_patch(arrow1)
        self.__ax[i].add_patch(arrow2)
        if self.__gene_df.iloc[i, 3] == '+':
            self.__ax[i].arrow(
                -self.__head_length,
                0,
                0.0001,
                0,
                fc='k',
                ec='k',
                head_width=self.__head_width,
                head_length=self.__head_length,
                zorder=4,
                )

            self.__ax[i].arrow(
                1,
                0,
                self.__head_length,
                0,
                fc='k',
                ec='k',
                head_width=self.__head_width,
                head_length=self.__head_length,
                length_includes_head=True,
                zorder=4,
                )
        else:
             # when its on the - strand
            self.__ax[i].arrow(
                0.0001,
                0,
                -0.04,
                0,
                fc='k',
                ec='k',
                head_width=self.__head_width,
                head_length=self.__head_length,
                length_includes_head=True,
                zorder=4,
                )
            self.__ax[i].arrow(
                1 + self.__head_length,
                0,
                -self.__head_length,
                0,
                fc='k',
                ec='k',
                head_width=self.__head_width,
                head_length=self.__head_length,
                length_includes_head=True,
                zorder=4,
                )

    def __plot_red_box(self, i):
        interval = (self.__coordinate.get_e() - self.__coordinate.get_b()) \
        / (self.__coordinate.get_e() - self.__coordinate.get_b() + 2 * self.__search_pad + 2 * self.__inputs.get_sequence_padding())                 # the size of the red box in the ax's coordinate
        red_box = patches.Rectangle(
            (0.5 - interval / 2, -0.05),
            interval,
            0.1,
            fill=False,
            color='r',
            lw=3,
            zorder=3,
            )
        self.__ax[i].add_patch(red_box)

    def __plot_ORF(self, gene_color, i, txStart, txEnd):
        self.__ax[i].hlines(
            0,
            cordinate_converter(
                self.__coordinate.get_b(), self.__coordinate.get_e(), 
                txStart, self.__search_pad, self.__inputs.get_sequence_padding()
                ),
            cordinate_converter(
                self.__coordinate.get_b(), self.__coordinate.get_e(), 
                txEnd, self.__search_pad, self.__inputs.get_sequence_padding()
                ),
            linestyle='-',
            lw=8,
            color=gene_color,
            zorder=1,
            )

    def __plot_genic(self, gene_color):
        # Plot red box size
        self.__ax[0].set_title(
            '%dbp' % self.__red_box_size, fontsize = size(self.__fs, 20)
            )
        #for each isoform
        for i in range(self.__gene_df.shape[0]):
            txStart = self.__gene_df.iloc[i, 4] 
            txEnd = self.__gene_df.iloc[i,   5]
            self.__plot_track(i)
            self.__plot_arrows(i)
            self.__plot_red_box(i)
            self.__plot_ORF(gene_color, i, txStart, txEnd)
            self.__plot_exons(gene_color, i)

    def __plot_track(self, i):
        self.__ax[i].hlines(0, 0, 1, linestyle='--', zorder=0)
        self.__ax[i].set_xlim([-self.__head_length, 1 + self.__head_length])
        self.__ax[i].annotate(
            '%s(ex count:%s)' % (
                self.__gene_df.iloc[i, 12],
                self.__gene_df.iloc[i, 8]), 
            xy=(0, 1),
            xycoords='axes fraction', fontsize=size(self.__fs, 25)
            )  # this is the fontsize for gene name and exon count

    def __plot_exons(self, gene_color, i):
        exStarts = re.findall(r'\d+', self.__gene_df.iloc[i, 9])
        exStarts = [int(i) for i in exStarts]
        exEnds = re.findall(r'\d+', self.__gene_df.iloc[i, 10])
        exEnds = [int(i) for i in exEnds]

        for j in range(np.size(exStarts)):

            ex = patches.Rectangle(
                (cordinate_converter(self.__coordinate.get_b(),
                 self.__coordinate.get_e(), exStarts[j], self.__search_pad, self.__inputs.get_sequence_padding()),
                 -(self.__head_width / 2)),
                (exEnds[j] - exStarts[j]) / (self.__coordinate.get_e()
                        - self.__coordinate.get_b() + 2 * self.__search_pad),
                self.__head_width,
                fill=True,
                color=gene_color,
                zorder=2,
                )
            self.__ax[i].add_patch(ex)
            self.__annotate_exons(exStarts, i, j)


    def __annotate_exons(self, exStarts, i, j):
        if cordinate_converter(self.__coordinate.get_b(),
                self.__coordinate.get_e(), exStarts[j],
                self.__search_pad, self.__inputs.get_sequence_padding()) > 0 \
            and cordinate_converter(self.__coordinate.get_b(),
                self.__coordinate.get_e(), exStarts[j],
                self.__search_pad, self.__inputs.get_sequence_padding()) < 1:
            if self.__gene_df.iloc[i, 3] == '+':
                exon_number = j + 1
            else:
                exon_number = np.size(exStarts) - j
            self.__ax[i].annotate(
            	'%d' % exon_number, 
            	xy=(
            		cordinate_converter(
            			self.__coordinate.get_b(), self.__coordinate.get_e(), 
                        exStarts[j], self.__search_pad, self.__inputs.get_sequence_padding()
        			),
    				self.__head_width / 2 - 0.01),
                xytext = (
                	cordinate_converter(
                		self.__coordinate.get_b(), self.__coordinate.get_e(), 
                        exStarts[j], self.__search_pad, self.__inputs.get_sequence_padding()
            		),
            		self.__head_width / 2 - 0.01 + 0.024), 
                fontsize=size(self.__fs,25)
            )