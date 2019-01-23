from helper import size
import matplotlib
matplotlib.use('Agg')  # force matplotlib not to use -Xwindow backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Grid:
    def __init__(self, inputs, gene_df):
        self.__inputs = inputs
        self.__gene_df = gene_df
        samp_count = self.__inputs.get_samp_count()
        self.__fs = samp_count * 4 / 10
        self.__figsize = (self.__fs * 10, self.__fs * 10)
        self.__f, self.__ax_dotplots = plt.subplots(samp_count, samp_count, figsize=self.__figsize)
        self.__box_dotplot = self.get_ax_dotplots()[0,0].get_position().bounds
        self.__ax_annotation = self.__get_ax_annotation()
        self.__ax_title = self.__get_ax_title()
        self.__ax_length = self.__get_ax_length()
        self.__motif_postion = self.__get_motif_postion()
        self.__ax_gc = self.__get_ax_gc()



    def get_fs(self):
        return self.__fs

    def get_figure(self):
        return self.__f

    def get_ax_dotplots(self):
        return self.__ax_dotplots

    def get_ax_annotation(self):
        return self.__ax_annotation

    def get_ax_title(self):
        return self.__ax_title

    def get_ax_length(self):
        return self.__ax_length

    def get_motif_position(self):
        return self.__motif_postion

    def get_ax_gc(self):
        return self.__ax_gc

    def __get_ax_gc(self):
        ax_gc = self.__f.add_axes((0.5, 0.1, 0.25, 0.25))
        ax_gc.axis('off')
        return ax_gc

    def __get_ax_annotation(self):
        box_dotplot5 = self.get_ax_dotplots()[5,5].get_position().bounds
        box_dotplot3 = self.get_ax_dotplots()[3,3].get_position().bounds
        is_genic = True
        if len(self.__gene_df) > 1:
            is_genic = False
        ax_annotaitons = []
        if not is_genic:  # this is true when there is no gene in the region
            pos = [self.__box_dotplot[0], 0.46538461538461534, 0.33561538461538465, 0.02132307692307693]
            ax_annotaitons.append(self.__f.add_axes(pos))
            ax_annotaitons[0].axis('off')
        else:
            for i in range(self.__gene_df[0].shape[0]):
                pos = [self.__box_dotplot[0], 0.46538461538461534 - i * 0.02632307692307693, 0.33561538461538465, 0.02132307692307693]
                ax_annotaitons.append(self.__f.add_axes(pos))
                ax_annotaitons[i].axis('off')
        return ax_annotaitons

    def __get_ax_title(self):
        pos = [0.5, self.__box_dotplot[1], 0, self.__box_dotplot[3] + size(self.__fs,0.025)]
        ax_title = self.__f.add_axes(pos)
        ax_title.axis('off')
        ax_title.set_position(pos)
        return ax_title

    def __get_ax_length(self):
        pos = [self.__box_dotplot[0] + 0.03, 0.5364615384615384, 0.19653846153846155, 0.05923076923076931]
        ax_len = self.__f.add_axes(pos)
        ax_len.spines['left'].set_color('w')
        ax_len.spines['right'].set_color('w')
        ax_len.spines['top'].set_color('w')
        return ax_len

    def __get_motif_postion(self):
        position = []
        position.append(self.__box_dotplot[0])
        position.append(0.09)
        return position

    def save_fig(self):
        self.__f.savefig('{}'.format(self.__inputs.get_outfile()), facecolor='w', bbox_inches='tight', dpi = 200, pad_inches = 0)