from helper import kmer_arr, kmer_iter
from helper import size
import numpy as np
import warnings

class Dotplot:
    def __init__(self, inputs, axes, fs):
        self.__inputs = inputs;
        self.__axes = axes
        self.__highlight = self.__inputs.get_highlight()
        self.__fs = fs
        self.__samp_count = self.__inputs.get_samp_count()
        self.__resolution = size(self.__fs, (72./300)**2*25)
    
    def plot_dotplots(self):
        for row in range(self.__samp_count):
            for col in range(self.__samp_count):
                self.__plot_ax(
                    self.__inputs.get_rec_list()[col],
                    self.__inputs.get_rec_list()[row],
                    self.__inputs.get_sliding_window_size(),
                    self.__axes[row, col],
                    self.__inputs.get_rec_list()[col].id,
                    self.__inputs.get_rec_list()[row].id,
                    row,
                    col
                    )

    def __turn_off_ax(self, ax):
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.tick_params(axis='x', colors='w')
        ax.tick_params(axis='y', colors='w')

    def __update(self, ax, len1, len2):
        ax.set_aspect('equal')
        ax.set_xlim([0, len1-self.__inputs.get_sliding_window_size()])
        ax.set_ylim([0, len2-self.__inputs.get_sliding_window_size()])
        ax.set_ylim(ax.get_ylim()[::-1])        # invert the axis
        ax.xaxis.tick_top()                     # and move the X-Axis  
        ax.tick_params(axis='x', pad=0)
        ax.tick_params(axis='y', pad=0)

    def __highlight_ax(self, row, col, ax):
        for color in self.__highlight:
            if row in color[1]:
                ax.spines['left'].set_color(color[0])
                ax.spines['left'].set_linewidth(size(self.__fs, 3))
                ax.tick_params(axis = 'y', colors = color[0])
            if col in color[1]:
                ax.spines['top'].set_color(color[0])
                ax.spines['top'].set_linewidth(size(self.__fs, 3))
                ax.tick_params(axis='x', colors=color[0])

    def __highlight_label(self, row, col, ax):
        for color in self.__highlight:
            if row in color[1]:
                ax.yaxis.label.set_color(color[0])
                ax.yaxis.label.set_weight(size(self.__fs, 20))
            if col in color[1]:
                ax.xaxis.label.set_color(color[0])
                ax.xaxis.label.set_weight(size(self.__fs, 20))

    def __plot_lable(self, row, col, ax, files_x, files_y):
        if row == 0:
            ax.set_xlabel(files_x)
            ax.xaxis.labelpad = size(self.__fs, 20)
            ax.xaxis.label.set_size(size(self.__fs, 15))
            ax.xaxis.set_label_position('top')
        if row == col:
            ax.set_ylabel(files_y)
            ax.yaxis.labelpad = size(self.__fs, 20)
            ax.yaxis.label.set_size(size(self.__fs, 15))

    def __plot_ax(self, seq_rec1, seq_rec2, sliding_window_size, ax, files_x, files_y, row = 0, col = 0):
        '''
        resolution: size of each dot
        '''
        if row > col:                # Turns the redundant self.__plot_axs off
            ax.axis('off')
        else:
            seq1 = str(seq_rec1.seq).upper()
            seq2 = str(seq_rec2.seq).upper()
            len1 = len(seq1)
            len2 = len(seq2)
            v_x, v_y = self.__super_string(seq1, seq2, sliding_window_size)
            if (v_x == []) | (v_y == []):
                self.__turn_off_ax(ax)
            else:
                ax.scatter(v_x,v_y, c = '#000000', s = self.__resolution, marker = '.', edgecolor= '', rasterized = True)
                self.__update(ax, len1, len2)
                self.__highlight_ax(row, col, ax)
            self.__plot_lable(row, col, ax, files_x, files_y)
            self.__highlight_label(row, col, ax)

    def __super_string(self, s1, s2, k):
        '''
        generate the two vectors that catalog the position of which two sequences have exact kmer
        input:
        s1:sequence1
        s2:sequence2
        k:size of kmer
        output:
        v_x:vector x
        v_y:vector y
        '''
        v_x = []
        v_y = []
        min_seqlen = min(len(s1), len(s2))  #check if the sequence is shorter than k
        if min_seqlen < k:                       #s1 would be on the x axis
            warnings.warn("seqeunce length smaller than sliding window size: %d, %d, %d"%(len(s1), len(s2), k))
            return v_x, v_y
        else:
            pass
        kmer = kmer_arr(s2, k)    #generate two kmer_itter from two strings and turn them into arrays
        for km, index in kmer_iter(s1, k):
            position = np.where(kmer == km)[0]
            for i in position:
                v_x.append(index-k)
                v_y.append(i)
        return v_x, v_y

