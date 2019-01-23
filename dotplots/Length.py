import numpy as np
from helper import size
class Length:
    def __init__(self, inputs, referance_length, fs):
        self.__inputs = inputs
        self.__referance_length = referance_length 
        self.__fs = fs
        self.__highlight = self.__inputs.get_highlight()


    def __get_seq_length(self):
        '''
        this function returns the legnthes of human and great apes        
        '''
        all_lengths = []
        for species in self.__highlight:
            lengths = []
            for sample in species[1]:
                lengths.append(max(len(self.__inputs.get_rec_list()[sample].seq) - 2
                              * self.__inputs.get_sequence_padding(), 0))
            all_lengths.append(lengths)
        return all_lengths

    # YOu need to fix this in the future so that the order and species are costomizable
    def graph_ave_length(self, ax):
        '''
        This function graphs the average length for human and great apes
        '''
        #get rid of all the zero length when calculating mean
        lengths = self.__get_seq_length()
        lengths_adj = []
        means = []
        mx = []
        mi = []
        stds = []
        for species in lengths:
            species_adj = [i for i in species if i != 0]
            if len(species_adj) == 0:
                species_adj.append(0);
            lengths_adj.append(species_adj)
            means.append(np.mean(species_adj))
            mx.append(np.max(species_adj))
            mi.append(np.min(species_adj))
            stds.append(np.std(species_adj))
        ax.plot(mx, list(range(len(self.__highlight))), "|", ms=20, mew=7, mec='k')
        ax.plot(mi, list(range(len(self.__highlight))), "|", ms=20, mew=7, mec='k')
        barlist = ax.barh(range(len(self.__highlight)), means, 0.8, xerr=stds, error_kw=dict(lw=10))
        for i in range(len(barlist)):
            barlist[i].set_color(self.__highlight[i][0])
        ax.set_yticks(list(range(len(self.__highlight))))
        labels = [i[2] for i in self.__highlight]
        ax.set_yticklabels(labels, fontsize = size(self.__fs, 30))
        # setting the ticksize and notations
        ax.xaxis.set_tick_params(labelsize=20)
        # setting bar colors
        ax.set_title('Average Length\nwith max, min and std', fontsize=size(self.__fs, 30), loc='left')
        ax.set_xlim(left=0)


