from helper import size
from Bio.SeqUtils import GC

class Gc:

    def __init__(self, inputs, search_pad, fs):
        self.__inputs = inputs
        self.__gc_list = self.__get_GC_content()
        self.__fs = fs
        self.__search_pad = search_pad

    def __get_GC_content(self):
        '''
        this function returns the GC content of the vntr region for each sample
        '''
        gc_list = []
        for sample in self.__inputs.get_rec_list():
            gc_list.append(GC(sample.seq[self.__inputs.get_sequence_padding():-self.__inputs.get_sequence_padding()]))
        return gc_list

    def plot(self, ax):
        gc_string = ''.join(str('%.2f\n' % float(e)) for e in
                            self.__gc_list)
        formations = (
            self.__search_pad,
            self.__inputs.get_sliding_window_size(),
            self.__inputs.get_outfile(),
            int(self.__inputs.get_sequence_padding()),
            gc_string,
            )
        ax.text(0, 0,

            '''
search padding size: %d
sliding window size: %d
file:%s
seqence_padding%d
GC contents\n%s

            '''% formations, fontsize = size(self.__fs, 20))


