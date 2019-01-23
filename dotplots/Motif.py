import matplotlib.pyplot as plt
from helper import size

MAX_MOTIF_LENGTH = 300

class Motif:

    def __init__(self, fig, motif_list, highlight, position, fs):
        self.__motif_list = motif_list
        for i in range(len(self.__motif_list)):
            for j in range(len(self.__motif_list[i])):
                self.__motif_list[i][j] = str(self.__motif_list[i][j])
        maxlen = 0
        for i in self.__motif_list:
            if len(i[0]) > maxlen:
                maxlen = len(i[0])
        for i in range(len(self.__motif_list)):
            self.__motif_list[i][0] = self.__motif_list[i][0] + ' ' * (maxlen - len(self.__motif_list[i][0]))
        self.__highlight = highlight
        self.__fig = fig
        self.__fs = fs
        self.position = position
    def plot(self):
        color = ''
        for i in range(len(self.__motif_list)):
            for sample in self.__highlight:
                if i in sample[1]:
                    color = sample[0]
            if len(self.__motif_list[i][1]) > MAX_MOTIF_LENGTH:
                self.__motif_list[i][1] = self.__motif_list[i][1][:MAX_MOTIF_LENGTH] + ">>>" + str(len(self.__motif_list[i][1]))
            self.__fig.text(
                self.position[0], self.position[1] - size(self.__fs, 0.01 * i), 
                ' '.join(self.__motif_list[i][:4]), color = color, size = size(self.__fs, 25), family = 'monospace')