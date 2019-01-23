from helper import size, Coordinate

class Title:
	def __init__(self, inputs, grid):
		self.__coordinate = inputs.get_coordinate()
		self.__grid = grid
		self.__fs = self.__grid.get_fs()
	def set_title(self):
	    '''
	    This program returns the title of the locus
	    '''
	    self.__grid.get_ax_title().set_title('%s' % self.__coordinate.get_info(),
	                               loc='center', fontsize = size(self.__grid.get_fs(), 50))