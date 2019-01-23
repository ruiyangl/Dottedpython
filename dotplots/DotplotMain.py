from Inputs import Inputs
from Region import Region
from Grid import Grid
from Dotplot import Dotplot
from Annotation import Annotation
from Title import Title
from Length import Length
from Motif import Motif
from Gc import Gc


def main():   ## Should I make is_plotting a feild or a function
    my_inputs = Inputs()
    my_region = Region(my_inputs, 5000)
    my_grid = Grid(my_inputs, my_region.get_gene_df())
    my_dotplot = Dotplot(my_inputs, my_grid.get_ax_dotplots(), my_grid.get_fs())

    print("plotting dotplots")
    my_dotplot.plot_dotplots()

    if my_inputs.is_plotting_annotation:
        print("plotting annotation")
        my_annotation = Annotation(my_inputs, my_grid.get_ax_annotation(), my_grid.get_fs(),
        my_region.get_gene_df(), my_region.get_red_box_size(), my_region.get_search_pad())
        my_annotation.plot_annotation()

    if my_inputs.is_plotting_length_stats:
        print("plotting length stats")
        my_length_stats = Length(my_inputs, my_region.get_red_box_size(), my_grid.get_fs())
        my_length_stats.graph_ave_length(my_grid.get_ax_length())

    if my_inputs.is_plotting_motif:
        print("plotting motif")
        my_motif = Motif(my_grid.get_figure(), my_inputs.get_motif(), my_inputs.get_highlight(), my_grid.get_motif_position(), my_grid.get_fs())
        my_motif.plot()

    if my_inputs.is_plotting_gc:
        print("plotting gc content")
        my_Gc = Gc(my_inputs, my_region.get_search_pad(), my_grid.get_fs())
        my_Gc.plot(my_grid.get_ax_gc())
        
    my_title = Title(my_inputs, my_grid)
    my_title.set_title()
    my_grid.save_fig()
if __name__ == '__main__':
    main()