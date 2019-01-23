class Region:
    def __init__(self, inputs, search_pad):
        self.__inputs = inputs
        self.__search_pad = search_pad
        self.__gene_df = self.__grab_gene()
        self.red_box_size = 0
    def __grab_gene(self):
        '''
        this method returns a data frame containing the gene(s) in the search area
        #-----------------------------------------------------------------------------#
        |there 6 possible relative position of our search area and a gene:            |
        |                                                                             |
        |     --------------------|      Search Area       |-------------------       |
        |       |gene1|         |gene2|     |gene3|      |gene4|        |gene5|       |
        |                   |                gene6                 |                  |
        |                                                                             |
        |gene 2-4 and 6 are considered in the search area                             |       
        |we let b and e to be the begning and end coordinate of our search area and   |
        |the coordiante for each gene is s-f, then this function discard all the gene |
        |that satisfy [(s<b)&(f<b)] or [(s>e)&(f>e)] i.e. gene 1 or gene2             |
        #-----------------------------------------------------------------------------#
        '''
        coordinate = self.__inputs.get_coordinate()
        b = coordinate.get_b()
        e = coordinate.get_e()
        self.__red_box_size = e - b
        chrom = coordinate.get_chrom()
        refgene = self.__inputs.get_refgene()
        sequence_padding = self.__inputs.get_sequence_padding()
        s = min(b, e) - self.__search_pad - sequence_padding #strat of the search area
        f = max(b, e) + self.__search_pad + sequence_padding #finish of the search area
        for (i, i_data) in refgene.groupby('chrom'):
            if i == chrom:
                temp = (s < i_data['txStart']) & (s < i_data['txEnd']) \
                    & (f < i_data['txStart']) & (f < i_data['txEnd']) \
                    | (s > i_data['txStart']) & (s > i_data['txEnd']) \
                    & (f > i_data['txStart']) & (f > i_data['txEnd'])
                genes = i_data[~temp]
            else:
                pass
        if genes.shape[0] == 0:
            df_txStart, df_txEnd = self.__make_sorted_df(refgene, chrom)
            return list([self.__binary_search(b, df_txEnd).iloc[0], self.__binary_search(e, df_txStart).iloc[1]])
        else:
            return list([genes])

    def get_gene_df(self):
        return self.__gene_df

    def get_red_box_size(self):
        return self.__red_box_size

    def get_search_pad(self):
        return self.__search_pad

    def __binary_search(self, a, df):
        '''
        do a binary search of a value on a sorted list
        
        inputs:
        a: value which should not exist in the list
        lis: asorted list to do the search on
        
        output:
        a tuple of the two value samller than and greater than a
        '''
        half = int(df.shape[0]/2)
        if df.shape[0] > 2:
            if df.iloc[half, 5] > a:
                return self.__binary_search(a, df[:half+1])
            else:
                return self.__binary_search(a, df[half:])
        else:
            return df

    def __make_sorted_df(self, refgene, chrom):
        '''
        make 2 lists out of refgene, one sorted by stStart, the other one sorted by txEnd
        '''
        rg_chrom = refgene['chrom'] == '{}'.format(chrom)
        df_txStart = refgene[rg_chrom].sort_values('txStart')
        df_txEnd = refgene[rg_chrom].sort_values('txEnd')
        return df_txStart, df_txEnd

