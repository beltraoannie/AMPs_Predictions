#importing libraries for data analysis and manipulation
import numpy as np
import pandas as pd
#importing Biopython library
from Bio.SeqUtils.ProtParam import ProteinAnalysis
#importing libraries for graph generation
import matplotlib.pyplot as plt
import seaborn as sns

class ProteinStats():
    """Preparing Sequence Lists and Dataframes
        Here were defined the functions that were used for sequence analysis by 
        biopython library (ProteinAnalysis), data analysis by  pandas library (pd)
    """

    def fasta(file):
        """This function receives a multi-fasta file and stores only the lines containing sequences in a list
        Args:
            file ([type]): [description]

        Returns:
            [type]: [description]
        """        
        sequencelst = []
        for line in file:
            try:
                if(">" not in line):
                    sequence = str(line).replace("\n", "")
                    sequencelst.append(sequence)
                return sequencelst
            except Exception as e:
                print(e)

    def seq_to_df(list):
        """This function iterates over all the sequences in the previously generated list and saves in a dataframe (table):
           Amino acid counting, size, molecular weight, gravy number and isoelectric point 
        Args:
            list ([type]): [description]

        Returns:
            [type]: [description]
        """        
        AA_Freq = []
        MW = []
        IP =[]
        GVY = []
        lenght =[]
        for j in list:
            try:
                a = ProteinAnalysis(str(j))
                AA_Freq.append(a.count_amino_acids())
                MW.append(a.molecular_weight())
                IP.append(a.isoelectric_point())
                GVY.append(a.gravy())
                lenght.append(len(j))

                df= pd.DataFrame.from_dict(AA_Freq)
                df['MW'] = MW
                df['lenght'] = lenght
                df['IP'] = IP
                df['GVY'] = GVY

                return df
            except Exception as e:
                print(e)

    """Graph Generation Fuctions"""

    def lenght_graph(frame, xlabel, ylabel):
        """
        Here were defined the functions that were used for graph generation using seaborn (sns) and matplotlib (plt) libraries
        For generating frequency analysis of sequence lenght, categories were defined using the cut function, 
        where the bins were defined at intervals ranging from 0 to the highest value (372 aa) found in the set 
        and the data was grouped according to the value of lenght column, then the graph was constructed using sns catplot function
        Args:
            frame ([type]): [description]
            xlabel ([type]): [description]
            ylabel ([type]): [description]

        Returns:
            [type]: [description]
        """        
        labels = [1,2,3,4,5,6,7,8,9]
        bins= [0, 30, 40, 50, 60, 70, 80, 90, 100, 372]
        try:
            frame = frame.copy()
            frame['lenght'],labels =pd.cut(frame['lenght'], bins=bins, labels=labels, retbins=True)
            print(frame.lenght.value_counts(normalize=True)*100)
            l = sns.catplot(x="lenght", kind="count",
                        color='black',data=frame);
            l.set_xticklabels(['≤ 30','≤ 40', '≤ 50', '≤ 60', '≤ 70', '≤ 80', '≤ 90', '≤ 100', '>100' ], rotation= 45)
            l.set_xlabels(xlabel)
            l.set_ylabels(ylabel);
            return l
        except Exception as e:
            print(e)
    
    def ip_graph(frame, xlabel, ylabel):
        """For generating frequency analysis of isoeletric point (IP), 
        categories (bins) were defined according to pH scale, ranging 1 to 14, 
        the data was grouped according to the value of the IP column, then the graph was constructed using sns catplot function

        Args:
            frame ([type]): [description]
            xlabel ([type]): [description]
            ylabel ([type]): [description]

        Returns:
            [type]: [description]
        """        
        labels = [1,2,3,4,5,6,7,8,9,10,11,12,13]
        bins= [1,2,3,4,5,6,7,8,9,10,11, 12,13,14]
        frame = frame.copy()
        frame['IP'],labels =pd.cut(frame['IP'], bins=bins, labels=labels, retbins=True)
        i = sns.catplot(x="IP", kind="count",
                    color='black', data=frame);
        i.set_xticklabels(['1','≤ 2','≤ 3','≤ 4','≤ 5', '≤ 6', '≤ 7', '≤ 8', '≤ 9', '≤ 10', '≤ 11', '≤ 12','≤ 13','≤ 14' ], rotation= 45)
        i.set_xlabels(xlabel)
        i.set_ylabels(ylabel)
        return i

    def mw_graph(frame, xlabel, ylabel):
        """*For generating frequency analysis of molecular weight (WP), 
        categories (bins) were defined at intervals ranging from 0 to the highest value (50.000 kDa) 
        found in the set the data was grouped according to the value of the WP column, 
        then the graph was constructed using sns catplot function*

        Args:
            frame ([type]): [description]
            xlabel ([type]): [description]
            ylabel ([type]): [description]

        Returns:
            [type]: [description]
        """        
        labels = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
        bins= [0,1000,2000,4000,6000,8000,10000,12000,14000,16000,18000,20000,30000,40000,50000]
        frame=frame.copy()
        frame['MW'],labels =pd.cut(frame['MW'], bins=bins, labels=labels, retbins=True)
        print(frame.MW.value_counts(normalize=True)*100)
        m = sns.catplot(x="MW", kind="count",
                    color='black',data=frame);
        m.set_xticklabels(['<1000','≤ 2000','≤ 4000','≤ 6000', '≤ 8000', '≤ 10.000', '≤ 12.000', '≤ 14.000', '≤ 16.000', '≤ 18.000', '≤ 20.000','≤ 30.000','≤ 40.000','≤ 50.000'], rotation= 45)
        m.set_xlabels(xlabel)
        m.set_ylabels(ylabel)
        return m
    def gvy_graph(frame, xlabel, ylabel):
        """For generating frequency analysis of Gravy Number (GVY), 
        categories (bins) were defined according to Gravy scale, 
        ranging -2 to 2, the data was grouped according to the value of the GVY column in 13 intervals, 
        then the graph was constructed using sns catplot function

        Args:
            frame ([type]): [description]
            xlabel ([type]): [description]
            ylabel ([type]): [description]

        Returns:
            [type]: [description]
        """        
        labels = [1,2,3,4,5,6,7,8,9,10,11,12,13]
        bins= [-2,-1.8,-1.6,-1.4,-1.2,-1,-0.5,0.5,1,1.2,1.4,1.6,1.8,2]
        frame = frame.copy()
        #print(pd.cut(GVY['GVY'], bins=bins, retbins=True))
        frame['GVY'],labels =pd.cut(frame['GVY'], bins=bins, labels=labels, retbins=True)
        print(frame.GVY.value_counts(normalize=True)*100)
        g = sns.catplot(x="GVY", kind="count",
                        color='black',data=frame);
        g.set_xticklabels(['≤ -2','≤ -1.8','≤ -1.6','≤ -1.4', '≤ -1.2','≤ -1','≤ 0,5','≤ 1', '≤ 1.2', '≤ 1.4', '≤ 1.6', '≤ 1.8', '≤ 2' ], rotation= 45)
        g.set_xlabels(xlabel)
        g.set_ylabels(ylabel)
        return g