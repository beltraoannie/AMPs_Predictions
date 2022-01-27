#importing libraries for data analysis and manipulation
import numpy as np
import pandas as pd
from proteinstats import ProteinStats

ps = ProteinStats()



details = {
            'Characteristic' : ['Number of sequences','Number of cysteines','Main AA Composition (Mean)','Main AA Composition (Mode)','Main AA Variation (Standard Deviation)','Length','Max Lenght','Min Lenght','Molecular Weight'], 
        }
        
# creating a Dataframe object 
df = pd.DataFrame(details)
df['Heveínas'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/phytAMP/hevein_phytamp.fasta')
df['Snakinas'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/phytAMP/snakin_phytamp.fasta')
df['Thioninas'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/phytAMP/thionin_phytamp.fasta')
df['LTP'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/phytAMP/ltp_phytamp.fasta')
df['Defesinas'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/phytAMP/defensin_phytamp.fasta')
df['Ciclotídeos'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/phytAMP/ciclotide_phytamp.fasta')


df.to_csv('stast.csv')
