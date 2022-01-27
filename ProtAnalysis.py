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
df['Heveínas'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/Originais/hev_all_blastdbcmd_1kp_cluster100')
df['Snakinas'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/Originais/snak_all_blastdbcmd_1kp_cluster100')
df['Thioninas'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/Originais/thion_all_blastdbcmd_1kp_cluster100')
df['LTP'] = ps.main('C:/Users/beltr/Documents/AMPs_Predictions/Originais/LTP_all_blastdbcmd_1kp_cluster100')


df.to_csv('stast.csv')
