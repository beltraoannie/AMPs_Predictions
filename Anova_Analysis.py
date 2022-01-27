from scipy.stats import f_oneway 
from os import listdir
import pandas as pd
import numpy as np

 
path = 'C:/Users/beltr/Documents/AMPs_Predictions/Graficos_Isadora'

primer_9_hevein_Lus10028377_ROOT = pd.read_csv(f'{path}/primer12B.csv') 
primer_9_hevein_Lus10028377_ABOVE = pd.read_csv(f'{path}/primer12AB.csv')
Hevein_Lus10006552_B11_BELLOW =  pd.read_csv(f'{path}/primer2B.csv')
Hevein_Lus10006552_B11_ABOVE = pd.read_csv(f'{path}/primer2AB.csv') 
LTP_Lus10026418_F9_ABOVE = pd.read_csv(f'{path}/primer3AB.csv')
LTP_Lus10026418_F9_BELLOW =  pd.read_csv(f'{path}/primer3B.csv')
Snakin_Lus10017212_A11_ABOVE = pd.read_csv(f'{path}/primer4AB.csv') 
Snakin_Lus10017212_A11_BELLOW = pd.read_csv(f'{path}/primer4B.csv')
Snakin_Lus10042203_3_BELLOW =  pd.read_csv(f'{path}/primer8B.csv')
Snakin_Lus10042203_3_ABOVE =  pd.read_csv(f'{path}/primer8AB.csv')
Hevein_Lus100000453_7_BELLOW =  pd.read_csv(f'{path}/primer9B.csv')
Hevein_Lus100000453_7_ABOVE =  pd.read_csv(f'{path}/primer9AB.csv')

frames  = [primer_9_hevein_Lus10028377_ROOT,primer_9_hevein_Lus10028377_ABOVE, 
Hevein_Lus10006552_B11_BELLOW,Hevein_Lus10006552_B11_ABOVE, 
LTP_Lus10026418_F9_ABOVE,LTP_Lus10026418_F9_BELLOW,
Snakin_Lus10017212_A11_ABOVE, Snakin_Lus10017212_A11_BELLOW,
Snakin_Lus10042203_3_BELLOW,Snakin_Lus10042203_3_ABOVE,
Hevein_Lus100000453_7_BELLOW,Hevein_Lus100000453_7_ABOVE]

primers = pd.concat(frames) 
primers['LOG'] = np.log(primers['RE'])
primers['Days'] = primers['Days'].astype(str)

def anova_amostras(primer,days,tecido,amostra,tipo):
    df = primers.loc[(primers['PS']==f'{tecido}') & (primers['Days']==f'{days}') & (primers['Primer'] == f'{primer}') & (primers['TRT'].isin([f'{amostra}','Mock']))][['TRT','LOG','PS']]
    stats = f_oneway(df.loc[df['TRT']==f'{amostra}']['LOG'].tolist(),df.loc[df['TRT']=='Mock']['LOG'].tolist())
    if tipo == 0:
        return stats[0]
    else:
        return stats[1]

Bellow = 'Above'
Diff = 'Above'

details = {
            'Primers' : ['primer_9_hevein_Lus10028377','Hevein_Lus10006552_B11', 
'LTP_Lus10026418_F9','Snakin_Lus10017212_A11','Snakin_Lus10042203_3','Hevein_Lus100000453_7'],

        'FUS_9': [anova_amostras('primer_9_hevein_Lus10028377','9',f'{Diff}','Fusarium',0),
                  anova_amostras('Hevein_Lus10006552_B11','9',f'{Bellow}','Fusarium',0),
                  anova_amostras('Hevein_Lus100000453_7','9',f'{Bellow}','Fusarium',0),
                  anova_amostras('LTP_Lus10026418_F9','9',f'{Bellow}','Fusarium',0),
                  anova_amostras('Snakin_Lus10017212_A11','9',f'{Bellow}','Fusarium',0),
                  anova_amostras('Snakin_Lus10042203_3','9',f'{Bellow}','Fusarium',0)],

        'FUS_14': [anova_amostras('primer_9_hevein_Lus10028377','14',f'{Diff}','Fusarium',0),
                  anova_amostras('Hevein_Lus10006552_B11','14',f'{Bellow}','Fusarium',0),
                  anova_amostras('Hevein_Lus100000453_7','14',f'{Bellow}','Fusarium',0),
                  anova_amostras('LTP_Lus10026418_F9','14',f'{Bellow}','Fusarium',0),
                  anova_amostras('Snakin_Lus10017212_A11','14',f'{Bellow}','Fusarium',0),
                  anova_amostras('Snakin_Lus10042203_3','14',f'{Bellow}','Fusarium',0)],

        'RHIZO_9': [anova_amostras('primer_9_hevein_Lus10028377','9',f'{Diff}','Rhizo',0),
                  anova_amostras('Hevein_Lus10006552_B11','9',f'{Bellow}','Rhizo',0),
                  anova_amostras('Hevein_Lus100000453_7','9',f'{Bellow}','Rhizo',0),
                  anova_amostras('LTP_Lus10026418_F9','9',f'{Bellow}','Rhizo',0),
                  anova_amostras('Snakin_Lus10017212_A11','9',f'{Bellow}','Rhizo',0),
                  anova_amostras('Snakin_Lus10042203_3','9',f'{Bellow}','Rhizo',0)],

        'RHIZO_14': [anova_amostras('primer_9_hevein_Lus10028377','14',f'{Diff}','Rhizo',0),
                  anova_amostras('Hevein_Lus10006552_B11','14',f'{Bellow}','Rhizo',0),
                  anova_amostras('Hevein_Lus100000453_7','14',f'{Bellow}','Rhizo',0),
                  anova_amostras('LTP_Lus10026418_F9','14',f'{Bellow}','Rhizo',0),
                  anova_amostras('Snakin_Lus10017212_A11','14',f'{Bellow}','Rhizo',0),
                  anova_amostras('Snakin_Lus10042203_3','14',f'{Bellow}','Rhizo',0)]
        }

details_stats = {
            'Primers' : ['primer_9_hevein_Lus10028377','Hevein_Lus10006552_B11', 
'LTP_Lus10026418_F9','Snakin_Lus10017212_A11','Snakin_Lus10042203_3','Hevein_Lus100000453_7'],

        'FUS_9': [anova_amostras('primer_9_hevein_Lus10028377','9',f'{Diff}','Fusarium',1),
                  anova_amostras('Hevein_Lus10006552_B11','9',f'{Bellow}','Fusarium',1),
                  anova_amostras('Hevein_Lus100000453_7','9',f'{Bellow}','Fusarium',1),
                  anova_amostras('LTP_Lus10026418_F9','9',f'{Bellow}','Fusarium',1),
                  anova_amostras('Snakin_Lus10017212_A11','9',f'{Bellow}','Fusarium',1),
                  anova_amostras('Snakin_Lus10042203_3','9',f'{Bellow}','Fusarium',1)],

        'FUS_14': [anova_amostras('primer_9_hevein_Lus10028377','14',f'{Diff}','Fusarium',1),
                  anova_amostras('Hevein_Lus10006552_B11','14',f'{Bellow}','Fusarium',0),
                  anova_amostras('Hevein_Lus100000453_7','14',f'{Bellow}','Fusarium',1),
                  anova_amostras('LTP_Lus10026418_F9','14',f'{Bellow}','Fusarium',1),
                  anova_amostras('Snakin_Lus10017212_A11','14',f'{Bellow}','Fusarium',1),
                  anova_amostras('Snakin_Lus10042203_3','14',f'{Bellow}','Fusarium',1)],

        'RHIZO_9': [anova_amostras('primer_9_hevein_Lus10028377','9',f'{Diff}','Rhizo',1),
                  anova_amostras('Hevein_Lus10006552_B11','9',f'{Bellow}','Rhizo',1),
                  anova_amostras('Hevein_Lus100000453_7','9',f'{Bellow}','Rhizo',1),
                  anova_amostras('LTP_Lus10026418_F9','9',f'{Bellow}','Rhizo',1),
                  anova_amostras('Snakin_Lus10017212_A11','9',f'{Bellow}','Rhizo',1),
                  anova_amostras('Snakin_Lus10042203_3','9',f'{Bellow}','Rhizo',1)],

        'RHIZO_14': [anova_amostras('primer_9_hevein_Lus10028377','14',f'{Diff}','Rhizo',1),
                  anova_amostras('Hevein_Lus10006552_B11','14',f'{Bellow}','Rhizo',1),
                  anova_amostras('Hevein_Lus100000453_7','14',f'{Bellow}','Rhizo',1),
                  anova_amostras('LTP_Lus10026418_F9','14',f'{Bellow}','Rhizo',1),
                  anova_amostras('Snakin_Lus10017212_A11','14',f'{Bellow}','Rhizo',1),
                  anova_amostras('Snakin_Lus10042203_3','14',f'{Bellow}','Rhizo',1)]
        }
        
# creating a Dataframe object 
df = pd.DataFrame(details_stats)
df1 = pd.DataFrame(details)

df.to_csv('Above_stats.csv')
df.to_csv('Above_f1.csv')
  