import pandas as pd
from tqdm import tqdm
import operator

'''
Documentation
This file is used to make a result3.csv file that matches WHO clades and tips from the subtree model
WHO_max finds the WHO clades related to each tip and by majority rule assigns that WHO clade to a subtree
percent_common gives the percentage of tips with the WHO_max out of the total number of tips from a subtree
pred is the prediction from the model's output
This is the simple version to be used in other languages such as R for easier loading of the csv file

INPUT:
Clade assignments tsv file relating the clades from the model to the WHO clades
Predictions csv file from the model with the subtree clades from the model and the labels and forecasts(predictions)
Metadata csv file that contains the subtree clades from the model and the tips contained in it separated by spaces
'''

who_df = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/clade_assignments/assignments/H3N2_h3n2.tsv', sep='\t')
predictions = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/RecentPrediction_Final.csv')
meta_df = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/output.txt', names=['name', 'tips'])

def add_tip(name, tips_dict):
    if name not in tips_dict.keys():
        tips_dict[name] = 1
    else:
        tips_dict[name] += 1

list_representatives = []
for pred_row in tqdm(predictions.iterrows(),total=len(predictions)):
    total_tips = {}
    row = pred_row[1]
    clade = row['Clade']
    true_label = row['Labels']
    pred_label = row['Forcasting_model']
    clade_tips = meta_df[meta_df['name'] == clade]['tips']
    clade_tips_list = clade_tips.values[0].split(' ')
    for tip in tqdm(clade_tips_list):
        who_clade = who_df[who_df['name'].values == tip]
        if(who_clade.empty):
            continue
        specific_clades = who_clade['cladesSpec'].real[0].split(',')
        for specific_clade in specific_clades:
            add_tip(specific_clade, total_tips)
    max_who_clade = max(total_tips.items(), key=operator.itemgetter(1))
    sum_who_clade = sum(total_tips.values())
    percent_common = max_who_clade[1] / sum_who_clade
    representative_clade = max_who_clade[0] if percent_common >= 0.4 else 'UNASSIGNED' # THRESHOLD
    list_representatives.append([clade, representative_clade, percent_common, pred_label])

out_df = pd.DataFrame(list_representatives, columns=['subtree', 'WHO_max', 'percent_common', 'pred'])
out_df['subtree'] = out_df['subtree'].astype(int)
out_df.to_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/result3.csv')