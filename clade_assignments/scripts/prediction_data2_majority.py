import pandas as pd
from tqdm import tqdm
import operator

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
    max_who_clades = [x for i, x in enumerate(total_tips.items()) if x == max(total_tips.items(), key=operator.itemgetter(1))]
    sum_who_clade = sum(total_tips.values())
    percent_common = [[i, x / sum_who_clade] for i, x in total_tips.items()]
    representative_clade = [x[0] for x in percent_common if x[1] >= 0.4] # THRESHOLD
    list_representatives.append([clade, representative_clade, percent_common, pred_label])

out_df = pd.DataFrame(list_representatives, columns=['subtree', 'WHO_max', 'percent_common', 'pred'])
out_df['subtree'] = out_df['subtree'].astype(int)
out_df.to_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/result4.csv')