import pandas as pd
from tqdm import tqdm

who_df = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/clade_assignments/assignments/H3N2_h3n2.tsv', sep='\t')
predictions = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/RecentPrediction_Final.csv')
meta_df = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/output.txt', names=['name', 'tips'])

true_tips = {}
total_tips = {}
def add_tip(name, tips_dict):
    if name not in tips_dict.keys():
        tips_dict[name] = 1
    else:
        tips_dict[name] += 1

for pred_row in tqdm(predictions.iterrows(),total=len(predictions)):
    row = pred_row[1]
    clade = row['Clade']
    true_label = row['Labels']
    pred_label = row['Forcasting_model']
    clade_tips = meta_df[meta_df['name'] == clade]['tips']
    clade_tips_list = clade_tips.values[0].split(' ')
    for tip in tqdm(clade_tips_list):
        who_clade = who_df[who_df['name'].str.contains(tip)]
        if(who_clade.empty):
            continue
        specific_clades = who_clade['cladesSpec'].real[0].split(',')
        for specific_clade in specific_clades:
            add_tip(specific_clade, total_tips)
            if(true_label == pred_label):
                add_tip(specific_clade, true_tips)

out_frame =[]
for pred_row in predictions.iterrows():
    row = pred_row[1]
    clade = row['Clade']
    clade_tips = meta_df[meta_df['name'] == clade]['tips']
    clade_tips_list = clade_tips.values[0].split(' ')
    for tip in clade_tips_list:
        try:
            num_trues = true_tips[tip]
        except:
            num_trues = 0
        ratio = num_trues / total_tips[tip]
        out_frame.append([tip, ratio, total_tips[tip]])
out_df = pd.DataFrame(out_frame, columns=['tip_name', 'ratio', 'size'])
out_df.to_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/result1.csv')