import pandas as pd
from tqdm import tqdm

who_df = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/clade_assignments/assignments/H3N2_h3n2.tsv', sep='\t')
predictions = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/RecentPrediction_Final.csv')
meta_df = pd.read_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/output.txt', names=['name', 'tips'])

true_tips = {}
total_tips = {}
def add_tip(name, tips_dict, value=1):
    if name not in tips_dict.keys():
        tips_dict[name] = value
    else:
        tips_dict[name] += value

for pred_row in tqdm(predictions.iterrows(),total=len(predictions)):
    row = pred_row[1]
    clade = row['Clade']
    true_label = row['Labels']
    pred_label = row['Forcasting_model']
    clade_tips = meta_df[meta_df['name'] == clade]['tips']
    clade_tips_list = clade_tips.values[0].split(' ')
    for tip in tqdm(clade_tips_list):
        add_tip(tip, total_tips)
        if(true_label == pred_label):
            add_tip(tip, true_tips)

who_tips_sum = {}
who_tips_nonzero = {}
who_tips_total = {}
tips_in_who_clade = {}
for pred_row in tqdm(predictions.iterrows(),total=len(predictions)):
    row = pred_row[1]
    clade = row['Clade']
    clade_tips = meta_df[meta_df['name'] == clade]['tips']
    clade_tips_list = clade_tips.values[0].split(' ')
    for tip in tqdm(clade_tips_list):
        who_clade = who_df[who_df['name'].values == tip]
        if(who_clade.empty):
            continue
        try:
            num_trues = true_tips[tip]
        except:
            num_trues = 0
        ratio = num_trues / total_tips[tip]
        specific_clades = who_clade['cladesSpec'].real[0].split(',')
        for specific_clade in specific_clades:
            add_tip(specific_clade, tips_in_who_clade, [[tip, total_tips[tip], ratio]])
            add_tip(specific_clade, who_tips_total)
            add_tip(specific_clade, who_tips_sum, ratio)
            if(ratio > 0):
                add_tip(specific_clade, who_tips_nonzero)

for key, value in who_tips_total.items():
    try:
        who_tips_nonzero[key] = who_tips_nonzero[key] / who_tips_total[key]
    except KeyError:
        who_tips_nonzero[key] = 0

# Create csv
csv_list = []
for key, value in tips_in_who_clade.items():
    for tip, size, ratio in value:
        csv_list.append([key, tip, size, ratio])

out_df = pd.DataFrame(csv_list, columns=['WHO', 'tip', 'size', 'ratio'])
out_df.to_csv('/Users/ivan/Documents/repos/hs19-flu/visualization/result2.csv')