import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.ioff()

sns.set_theme(style="whitegrid", palette="pastel", font_scale=2)
# SMALL_SIZE = 15
# MEDIUM_SIZE = 25
# BIGGER_SIZE = 36

# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def convertplif2png(plif_out_file, occupancy=0.1, plot_width=15, plot_height=5):

    plt.figure(figsize=(plot_width, plot_height))
    plt.xticks(rotation=90)
    new_names = {"hbacceptor": "A", "anionic": "N", "hydrophobic": "H",
                 "pistacking": "pi-s", "hbdonor": "D", "pication": "pi+", "cationic": "P"}

    label_colors = {"A": "red", "D": "forestgreen", "N": "blue", "P": "magenta",
                    "H": "orange", "pi+": "black", "pi-s": 'black'}

    df = pd.read_csv(plif_out_file, sep='\t')
    subdf = pd.melt(df, id_vars=['Frame'])
    subdf['residue'] = subdf['variable'].apply(lambda x: '.'.join(x.split('.')[:-1]))
    subdf['interaction'] = subdf['variable'].apply(lambda x: x.split('.')[-1])

    df_occupancy = subdf.groupby('variable').apply(lambda x: round(x['value'].sum() / len(x), 1))
    contact_list_within_occupancy =  df_occupancy[df_occupancy >= occupancy].index.to_list()

    subdf_occupancy = subdf[(subdf['variable'].isin(contact_list_within_occupancy)) & (subdf['value'])]

    c = [label_colors.get(i.split('.')[-1], 'grey') for i in subdf['variable']]
    plot = plt.scatter(subdf['variable'], subdf['fname'], c=c, s=85)
    # plt.xticks(subdf['variable'], labels=subdf['variable'], rotation=90)
    # plot = plt.plot(subdf['variable'], subdf['fname'], color='black', marker='o', linewidth=2, markersize=5)

    # for ticklabel, color in zip(plt.gca().get_xticklabels(), c):
    #     ticklabel.set_color(color)
    #     ticklabel.set_rotation(90)

    output_name = os.path.join(os.path.dirname(plif_out_file), f"{os.path.basename(plif_out_file).strip('.csv')}_occupancy{occupancy}.png")
    plot.figure.savefig(output_name, bbox_inches="tight")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Returns the formal charge for the molecule using RDKiT''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.csv')
    args = parser.parse_args()

    convertplif2png(args.input)