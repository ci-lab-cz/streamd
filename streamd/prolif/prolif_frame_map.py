import argparse
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from plotnine import (ggplot, geom_point, aes, theme, element_text, element_blank,
                      theme_bw, scale_color_manual, element_rect, facet_wrap, labs, scale_x_continuous, element_line,facet_grid )
plt.ioff()

def convertplifbyframe2png(plif_out_file, plot_width=15, plot_height=10, occupancy=0, filter_only_hydrophobic=False, base_size=12):

    label_colors = {"hbacceptor": "red", "hbdonor": "forestgreen", "anionic": "blue", "cationic": "magenta",
                    "hydrophobic": "orange", "pication": "black", "cationpi": "darkblue",
                    "pistacking": 'darkslategray', 'metalacceptor': 'cyan'}

    df = pd.read_csv(plif_out_file, sep='\t')
    subdf = pd.melt(df, id_vars=['Frame'])

    subdf['resi'] = subdf['variable'].apply(lambda x: int(re.findall('[0-9]+', x)[0]))
    subdf['chain'] = subdf['variable'].apply(lambda x: x.split('.')[1])

    subdf['residue'] = subdf['variable'].apply(lambda x: '.'.join(x.split('.')[:-1]).upper())
    subdf['interaction'] = subdf['variable'].apply(lambda x: x.split('.')[-1])  # .replace(new_names, regex=True)
    subdf['color'] = subdf['interaction'].apply(lambda x: label_colors.get(x, 'grey'))

    occupancy_df = subdf.groupby('variable').apply(lambda x: round(x['value'].sum() / len(x), 1))
    contact_list_within_occupancy = occupancy_df[occupancy_df >= occupancy].index.to_list()
    subdf_occupancy = subdf.loc[(subdf['variable'].isin(contact_list_within_occupancy)) & (subdf['value']), :]
#    subdf = subdf[subdf['value']]
    if filter_only_hydrophobic:
        not_only_H_residue_mask = subdf_occupancy.groupby('residue').apply(lambda x: True if x['interaction'].unique().tolist() != ['hydrophobic'] else False)
        not_only_H_residue_list = not_only_H_residue_mask[not_only_H_residue_mask].index
        subdf_occupancy = subdf_occupancy.loc[subdf_occupancy['residue'].isin(not_only_H_residue_list),:]

    subdf_occupancy.loc[:,'Time, ns'] = subdf_occupancy.loc[:, 'Frame'].apply(lambda x: x/100)

    subdf_occupancy = subdf_occupancy.sort_values(by=['chain', 'resi', 'interaction']).reset_index(drop=True)
    subdf_occupancy['residue'] = pd.Categorical(subdf_occupancy.residue, categories=pd.unique(subdf_occupancy.residue))

    plot = (ggplot(subdf_occupancy, aes(y="interaction",x="Time, ns", color="interaction"))+ theme_bw(base_size)+
            geom_point(alpha=0.9, shape='|', size=5) +
    theme(panel_background = element_rect(fill="white", color="white"),
          panel_grid = element_blank(),
          panel_grid_major = element_line(linetype = 'solid', colour = "white"),
          strip_background = element_rect(fill="white", colour="white"),
          strip_text_y=element_text(angle = 0),
          axis_text_y=element_blank(),
          #axis_ticks_y=element_blank(),
          legend_title = element_text(), legend_text = element_text(),
          legend_position="bottom",
          legend_key=element_rect(color="white"),
          legend_box_spacing=0,
          )+
        # scale_x_continuous(breaks = range(min(subdf_occupancy['Frame']),
        #                max(subdf_occupancy['Frame']), 1000))
            facet_grid('residue', scales = 'free_y')+
        labs(y='', title='', x= '\nTime, ns')+ scale_color_manual(values = label_colors, na_value="white"))

    output_name = os.path.join(os.path.dirname(plif_out_file),
                               f"{os.path.basename(plif_out_file).strip('.csv')}_framemap")
    if occupancy > 0:
        output_name = f'{output_name}_occupancy{occupancy}'
    if filter_only_hydrophobic:
        output_name = f'{output_name}_filterNotHydrophobicyOnly'

    output_name = f'{output_name}.png'

    if plot_width and plot_height:
        plot.save(output_name, width=plot_width, height=plot_height, dpi=300, verbose=False)
    else:
        plot.save(output_name, dpi=300, verbose=False)

def main():
    parser = argparse.ArgumentParser(description='''Draw prolif plot for analysis of contacts by each frame of the unique ligand''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, nargs='*',
                        help='input file with prolif output for the unique molecule. Supported formats: *.csv.'
                             ' Ex: plifs.csv')
    parser.add_argument('--occupancy', metavar='float', default=0, type=float,
                        help='minimum occupancy of the unique contacts to show. Show all contacts by default')
    parser.add_argument('--filt_only_H', action='store_true', default=False,
                        help='filt residues where only hydrophobic contacts occur')
    parser.add_argument('--width', metavar='int', default=15, type=int,
                        help='width of the output picture')
    parser.add_argument('--height', metavar='int', default=10, type=int,
                        help='height of the output picture')
    parser.add_argument('--base_size', metavar='int', default=12, type=int,
                        help='base size of the output picture')
    args = parser.parse_args()
    for input_file in args.input:
        convertplifbyframe2png(input_file, plot_width=args.width, plot_height=args.height,
                               occupancy=args.occupancy, filter_only_hydrophobic=args.filt_only_H, base_size=args.base_size)


if __name__ == '__main__':
    main()