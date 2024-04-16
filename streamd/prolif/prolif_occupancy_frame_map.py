import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from plotnine import (ggplot, geom_point, aes, theme, element_text, element_blank,
                      theme_bw, scale_color_manual, element_rect, facet_wrap, labs, scale_x_continuous, element_line,facet_grid )
plt.ioff()

def convertplif2png(plif_out_file, occupancy=0, plot_width=15, plot_height=5):

    label_colors = {"hbacceptor": "red", "hbdonor": "forestgreen", "anionic": "blue", "cationic": "magenta",
                    "hydrophobic": "orange", "pication": "black", "pistacking": 'black', 'metalacceptor': 'pink'}

    df = pd.read_csv(plif_out_file, sep='\t')
    subdf = pd.melt(df, id_vars=['Frame'])
    subdf['residue'] = subdf['variable'].apply(lambda x: '.'.join(x.split('.')[:-1]).upper())
    subdf['interaction'] = subdf['variable'].apply(lambda x: x.split('.')[-1])#.replace(new_names, regex=True)
    subdf['color'] = subdf['interaction'].apply(lambda x: label_colors.get(x, 'grey'))
    occupancy_df = subdf.groupby('variable').apply(lambda x: round(x['value'].sum() / len(x), 1))
    occupancy_30 = occupancy_df[occupancy_df >= occupancy]
    contact_list_within_occupancy = occupancy_30.index.to_list()
    subdf_occupancy = subdf[(subdf['variable'].isin(contact_list_within_occupancy)) & (subdf['value'])]
    # subdf_occupancy = subdf[(subdf['variable'].isin(contact_list_within_occupancy))]

    subdf_occupancy.loc[:,'Time (ns)'] = subdf_occupancy.loc[:, 'Frame'].apply(lambda x: x/100)
    plot = (ggplot(subdf_occupancy, aes(y="interaction",x="Time (ns)", color="interaction"))+ theme_bw(12)+
            geom_point(alpha=0.9, shape='|', size=5) +
    theme(panel_background = element_rect(fill="white", color="white"),
          panel_grid = element_blank(),
          panel_grid_major = element_line(linetype = 'solid', colour = "white"),
          strip_background = element_rect(fill="white", colour="white"),
          strip_text_y=element_text(angle = 0),
          axis_text_y=element_blank(),
          axis_ticks_y=element_blank(),
          legend_title = element_text(), legend_text = element_text(),
          legend_position="bottom",
          legend_key=element_rect(color="white"),
          legend_box_spacing=0,
          )+
        # scale_x_continuous(breaks = range(min(subdf_occupancy['Frame']),
        #                max(subdf_occupancy['Frame']), 1000))
            facet_grid('residue', scales = 'free_y')+
        labs(y='', title='')+ scale_color_manual(values = label_colors, na_value="white"))

    output_name = os.path.join(os.path.dirname(plif_out_file),
                               f"{os.path.basename(plif_out_file).strip('.csv')}_occupancy{occupancy}.png")

    if plot_width and plot_height:
        plot.save(output_name, width=plot_width, height=plot_height, dpi=300, verbose=False)
    else:
        plot.save(output_name, dpi=300, verbose=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Returns the formal charge for the molecule using RDKiT''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.csv')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.csv')
    args = parser.parse_args()

    convertplif2png(args.input)