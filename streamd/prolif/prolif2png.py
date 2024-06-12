import argparse
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
from plotnine import ggplot, geom_point, aes, theme, element_text, element_blank, theme_bw, scale_color_manual, element_rect, scale_x_discrete
plt.ioff()

# def calculate_figure_size(num_data_points_x, num_data_points_y):
#     # aspect_ratio = (num_data_points_x/num_data_points_y)/10
#     # # aspect_ratio = 16/12
#     # pixels_per_point_x = num_data_points_x*10
#     # pixels_per_point_y = num_data_points_y*10
#     # width = (pixels_per_point_x * aspect_ratio) ** 0.5
#     # height = (pixels_per_point_y / aspect_ratio) ** 0.5
#     pixels_per_point_x = 53
#     pixels_per_point_y = 40
#
#     width = num_data_points_x * pixels_per_point_x / 300
#     height = num_data_points_y * pixels_per_point_y / 300
#     return (width, height)

# def align_by_spaces(var):
#     split_row = var.split('.')
#     chain = f'.{split_row[1]}' if len(split_row) > 2 else ''
#     # resid_chain = f' {split_row[0]}{chain}'.ljust(10, '*')
#     #
#     # contact = split_row[-1].rjust(7,'_')
#     # return f'{resid_chain} {contact}'
#     return f'{split_row[0]}{chain} {split_row[-1]}'

def convertprolif2png(plif_out_file, occupancy=0.6, plot_width=None, plot_height=None, base_size=12):

    new_names = {"HBACCEPTOR": "A", "ANIONIC": "N", "HYDROPHOBIC": "H", 'METALACCEPTOR':'MeA',
                 "PISTACKING": "pi-s", "HBDONOR": "D", "PICATION": "pi+", "CATIONPI": "+pi", "CATIONIC": "P"}
    label_colors = {"hbacceptor": "red", "hbdonor": "forestgreen", "anionic": "blue", "cationic": "magenta",
                    "hydrophobic": "orange", "pication": "black",
                    "cationpi": "darkblue", "pistacking": 'darkslategray', 'metalacceptor': 'cyan'}

    df = pd.read_csv(plif_out_file, sep='\t')
    df_occup = df.drop('Frame', axis=1).groupby('Name').apply(lambda x: round(x.sum() / len(x), 1))

    subdf = pd.melt(df_occup.reset_index(), id_vars=['Name'])
    subdf = subdf[(subdf['value'] >= occupancy) & (subdf['value'] > 0)]

    subdf['resi'] = subdf['variable'].apply(lambda x: int(re.findall('[0-9]+', x)[0]))
    subdf['chain'] = subdf['variable'].apply(lambda x: x.split('.')[1])
    #subdf = subdf.sort_values(by=['chain','resi'])

    subdf['interaction'] = subdf['variable'].apply(lambda x: x.split('.')[-1])

    subdf['variable'] = subdf['variable'].str.upper().replace(new_names, regex=True)
    # subdf['variable_fill'] = subdf['variable'].apply(align_by_spaces)
    # reorder for plotnine
    subdf = subdf.sort_values(by=['chain','resi', 'interaction']).reset_index(drop=True)
    subdf['variable'] = pd.Categorical(subdf.variable, categories=pd.unique(subdf.variable))

    plot = (ggplot(subdf)+ geom_point(aes( x="variable", y="Name", color = "interaction"))
        + theme_bw(base_size=base_size) + scale_color_manual(values = label_colors)+
            theme(axis_text_x=element_text(rotation=90, hjust=0.5),
                axis_text_y=element_text(vjust=0.5),
                  axis_title_y=element_blank(),
                axis_title_x=element_blank(),
                legend_position = "bottom",
                panel_border = element_blank(),
                panel_background = element_blank(),
                legend_key=element_rect(color = "white"),
                  legend_box_spacing=0,
    )
            )

    output_name = os.path.join(os.path.dirname(plif_out_file),
                               f"{os.path.basename(plif_out_file).rstrip('.csv')}_occupancy{occupancy}.png")

    if plot_width and plot_height:
        plot.save(output_name, width=plot_width, height=plot_height, dpi=300, verbose=False)
    else:
        plot.save(output_name, dpi=300, verbose=False)

def main():
    parser = argparse.ArgumentParser(description='''Draw prolif plot for analysis binding mode of multiple ligands''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, nargs='+',
                        help='input file with compound. Supported formats: *.csv')
    parser.add_argument('--occupancy', metavar='float', default=0.6, type=float,
                        help='minimum occupancy of the unique contacts to show')
    parser.add_argument('--width', metavar='int', default=15, type=int,
                        help='width of the output picture')
    parser.add_argument('--height', metavar='int', default=10, type=int,
                        help='height of the output picture')
    parser.add_argument('--base_size', metavar='int', default=12, type=int,
                        help='base size of the output picture')

    args = parser.parse_args()

    for input_file in args.input:
        convertprolif2png(input_file, occupancy=args.occupancy,
                          plot_width=args.width, plot_height=args.height,
                          base_size=args.base_size)


if __name__ == '__main__':
    main()
