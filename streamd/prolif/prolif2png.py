"""Convert ProLIF aggregated outputs to occupancy heatmap plots."""

import argparse
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
from plotnine import ggplot, geom_point, geom_text, aes, theme, element_text, element_blank, theme_bw, scale_color_manual, element_rect, scale_x_discrete
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

def convertprolif2png(plif_out_file, output=None, occupancy=0.6,
                      plot_width=None, plot_height=None,
                      point_size=3, base_size=12, show_percentage=True):
    '''

    :param plif_out_file:
    :param occupancy:
    :param plot_width:
    :param plot_height:
    :param point_size:
    :param base_size:
    :param show_percentage: if True, draw the occupancy percentage as a label above each dot
    :return:
    '''

    new_names = {"HBACCEPTOR": "A", "ANIONIC": "N", "HYDROPHOBIC": "H", 'METALACCEPTOR':'MeA',
                 "PISTACKING": "pi-s", "HBDONOR": "D", "PICATION": "pi+", "CATIONPI": "+pi", "CATIONIC": "P",
                 "XBDONOR": "XD", "XBACCEPTOR": "XA", "WATERBRIDGE": "W"}
    label_colors = {"hbacceptor": "red", "hbdonor": "forestgreen", "anionic": "blue", "cationic": "magenta",
                    "hydrophobic": "orange", "pication": "black",
                    "cationpi": "darkblue", "pistacking": 'darkslategray', 'metalacceptor': 'cyan',
                    "xbdonor": "blueviolet", "xbacceptor": "saddlebrown", "waterbridge": "teal"}

    df = pd.read_csv(plif_out_file, sep='\t')
    id_columns = ['Name','directory'] if 'directory' in df.columns else ['Name']
    # keep the exact occupancy fraction; rounding here to 1 decimal (0.1) collapsed all
    # contacts below ~5% occupancy to 0.0 (dropped by the value > 0 filter) and forced
    # every percentage label to a multiple of 10
    df_occup = df.drop(['Frame'], axis=1).groupby(id_columns).apply(lambda x: x.sum() / len(x))

    subdf = pd.melt(df_occup.reset_index(), id_vars=id_columns)
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

    if show_percentage:
        def _format_percentage(value):
            pct = value * 100
            # normally show a whole-number percentage; when no occupancy cutoff is applied,
            # keep 2 decimals for sub-1% contacts so rare interactions are not shown as "0%"
            if occupancy == 0 and round(pct) == 0:
                return f'{pct:.2f}%'
            return f'{round(pct)}%'

        subdf['percentage'] = subdf['value'].apply(_format_percentage)

    plot = (ggplot(subdf)+ geom_point(aes( x="variable", y="Name",
                                           color = "interaction",
                                          ), size=point_size)
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

    if show_percentage:
        plot = plot + geom_text(aes(x="variable", y="Name", label="percentage"),
                                size=base_size * 0.55, nudge_y=0.18, color="black",
                                show_legend=False)

    if output is None:
        output_name = os.path.join(os.path.dirname(plif_out_file),
                               f"{os.path.basename(plif_out_file).rstrip('.csv')}_occupancy{occupancy}.png")
    else:
        output_name = os.path.splitext(output)[0]

    if plot_width and plot_height:
        plot.save(output_name, width=plot_width, height=plot_height, dpi=300, verbose=False)
    else:
        plot.save(output_name, dpi=300, verbose=False)

def main():
    """CLI for converting ProLIF CSV outputs into PNG plots."""
    parser = argparse.ArgumentParser(description='''Draw prolif plot for analysis binding mode of multiple ligands''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, nargs='+',
                        help='input file with compound. Supported formats: *.csv')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output prefix name')
    parser.add_argument('--occupancy', metavar='float', default=0.6, type=float,
                        help='minimum occupancy of the unique contacts to show')
    parser.add_argument('--width', metavar='int', default=15, type=int,
                        help='width of the output picture')
    parser.add_argument('--height', metavar='int', default=10, type=int,
                        help='height of the output picture')
    parser.add_argument('--base_size', metavar='int', default=12, type=int,
                        help='base size of the output picture')
    parser.add_argument('--point_size', metavar='int', default=3, type=int,
                        help='dots size of the output picture')
    parser.add_argument('--no-show_percentage', default=False, action='store_true',
                        help='do not show the occupancy percentage label above each dot '
                             '(percentages are shown by default)')

    args = parser.parse_args()

    for input_file in args.input:
        convertprolif2png(input_file, output=args.output, occupancy=args.occupancy,
                          plot_width=args.width, plot_height=args.height,
                          point_size=args.point_size, base_size=args.base_size,
                          show_percentage=not args.no_show_percentage)


if __name__ == '__main__':
    main()
