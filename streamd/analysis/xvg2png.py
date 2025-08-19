"""Convert GROMACS XVG output to CSV and PNG visualizations."""

import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd

from streamd.utils.utils import backup_prev_files


def convertxvg2png(xvg_file, system_name=None, transform_nm_to_A=False):
    """Parse an XVG file and save corresponding CSV and PNG plots.

    Parameters
    ----------
    xvg_file : str | os.PathLike
        Path to the GROMACS ``.xvg`` file.
    system_name : str | None, optional
        System name appended to the plot title.
    transform_nm_to_A : bool, optional
        If ``True`` convert values from nanometres to Ångström.
    """

    def check_if_value_found(value):
        """Return the first element of ``value`` or an empty string.

        Parameters
        ----------
        value : list[str]
            Regex match result.

        Returns
        -------
        str
            First matched value or an empty string.
        """
        if value:
            return value[0]
        else:
            return ''

    csv_file = xvg_file.replace('.xvg', '.csv')
    png_file = xvg_file.replace('.xvg', '.png')

    if os.path.isfile(csv_file):
        backup_prev_files(file_to_backup=csv_file)
    if os.path.isfile(png_file):
        backup_prev_files(file_to_backup=png_file)

    with open(xvg_file) as inp:
        data = inp.readlines()

    data_coord = [i for i in data if not (i.startswith('#') or i.startswith('@'))]
    title = check_if_value_found([i for i in data if ('@' in i and 'title' in i)]).replace('@    title ','').replace('"','')
    if system_name is not None:
        title = f'{title} {system_name} complex'
    subtitle = (check_if_value_found([i for i in data if ('@' in i and 'subtitle' in i)]).
                replace('@ subtitle ','').replace('"',''))
    legend_list = [i.replace(r'\s','').replace(r'\N','').split('legend')[1] for i in data if ('@ s' in i and 'legend' in i)]

    xaxis = [i for i in data if (i.startswith('@') and 'xaxis' in i)]
    yaxis = [i for i in data if (i.startswith('@') and 'yaxis' in i)]

    coords = [[float(j) for j in i.strip().split(' ') if j] for i in data_coord]

    xaxis = xaxis[0].strip().replace('@    xaxis  label ', '').replace('"', '') if xaxis else 'OX'
    yaxis = yaxis[0].strip().replace('@    yaxis  label ', '').replace('"', '') if yaxis else 'OY'

    plot1 = None
    plt.ioff()
    plt.rcParams.update({'font.size': 15})
    plt.figure(figsize=(15, 12))
    plt.title(title)
    plt.xlabel(xaxis)
    if transform_nm_to_A:
        yaxis_A = yaxis.replace('nm', 'Å')
        plt.ylabel(yaxis_A)
    else:
        plt.ylabel(yaxis)

    if coords and len(coords[0]) == 2:
        d = pd.DataFrame(coords, columns=[xaxis, yaxis])
        d.to_csv(csv_file, sep='\t', index=False)
        if transform_nm_to_A and 'nm' in yaxis.lower():
            #logging.warning(f'INFO: {xvg_file} nm ({yaxis}) values are converted in Angstrom ({yaxis_A})')
            plot1 = plt.plot(d[xaxis], d[yaxis].apply(lambda x: round(x*10,3)), marker='o', linewidth=2, markersize=5)
        else:
            plot1 = plt.plot(d[xaxis], d[yaxis], marker='o', linewidth=2, markersize=5)

    elif coords and len(coords[0]) > 2 and legend_list:
        d = pd.DataFrame(coords, columns=[xaxis]+legend_list)
        d.to_csv(csv_file, sep='\t', index=False)
        plt.legend(legend_list, loc='upper right', bbox_to_anchor=(0, 0), borderaxespad=-1)
        if transform_nm_to_A and 'nm' in yaxis.lower():
            #logging.warning(f'INFO: {xvg_file} nm ({yaxis}) values are converted in Angstrom ({yaxis_A})')
            plot1 = plt.plot(d[xaxis], d[legend_list].apply(lambda x: round(x*10,3)), marker='o', linewidth=2, markersize=3)
        else:
            plot1 = plt.plot(d[xaxis], d[legend_list], marker='o', linewidth=2, markersize=3)
    if plot1:
        plot1 = plot1[0]
        plt.legend(legend_list)
        plot1 = plot1.figure.suptitle(subtitle)
        plot1.figure.savefig(png_file, bbox_inches="tight")
    plt.clf()
    plt.close('all')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Draw png plot for xvg outputs of md analysis''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input file with compound. Supported formats: *.csv')
    parser.add_argument('--convert_to_A', action='store_true', default=False,
                        help='convert value from nm to A')
    args = parser.parse_args()

    convertxvg2png(args.input, transform_nm_to_A=args.convert_to_A)