"""Utilities for aggregating and plotting RMSD analysis results."""

import argparse
from datetime import datetime
from functools import partial
import os
import pandas as pd
import logging
import re

from streamd.analysis.plot_build import plot_rmsd_mean_std
from streamd.utils.utils import filepath_type


def _parse_system_name(system: str) -> tuple[str, int]:
    match = re.match(r"(.*)_replica(\d+)$", system)
    if match:
        return match.group(1), int(match.group(2))
    return system, 1


def _ensure_replica_cols(df: pd.DataFrame) -> pd.DataFrame:
    if 'protein_name' not in df.columns or 'replica' not in df.columns:
        protein_name, replica = _parse_system_name(str(df['system'].iloc[0]))
        df.loc[:, 'protein_name'] = protein_name
        df.loc[:, 'replica'] = replica
    return df


def merge_rmsd_csv(csv_files, out):
    """Concatenate multiple RMSD CSV files into one DataFrame.

    Parameters
    ----------
    csv_files : list[str]
        Paths to the CSV files with RMSD data to merge.
    out : str | os.PathLike
        Destination path where the merged table will be written.

    Returns
    -------
    pandas.DataFrame
        Combined DataFrame containing data from all ``csv_files``.
    """
    all_data_list = []
    csv_files.sort()
    for i in csv_files:
        data = pd.read_csv(i, sep='\t')
        data = _ensure_replica_cols(data)
        all_data_list.append(data)

    merged_data = pd.concat(all_data_list)
    merged_data.to_csv(out, sep='\t', index=False)
    return merged_data


def calc_mean_std_by_ranges_time(rmsd_data, time_ranges, rmsd_system='backbone', system_cols=['ligand_name','system']):
    """Return mean and standard deviation of RMSD per time range.

    Parameters
    ----------
    rmsd_data : pandas.DataFrame
        Table with at least ``time(ns)`` and RMSD columns.
    time_ranges : list[tuple[float, float]]
        Time intervals (in ns) over which to compute statistics.
    rmsd_system : str, optional
        Column name containing RMSD values to analyse.
    system_cols : list[str], optional
        Columns identifying unique systems in ``rmsd_data``.

    Returns
    -------
    pandas.DataFrame
        Mean and standard deviation of RMSD for each time range.
    """
    res_list = []
    for start, end in time_ranges:
        key = f'{start}-{end}ns'
        mean = rmsd_data.loc[(start <= rmsd_data['time(ns)']) & (rmsd_data['time(ns)'] <= end), system_cols+[rmsd_system]].groupby(
            system_cols).mean().reset_index().rename({rmsd_system: 'RMSD_mean'}, axis='columns').round(2)
        std = rmsd_data.loc[
            (start <= rmsd_data['time(ns)']) & (rmsd_data['time(ns)'] <= end), system_cols+ [rmsd_system]].groupby(
            system_cols).std().reset_index().rename({rmsd_system: 'RMSD_std'}, axis='columns').round(2)
        res_tmp = pd.merge(mean, std, on=system_cols)
        res_tmp.loc[:, 'rmsd_system'] = rmsd_system
        res_tmp.loc[:, 'time_range'] = key
        res_list.append(res_tmp)

    res = pd.concat(res_list)
    return res


def make_lower_case(df, cols):
    """Lower-case specified columns and warn on missing values.

    Parameters
    ----------
    df : pandas.DataFrame
        Input table whose columns will be normalised.
    cols : list[str]
        Column names to convert to lower case.

    Returns
    -------
    pandas.DataFrame
        The modified DataFrame.
    """
    for col in cols:
        if df[col].isna().any():
            logging.warning(f'The RMSD DataFrame column {col} contains the None value, '
                            f'which could potentially cause an issue for MD system identification. '
                            f'Check the RMSD Input/Output data for inconsistencies.')
        df[col] = df[col].astype(str).str.lower()
    return df


def run_rmsd_analysis(rmsd_files, wdir, unique_id, time_ranges=None,
                      rmsd_type_list=['backbone', 'ligand'], paint_by_fname=None,
                      title=None):
    """Execute RMSD analysis workflow and write summary files.

    Parameters
    ----------
    rmsd_files : list[str]
        Paths to input RMSD tables.
    wdir : str | os.PathLike
        Output directory for analysis results.
    unique_id : str
        Identifier appended to generated file names.
    time_ranges : list[tuple[float, float]] | None, optional
        Optional list of time intervals in nanoseconds.
    rmsd_type_list : list[str], optional
        RMSD column names to process from the input files.
    paint_by_fname : str | None, optional
        Optional CSV file defining colouring for the interactive plot.
    title : str | None, optional
        Custom title for the generated HTML plot.

    Returns
    -------
    None
        The function writes output files to ``wdir`` and returns ``None``.
    """
    if len(rmsd_files) > 1:
        rmsd_merged_data = merge_rmsd_csv(rmsd_files, os.path.join(wdir, f'rmsd_all_systems_{unique_id}.csv'))
    else:
        rmsd_merged_data = _ensure_replica_cols(pd.read_csv(rmsd_files[0], sep='\t'))

    base_cols = ['system', 'protein_name', 'replica']
    no_ligand = all(rmsd_merged_data['ligand_name'].isna())
    system_cols = base_cols if no_ligand else base_cols + ['ligand_name']

    #Use directory column as system column in case of replicate runs where ligand_name and system columns are the same,
    # but working_directories are different
    # important for paint_by file where all system columns should be presented and used
    if 'directory' in rmsd_merged_data.columns:
        max_num_unique_dirs = rmsd_merged_data.groupby(system_cols).apply(
            lambda x: len(x['directory'].unique()), include_groups=False).reset_index(drop=True).max()
        if max_num_unique_dirs > 1:
            system_cols.append('directory')

    lower_cols = [c for c in system_cols if c != 'replica']
    rmsd_merged_data = make_lower_case(rmsd_merged_data, cols=lower_cols)

    if time_ranges is None:
        start = rmsd_merged_data['time(ns)'].min()
        end = rmsd_merged_data['time(ns)'].max()
        mean = end // 2
        duration = end - start
        if duration > 1:
            #time_ranges = [(start, 1), (start, mean), (mean, end),  (end - 1, end), (start, end)]
            time_ranges = [ (start, end), (mean, end),  (end - 1, end)]
        else:
            time_ranges = [(start, end)]

        # df_mean_std = pd.DataFrame()
    mean_std_list = []
    for rmsd_system in rmsd_type_list:
        mean_std_data = calc_mean_std_by_ranges_time(rmsd_data=rmsd_merged_data,
                                                     time_ranges=time_ranges,
                                                     rmsd_system=rmsd_system,
                                                     system_cols=system_cols)
        mean_std_list.append(mean_std_data)

    df_mean_std = pd.concat(mean_std_list)
    df_mean_std.to_csv(os.path.join(wdir, f'rmsd_mean_std_time-ranges_{unique_id}.csv'), index=False, sep='\t')
    if paint_by_fname:
        if os.path.isfile(paint_by_fname):
            paint_by_data = pd.read_csv(paint_by_fname, sep='\t')
            merge_cols = ['protein_name'] if no_ligand else ['protein_name', 'ligand_name']
            if not all([i in paint_by_data.columns for i in merge_cols]):
                logging.warning(f'Cannot paint by custom column in html rmsd analysis. Missing columns: {merge_cols} in {paint_by_fname}')
                paint_by_col = 'rmsd_system'
                show_legend = False
            else:
                paint_by_data = make_lower_case(paint_by_data, cols=merge_cols)
                paint_by_col = [i for i in paint_by_data.columns if i not in merge_cols][0]
                df_mean_std = pd.merge(df_mean_std, paint_by_data.loc[:,
                                                    merge_cols+[paint_by_col]], on=merge_cols)
                show_legend = True

        else:
            raise FileNotFoundError(f'Cannot find file {paint_by_fname} for html painting')
    else:
        paint_by_col = 'replica' if 'replica' in df_mean_std.columns else 'rmsd_system'
        show_legend = False if paint_by_col == 'rmsd_system' else True
        df_mean_std[paint_by_col] = df_mean_std[paint_by_col].astype(str)

    plot_rmsd_mean_std(data=df_mean_std, paint_by_col=paint_by_col, show_legend=show_legend,
                        out_name=os.path.join(wdir, f'rmsd_mean_std_time-ranges_{unique_id}.html'),
                       title=title)


# def create_html_rmsd_mean_std(data, paint_by_col,show_legend, out_name):
#     plot_rmsd_mean_std(data=data, paint_by_col=paint_by_col,
#                        show_legend=show_legend,
#                        out_name=out_name)

def main():
    """CLI entry point for RMSD analysis.

    Parameters
    ----------
    None
    """
    parser = argparse.ArgumentParser(description='''Run rmsd analysis for StreaMD output files''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, nargs='+',
                        help='input file(s) with rmsd. Supported formats: *.csv. Required columns: '
                             'time(ns)\tbackbone\tligand\tligand_name\tsystem')
    parser.add_argument('--rmsd_type', metavar='backbone', required=False, nargs='+',
                        help='Column names in the input file to use for the analysis', default=['backbone', 'ligand'])
    parser.add_argument('--time_ranges', metavar='0-1 5-10 9-10', required=False, nargs='+',
                        help='Time ranges in ns. Default: start-end, middle-end, end-1ns', default=None)
    parser.add_argument('-d','--wdir', metavar='dirname', required=False,
                        type=partial(filepath_type, check_exist=False, create_dir=True),
                        help='Output files directory', default=None)
    parser.add_argument('--paint_by', default='',
                        help='File to paint by additional column. '
                             'Required columns: '
                             '- Protein-ligand simulations: protein_name\tligand_name.'
                             '- Protein only in water simulations: protein_name'
                             ' Sep: /\t. '
                             'The plot will be painted by any other than system and ligand_name column.')
    parser.add_argument('-o','--out_suffix', default=None,
                        help='Suffix for output files')
    parser.add_argument('--title', default=None,
                        metavar='RMSD Mean vs RMSD Std',
                        help='Title for html plot. Default: RMSD Mean vs RMSD Std')

    # add description about painting
    args = parser.parse_args()

    out_time = f'{datetime.now().strftime("%d-%m-%Y-%H-%M-%S")}'
    if args.time_ranges is not None:
        time_ranges = [(float(i.split('-')[0]),float(i.split('-')[1])) for i in args.time_ranges]
        logging.info(f"Will use the user's time ranges: {time_ranges}")
    else:
        time_ranges = None

    if args.wdir is None:
        wdir = os.getcwd()
    else:
        wdir = args.wdir

    if args.out_suffix:
        unique_id = f'{out_time}_{args.out_suffix}'
    else:
        unique_id = out_time

    run_rmsd_analysis(rmsd_files=args.input,
                      paint_by_fname=args.paint_by,
                      wdir=wdir, unique_id=unique_id, time_ranges=time_ranges,
                      rmsd_type_list=args.rmsd_type,
                      title=args.title)


if __name__ == '__main__':
    main()





