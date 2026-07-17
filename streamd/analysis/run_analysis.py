"""Utilities for aggregating and plotting RMSD analysis results."""

import argparse
from datetime import datetime
from functools import partial
import math
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


def _resolve_rmsd_columns(rmsd_files, rmsd_type_list):
    """Validate base columns and return the requested metrics that are usable.

    The base columns ``time(ns)``, ``system`` and ``ligand_name`` are structurally
    required and must be present in every input file. Requested RMSD metrics are
    handled leniently so a single incomplete system cannot abort an entire batch:
    a metric missing from *some* files is kept (those systems contribute empty
    values) with a warning, a metric missing from *all* files is dropped with a
    warning, and a ``ValueError`` is raised only when no requested metric survives.
    """
    if not rmsd_files:
        raise ValueError("At least one RMSD input file is required")

    if not rmsd_type_list:
        raise ValueError("At least one RMSD metric must be requested")

    required_base_columns = {'time(ns)', 'system', 'ligand_name'}
    per_file_columns = {}
    for rmsd_file in rmsd_files:
        columns = set(pd.read_csv(rmsd_file, sep='\t', nrows=0).columns)
        per_file_columns[rmsd_file] = columns
        missing_base = sorted(required_base_columns - columns)
        if missing_base:
            raise ValueError(
                f"{rmsd_file} is missing required base column(s): {', '.join(missing_base)}"
            )

    all_columns = set().union(*per_file_columns.values())
    usable_metrics = []
    for metric in rmsd_type_list:
        if metric not in all_columns:
            logging.warning(
                "Skipping RMSD metric '%s': it is not present in any input file.", metric)
            continue
        files_without_metric = sorted(
            rmsd_file for rmsd_file, columns in per_file_columns.items()
            if metric not in columns)
        if files_without_metric:
            logging.warning(
                "RMSD metric '%s' is missing from %d input file(s); the affected "
                "systems will have empty '%s' values: %s",
                metric, len(files_without_metric), metric, ', '.join(files_without_metric))
        usable_metrics.append(metric)

    if not usable_metrics:
        raise ValueError(
            "None of the requested RMSD metrics "
            f"({', '.join(rmsd_type_list)}) are present in the input file(s)"
        )
    return usable_metrics


def _default_time_ranges(start, end):
    """Build normalized default convergence windows from available trajectory bounds."""
    duration = end - start

    if duration <= 1:
        return [
            (
                round(float(start), 10),
                round(float(end), 10),
            )
        ]

    exact_midpoint = (start + end) / 2
    midpoint = math.ceil(exact_midpoint)

    if midpoint >= end:
        midpoint = math.floor(exact_midpoint)

    candidates = [
        (start, end),
        (midpoint, end),
        (end - 1.0, end),
    ]

    time_ranges = []
    for range_start, range_end in candidates:
        normalized_range = (
            round(float(range_start), 10),
            round(float(range_end), 10),
        )
        if normalized_range not in time_ranges:
            time_ranges.append(normalized_range)

    return time_ranges


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
    rmsd_type_list = _resolve_rmsd_columns(rmsd_files, rmsd_type_list)

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
        max_num_unique_dirs = rmsd_merged_data.groupby(system_cols)['directory'].nunique().max()
        if max_num_unique_dirs > 1:
            system_cols.append('directory')

    lower_cols = [c for c in system_cols if c != 'replica']
    rmsd_merged_data = make_lower_case(rmsd_merged_data, cols=lower_cols)

    if time_ranges is None:
        start = rmsd_merged_data['time(ns)'].min()
        end = rmsd_merged_data['time(ns)'].max()
        time_ranges = _default_time_ranges(start, end)

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
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, nargs='+',
                        help='Input RMSD TSV file(s). Required base columns: '
                             'time(ns), system, ligand_name. A column requested through '
                             '--rmsd_type that is missing from some files is skipped for '
                             'those systems (with a warning); one missing from all files '
                             'is dropped.')
    parser.add_argument('--rmsd_type', metavar='COLUMN', required=False, nargs='+',
                        help='RMSD column names to aggregate, for example '
                             'backbone, ligand, ligand_local, or ActiveSite5.0A',
                        default=['backbone', 'ligand', 'ligand_local', 'ActiveSite5.0A'])
    parser.add_argument('--time_ranges', metavar='0-1 5-10 9-10', required=False, nargs='+',
                        help='Time ranges in nanoseconds. Default: start-end, middle-end, and the final 1 ns.',
                        default=None)
    parser.add_argument('-d','--wdir', metavar='dirname', required=False,
                        type=partial(filepath_type, check_exist=False, create_dir=True),
                        help='Output files directory', default=None)
    parser.add_argument('--paint_by', default='',
                        help='File to paint by additional column. '
                             'Required columns: '
                             'protein-ligand simulations: protein_name, ligand_name; '
                             'protein-only simulations: protein_name. '
                             'The plot is coloured using the first column other than '
                             'protein_name and ligand_name.')
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



