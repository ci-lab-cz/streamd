import logging
import math
import subprocess
import time

from dask.distributed import Client
from rdkit import Chem


def init_dask_cluster(n_tasks_per_node, ncpu, hostfile=None):
    '''

    :param n_tasks_per_node: number of task on a single server
    :param ncpu: number of cpu on a single server
    :param hostfile:
    :return:
    '''
    if hostfile:
        with open(hostfile) as f:
            hosts = [line.strip() for line in f]
            n_servers = sum(1 if line.strip() else 0 for line in f)
    else:
        n_servers = 1

    n_workers = n_servers * n_tasks_per_node
    n_threads = math.ceil(ncpu / n_tasks_per_node)
    if hostfile is not None:
        cmd = f'dask ssh --hostfile {hostfile} --nworkers {n_workers} --nthreads {n_threads} &'
        subprocess.check_output(cmd, shell=True)
        time.sleep(10)
        dask_client = Client(hosts[0] + ':8786', connection_limit=2048)
    else:
        dask_client = Client(n_workers=n_workers, threads_per_worker=n_threads)  # to run dask on a single server

    dask_client.forward_logging(level=logging.INFO)
    return dask_client


def calc_dask(func, main_arg, dask_client, dask_report_fname=None, **kwargs):
    main_arg = iter(main_arg)
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    if dask_client is not None:
        from dask.distributed import as_completed, performance_report
        # https://stackoverflow.com/a/12168252/895544 - optional context manager
        from contextlib import contextmanager
        none_context = contextmanager(lambda: iter([None]))()
        with (performance_report(filename=dask_report_fname) if dask_report_fname is not None else none_context):
            nworkers = len(dask_client.scheduler_info()['workers'])
            futures = []
            for i, arg in enumerate(main_arg, 1):
                futures.append(dask_client.submit(func, arg, **kwargs))
                if i == nworkers:  # you may submit more tasks then workers (this is generally not necessary if you do not use priority for individual tasks)
                    break
            seq = as_completed(futures, with_results=True)
            for i, (future, results) in enumerate(seq, 1):
                yield results
                del future
                try:
                    arg = next(main_arg)
                    new_future = dask_client.submit(func, arg, **kwargs)
                    seq.add(new_future)
                except StopIteration:
                    continue
