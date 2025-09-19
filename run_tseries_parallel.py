"""
Orchestrates the time-series synthesis pipeline.
- Loads configuration from YAML.
- Prepares datasets.
- Allocates matrices.
- Computes time-series.
- Saves results to disk.

author: Agesinaldo Silva
date: June, 2024
"""

from __future__ import annotations

import sys
import logging
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any

import numpy as np

# Import existing helpers from the original module.
# These functions are expected to exist in your current codebase.
from tseries_parallel import (
    pre_config,
    set_empty_matrix,
    create_matrix,
    save_table,
)
from util.store import Store


# ----------------------------- Dataclasses -----------------------------------

@dataclass(frozen=True)
class DatasetBundle:
    """Container for a single dataset and its metadata required by the pipeline."""
    doma: np.ndarray
    resp: np.ndarray
    respc: np.ndarray
    scale: float
    skr_sdr: Tuple[float, float]  # (skr, sdr)


@dataclass(frozen=True)
class Limits:
    """Holds vlim/mlim/rlim arrays produced by set_empty_matrix (for clarity)."""
    vlim: np.ndarray
    mlim: np.ndarray
    rlim: np.ndarray


@dataclass(frozen=True)
class RunConfig:
    """All inputs produced by pre_config, consolidated for readability."""
    spectrum: Any
    offset: float
    x_size: int
    central_range: np.ndarray
    data_set_list: List[Dict[str, np.ndarray]]
    respc_list: List[np.ndarray]
    scale_list: List[float]
    par_map_list: List[Tuple[float, float]]
    x_spec_full: np.ndarray
    conve_mod: Any
    tranx_num: int


@dataclass(frozen=True)
class CreateMatrixResult:
    """Result object for create_matrix so we don't return a bare int."""
    center_freq_pct: int


# ------------------------------ Adapters -------------------------------------

def _load_run_config(
    *, task_name: str, store: Store, config_file: str, output_path: str
) -> RunConfig:
    """
    Adapter that calls the legacy pre_config and packages its many outputs into
    a single RunConfig object.
    """
    (spectrum, offset, x_size, central_range,
     data_set_list, respc_list, scale_list, par_map_list,
     x_spec_full, conve_mod, tranx_num) = pre_config(
        task_name=task_name,
        store=store,
        config_file=config_file,
        output_path=output_path,
    )

    return RunConfig(
        spectrum=spectrum,
        offset=offset,
        x_size=x_size,
        central_range=central_range,
        data_set_list=data_set_list,
        respc_list=respc_list,
        scale_list=scale_list,
        par_map_list=par_map_list,
        x_spec_full=x_spec_full,
        conve_mod=conve_mod,
        tranx_num=tranx_num,
    )


def _call_create_matrix(
    *,
    spectrum: Any,
    x_size: int,
    central_range: np.ndarray,
    data_set: Dict[str, np.ndarray],
    x_spec_full: np.ndarray,
    mat_tseries: np.ndarray,
    offset: float,
    limits: Limits,
    vlimmax: float,
) -> CreateMatrixResult:
    """
    Wrapper around legacy create_matrix. Returns a small typed result object.
    NOTE: We are not altering the internal algorithm of the legacy function
    here (e.g., center-bin selection), only wrapping its return value so that
    upstream code is more readable.
    """
    center_freq_pct = create_matrix(
        SP=spectrum,
        x_size=x_size,
        central_range=central_range,
        data_set=data_set,
        x_spec_full=x_spec_full,
        mat_tseries=mat_tseries,
        offset=offset,
        lim=(limits.vlim, limits.mlim, limits.rlim),
        vlimmax=vlimmax,
    )
    return CreateMatrixResult(center_freq_pct=center_freq_pct)


# --------------------------- Refactored entrypoint ----------------------------

def run_tseries_parallel(
    *,
    store: Store,
    config_file: str,
    output_path: str,
    task_name: str = "tseries",
    logger: Optional[logging.Logger] = None,
) -> None:
    """
    Orchestrate the time-series synthesis pipeline for one or more datasets.

    Parameters
    ----------
    store : Store
        Backing store with configuration and (optionally) cached datasets.
    config_file : str
        YAML path with spectral and run parameters.
    output_path : str
        Directory to write outputs (per save_table semantics).
    task_name : str, default "tseries"
        When "tseries", loads h5 files from disk; otherwise uses `store`.
    logger : logging.Logger, optional
        If provided, progress goes to this logger instead of print().
    """
    log = logger or logging.getLogger(__name__)

    cfg = _load_run_config(
        task_name=task_name,
        store=store,
        config_file=config_file,
        output_path=output_path,
    )

    # Normalize inputs into explicit bundles so callsites are self-documenting.
    bundles: List[DatasetBundle] = [
        DatasetBundle(
            doma=ds["doma"],
            resp=ds["resp"],
            respc=rc,
            scale=sc,
            skr_sdr=pm,  # tuple: (skr, sdr)
        )
        for ds, rc, sc, pm in zip(
            cfg.data_set_list, cfg.respc_list, cfg.scale_list, cfg.par_map_list
        )
    ]

    for b in bundles:
        log.info("Loaded dataset | doma=%s resp=%s", b.doma.shape, b.resp.shape)

        # Allocate matrices/limits for this dataset.
        mat_tseries, (vlim, mlim, rlim), vlimmax, (probv, probm) = set_empty_matrix(
            SP=cfg.spectrum,
            offset=cfg.offset,
            doma_size=b.doma.shape,
            resp_size=b.resp.shape,
            respc=b.respc,
            scale=b.scale,
            tranx_num=cfg.tranx_num,
        )

        limits = Limits(vlim=vlim, mlim=mlim, rlim=rlim)

        # Compute the time-series matrices for the current dataset.
        cm_res = _call_create_matrix(
            spectrum=cfg.spectrum,
            x_size=cfg.x_size,
            central_range=cfg.central_range,
            data_set={"doma": b.doma, "resp": b.resp},
            x_spec_full=cfg.x_spec_full,
            mat_tseries=mat_tseries,
            offset=cfg.offset,
            limits=limits,
            vlimmax=vlimmax,
        )

        # Persist all outputs for this dataset.
        save_table(
            conve_mod=cfg.conve_mod,
            x_spec_full=cfg.x_spec_full,
            freq=cm_res.center_freq_pct,   # clearer name locally
            par=b.skr_sdr,                  # (skr, sdr)
            lim=(limits.vlim, limits.mlim, limits.rlim),
            vlimmax=vlimmax,
            prob=(probv, probm),
            offset=cfg.offset,
            scale=b.scale,
            output_path=output_path,
        )

if __name__ == '__main__':
   
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} /path/to/output  config.yaml'.format(sys.argv[0]))

    output_path = sys.argv[1]
    config_file = sys.argv[2]

    store = Store('','')
    run_tseries_parallel(
        store=store,
        config_file=config_file,
        output_path=output_path,
        task_name ="tseries"
    )