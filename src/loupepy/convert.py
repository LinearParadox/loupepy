from anndata import AnnData # type: ignore
import pandas as pd
import numpy as np
import os
import typing
from loupepy.utils import _validate_anndata, _validate_obs
import h5py # type: ignore

def create_loupe(anndata: AnnData, output_file:str, layer: str | None = None, tmp_file: str="tmp.h5ad", 
                 loupe_converter_path: str | None = None, dims: list[str] | None = None, 
                 obs_keys: list[str] | None=None) -> None:
    ''''
    Creates a temp h5 file and calls the loupe converter executable for the conversion
    '''
    if not os.path.exists(os.path.basename(output_file)):
        raise ValueError('Output file does not exist')
    if layer is None:
        _validate_anndata(anndata)
    else:
        _validate_anndata(anndata, layer)
    if obs_keys:
        obs = anndata.obs.loc[:, obs_keys]
        obs = _validate_obs(obs)
    else:
        obs = anndata.obs
        obs = _validate_obs(obs)

    with h5py.File(tmp_file, 'w') as f:
        pass