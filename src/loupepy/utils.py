from pandas import Series, DataFrame
from typing import Union
from numpy.typing import ArrayLike
from anndata import AnnData  # type: ignore
import numpy as np
import scipy.sparse as sp
import os
import platform
import warnings

def _validate_barcodes(barcodes: Series) -> bool:
    barcode_regex = "^(.*[:_])?([ACGT]{14,})([:_].*)?$"
    barcode_gem_regex = "^(.*[:_])?([ACGT]{14,})-(\\d+)([:_].*)?$"
    visium_hd_regex = "^(.*[:_])?(s_\\d{3}um_\\d{5}_\\d{5})([:_].*)?$"
    visium_hd_gem_regex = "^(.*[:_])?(s_\\d{3}um_\\d{5}_\\d{5})-(\\d+)([:_].*)?$"
    xenium_cell_id_regex = "^(.*[:_])?([a-p]{1,8})-(\\d+)([:_].*)?$"
    for n in [barcode_gem_regex, barcode_regex, visium_hd_regex, visium_hd_gem_regex, xenium_cell_id_regex]:
        if not barcodes.str.fullmatch(n).all():
            return False
    return True

def _validate_counts(mat: Union[ArrayLike, sp.spmatrix])  -> None:
    if sp.issparse(mat) and not bool(np.isnan(mat.data).any() or np.isinf(mat.data).any()):
        raise ValueError('Counts matrix contains NaN! This is not compatible with loupe converter!')
    if bool(np.isnan(mat).any() or np.isinf(mat).any()): # type: ignore
        #ignored due to weird ufunc mypy error
        raise ValueError('Counts matrix contains inf! This is not compatible with loupe converter!')
    
def _validate_anndata(anndata: AnnData, layer: str | None = None) -> None:
    if not isinstance(anndata, AnnData):
        raise ValueError('Input is not an AnnData object!')
    if not _validate_barcodes(anndata.obs.index):
        raise ValueError('Barcodes do not match the format required for loupeconverter!')
    if anndata.obs.n_obs == 0:
        raise ValueError('No observations found in the anndata object!')
    if anndata.var.n_vars == 0:
        raise ValueError('No vars found in the anndata object!')
    if layer is None:
        _validate_counts(anndata.X)
    else:
        _validate_counts(anndata.layers[layer])
    if anndata.var.index.equals("").any():
        raise ValueError('Empty vars found in the anndata object!')
    if anndata.obs.index.duplicated().any():
        raise ValueError('Duplicate barcodes found in the anndata object!')
    if anndata.var.index.duplicated().any():
        raise ValueError('Duplicate vars found in the anndata object!')
    
def _get_loupe_path() -> str:
    '''
    Returns the path to the loupeR binary
    '''
    if platform.system().startswith('linux'):
        path = os.environ['HOME'] + '/.local/bin/loupe/louper'
    elif platform.system().startswith('win'):
        path = os.environ['LOCALAPPDATA'] + '/Programs/louper/louper.exe'
    elif platform.system().startswith('darwin') or platform.system().startswith('apple'):
        path = os.environ['HOME'] + '/Applications/loupe/louper.app'
    else:
        raise OSError('Operating system not supported')
    return path

def _validate_obs(obs: DataFrame) -> DataFrame:
    for col in obs.columns:
        if not obs[col].dtype == 'category':
            obs[col] = obs[col].astype('category')
            warnings.warn(f'Column {col} was not a category and has been converted to a category')
        if len(obs[col].cat.categories) > 32768:
            warnings.warn(f'Column {col} has more than 32768 categories, skipping')
            obs.drop(col, axis=1, inplace=True)
    return obs
