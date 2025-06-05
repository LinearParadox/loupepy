from pandas import Series, DataFrame
from typing import Union
from numpy.typing import ArrayLike
from anndata import AnnData  # type: ignore
import numpy as np
from pathlib import Path
import scipy.sparse as sp
from .setup import _get_install_path
import os
import platform
import warnings

def _validate_barcodes(barcodes: Series) -> bool:
    """
    Validate the barcodes in the anndata object.
    Args:
        barcodes (Series): Barcodes to validate.
    Returns:
        bool: True if the barcodes are valid, False otherwise.
    """
    barcode_regex = r"^(.*[:_])?([ACGT]{14,})([:_].*)?$"
    barcode_gem_regex = r"^(.*[:_])?([ACGT]{14,})-(\d+)([:_].*)?$"
    visium_hd_regex = r"^(.*[:_])?(s_\d{3}um_\d{5}_\d{5})([:_].*)?$"
    visium_hd_gem_regex = r"^(.*[:_])?(s_\d{3}um_\d{5}_\d{5})-(\d+)([:_].*)?$"
    xenium_cell_id_regex = r"^(.*[:_])?([a-p]{1,8})-(\d+)([:_].*)?$"
    for n in [barcode_gem_regex, barcode_regex, visium_hd_regex, visium_hd_gem_regex, xenium_cell_id_regex]:
        if barcodes.str.fullmatch(n).all():
            return True
    return False

def _validate_counts(mat: Union[ArrayLike, sp.spmatrix])  -> None:
    '''
    Validate the counts matrix.
    Args:
        mat (Union[ArrayLike, sp.spmatrix]): Counts matrix to validate.
        Returns:
        None
            
        Raises:
        ValueError: If the counts matrix is not valid.
        '''
    if sp.issparse(mat) and bool(np.isnan(mat.data).any() or np.isinf(mat.data).any()):
        raise ValueError('Counts matrix contains NaN! This is not compatible with loupe converter!')
    elif not sp.issparse(mat) and bool(np.isnan(mat).any() or np.isinf(mat).any()): # type: ignore
        #ignored due to weird ufunc mypy error
        raise ValueError('Counts matrix contains inf! This is not compatible with loupe converter!')
    
def _validate_anndata(anndata: AnnData, layer: str | None = None) -> None:
    if not isinstance(anndata, AnnData):
        raise ValueError('Input is not an AnnData object!')
    if not _validate_barcodes(anndata.obs.index):
        raise ValueError('Barcodes do not match the format required for loupeconverter!')
    if anndata.n_obs == 0:
        raise ValueError('No observations found in the anndata object!')
    if anndata.n_vars == 0:
        raise ValueError('No vars found in the anndata object!')
    if layer is None:
        _validate_counts(anndata.X)
    else:
        _validate_counts(anndata.layers[layer])
    if (anndata.var.index == "").any().any():
        raise ValueError('Empty vars found in the anndata object!')
    if anndata.obs.index.duplicated().any():
        raise ValueError('Duplicate barcodes found in the anndata object!')
    if anndata.var.index.duplicated().any():
        raise ValueError('Duplicate vars found in the anndata object!')
    
def _get_loupe_path() -> Path:
    '''
    Returns the path to the default loupe-converter install location
    '''
    path = _get_install_path()
    path = path / 'loupe_converter'
    if not os.path.exists(path):
        raise ValueError('Loupe converter path does not exist')
    return path


def _validate_obs(obs: DataFrame) -> DataFrame:
    """
    Validate the obs dataframe.
    Args:
        obs (DataFrame): obs dataframe to validate.
    Returns:
        DataFrame: Validated obs dataframe with invalid columns dropped, and category columns converted to category dtype.
    """

    for col in obs.columns:
        if not obs[col].dtype == 'category':
            obs.drop(col, axis=1, inplace=True)
        elif len(obs[col].cat.categories) > 32768:
            warnings.warn(f'Column {col} has more than 32768 categories, skipping')
            obs.drop(col, axis=1, inplace=True)
    return obs
