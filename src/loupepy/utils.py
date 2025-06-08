from pandas import Series, DataFrame
from typing import Union, List
from numpy.typing import ArrayLike
from anndata import AnnData  # type: ignore
import numpy as np
from pathlib import Path
import scipy.sparse as sp
import logging
from scipy.sparse import csc_matrix
from .setup import _get_install_path
import os
from numpy import ndarray
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


def _validate_obs(obs: DataFrame, strict=False) -> DataFrame:
    """
    Validate the obs dataframe.
    Args:
        obs (DataFrame): obs dataframe to validate.
    Returns:
        DataFrame: Validated obs dataframe with invalid columns dropped, and category columns converted to category dtype.
    """

    for col in obs.columns:
        if not obs[col].dtype == 'category':
            if strict:
                raise ValueError(f'Column {col} is not categorical, which is required for Loupe. '
                                 f'Please check that this is truly categorical data.')
            obs.drop(col, axis=1, inplace=True)
            logging.warning(f'Column {col} is not categorical, dropping from final obs dataframe.')
        elif len(obs[col].cat.categories) > 32768:
            if strict:
                raise ValueError(f'Column {col} has more than 32768 categories, which '
                                 f'is not supported by Loupe. Please check that this is truly categorical data.')
            warnings.warn(f'Column {col} has more than 32768 categories, skipping')
            obs.drop(col, axis=1, inplace=True)

def _validate_obsm(obsm: dict[str, ndarray], obsm_keys: list[str]|None = None, strict: bool = False) -> list[str]:
    """
    Validate the obsm dictionary.
    Args:
        obsm (dict[str, ndarray]): obsm dictionary to validate.
        strict (bool): If True, will raise an error if any of the arrays are not 2D.
    Returns:
        list[str] : List of valid keys in the obsm dictionary.
    """
    valid_keys = []
    if obsm_keys is None:
        obsm_keys = list(obsm.keys())
    for key in obsm_keys:
        if not isinstance(obsm[key], np.ndarray):
            if strict:
                raise ValueError(f'Obsm key {key} has invalid type {type(obsm[key])}. Must be a numpy array.')
            warnings.warn(f'Obsm key {key} has invalid type {type(obsm[key])}. Dropping from output.')
        elif obsm[key].shape[1] != 2:
            if strict:
                raise ValueError(f'Obsm key {key} has invalid shape {obsm[key].shape}. '
                                 f'Must be an array with shape (n_cells, 2).')
            warnings.warn(f'Obsm key {key} has invalid shape {obsm[key].shape}. Dropping from output.')
        else:
            valid_keys.append(key)
    return valid_keys

def get_count_matrix(anndata: AnnData, layer: str | None = None) -> csc_matrix:
    """
    Get the counts matrix from an AnnData object in the format for loupe converter.
    Args:
        anndata (AnnData): AnnData object to get the counts matrix from.
        layer (str | None): Layer to get the counts matrix from. If None, will use the X attribute.
    Returns:
        csc_, sparse matrix: Counts matrix in the format for loupe converter.
    """
    if layer is None:
        return csc_matrix(anndata.X.T)
    else:
        return csc_matrix(anndata.X)

def get_obs(anndata: AnnData, obs_keys: List[str]|None = None, strict: bool = False) -> DataFrame:
    """
    Get the obs dataframe from an AnnData object in the format for loupe converter.
    Args:
        anndata (AnnData): AnnData object to get the obs dataframe from.
        obs_keys (str | None): Keys to subset the obs dataframe. If None, will use all valid keys.
        strict (bool): If True, will raise an error if any of the columns are not categorical.
    Returns:
        DataFrame: Obs dataframe in the format for loupe converter.
    """
    obs = anndata.obs.copy()
    if obs_keys:
        obs = obs.loc[:,obs_keys]
    _validate_obs(obs, strict)
    return obs

def get_obsm(anndata: AnnData, obsm_keys: List[str] | None = None, strict: bool = False) -> dict[str, ndarray]:
    """
    Get the obsm dictionary from an AnnData object in the format for loupe converter.
    Args:
        anndata (AnnData): AnnData object to get the obsm dictionary from.
        obsm_keys (str | None): Keys to subset the obsm dictionary. If None, will use all valid keys.
        strict (bool): If True, will raise an error if any of the arrays are not 2D. If false, will drop invalid keys.
    Returns:
        dict[str, ndarray]: Obsm dictionary in the format for loupe converter.
    """
    obsm = anndata.obsm.copy()
    valid_keys = _validate_obsm(obsm, obsm_keys, strict)
    return {key: obsm[key] for key in valid_keys}
