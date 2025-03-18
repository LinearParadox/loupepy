import re
from pandas import Series
from typing import List, Union
from numpy.typing import ArrayLike
from anndata import AnnData  # type: ignore
import numpy as np
import scipy.sparse as sp

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

def _validate_counts(mat: Union[ArrayLike, sp.spmatrix]) -> bool:
    if sp.issparse(mat):
        return bool(np.isnan(mat.data).any() or np.isinf(mat.data).any())
    else:
        return bool(np.isnan(mat).any() or np.isinf(mat).any()) # type: ignore
    
def _validate_anndata(anndata: AnnData):
    if not isinstance(anndata, AnnData):
        raise ValueError('Input is not an AnnData object!')
    if not _validate_barcodes(anndata.obs.index):
        raise ValueError('Barcodes do not match the format required for loupeconverter!')
    if anndata.obs.n_obs == 0:
        raise ValueError('No observations found in the anndata object!')
    if anndata.var.n_vars == 0:
        raise ValueError('No vars found in the anndata object!')
    if not _validate_counts(anndata.X):
        raise ValueError('Counts matrix contains NaN or Inf values')
    if anndata.obs.index.equals("").any():
        raise ValueError('Empty barcodes found in the anndata object!')
    if anndata.var.index.equals("").any():
        raise ValueError('Empty vars found in the anndata object!')
    if anndata.obs.index.duplicated().any():
        raise ValueError('Duplicate barcodes found in the anndata object!')
    if anndata.var.index.duplicated().any():
        raise ValueError('Duplicate vars found in the anndata object!')