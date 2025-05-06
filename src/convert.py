import os
import sys
from anndata import AnnData # type: ignore
import pandas as pd
from scipy.sparse import csr_matrix
import h5py # type: ignore
from typing import Any, Union, List
from array import array
from utils import _validate_anndata, _validate_obs, _get_loupe_path


def _create_string_dataset(obj: h5py.Group, key: str, strings: Union[List[str], pd.Series[str], str]) -> None:
    '''
    Creates a dataset in the h5 file with the given strings
    '''
    if len(strings) == 0:
        max_len = 1
    elif isinstance(strings, pd.Series):
        strings = strings.tolist()
        max_len = max(len(s) for s in strings)
    elif isinstance(strings, str):
        max_len = len(strings)
    else:
        max_len = max(len(s) for s in strings)
    dtype = h5py.string_dtype(encoding='ascii', length=max_len)
    d=obj.create_dataset(key, strings, dtype=dtype)
    d.close()

def _write_matrix(f: h5py.File, matrix: csr_matrix, features: pd.Series, barcodes: pd.Series, feature_ids: list[str]|pd.Series|None = None) -> None:
    '''
    Writes the matrix to the h5 file
    '''
    matrix_group = f.create_group('matrix')
    features_group = matrix_group.create_group('features')
    _create_string_dataset(matrix_group, 'barcodes', barcodes)
    matrix_group.create_dataset('data', data=matrix.data, dtype=int)
    matrix_group.create_dataset('indices', data=matrix.indices, dtype=int)
    matrix_group.create_dataset('indptr', data=matrix.indptr, dtype=int)
    matrix_group.create_dataset('shape', data=matrix.shape, dtype=int)
    matrix_group.close()
    _create_string_dataset(features, 'id', features)
    if feature_ids is None:
        feature_ids = [f"feature_{i}" for i in range(1, len(features) + 1)]
    if len(feature_ids) != len(features):
        raise ValueError('Length of feature ids does not match the length of features')
    _create_string_dataset(features_group, 'name', features)
    _create_string_dataset(features_group, 'id', feature_ids)
    _create_string_dataset(features_group, 'feature_type', ["Gene Expression"] * len(features))
    _create_string_dataset(features_group, "_all_tag_keys", "")
    features_group.close()

def _write_clusters(f: h5py.File, obs: pd.DataFrame) -> None:
    '''
    Writes the clusters to the h5 file
    '''
    cluster_group = f.create_group('clusters')
    for i in obs.columns:
        name = i
        cluster = obs.loc[:, i]
        group = cluster_group.create_group(name)
        _create_string_dataset(group, "name", [i])
        _create_string_dataset(group, "group_names", cluster.cat.categories.tolist())
        group.create_dataset("assignments", data=cluster.cat.codes, dtype=int)
        group.create_dataset("score", 0.0)
        _create_string_dataset(group, "clustering_type", ["unknown"])
        group.close()
    cluster_group.close()

def _write_projection(f: h5py.Group, dim: array, name: str) -> None:
    '''
    Writes the projections to the h5 file

    Args:
        f (h5py.Group): h5py group to write to
        dim (array): projection data
        name (str): name of the projection
    '''
    projection_group = f.create_group(name)
    _create_string_dataset(projection_group, "name", name)
    _create_string_dataset(projection_group, "method", name)
    projection_group.create_dataset("data", data=dim)
    projection_group.close()

def create_loupe(anndata: AnnData, output_file:str, layer: str | None = None, tmp_file: str="tmp.h5ad",
                 loupe_converter_path: str | None = None, dims: list[str] | None = None,
                 obs_keys: list[str] | None=None, 
                 feature_ids: list["str"]|pd.Series|None = None) -> None:
    ''''
    Creates a temp h5 file and calls the loupe converter executable for the conversion
    Args:
        anndata (AnnData): AnnData object to convert.
        output_file (str): Path to the output file.
        layer (str | None, optional): Layer to use. Defaults to None.
        tmp_file (str, optional): Path to the temp file. Defaults to "tmp.h5ad".
        loupe_converter_path (str | None, optional): Path to the loupe converter executable. Defaults to None.
        dims (list[str] | None, optional): Dimensions to use. Defaults to None.
        obs_keys (list[str] | None, optional): Keys of obs to use. Defaults to None.
        feature_ids (list["str"]|pd.Series|None, optional): Feature ids. Defaults to None.
    Raises:
        ValueError: If the output file does not exist.
        ValueError: If the layer is not valid.
        ValueError: If the obs keys are not valid.
        ValueError: If the feature ids are not valid.
    '''
    if not os.path.exists(os.path.basename(output_file)):
        raise ValueError('Output file does not exist')
    if loupe_converter_path is None:
        loupe_converter_path = _get_loupe_path()
    if not os.path.exists(loupe_converter_path):
        raise ValueError('Loupe converter path does not exist')
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
        features = anndata.var_names
        barcodes = anndata.obs_names
        if layer is None:
            _write_matrix(f, anndata.X.T, features, barcodes, feature_ids)
        else:
            _write_matrix(f, anndata.layers[layer].T, features, barcodes, feature_ids)
        _write_clusters(f, obs)
        projections = f.create_group('projections')
        if dims is None:
            dims = anndata.obsm.keys()
        for n in dims:
            if n not in anndata.obsm.keys():
                raise ValueError(f'{n} is not a valid projection')
            dim = anndata.obsm[n]
            _write_projection(projections, dim, n)
        projections.close()
        f.close()

        