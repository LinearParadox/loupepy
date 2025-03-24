import os
from anndata import AnnData # type: ignore
import pandas as pd
from scipy.sparse import csr_matrix
import h5py # type: ignore
from loupepy.utils import _validate_anndata, _validate_obs


def _create_string_dataset(obj: h5py.Group, key: str, strings: list[str] | pd.Series) -> None:
    '''
    Creates a dataset in the h5 file with the given strings
    '''
    if len(strings) == 0:
        max_len = 1
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
    _create_string_dataset(features_group, "_all_tag_keys", [""])
    features_group.close()

def _write_clusters(f: h5py.File, obs: pd.DataFrame) -> None:
    '''
    Writes the clusters to the h5 file
    '''
    cluster_group = f.create_group('clusters')
    for i in obs.columns:
        name = i
        cluster = df.loc[:, i]
        group=cluster_group.create_group(name)
        _create_string_dataset(group, "name", name)
        _create_string_dataset(group, "group_names", cluster.cat.categories)
        group.create_dataset("assignments", data=cluster.cat.codes, dtype=int)
        group.create_dataset("score", 0.0)
        _create_string_dataset(group, "clustering_type", "unknown")
        group.close()
    cluster_group.close()


        

def create_loupe(anndata: AnnData, output_file:str, layer: str | None = None, tmp_file: str="tmp.h5ad",
                 loupe_converter_path: str | None = None, dims: list[str] | None = None,
                 obs_keys: list[str] | None=None, feature_ids: list["str"]|pd.Series|None = None) -> None:
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
        features = anndata.var_names
        barcodes = anndata.obs_names
        _write_matrix(f, anndata.X.T, features, barcodes, feature_ids)
        _write_clusters(f, obs)