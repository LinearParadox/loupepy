import os
from os import PathLike

from anndata import AnnData # type: ignore
from anndata import __version__ as ad_version
import pandas as pd
from scipy.sparse import csc_matrix
import h5py # type: ignore
from typing import List
import logging
from numpy import ndarray
import platform
from .utils import _validate_anndata, _get_loupe_path, get_obs, get_obsm, get_count_matrix
from . import __version__



def _create_string_dataset(obj: h5py.Group, key: str, strings: List[str]|pd.Series|str|pd.Index) -> None:
    '''
    Creates a dataset in the h5 file with the given strings
    '''
    if len(strings) == 0:
        max_len = 1
    elif isinstance(strings, str):
        max_len = len(strings)
        strings = [strings]
    elif isinstance(strings, pd.Series) or isinstance(strings, pd.Index):
        strings = strings.astype(str)
        strings = strings.to_list()
        max_len = max(len(s) for s in strings)
    else:
        max_len = max(len(str(s)) for s in strings)
    dtype = h5py.string_dtype(encoding='ascii', length=max_len)
    if strings == "":
        #for the special case of an empty string
        #matches r behavior
        obj.create_dataset(name=key, dtype=dtype, shape=(0,))
    else:
        obj.create_dataset(name=key, data=strings, dtype=dtype)

def _write_metadata(f: h5py.File) -> None:
    '''
    Writes the metadata to the h5 file
    '''
    metadata=f.create_group('metadata')
    meta = {}
    h5_version = h5py.h5.HDF5_VERSION_COMPILED_AGAINST
    h5_version = f"{h5_version[0]}.{h5_version[1]}.{h5_version[2]}"
    meta["tool"] = "LoupePy"
    meta["tool_version"] = __version__
    meta["os"] = platform.system()
    meta["platform"] = platform.platform()
    meta["language"] = f"Python -- {platform.python_version()}"
    meta["h5py_version"] = h5py.__version__
    meta["anndata_version"] = ad_version
    meta["hdf5_version"] = h5_version
    for k,v in meta.items():
        v = str(v)
        _create_string_dataset(metadata, k, v)

def _write_matrix(f: h5py.File, matrix: csc_matrix,
                  features: pd.Series|pd.Index, barcodes: pd.Index|pd.Series,
                  feature_ids: list[str]|pd.Series|None = None) -> None:
    '''
    Writes the matrix to the h5 file
    '''
    if not isinstance(matrix, csc_matrix):
        f.close()
        raise ValueError('Matrix must be a csc_matrix')
    matrix_group = f.create_group('matrix')
    features_group = matrix_group.create_group('features')
    _create_string_dataset(matrix_group, 'barcodes', barcodes)
    matrix_group.create_dataset('data', data=matrix.data, dtype='i4') # loupe expects int32
    matrix_group.create_dataset('indices', data=matrix.indices, dtype=int)  # type: ignore
    matrix_group.create_dataset('indptr', data=matrix.indptr, dtype=int)  # type: ignore
    matrix_group.create_dataset('shape', data=matrix.shape, dtype=int)
    if feature_ids is None:
        feature_ids = [f"feature_{i}" for i in range(1, len(features) + 1)]
    if len(feature_ids) != len(features):
        raise ValueError('Length of feature ids does not match the length of features')
    _create_string_dataset(features_group, 'name', features)
    _create_string_dataset(features_group, 'id', feature_ids)
    _create_string_dataset(features_group, 'feature_type', ["Gene Expression"] * len(features))
    _create_string_dataset(features_group, "_all_tag_keys", "")

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
        group.create_dataset(name="score", shape=(1,), data=[0.0])
        _create_string_dataset(group, "clustering_type", "unknown")

def _write_projection(f: h5py.Group, dim: ndarray, name: str) -> None:
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
    projection_group.create_dataset("data", data=dim.T)

def create_loupe_from_anndata(anndata: AnnData, output_cloupe: str | PathLike = "cloupe.cloupe",
                              layer: str | None = None, tmp_file: str|PathLike ="tmp.h5",
                              loupe_converter_path: str | None | PathLike = None, dims: list[str] | None = None,
                              obs_keys: list[str] | None=None, feature_ids: list["str"]|pd.Series|None = None,
                              strict_checking: bool = False, clean_tmp_file: bool=True, force: bool = False,
                              test_mode=False) -> None:
    ''''
    Creates a temp h5 file and calls the loupe converter executable for the conversion
    Args:

        anndata (AnnData): AnnData object to convert.
        output_cloupe (str): Path to the output file.
        layer (str | None, optional): Layer to use. Defaults to None.
        tmp_file (str, optional): Path to the temp file. Defaults to "tmp.h5ad".
        loupe_converter_path (str | None, optional): Path to the loupe converter executable. Defaults to None.
        dims (list[str] | None, optional): Dimensions to use. Defaults to None.
        obs_keys (list[str] | None, optional): Keys of obs to subset to. Defaults to None.
        feature_ids (list["str"]|pd.Series|None, optional): Feature ids. Defaults to None.
        strict_checking (bool, optional): Whether to perform strict checking of the AnnData object. Defaults to False.
        If strict_checking is True, the function will raise an error if projections are not of the (ncells, 2) or if
        categoricals have more than 32768 categories. Default behavior is to drop these without an error.
        clean_tmp_file: whether to delete the temporary file after conversion. Defaults to True.
        force: whether to overwrite the cloupe file if it already exists. Defaults to False.
        test_mode: If True, will not run the loupe converter and will only write the h5 file.
    Raises:
        ValueError: If the output file does not exist.
        ValueError: If the layer is not valid.
        ValueError: If the obs keys are not valid.
        ValueError: If the feature ids are not valid.
    '''
    if not test_mode:
        if loupe_converter_path is None:
            loupe_converter_path = _get_loupe_path()
        if not os.path.exists(loupe_converter_path) and not test_mode:
            if loupe_converter_path is None:
                raise ValueError('Loupe converter Not found at default install location.'
                                 'Please run loupepy.setup() to install it or provide a path to the executable.')
            else:
                raise ValueError(f'Loupe converter not found at {loupe_converter_path}. Please provide a valid path.')
    if test_mode:
        logging.warning("Test mode is enabled. Loupe file will not be created.")
        clean_tmp_file = False
    _validate_anndata(anndata, layer)
    obs = get_obs(anndata, obs_keys=obs_keys, strict=strict_checking)
    mat = get_count_matrix(anndata, layer=layer)
    projections = get_obsm(anndata, obsm_keys=dims, strict=strict_checking)
    if len(projections) == 0:
        raise ValueError("No valid projections!")
    create_loupe(mat, obs, anndata.var, projections, tmp_file,
                 loupe_converter_path, output_path=output_cloupe, clean_tmp_file=clean_tmp_file,
                 feature_ids=feature_ids, force=force, test_mode=test_mode)


def create_loupe(mat: csc_matrix,
                 obs: pd.DataFrame,
                 var: pd.DataFrame,
                 obsm: dict[str, ndarray],
                 tmp_file: PathLike | str,
                 loupe_converter_path: str|PathLike|None = None,
                 output_path: PathLike|str = "cloupe.cloupe",
                 clean_tmp_file: bool = True,
                 force: bool = False,
                 feature_ids: list[str]|pd.Series|None = None,
                 test_mode=False) -> None:
    '''
    Creates a loupe file from a matrix, obs, var and obsm.
    Args:
        mat: csc matrix of shape (n_features, n_cells)
        obs: obs dataframe of shape (n_cells, n_obs)
        var: var (genes) dataframe of shape (n_features, n_vars)
        obsm: dimension reductions or other projections
        tmp_file: where to write the temporary h5 file
        loupe_converter_path: path to the loupe converter executable. If None, will use the default path.
        output_path: path to the output loupe file.
        clean_tmp_file: whether to delete the temporary file after conversion
        force: whether to overwrite the cloupe file if it already exists
        feature_ids: Feature ids to use. If None, will use the default feature ids.
        Seems to be not used in loupe browser
        test_mode: will not write a cloupe file. mainly used for testing
    Returns:
        None
    '''
    _write_hdf5(mat, obs, var, obsm, tmp_file, feature_ids=feature_ids)
    if test_mode:
        return
    if not test_mode:
        if loupe_converter_path is None:
            loupe_converter_path = _get_loupe_path()
        if not os.path.exists(loupe_converter_path):
            raise ValueError('Loupe converter path does not exist')
    cmd = f"{loupe_converter_path} create --input={tmp_file} --output={output_path}"
    if force:
        cmd += " --force"
    os.system(cmd)
    if clean_tmp_file:
        os.remove(tmp_file)



def _write_hdf5(mat: csc_matrix,
                obs: pd.DataFrame,
                var: pd.DataFrame,
                obsm: dict[str, ndarray],
                file_path: PathLike|str,
                feature_ids: list[str]|pd.Series|None = None,) -> None:
    try:
        with h5py.File(file_path, 'w') as f:
            features = var.index
            barcodes = obs.index
            _write_matrix(f, mat, features, barcodes, feature_ids)
            _write_clusters(f, obs)
            projections = f.create_group('projections')
            for n in obsm.keys():
                _write_projection(projections, obsm[n], n)
            _write_metadata(f)
            f.close()
    except ValueError:
        logging.error("Something went wrong while writing the h5 file. Please check the input data.")
        os.remove(file_path)


