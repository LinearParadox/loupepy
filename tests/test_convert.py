import pytest
import h5py  # type: ignore
import os
import numpy as np
import pandas as pd
import scanpy as sc  # type: ignore
from scipy.sparse import diags
from scipy.sparse import csc_matrix
from loupepy.convert import create_loupe_from_anndata, create_loupe  # type: ignore
from loupepy.utils import get_obs, get_obsm, get_count_matrix  # type: ignore

def reverse_engineer_counts(adata, n_counts_column="n_counts"):
    n_counts = adata.obs[n_counts_column].values
    counts_matrix = adata.X.expm1()
    counts = diags(n_counts) * counts_matrix
    adata.X = counts
    return adata

def cluster_dict(levels):
    """
    Create a dictionary mapping cluster levels to their assignments.
    """
    mapping_dict = {}
    for level, type in enumerate(levels):
        mapping_dict[level] = str(type.astype(str))
    return mapping_dict

def get_equivalent_matrix(file):
    """
    Extracts the equivalent matrix from the HDF5 file.
    """
    mat = csc_matrix((file['matrix/data'], file['matrix/indices'], file['matrix/indptr']),
                     shape=file['matrix/shape']).toarray()
    valid_barcodes = np.array(file['matrix/barcodes']).astype(str)
    valid_features=pd.Series(np.array(file['matrix/features/name']).astype(str))
    valid_features = np.array(valid_features.str.replace("_", "-")) # some genes have _ replaced with - in R.
    #we replace to ensure consistency
    #we sort the matrix to make sure the order of barcodes and features is consistent
    sort_col = np.argsort(valid_barcodes)
    sort_row = np.argsort(valid_features)
    mat=mat[:, sort_col]
    mat=mat[sort_row, :]
    return [mat, valid_barcodes[sort_col], valid_features[sort_row]]


@pytest.fixture(scope="module")
def adata_for_loupe():
    """Provides the AnnData object for Loupe conversion."""
    file = os.path.join(os.path.dirname(__file__), "data", "pbmc3k_processed.h5ad")
    adata_raw = sc.read_h5ad(file).raw.to_adata()
    return reverse_engineer_counts(adata_raw)

@pytest.fixture(scope="module")
def valid_h5():
    """Provides the path to the expected HDF5 file."""
    return os.path.join(os.path.dirname(__file__), "data", "loupeR_output.h5")

@pytest.fixture(scope="module")
def generate_h5_file(adata_for_loupe, tmp_path_factory):
    """
    Generates the Loupe HDF5 file and returns its path.
    Uses tmp_path_factory for a module-scoped temporary directory.
    """
    output_dir = tmp_path_factory.mktemp("generated_h5_data")
    generated_file_path = output_dir / "loupepy.h5"
    create_loupe_from_anndata(adata_for_loupe, tmp_file=generated_file_path, test_mode=True)
    return str(generated_file_path)

@pytest.fixture(scope="module")
def generate_subset(adata_for_loupe, tmp_path_factory):
    """
    Generates the Loupe HDF5 file and returns its path.
    Uses tmp_path_factory for a module-scoped temporary directory.
    """

    output_dir = tmp_path_factory.mktemp("generated_subset")
    generated_subset = output_dir / "loupepy.h5"
    adata_for_loupe.obs["some_cat"] = np.random.randint(0, 2, adata_for_loupe.n_obs)
    adata_for_loupe.obs["some_cat"] = adata_for_loupe.obs["some_cat"].astype("category")
    create_loupe_from_anndata(adata_for_loupe, tmp_file=generated_subset, test_mode=True,
                              dims=["X_umap"], obs_keys=["some_cat"])
    return str(generated_subset)

@pytest.fixture(scope="module")
def generate_manually(adata_for_loupe, tmp_path_factory):
    """
    Generates the Loupe HDF5 file and returns its path.
    Uses tmp_path_factory for a module-scoped temporary directory.
    Generates the file manually without using the create_loupe_from_anndata function.
    """
    output_dir = tmp_path_factory.mktemp("manually_generated")
    generated_subset = output_dir / "loupepy.h5"
    mat = get_count_matrix(adata_for_loupe)
    obs = get_obs(adata_for_loupe)
    obsm = get_obsm(adata_for_loupe)
    create_loupe(mat, obs, adata_for_loupe.var, obsm, tmp_file=generated_subset, test_mode=True)
    return str(generated_subset)

@pytest.fixture(scope="module")
def yield_tests():
    """
    Fixture to yield the test functions.
    This is useful if you want to run tests dynamically or in a specific order.
    """


@pytest.mark.parametrize("generate_h5_file", ["generate_h5_file", "generate_manually"], indirect=True)
def test_obs(generate_h5_file, valid_h5):
    """Test if the obs data is correctly written to the HDF5 file."""
    with h5py.File(valid_h5, "r") as valid:
        with h5py.File(generate_h5_file, "r") as generated:
            valid_levels=cluster_dict(valid["clusters/louvain/group_names"])
            generated_levels=cluster_dict(generated["clusters/louvain/group_names"])
            valid_types = pd.Series(valid["clusters/louvain/assignments"]).map(valid_levels).value_counts()
            generated_types = pd.Series(generated["clusters/louvain/assignments"]).map(generated_levels).value_counts()
            assert (valid_types == generated_types).all()

@pytest.mark.parametrize("generate_h5_file", ["generate_h5_file", "generate_manually"], indirect=True)
def test_projections(generate_h5_file, valid_h5):
    """Test if the projections data is correctly written to the HDF5 file."""
    with h5py.File(valid_h5, "r") as valid:
        with h5py.File(generate_h5_file, "r") as generated:
            generated_projections = pd.DataFrame(data=np.array(generated['projections/X_umap/data']).T,
                                             index=np.array(generated['matrix/barcodes']).astype(str))
            valid_projections = pd.DataFrame(data=np.array(valid['projections/Xumap_/data']).T,
                                                                index=np.array(valid['matrix/barcodes']).astype(str))
            generated_projections.sort_index(inplace=True)
            valid_projections.sort_index(inplace=True)
            assert valid_projections.equals(generated_projections)

@pytest.mark.parametrize("generate_h5_file", ["generate_h5_file", "generate_manually"], indirect=True)
def test_matrix(generate_h5_file, valid_h5):
    """Test if the matrix data is correctly written to the HDF5 file."""
    with h5py.File(valid_h5, "r") as valid:
        with h5py.File(generate_h5_file, "r") as generated:
            valid_matrix = get_equivalent_matrix(valid)
            generated_matrix = get_equivalent_matrix(generated)
            for x,y in zip(valid_matrix, generated_matrix):
                assert np.array_equal(x, y)

def test_subsetting(generate_subset):
    with h5py.File(generate_subset, "r") as f:
        assert set(f["projections"].keys()) == {"X_umap"}
        assert set(f['clusters'].keys()) == {"some_cat"}
