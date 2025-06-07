import pandas as pd
import pytest
import scanpy as sc  # type: ignore
from scipy.sparse import diags
import numpy as np
from loupepy.utils import (_validate_anndata, _validate_obs, _validate_counts, _validate_barcodes)  # type: ignore
import os

def reverse_engineer_counts(adata, n_counts_column="n_counts"):
    """
    The generic scanpy pipeline converts the contents of `adata.X` into log-transformed
    normalized counts and stores it in a`data.raw.X`.

    This function effectively returns an AnnData object with the same `obs` and `var` dataframes
    that contain the original counts in `X`

    Copied from ttps://github.com/yelabucsf/scrna-parameter-estimation
    """
    n_counts = adata.obs[n_counts_column].values
    counts_matrix = adata.X.expm1()
    counts = diags(n_counts) * counts_matrix
    adata.X = counts
    return adata


@pytest.fixture
def mock_data():
    file = os.path.join(os.path.dirname(__file__), "data", "pbmc3k_processed.h5ad")
    adata_raw = sc.read_h5ad(file).raw.to_adata()
    return reverse_engineer_counts(adata_raw)

@pytest.fixture
def generate_long_df():
    df = pd.DataFrame(np.zeros((40000, 3)))
    df['long'] = np.arange(40000)
    return df.astype('category')

def test_validate_counts(mock_data):
    assert _validate_counts(mock_data.X) is None

def test_validate_obs(mock_data):
    _validate_obs(mock_data.obs)
    assert len(mock_data.obs.columns) == 1

def test_validate_barcodes(mock_data):
    assert _validate_barcodes(mock_data.obs.index) is True

def test_long_category(generate_long_df):
    _validate_obs(generate_long_df)
    assert len(generate_long_df.columns) == 3

def test_valid_anndata(mock_data):
    assert _validate_anndata(mock_data) is None

def test_invalid_index(mock_data):
    mock_data.obs = mock_data.obs.add_suffix("-abc", axis='index')
    assert not _validate_barcodes(mock_data.obs.index)

