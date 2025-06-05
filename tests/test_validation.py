import pandas as pd
import pytest
import scanpy as sc  # type: ignore
from scipy.sparse import diags
import numpy as np
from loupepy.utils import (_validate_anndata, _validate_obs, _validate_counts, _validate_barcodes)  # type: ignore

def reverse_engineer_counts(adata, n_counts_column="n_counts"):
    """
    The generic scanpy pipeline converts the contents of `adata.X` into log-transformed
    normalized counts and stores it in a`data.raw.X`.

    This function effectively returns an AnnData object with the same `obs` and `var` dataframes
    that contain the original counts in `X`

    Copied from ttps://github.com/yelabucsf/scrna-parameter-estimation
    """
    n_counts = adata.obs[n_counts_column].values

    counts = diags(n_counts) * adata.raw.X.expm1()

    return sc.AnnData(X=counts, obs=adata.obs, var=adata.raw.var)


@pytest.fixture
def mock_data():
    return reverse_engineer_counts(sc.datasets.pbmc3k_processed())

@pytest.fixture
def generate_long_df():
    df = pd.DataFrame(np.zeros((40000, 3)))
    df['long'] = np.arange(40000)
    return df.astype('category')

def test_validate_counts(mock_data):
    assert _validate_counts(mock_data.X) is None

def test_validate_obs(mock_data):
    assert len(_validate_obs(mock_data.obs).columns) == 1

def test_validate_barcodes(mock_data):
    assert _validate_barcodes(mock_data.obs.index) is True

def test_long_category(generate_long_df):
    assert len(_validate_obs(generate_long_df).columns) == 3

def test_valid_anndata(mock_data):
    assert _validate_anndata(mock_data) is None

def test_invalid_index(mock_data):
    mock_data.obs = mock_data.obs.add_suffix("-abc", axis='index')
    assert _validate_barcodes(mock_data.obs.index) == False

