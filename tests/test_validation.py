import pytest
import scanpy as sc
from scipy.sparse import diags

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
    return reverse_engineer_counts(sc.datasets.pbmc3k_processed)



def test_validate_counts():


