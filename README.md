While this tool is independently maintained and not officially supported by 10x Genomics, it relies on code from [LoupeR](https://github.com/10xGenomics/loupeR), an R package developed by 10x Genomics. LoupeR is licensed for use only in connection with data generated from 10x Genomics products. For more information, please refer to the [10x End User Software License](https://www.10xgenomics.com/legal/end-user-software-license-agreement). By using the setup function of this tool, you will need to agree to this license with an interactive prompt. 

## Installation

This can be simply installed from pip, using:  
`pip install loupepy`

You must run the setup function and agree to the eula. This only needs to be done once.

```{py}
from loupepy import setup

setup()

> This tool is independently maintained, but utilizes 10x genomics tools to perform conversions
> By using this tool, you agree to the 10X Genomics End User License Agreement (https://www.10xgenomics.com/legal/end-user-software-license-agreement)
> Do you agree to the terms of the EULA? (yes/y or no/n)

y
```
Will install the tool to the default program install directory in your OS. A custom path for the install can also be provided, however you will need to include it in subsequent commands using the:
`loupe_converter_path` argument.

## Main Functions

### Notes

There are some differences between this one and the R version. Mainly, the R version will throw an error if a categorical has too many categories, or a dimensionality reduction isn't 2 x n_cells.  Here, these are dropped with a warning. Enabling the strict parameter will mimic R's behavior.

### `create_loupe_from_anndata(anndata, output_cloupe="cloupe.cloupe", layer=None, ...)`

Creates a Loupe Browser compatible `.cloupe` file from an `AnnData` object.

**Parameters:**
- `anndata`: AnnData object containing single-cell data
- `output_cloupe`: Path to save the output file (default: "cloupe.cloupe")
- `layer`: Layer in AnnData to use for counts (default: None, uses `.X`)
- `dims`: List of dimension reduction keys from `adata.obsm` (default: None, uses all valid projections)
- `obs_keys`: List of categorical annotations from `adata.obs` to include (default: None, uses all valid categories)
- `strict_checking`: Whether to strictly validate inputs (default: False)
- `tmp_file`: Temporary file path for intermediate data (default: 'tmp.h5')
- `dims`: List of dimension reduction keys from `adata.obsm` to include (default: None, uses all valid projections)
- `obs_keys`: List of categorical annotations from `adata.obs` to include (default: None, uses all valid categories)
- `loupe_converter_path`: Path to the Loupe converter binary (default: None, uses default installation path)
- `clean_tmp_file`: Whether to delete the temporary file after creation (default: True)
- `force`: Whether to overwrite existing cloupe file (default: False)
- `test_mode`: Whether to run in test mode. Will not create a cloupe file.

**Example:**
```python
import scanpy as sc
import loupepy
from scipy.sparse import diags

# Load data
adata = sc.datasets.pbmc3k_processed().raw.to_adata()
# Revert to raw counts, as loupee converter only takes raw counts
n_counts = adata.obs["n_counts"].values
counts_matrix = adata.X.expm1()
counts = diags(n_counts) * counts_matrix
adata.X = counts

# Create Loupe file with default settings
loupepy.create_loupe_from_anndata(adata, "my_results.cloupe")

# Create with specific dimensions and annotations
loupepy.create_loupe_from_anndata(
    adata,
    "my_custom_results.cloupe",
    dims=["X_umap"],
    obs_keys=["leiden"]
)
```
## Utility Functions

### `get_obs(anndata, obs_keys=None, strict=False)`

Extracts categorical observation data from an AnnData object.

**Parameters:**
- `anndata`: AnnData object containing the data
- `obs_keys`: List of keys to extract from `adata.obs` (default: None, uses all valid categories)
- `strict`: Whether to strictly validate the keys. If true will throw errors on invalid keys instead of dropping
(default: False)

### `get_obsm(anndata, obsm_keys=None, strict=False)`

Extracts dimensional reduction data from an AnnData object.
**Parameters:**
- `anndata`: AnnData object containing the data
- `obsm_keys`: List of keys to extract from `adata.obsm` (default: None, uses all valid dimension reductions)
- `strict`: Whether to strictly validate the keys. If true will throw errors on invalid keys instead of dropping
(default: False)
### `get_count_matrix(anndata, layer=None)`

Gets the count matrix from an AnnData object in the format required for Loupe.

**Parameters:**
- `anndata`: AnnData object containing the data
- `layer`: Layer in AnnData to use for counts (default: None, uses `.X`)

---

### `create_loupe(mat, obs, var, obsm, tmp_file, ...)`

Lower-level function to create a Loupe file from raw data components.

**Parameters:**
- `mat`: CSC matrix of counts (shape: n_features × n_cells)
- `obs`: DataFrame of categorical cell annotations
- `var`: DataFrame of gene/feature information
- `obsm`: Dictionary of dimension reduction embeddings
- `tmp_file`: Path for temporary file
- `output_path`: Path to save the output file
- `strict_checking`: Whether to strictly validate inputs (default: False)
- `loupe_converter_path`: Path to the Loupe converter binary (default: None, uses default installation path)'
- `clean_tmp_file`: Whether to delete the temporary file after creation (default: True)
- `force`: Whether to overwrite existing cloupe file (default: False)
- `test_mode`: Whether to run in test mode. Will not create a cloupe file.

**Example:**
```python
import scanpy as sc
import loupepy
from scipy.sparse import diags

# Load data
adata = sc.datasets.pbmc3k_processed().raw.to_adata()
# Revert to raw counts, as loupee converter only takes raw counts
n_counts = adata.obs["n_counts"].values
counts_matrix = adata.X.expm1()
counts = diags(n_counts) * counts_matrix
adata.X = counts

mat = loupepy.get_count_matrix(adata)
obs = loupepy.get_obs(adata, obs_keys=["leiden"])
var = adata.var
obsm = loupepy.get_obsm(adata, obsm_keys=["X_umap"])
loupepy.create_loupe(mat,obs,var,obsm,"my_results.cloupe")
```

## Setup Functions

### `setup(path=None)`

Downloads and installs the Loupe converter binary.

**Parameters:**
- `path`: Installation directory path (default: None, uses OS-specific default location)

**Example:**
```python
import loupepy

# Install the Loupe converter to the default location
loupepy.setup()
```

### `eula(path=None)`

**Parameters:**
- `path`: Installation path (default: None, uses OS-specific default install location. Usually same path as setup)

### `eula_reset(path=None)`

Resets the EULA agreement and removes the installed Loupe converter binary.

## Acknowledgments

A special thank you to the Lowy Lab for their support.  
Thank you to 10x Genomics for supporting this and allowing the use of their tool.

