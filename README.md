A utility to convert anndata files into loupe files visualizable in 10xs loupe browser.

This package widely mirrors the behavior of the R package, with some minor differences:  

- Cell based annotations are automatically filtered out if they contain more than the max number of categories.
  (This was an issue in R when converting from anndata).
- Easier to pass specific embeddings to the function, which is especially useful when using scVI.
  
Note this package requires an agreement to the 10x loupe converter license, which can be agreed to and downloaded by  
running the setup function.
** I am not an employee of 10x, nor associated with them **. This package simply writes a file which the converter can read  
and calls it.

# This is still not done and a work in progress #
