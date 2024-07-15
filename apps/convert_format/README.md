# Format Conversion for ScRNA-seq Object

## From `zarr` to `h5ad`: `zarr2h5ad.py`

Load a zarr/zarr.zip file (e.g., `pegasus` `MultimodalData` object) as input, and convert it to `anndata` in h5ad format (e.g., `scanpy` object).

Usage:
```
singularity exec -B ... robust-sc.sif python3 zarr2h5ad.py /path/to/input/file/name /path/to/output/file/name
```

The output `h5ad` file contains an anndata with one or two associated matrices, according to the current matrix in `pegasus MultimodalData` object:
- the current matrix (e.g., `counts.log_norm`) is binded to `.X`;
- the precedent matrix (e.g., `counts`) is binded to `.raw.X`.

The precedent matrix is determined by the delimiter `.`, see [PegasusIO GitHub Page](https://github.com/lilab-bcb/pegasusio/blob/master/pegasusio/unimodal_data.py) line 597-637.
