# Format Conversion for Single-cell RNA-seq Object

## From *MultimodalData* to *anndata*: `pegasus2anndata.py`

Load a `pegasus` *MultimodalData* object as input, and convert it to *anndata* in `h5ad` format.

Usage:
```
singularity exec -B ... robust-sc.sif python3 zarr2h5ad.py /path/to/input/file/name /path/to/output/file/name
```

The output `h5ad` file contains an *anndata* with one or two associated matrices, according to the *current matrix* in pegasus *MultimodalData* object:
- the current matrix (e.g., `counts.log_norm`) is bound to `.X`;
- the precedent matrix (e.g., `counts`) is bound to `.raw.X`;
- if the precedent matrix does not exist, put it as `None`.

The precedent matrix is determined by the delimiter `.`, see [PegasusIO GitHub Page](https://github.com/lilab-bcb/pegasusio/blob/master/pegasusio/unimodal_data.py) lines 597-637.
