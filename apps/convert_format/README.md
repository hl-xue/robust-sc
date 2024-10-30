# Format Conversion for Single-cell RNA-seq Object

Load a `pegasus` *MultimodalData* object as input, and convert it to:
- sparse count matrix in the given folder,
- `.h5` file, 
- *AnnData* object stored in `.h5ad` format.

Usage:
```
singularity exec -B ... robust-sc.sif python3 convert_format.py /path/to/input/file/name /path/to/output/
```

If the output is a folder (for sparse count matrix) or `.h5` file, it contains:
- the `raw.X` matrix if it exists;
- the `counts` matrix if it exists and the `raw.X` matrix does not exist;
- the current matrix if neither `raw.X` nor `counts` exists, with a warning.

If the output is a `.h5ad` file, it contains an *AnnData* with one or two associated matrices:
- the `X` matrix or the current matrix (e.g., `counts.log_norm`) is bound to `.X`;
- if a `raw.X` matrix or `counts` matrix exists, it is bound to `raw.X` and used to construct `adata.raw`;
- if `raw.X` matrix or `counts` matrix does not exist, it is assigned as `None`.

The part of script for convertion to *AnnData* was adapted from lines 597-637 in [PegasusIO GitHub Page](https://github.com/lilab-bcb/pegasusio/blob/master/pegasusio/unimodal_data.py).
