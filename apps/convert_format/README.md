# Format Conversion for Single-cell RNA-seq Object

Load a `pegasus` *MultimodalData* object as input, and convert it to:
- matrix count tables in the given folder,
- `.h5` file, 
- *AnnData* object stored in `.h5ad` format.

Usage:
```
singularity exec -B ... robust-sc.sif python3 convert_format.py /path/to/input/file/name /path/to/output/
```

The input file name only allows `.zarr` or `.zarr.zip`.

As for the output:
- If the output is a folder or `.h5` file, it contains the raw count matrix.
- If the output is a `.h5ad` file, it contains an *AnnData* with one or two associated matrices, according to the *current matrix* in pegasus *MultimodalData* object:
    - the `X` matrix or the current matrix (e.g., `counts.log_norm`) is bound to `.X`;
    - the `raw.X` matrix or the `counts` matrix is bound to `.raw.X`;
    - if the `raw.X` matrix or precedent matrix does not exist, put it as `None`.

The part of script for convertion to *AnnData* was adapted from lines 597-637 in [PegasusIO GitHub Page](https://github.com/lilab-bcb/pegasusio/blob/master/pegasusio/unimodal_data.py).
