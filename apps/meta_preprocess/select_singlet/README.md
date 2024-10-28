# Singlet Selection

This App practises a singlet selection task with two versions of Scrublet implementation: incorporated in `Pegasus` and `Scanpy` packages. The App also reports a simple comparison of identified non-singlets between the two packages with a contingency table for reference.

Four modes of filtration are supported:
- pegasus: remove non-singlets predicted by `Pegasus` implementation;
- scanpy: remove non-singlets predicted by `Scanpy` implementation;
- all: only remove non-singlets predicted by *both* implementations, as a *permissive* filtration;
- any: remove non-singlets predicted by *any* implementations, as a *strict* filtration.

Please cite the relevant works of [Scrublet](https://doi.org/10.1016/j.cels.2018.11.005), [Pegasus](https://doi.org/10.1038/s41592-020-0905-x) and [Scanpy](https://doi.org/10.1186/s13059-017-1382-0) accordingly to the usage.

To run the analysis, please do:
```bash
python3 select_singlet.py [-h] [-m METHOD] [-s OUT_SAMPLE_LABEL] input output
```
```
positional arguments:
  input                 input file (.h5, .zarr, .zarr.zip, .h5ad, .h5ad.zip, matrix folder)
  output                output folder

options:
  -h, --help                            show this help message and exit
  -m STR, --method STR                  singlet selection method [STR]
  -s STR, --out_sample_label STR        sample name as output file name prefix [STR]
```

Four output files will be generated in the indicated output folder:
- `{SAMPLE_NAME}_doublet_score_density.pdf`: plots summarising `Pegasus` prediction of doublet scores.
- `{SAMPLE_NAME}_prediction_contingency.csv`: contingency table of predicted doublets between `Pegasus` and `Scanpy`.
- `{SAMPLE_NAME}_full.zarr.zip`: the complete object with non-singlet barcodes being marked.
- `{SAMPLE_NAME}_singlets.zarr.zip`: the filtered object, only with singlet barcodes.
