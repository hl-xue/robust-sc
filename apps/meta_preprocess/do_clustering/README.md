# Preprocessing Workflow till Clustering Step

The workflow is built with the `Pegasus` package which is reported to outperform `Scanpy` on omputational effeciency, and provides a submodule to optimise the selection of top components after PCA.

The publication of the `Pegasus` package can be found [here](https://doi.org/10.1038/s41592-020-0905-x), please also cite their work if using this module.

To run the analysis, please do:
```bash
python3 run_till_clustering.py [-h] [-p INITIAL_PC_NUM] [-r LEIDEN_RESOLUTION] [-c CLUSTERING_COLNAME] [-s OUT_SAMPLE_LABEL] [-g] input output
```
```
positional arguments:
  input                 input file (.h5, .zarr, .zarr.zip, .h5ad, .h5ad.zip, matrix folder)
  output                output folder

options:
  -h, --help                                                        show this help message and exit.
  -p INITIAL_PC_NUM, --initial_pc_num INITIAL_PC_NUM                initial PC number to reduce [INT].
  -r LEIDEN_RESOLUTION, --leiden_resolution LEIDEN_RESOLUTION       resolution for leiden clustering [FLOAT].
  -c CLUSTERING_COLNAME, --clustering_colname CLUSTERING_COLNAME    column name of clustering results [STR].
  -s OUT_SAMPLE_LABEL, --out_sample_label OUT_SAMPLE_LABEL          sample name as output file name prefix [STR].
  -g, --regress_cellcycle                                           whether to regress cell cycle.
```

Three output files will be generated in the indicated output folder:
- `{SAMPLE_NAME}_HVFplot.pdf`: a variance-mean scatter plot of features with the highly variable ones being marked.
- `{SAMPLE_NAME}_elbowPlot.pdf`: elbowplot with suggested number of top PCs based on random matrix theory.
- `{SAMPLE_NAME}_tillClustering.zarr.zip`: processed object till the clustering step.
