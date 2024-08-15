# Adaptive QC with `ddqc`

`ddqc` is a package that apply quality control for single cell/nucleus RNA-seq data.
Its publication can be found [here](https://doi.org/10.1186/s13059-022-02820-w).

In addition, users can apply a manual filtration with hard threshold setting of:
- minimum detected gene number per cell;
- minimum UMI count per cell;
- maximum mitochondrial gene UMI count percentage per cell.

To run the analysis, please do:
```bash
run_ddqc.py [-h] [-u MIN_UMI_PER_CELL] [-g MIN_GENE_PER_CELL] [-M MAX_MITO_PER_CELL] [-s SAMPLE_NAME] input output
```
```
positional arguments:
  input                 input file (.h5, .zarr, .zarr.zip, .h5ad, .h5ad.zip, matrix folder).
  output                output folder.

options:
  -h, --help                                                        show this help message and exit
  -u MIN_UMI_PER_CELL, --min_umi_per_cell MIN_UMI_PER_CELL          minimum UMI count per cell [INT].
  -g MIN_GENE_PER_CELL, --min_gene_per_cell MIN_GENE_PER_CELL       minimum expressed gene number per cell [INT].
  -M MAX_MITO_PER_CELL, --max_mito_per_cell MAX_MITO_PER_CELL       maximum mitochondrial gene UMI percentage per cell [INT].
  -s SAMPLE_NAME,       --sample_name SAMPLE_NAME                   sample name to write in output folder [STR].
```

Three output files will be generated:
- `{SAMPLE_NAME}_dfqc_table.csv`: the `ddqc` filtering results in table format.
- `{SAMPLE_NAME}_full.zarr.zip`: the complete object with removed barcodes being marked.
- `{SAMPLE_NAME}_filtered.zarr.zip`: the filtered object, only with barcodes passing `ddqc` and manual removal steps.