# Metacell Construction

The workflow builds metacells with the `SuperCell` package. It takes an `h5` file as input (with unprocessed counts), and outputs a `csv` table associating each single cells to its metacell.

The publication of the `SuperCell` package can be found [here](https://doi.org/10.1186/s12859-022-04861-1).

A nice review article on *metacell analysis* can be found [here](https://doi.org/10.1038/s44320-024-00045-6).

To run the analysis, please do:
```bash
Rscript SuperCell.R [-h] [-n NUM.HVG] [-p NUM.PC] [-g GAMMA] input output
```
```
positional arguments:
  input                 input file in .h5 format
  output                output table in .csv format

optional arguments:
  -h, --help                        show this help message and exit
  -n NUM.HVG, --num.hvg NUM.HVG     number of HVGs to be used [1000]
  -p NUM.PC, --num.pc NUM.PC        number of PCs to be used [50]
  -g GAMMA, --gamma GAMMA           graining level [15]
```

A `csv` table with two columns will be output.
The first column contains cell barcodes, and the second one indicates associated metacells
