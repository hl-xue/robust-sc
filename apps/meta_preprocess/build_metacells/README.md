# Metacell Construction

The workflow builds metacells with the `SuperCell` package. It takes an `h5` file as input (with unprocessed counts), and outputs a `csv` table associating each single cells to its metacell.

The publication of the `SuperCell` package can be found [here](https://doi.org/10.1186/s12859-022-04861-1).
Please cite their article if using this script.

A nice review article on general *metacell analysis* can be found [here](https://doi.org/10.1038/s44320-024-00045-6).

To run the analysis, please do:
```bash
SuperCell.R [-h] [-f FEATURE.NUM] [-c PC.NUM] [-k NEIGHBOUR.NUM] [-g GAMMA] input out.prefix
```
```
positional arguments:
  input                 input file in .h5 format
  output.prefix         output prefix for .csv and .RData files

optional arguments:
  -h, --help                            show this help message and exit
  -f, --feature.num FEATURE.NUM         number of variable features to be used [1000]
  -c, --pc.num PC.NUM                   number of principal components to be used [50]
  -k, --neighbour.num NEIGHBOUR.NUM     number of neighbours for k-NN [5]
  -g, --gamma GAMMA                     graining level [15]
```

A `csv` table with two columns will be output.
The first column contains cell barcodes, and the second one indicates associated metacells
