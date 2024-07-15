"""
Load a Pegasus MultimodalData file as input, and convert it to anndata in h5ad format.
Usage: python3 pegasus2anndata.py /path/to/input/file/name /path/to/output/file/name
Author: Haoliang Xue
"""

import sys
import os
import gc
import anndata
import pegasus as pg


def Convert(data_src):
    # Rewritten from the to_anndata() method on https://github.com/lilab-bcb/pegasusio/blob/master/pegasusio/unimodal_data.py
    matrix_keys = data_src.list_keys("matrix")
    if "X" in matrix_keys:
        X_key = "X"
        if "raw.X" in matrix_keys:
            raw_key = "raw.X"
    else:
        X_key = data_src.current_matrix()
        components = X_key.split(".")
        if len(components) > 1:
            raw_key = '.'.join(components[:-1])
            if raw_key not in matrix_keys:
                raw_key = None
    layers = {}
    for matkey in matrix_keys:
        layers[matkey] = data_src.get_matrix(matkey)
    if raw_key != None:
        raw = anndata.AnnData(X = data_src.get_matrix(raw_key), dtype = data_src.get_matrix(raw_key).dtype, var = data_src.var[["featureid"]])
    else:
        raw = None
    return anndata.AnnData(X = data_src.get_matrix(X_key), dtype = data_src.get_matrix(X_key).dtype,
                           obs = data_src.obs, var = data_src.var, uns = data_src.uns,
                           obsm = data_src.obsm, varm = data_src.varm, obsp = data_src.obsp, varp = data_src.varp,
                           layers = layers, raw = raw)


if __name__ == "__main__":
    # Parse command-line arguments
    assert len(sys.argv) == 3
    infile, outfile = sys.argv[1:]
    print(f"* Input file: {infile}.", sep="", file=sys.stderr)
    print(f"* Output file: {outfile}.\n", sep="", file=sys.stderr)
    if not os.path.isfile(infile):
        print(f"\033[91m* ERROR: Input file {infile} does not exist.\033[0m", file=sys.stderr)
        assert False
    if not os.path.isdir(os.path.dirname(outfile)):
        print(f"\033[91m* ERROR: Output folder {os.path.dirname(outfile)} does not exist.\033[0m", file=sys.stderr)
        assert False
    if os.path.exists(outfile):
        print(f"\033[91m* ERROR: Output file already exists. Please manually remove it to avoid accidental overwriting.\033[0m", file=sys.stderr)
        assert False

    # Load input file
    data = pg.read_input(infile)
    print("* Loaded object:\n", data, sep="", file=sys.stderr)
    print("* Example counts in current matrix:", data.X[data.X > 10], file=sys.stderr)
    gc.collect()

    # Write output file
    Convert(data).write_h5ad(outfile)
    print("\033[92m* Conversion complete!", file=sys.stderr)
