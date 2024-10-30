import sys
import os
import argparse


def Convert2Anndata(data_src):
    # Rewritten from the to_anndata() method on https://github.com/lilab-bcb/pegasusio/blob/master/pegasusio/unimodal_data.py
    matrix_keys = data_src.list_keys("matrix")
    if "X" in matrix_keys:
        X_key = "X"
        if "raw.X" in matrix_keys:
            raw_key = "raw.X"
    else: # if `X` is not a matrix key
        X_key = data_src.current_matrix() # set current matrix as .X
        if (X_key == "counts") or ("counts" not in matrix_keys): # if current matrix is "counts" or if "counts" does not exist
            raw_key = None # set .raw.X as None
        else:
            raw_key = "counts" # set .raw.X as "counts"
    if raw_key != None:
        raw = anndata.AnnData(X = data_src.get_matrix(raw_key),
                              dtype = data_src.get_matrix(raw_key).dtype,
                              var = data_src.var[["featureid"]])
    else:
        raw = None
    layers = {matkey: data_src.get_matrix(matkey) for matkey in matrix_keys}
    return anndata.AnnData(X = layers[X_key], dtype = layers[X_key].dtype,
                           obs = data_src.obs, var = data_src.var, uns = data_src.uns,
                           obsm = data_src.obsm, varm = data_src.varm,
                           obsp = data_src.obsp, varp = data_src.varp,
                           layers = layers, raw = raw)


def Test(converted_data, src_data):
    ref_data = src_data.to_anndata()
    assert (converted_data.X != ref_data.X).sum() == 0
    assert (converted_data.raw.X != ref_data.raw.X).sum() == 0
    assert (converted_data.obs != ref_data.obs).sum().sum() == 0
    assert (converted_data.var != ref_data.var).sum().sum() == 0


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Convert MultimodalData data structure (used in Pegasus) to: " +
                                                 "matrix folder, general h5 file or anndata in h5ad.")
    parser.add_argument("input", help="input file")
    parser.add_argument("output", help="output file (directory, .h5, .h5ad)")
    args = vars(parser.parse_args())
    _, format = os.path.splitext(args["output"])
    if format == "":
        format = "matrix folder"
    for k, v in args.items():
        print(f"{k}: {v}", file=sys.stderr)
    print("", file=sys.stderr)

    # Check input folder and create output folder if needed
    if format not in ["matrix folder", ".h5", ".h5ad"]:
        print(f"\033[91m* ERROR: unknown format {format} to convert, only allows h5ad, h5 or directory.\033[0m",
              file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args["input"]):
        print(f"\033[91m* ERROR: input file {args['input']} does not exist.\033[0m",
              file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(os.path.dirname(args["output"])):
        print(f"\033[91m* ERROR: output folder {os.path.dirname(args['output'])} does not exist.\033[0m",
              file=sys.stderr)
        sys.exit(1)
    if os.path.exists(args["output"]):
        print("\033[91m* ERROR: output file already exists.\n" + 
              "Please manually remove it to avoid accidental overwriting.\033[0m",
              file=sys.stderr)
        sys.exit(1)
    print(f"* Converting to {format}.", file=sys.stderr)

    import warnings
    from numba.core.errors import NumbaDeprecationWarning
    warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
    warnings.simplefilter(action="ignore", category=FutureWarning)

    import pegasus as pg
    import anndata

    # Load input file
    data = pg.read_input(args["input"])
    print("* Loaded object:\n", data, sep="", file=sys.stderr)

    # Write output file and test
    if format == ".h5ad":
        adata = Convert2Anndata(data)
        adata.write_h5ad(args["output"])
        Test(adata, data)
        print("* Matrix value extract for .X:", adata.X[adata.X > 5], file=sys.stderr)
        if adata.raw is not None:
            print("* Matrix value extract for .raw.X:", adata.raw.X[adata.raw.X > 5], file=sys.stderr)
            print("* Matrix value extract for .layers['raw.X']:", adata.layers["raw.X"][adata.layers["raw.X"] > 5],
                  file=sys.stderr)
    else: # .h5 or directory
        matrix_keys = data.list_keys("matrix")
        if "raw.X" in matrix_keys:
            data.select_matrix("raw.X")
        elif "counts" in matrix_keys:
            data.select_matrix("counts")
        else:
            print(f"\033[93m* WARNING: neither raw.X nor counts exists as matrix key, keeping current matrix.\033[0m",
                  file=sys.stderr)
        pg.write_output(data, args["output"])
        print("* Matrix value extract:", data.X[data.X > 5], file=sys.stderr)

    print("\033[92m* Conversion complete!\033[0m", file=sys.stderr)
