import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

import gc
import os
import sys
import argparse
import pegasus as pg
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Identify and select singlet barcodes.")
    parser.add_argument("input", help="input file (.h5, .zarr, .zarr.zip, .h5ad, .h5ad.zip, matrix folder)")
    parser.add_argument("output", help="output folder")
    parser.add_argument("-m", "--method", default="pegasus", type=str, help="singlet selection method [STR]")
    parser.add_argument("-s", "--out_sample_label", default="sample", type=str, help="sample name as output file name prefix [STR]")
    args = vars(parser.parse_args())
    for k, v in args.items():
        print(f"{k}: {v}", file=sys.stderr)
    print("", file=sys.stderr)
    args["method"] = args["method"].lower()
    if args["method"] not in ["pegasus", "scanpy", "any", "all"]:
        print(f"\033[91m* ERROR: unknown method ({args['method']}) for singlet identification.\n" +
              "         Please choose one from ['pegasus', 'scanpy', 'all', 'any'].\033[0m",
              file=sys.stderr)

    # Check input folder and create output folder if needed
    if not os.path.exists(args["input"]):
        print(f"\033[91m* ERROR: Input file {args['input']} does not exist.\033[0m", file=sys.stderr)
        assert False
    if not args["output"].endswith("/"):
        args["output"] = f"{args['output']}/"
    if not os.path.exists(args["output"]):
        print(f"* Creating output folder {args['output']}...", file=sys.stderr)
        os.mkdir(args["output"])

    # Stop if output file already exists, to avoid accidental overwriting
    if os.path.exists(f"{args['output']}{args['sample_name']}_singlets.zarr.zip"):
        print(f"\033[91m* ERROR: Output file already exists." + 
              "         Please manually remove it to avoid accidental overwriting.\033[0m",
              file=sys.stderr)
        assert False

    # Load data
    data = pg.read_input(args["input"])
    
    # Run with Pegasus
    data_test = data.copy()
    data_test.obs["Sample"] = args["out_sample_label"]
    pg.infer_doublets(data_test, channel_attr="Sample", random_state=0)
    plt.savefig(f"{args['output']}{args['out_sample_label']}_doublet_score_density.pdf",
                bbox_inches="tight", facecolor="white")
    plt.close()
    data_test.obs.rename(columns={"pred_dbl": "pred_dbl.pegasus", "doublet_score": 
                                  "pred_dbl_score.pegasus"}, inplace=True)
    
    # Run with Scanpy
    adata_test = data.to_anndata()
    sc.external.pp.scrublet(adata_test, random_state=0)
    adata_test.obs.rename(columns={"predicted_doublet": "pred_dbl.scanpy",
                                   "doublet_score": "pred_dbl_score.scanpy"}, inplace=True)
    
    # Summarise and compare prediction
    assert (data_test.obs.index == adata_test.obs.index).all()
    pred_tab = pd.concat([data_test.obs[["pred_dbl_score.pegasus", "pred_dbl.pegasus"]],
                          adata_test.obs[["pred_dbl_score.scanpy", "pred_dbl.scanpy"]]],
                         axis=1, verify_integrity=True)
    assert (pred_tab.index == data.obs.index).all()
    del data_test, adata_test
    gc.collect()
    cmp_tab = pd.crosstab(pred_tab["pred_dbl.pegasus"], pred_tab["pred_dbl.scanpy"])
    cmp_tab.to_csv(f"{args['output']}{args['out_sample_label']}_prediction_contingency.csv")
    
    # Assign prediction results
    if args["method"] in ["pegasus", "scanpy"] :
        data.obs["pred_dbl_score"] = pred_tab[f"pred_dbl_score.{args['method']}"].tolist()
        data.obs["pred_dbl"] = pred_tab[f"pred_dbl.{args['method']}"].tolist()
    elif args["method"] == "all": # marking as doublet if identified by all methods
        data.obs["pred_dbl_score"] = pred_tab[["pred_dbl_score.pegasus", "pred_dbl_score.scanpy"]].min(axis=1)
        data.obs["pred_dbl"] = (pred_tab["pred_dbl.pegasus"] & pred_tab["pred_dbl.scanpy"])
    elif args["method"] == "any": # marking as doublet if identified by any method
        data.obs["pred_dbl_score"] = pred_tab[["pred_dbl_score.pegasus", "pred_dbl_score.scanpy"]].max(axis=1)
        data.obs["pred_dbl"] = (pred_tab["pred_dbl.pegasus"] | pred_tab["pred_dbl.scanpy"])
    else:
        assert False
           
    # Write outputs
    ## All cells with filtration marks
    pg.write_output(data, f"{args['output']}{args['sample_name']}_full.zarr.zip")
    ## Filtered cells
    pg.write_output(data[~data.obs["pred_dbl"], :].copy(),
                    f"{args['output']}{args['sample_name']}_singlets.zarr.zip")

    print(f"=====> Singlet selection of {args['sample_name']} has been completed!\n", file=sys.stderr)
