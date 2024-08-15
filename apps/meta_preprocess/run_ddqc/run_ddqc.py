import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

import sys
import os
import pegasus as pg
import ddqc
import argparse


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Sample quality control using ddqc.")
    parser.add_argument("input", help="input file (.h5, .zarr, .zarr.zip, .h5ad, .h5ad.zip, matrix folder).")
    parser.add_argument("output", help="output folder.")
    parser.add_argument("-u", "--min_umi_per_cell", default=0, type=int, help="minimum UMI count per cell [INT].")
    parser.add_argument("-g", "--min_gene_per_cell", default=0, type=int, help="minimum expressed gene number per cell [INT].")
    parser.add_argument("-M", "--max_mito_per_cell", default=100, type=int, help="maximum mitochondrial gene UMI percentage per cell [INT].")
    parser.add_argument("-s", "--sample_name", default="sample", type=str, help="sample name to write in output folder [STR].")
    args = vars(parser.parse_args())
    for k, v in args.items():
        print(k, v, file=sys.stderr)
    print("", file=sys.stderr)

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
    if os.path.exists(f"{args['output']}{args['sample_name']}_filtered.zarr.zip"):
        print(f"\033[91m* ERROR: Output file already exists. Please manually remove it to avoid accidental overwriting.\033[0m",
              file=sys.stderr)
        assert False

    # Processing by ddqc
    ## Run ddqc
    data = pg.read_input(args["input"])
    df_qc = ddqc.ddqc_metrics(data, clustering_method="leiden", random_state=0,
                              ribo_prefix="^RP[SL]|^RPLP|^RPSA", return_df_qc=True,
                              display_plots=False)
    df_qc.to_csv(f"{args['output']}/{args['sample_name']}_dfqc_table.csv")
    ## Checking df_qc with data object and integrate df_qc into data.obs
    for col in set(data.obs.columns).intersection(set(df_qc.columns)):
        if col in ["percent_mito", "percent_ribo"]:
            assert (data.obs[col] - df_qc[col]).max() < 1e-5
        else:
            assert (data.obs[col] == df_qc[col]).all()
    data.obs.rename(columns={"passed_qc": "passed_ddqc"}, inplace=True)

    # Post-ddqc manual filtration
    data.obs["manual_rm"] = (data.obs["n_counts"] < args["min_umi_per_cell"]) | \
        (data.obs["n_genes"] < args["min_gene_per_cell"]) | \
        (data.obs["percent_mito"] > args["max_mito_per_cell"])
    data.obs["passed_qc"] = data.obs["passed_ddqc"] & (~data.obs["manual_rm"])
    
    # Write outputs
    ## All cells with filtration marks
    pg.write_output(data, f"{args['output']}{args['sample_name']}_full.zarr.zip")
    ## Filtered cells
    pg.write_output(data[data.obs["passed_qc"], :].copy(), f"{args['output']}{args['sample_name']}_filtered.zarr.zip")
    
    print(f"=====> Adaptive QC of {args['sample_name']} has been completed!\n", file=sys.stderr)
