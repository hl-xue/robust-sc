import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

import os
import sys
import argparse
import pegasus as pg
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Sample preprocessing till the clustering step.")
    parser.add_argument("input", help="input file (.h5, .zarr, .zarr.zip, .h5ad, .h5ad.zip, matrix folder)")
    parser.add_argument("output", help="output folder")
    parser.add_argument("-p", "--initial_pc_num", default=100, type=int, help="initial PC number to reduce [INT]")
    parser.add_argument("-r", "--leiden_resolution", default=2.0, type=float, help="resolution for leiden clustering [FLOAT]")
    parser.add_argument("-c", "--clustering_colname", default="leiden_v00", type=str, help="column name of clustering results [STR]")
    parser.add_argument("-s", "--out_sample_label", default="sample", type=str, help="sample name as output file name prefix [STR]")
    parser.add_argument("-g", "--regress_cellcycle", action="store_true", help="whether to regress cell cycle")
    args = vars(parser.parse_args())
    for k, v in args.items():
        print(f"{k}: {v}", file=sys.stderr)
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
    if os.path.exists(f"{args['output']}{args['out_sample_label']}_tillClustering.zarr.zip"):
        print(f"\033[91m* ERROR: Output file already exists.\n" + 
              "         Please manually remove it to avoid accidental overwriting.\033[0m",
              file=sys.stderr)
        assert False

    # Load data and print to check
    data = pg.read_input(args["input"])
    print("* Filtered data object:", file=sys.stderr)
    print(data, file=sys.stderr)
    
    # Preprocessing till clustering
    ## Count normalisation, log1p transformation and HVG detection
    pg.identify_robust_genes(data)
    pg.log_norm(data)
    pg.highly_variable_features(data)
    pg.hvfplot(data, top_n=50, panel_size=(9, 5), dpi=600)
    plt.savefig(f"{args['output']}{args['out_sample_label']}_HVFplot.pdf", bbox_inches="tight", facecolor="white")
    plt.close()
    ## PCA
    pg.pca(data, n_components=args["initial_pc_num"], random_state=0)
    pg.elbowplot(data, rep="pca", panel_size=(5, 3), dpi=600)
    plt.savefig(f"{args['output']}{args['out_sample_label']}_elbowPlot.pdf", bbox_inches="tight", facecolor="white")
    plt.close()
    ncomps = data.uns["pca_ncomps"]
    if ncomps == args["initial_pc_num"]:
        print(f"\033[93m* WARNING: optimised PC number equals to indicated initial number.\n" +
              "           Please increase the --initial_pc_num [INT] argument to ensure the best optimisation.\033[0m",
              file=sys.stderr)
    print(f"* Recommanded to use {ncomps} PCs.", file=sys.stderr)
    ## Compute cell cycle score and regress it out if indicated
    pg.calc_signature_score(data, "cell_cycle_human")
    if args["regress_cellcycle"]:
        print("* Regressing out cell cycle score...", file=sys.stderr)
        pca_key = pg.regress_out(data, attrs=["G1/S", "G2/M"])
    else:
        pca_key = "pca"
    ## Neighbourhood graph, clustering and UMAP
    pg.neighbors(data, K=100, rep=pca_key, n_comps=ncomps, random_state=0)
    pg.leiden(data, rep=pca_key, resolution=args["leiden_resolution"], random_state=0, class_label=args["clustering_colname"])

    # Write output data
    pg.write_output(data, f"{args['output']}{args['out_sample_label']}_tillClustering.zarr.zip")

    print("=====> Preprocessing completed!\n", file=sys.stderr)
