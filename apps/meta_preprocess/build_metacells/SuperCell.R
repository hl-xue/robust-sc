rm (list = ls())

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description="Metacell construction.")
parser$add_argument("input", nargs=1, help="input file in .h5 format")
parser$add_argument("output", nargs=1, help="output table in .csv format")
parser$add_argument("-n", "--num.hvg", default=1000, type="integer", help="number of HVGs to be used [1000]")
parser$add_argument("-p", "--num.pc", default=50, type="integer", help="number of PCs to be used [50]")
parser$add_argument("-g", "--gamma", default=15, type="double", help="graining level [15]")
args <- parser$parse_args()

for (arg in c("input", "output", "num.hvg", "num.pc", "gamma")) {
    cat(arg, ": ", args[[arg]], "\n", sep = "", file = stderr())
}
cat("\n", file = stderr())

if (!file_test("-f", args$input)) {
    stop("\033[91minput file ", args$input, " does not exist.\033[0m")
}
if (file_test("-f", args$output)) {
    stop("\033[91moutput file already exists.\nPlease manually remove it to avoid accidental overwriting.\033[0m")
}
if (!stringr::str_ends(args$output, ".csv")) {
    cat("\033[93mWarning: output file", args$output, "does not have .csv extension.\033[0m\n")
}

suppressPackageStartupMessages(library(SuperCell))
suppressPackageStartupMessages(library(Seurat))

sc <- Read10X_h5(h5.in) %>%
    CreateSeuratObject(names.delim = "-") %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "disp", nfeatures = args$num.hvg)
hvgs <- VariableFeatures(sc)

mc <- SCimplify(X = GetAssayData(sc), genes.use = hvgs, gamma = args$gamma, n.pc = args$num.pc)
df.membership <- data.frame(mc$membership) %>%
    dplyr::rename(metacell = mc.membership) %>%
    tibble::rownames_to_column("barcodekey")
write.csv(df.membership, file = csv.out, quote = FALSE, row.names = FALSE)
