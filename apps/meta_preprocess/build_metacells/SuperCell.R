rm (list = ls())

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description="Metacell construction.")
parser$add_argument("input", nargs=1, help="input file in .h5 format")
parser$add_argument("out.prefix", nargs=1, help="output file prefix")
parser$add_argument("-f", "--feature.num", default=1000, type="integer", help="number of variable features to be used [1000]")
parser$add_argument("-c", "--pc.num", default=50, type="integer", help="number of principal components to be used [50]")
parser$add_argument("-k", "--neighbour.num", default=5, type="integer", help="number of neighbours for k-NN [5]")
parser$add_argument("-g", "--gamma", default=15, type="double", help="graining level [15]")
arg_lst <- parser$parse_args()

for (arg in c("input", "out.prefix", "feature.num", "pc.num", "neighbour.num", "gamma")) {
    cat(arg, ": ", arg_lst[[arg]], "\n", sep = "", file = stderr())
}
cat("\n", file = stderr())

if (!file_test("-f", arg_lst$input)) {
    stop("\033[91minput file ", arg_lst$input, " does not exist.\033[0m")
}
if (file_test("-f", paste0(arg_lst$out.prefix, ".csv")) || file_test("-f", paste0(arg_lst$out.prefix, ".RData"))) {
    stop("\033[91moutput file already exists.\n       Please manually remove it to avoid accidental overwriting.\033[0m")
}

suppressPackageStartupMessages(library(SuperCell))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(magrittr))

mc.smp <- Read10X_h5(arg_lst$input, use.names = FALSE) %>%
    CreateSeuratObject(names.delim = "-") %>%
    NormalizeData() %>%
    GetAssayData() %>%
    SCimplify(n.var.genes = arg_lst$hvg.num, gamma = arg_lst$gamma,
              k.knn = arg_lst$neighbour.num, n.pc = arg_lst$pc.num,
              do.approx = FALSE, seed = 0)
mc.smp$sample <- supercell_assign(rep(stringr::str_remove(basename(arg_lst$input), "\\.h5"),
                                      length(mc.smp$membership)),
                                  mc.smp$membership)

df.membership <- data.frame(mc.smp$membership) %>%
    dplyr::rename(metacell = mc.smp.membership) %>%
    tibble::rownames_to_column("barcodekey")

write.csv(df.membership,file = paste0(arg_lst$out.prefix, ".csv"),
          quote = FALSE, row.names = FALSE)
save(mc.smp, file = paste0(arg_lst$out.prefix, ".RData"))

cat("\033[92m* Processing complete!\n\033[0m", file = stderr())
