"""
Load a zarr/zarr.zip file as input, and convert it to h5ad format.
Usage: python3 zarr2h5ad.py /path/to/input/file/name /path/to/output/file/name
Author: Haoliang Xue
"""

import sys
import os
import pegasus as pg

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
    print("* Example counts:", data.X[data.X > 10], file=sys.stderr)

    # Write output file
    data.to_anndata().write_h5ad(outfile)
    print("\033[92m* Conversion complete!", file=sys.stderr)
