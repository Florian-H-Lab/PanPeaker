
import os
import argparse
import logging
import with_control

####################
##   ARGS INPUT   ##
####################

tool_description = """
The tool needs bam files which are correctly sorted.
Provide a list of signal (CLIP) files and a list of 
control (Input, IgG, GFP, or others) files. The tool
generates a list of roust peaks together with a Venn diagram 
shared between all differen peak calling algorithms, 
located in the folder robust_peaks. All peaks found by 
a pairwise comparison of the supported peakcaller are 
found in pairwise_peaks. PanPeaker generates further a 
list of peaks for each algorithm, found under the 
individual name (PEAKachu, Piranha, and PureCLIP).

You can provide for PEAKachu, Piranha and PureCLIP
parameters. If not, then PanPeaker will use the default
settings of each algorithm. 
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description, usage='%(prog)s [options]',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

# version
parser.add_argument(
    "-v", "--version", action="version", version="%(prog)s 0.0")

# positional arguments
parser.add_argument(
    "-i", "--input_signal_bam",
    metavar='*.bam',
    required=True,
    nargs="+",
    help="List of paths to the signal bam files.")
parser.add_argument(
    "-b", "--input_signal_bai",
    metavar='*.bai',
    required=True,
    nargs="+",
    help="List of paths to the signal bai files.")
parser.add_argument(
    "-g", "--genome_file",
    metavar='*.fa',
    required=True,
    help="Path to genome file (fasta).")

# optional arguments
parser.add_argument(
    "-c", "--input_control_bam",
    metavar='*.bam',
    nargs="+",
    help="List of paths to the control bam files.")
parser.add_argument(
    "-k" "--input_control_bai",
    metavar='*.bai',
    nargs="+",
    help="List of paths to the control bai files.")
parser.add_argument(
    "-o", "--output_folder",
    metavar='path/',
    help="Write results to this path.")
parser.add_argument(
    "-s", "--output_name",
    metavar='prefix',
    help="Use this to name the output file.")
parser.add_argument(
    "--peakachu",
    metavar='...',
    help="Parameter list for PEAKachu")
parser.add_argument(
    "--piranha",
    metavar='...',
    help="Parameter list for Piranha")
parser.add_argument(
    "--pureclip",
    metavar='...',
    help="Parameter list for PureCLIP")
parser.add_argument(
    "--threads",
    metavar='int',
    default='1',
    help="Threads for PureCLIP")

######################
##   CHECKS INPUT   ##
######################

parser.add_argument(
    "-d", "--debug",
    help="Print lots of debugging information",
    action="store_true")

args = parser.parse_args()
if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
else:
    logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
logging.info("Parsed arguments:")
logging.info("interval_file: '{}'".format(args.input_signal_bam))
if args.output_folder:
    logging.info("outfile: enabled writing to file")
    logging.info("outfile: '{}'".format(args.output_folder))
logging.info("")

###################
##   CODE BODY   ##
###################

bool_control = 1
if( args.input_control_bam ):
    bool_control = 0
    print("[NOTE] Controls provided")
else:
    print("[NOTE] Controls not provided")

with_control.peakcalling(args.input_signal_bam, args.input_signal_bai, args.input_control_bam, args.input_control_bai,
                         args.genome_file, args.output_folder, args.output_name, args.peakachu, args.piranha,
                         args.pureclip, args.threads)