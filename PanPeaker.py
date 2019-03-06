
import argparse
import logging
import os
import peak_calling_with_control as call_ctl
import peak_calling_without_control as call
import peakcaller_comparison as peakcomparison

#TODO allow for different bam files for PureCLIP and Piranha

#########################
##   NECESSARY TOOLS   ##
#########################

# PEAKachu
# Piranha
# PureCLIP
# bedtools

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
parser.add_argument(
    "--chr_sizes",
    metavar='*.txt',
    required=True,
    help="Path to genomes chromosome sizes (txt).")

# optional arguments
parser.add_argument(
    "-c", "--input_control_bam",
    metavar='*.bam',
    nargs="+",
    help="List of paths to the control bam files.")
parser.add_argument(
    "-k", "--input_control_bai",
    metavar='*.bai',
    nargs="+",
    help="List of paths to the control bai files.")
parser.add_argument(
    "-o", "--output_folder",
    metavar='path/',
    default=os.getcwd(),
    help="Write results to this path.")
parser.add_argument(
    "-s", "--output_name",
    metavar='prefix',
    help="Use this to name the output file.")
parser.add_argument(
    "--peakachu",
    metavar='...',
    default="--max_insert_size 200 --mad_multiplier 2.0 --fc_cutoff 2.0 --padj_threshold 0.05",
    help="Parameter list for PEAKachu")
parser.add_argument(
    "--piranha",
    metavar='...',
    default="-p 0.05",
    help="Parameter list for Piranha")
parser.add_argument(
    "--pureclip",
    metavar='...',
    default="",
    help="Parameter list for PureCLIP")
parser.add_argument(
    "-nt", "--threads",
    metavar='int',
    default="1",
    help="Threads for PanPeaker")
parser.add_argument(
    "--signal_bam_peakachu",
    metavar='*.bam',
    nargs="+",
    help="List of paths to the signal bam files for PEAKachu.")
parser.add_argument(
    "--control_bam_peakachu",
    metavar='*.bam',
    nargs="+",
    help="List of paths to the signal bam files for PEAKachu.")

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

print("[START]")

# check if controls are provided
bool_control = 1
if( args.input_control_bam ):
    bool_control = 0
    print("[NOTE] Controls provided")
else:
    print("[NOTE] Controls not provided")

# define input for PEAKachu base don user choice
signal_bam_files_peakachu = args.input_signal_bam
control_bam_files_peakachu = args.input_control_bam
if( args.signal_bam_peakachu ):
    signal_bam_files_peakachu = args.signal_bam_peakachu
if( args.control_bam_peakachu ):
    control_bam_files_peakachu = args.control_bam_peakachu

if( bool_control == 0 ):
    print("[NOTE] Running PureCLIP")
    call_ctl.peakcalling_pureclip(args.input_signal_bam, args.input_signal_bai, args.input_control_bam,
                                  args.input_control_bai, args.genome_file, args.chr_sizes, args.output_folder,
                                  args.pureclip.strip("\""), args.threads)

    print("[NOTE] Running Piranha")
    call_ctl.peakcalling_piranha(args.input_signal_bam, args.input_control_bam,
                                 args.output_folder, args.piranha.strip("\""), args.chr_sizes)

    print("[NOTE] Running PEAKachu")
    call_ctl.peakcalling_peakachu(args.input_signal_bam, args.input_control_bam,
                                 args.output_folder, args.peakachu.strip("\""), args.chr_sizes)
else:
    print("[NOTE] Running PureCLIP")
    call.peakcalling_pureclip(args.input_signal_bam, args.input_signal_bai, args.genome_file,
                              args.chr_sizes, args.output_folder, args.pureclip.strip("\""), args.threads)

    print("[NOTE] Running Piranha")
    call.peakcalling_piranha(args.input_signal_bam, args.output_folder, args.piranha.strip("\""), args.chr_sizes)

    print("[NOTE] Running PEAKachu")
    call.peakcalling_peakachu(args.input_signal_bam, args.output_folder, args.peakachu.strip("\""), args.chr_sizes)

peakcomparison.peakcaller_comparison(args.output_folder)

#TODO parameter optimization

peakcomparison.idr(args.output_folder, 123, 1000, args.threads)

print("[END]")