
import argparse
import logging
import os

import numpy

import matplotlib.pyplot as plt
import peak_calling_with_control as call_ctl
import peak_calling_without_control as call
import peakcaller_comparison as peakcomparison
import panpeaker_refinement as refi

#########################
##   NECESSARY TOOLS   ##
#########################

# PEAKachu
# Piranha
# PureCLIP
# bedtools
# matplotlib_venn

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
parser.add_argument(
    "--para_sets",
    metavar='*.txt',
    help="Paths to the parameter set file.")
parser.add_argument(
    "--seed",
    metavar='int',
    default=0,
    help="Set seed for IDR calculation.")
parser.add_argument(
    "--refinement",
    action='store_true',
    help="Activate the refinement module. Panpeaker will merge significant peaks and creates an additional plot")
parser.add_argument(
    "--adj_pval",
    metavar='float',
    default=0.05,
    help="Set the adjusted p-value (BH) threshold for the refinement module, to filer peaks.")

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

##TODO CHECK if bam and bai array have the same length

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

final_parameterset_dict = dict()
final_parameterset_dict["PEAKachu"] = args.peakachu.strip("\"")
final_parameterset_dict["Piranha"] = args.piranha.strip("\"")
final_parameterset_dict["PureCLIP"] = args.pureclip.strip("\"")

if ( args.para_sets ):

    print("[NOTE] Generate list of paramter sets.")
    para_dict = dict()
    para_dict["PEAKachu"] = list()
    para_dict["Piranha"] = list()
    para_dict["PureCLIP"] = list()

    parameterset_file = open(args.para_sets, "r")
    headline = parameterset_file.readline()
    headline_list = headline.strip("\n").split("\t")

    print(headline_list)

    for line in parameterset_file:
        para_sets = line.strip("\n").split("\t")

        for i in range(0, len(para_sets)):
            if ( para_sets[i] == "none" ):
                if ( headline_list[i] == "PEAKachu" ):
                    para_dict[headline_list[i]].append(args.peakachu.strip("\""))
                if (headline_list[i] == "Piranha"):
                    para_dict[headline_list[i]].append(args.piranha.strip("\""))
                if (headline_list[i] == "PureCLIP"):
                    para_dict[headline_list[i]].append(args.pureclip.strip("\""))
            else:
                para_dict[headline_list[i]].append(para_sets[i])

    num_parametersets = len(para_dict["PEAKachu"])

    num_robust_peaks_list = [-1] * num_parametersets

    # Generate a log file for the parameter testing in case it breaks down for a set.
    log_parameter_finding_file = open(args.output_folder + "/log_parameter_finding_file.txt", "w")
    log_parameter_finding_file.write(headline.strip("\n") + "\tNumber of Robust Peaks\n")
    log_parameter_finding_file.close()

    print("[NOTE] Find best parameter set.")
    for i in range(0, num_parametersets):
        print("... testing paramter set " + str(i+1))

        if( bool_control == 0 ):
            print("[NOTE] Running PureCLIP")
            call_ctl.peakcalling_pureclip(args.input_signal_bam, args.input_signal_bai, args.input_control_bam,
                                          args.input_control_bai, args.genome_file, args.chr_sizes, args.output_folder,
                                          para_dict["PureCLIP"][i], args.threads)

            print("[NOTE] Running Piranha")
            call_ctl.peakcalling_piranha(args.input_signal_bam, args.input_control_bam,
                                         args.output_folder, para_dict["Piranha"][i], args.chr_sizes)

            print("[NOTE] Running PEAKachu")
            call_ctl.peakcalling_peakachu(signal_bam_files_peakachu, control_bam_files_peakachu,
                                         args.output_folder, para_dict["PEAKachu"][i], args.chr_sizes)
        else:
            print("[NOTE] Running PureCLIP")
            call.peakcalling_pureclip(args.input_signal_bam, args.input_signal_bai, args.genome_file,
                                      args.chr_sizes, args.output_folder, para_dict["PureCLIP"][i], args.threads)

            print("[NOTE] Running Piranha")
            call.peakcalling_piranha(args.input_signal_bam, args.output_folder,
                                     para_dict["Piranha"][i], args.chr_sizes)

            print("[NOTE] Running PEAKachu")
            call.peakcalling_peakachu(signal_bam_files_peakachu, control_bam_files_peakachu,
                                      para_dict["PEAKachu"][i], args.chr_sizes)

        num_robust_peaks_list[i] = peakcomparison.peakcaller_comparison(args.output_folder)

        # create line
        next_line_for_log_file = ""
        for j in range(0, len(headline_list)):
            next_line_for_log_file = next_line_for_log_file + para_dict[headline_list[j]][i] + "\t"

        log_parameter_finding_file = open(args.output_folder + "/log_parameter_finding_file.txt", "a")
        log_parameter_finding_file.write("{}\t{}\n".format(next_line_for_log_file, num_robust_peaks_list[i]))
        log_parameter_finding_file.close()

    print("... generate barplot of the number of robust peaks.")

    y_pos = numpy.arange(len(num_robust_peaks_list))
    bars = ["Set" + str(x+1) for x in range(0, len(num_robust_peaks_list))]
    f = plt.figure()
    plt.bar(y_pos, num_robust_peaks_list)
    plt.xticks(y_pos, bars)
    plt.xlabel('Parameter Sets')
    plt.ylabel('Number of Robust Peaks (In All Peakcallers)')
    plt.show()
    f.savefig(args.output_folder + "/parameter_set_plot.pdf", bbox_inches='tight')

    print("... save best parameter set.")
    best_set_index = numpy.argmax(num_robust_peaks_list)
    final_parameterset_dict["PEAKachu"] = para_dict["PEAKachu"][best_set_index]
    final_parameterset_dict["Piranha"] = para_dict["Piranha"][best_set_index]
    final_parameterset_dict["PureCLIP"] = para_dict["PureCLIP"][best_set_index]

    log_parameter_finding_file = open(args.output_folder + "/log_parameter_finding_file.txt", "a")
    log_parameter_finding_file.write("Using Set:\t{}\t{}\t{}\n".format(para_dict["PEAKachu"][best_set_index],
                                                                 para_dict["Piranha"][best_set_index],
                                                                 para_dict["PureCLIP"][best_set_index]))
    log_parameter_finding_file.close()

print("[NOTE] Execute peak calling with best or chosen parameter set.")
if( bool_control == 0 ):
    print("[NOTE] Running PureCLIP")
    call_ctl.peakcalling_pureclip(args.input_signal_bam, args.input_signal_bai, args.input_control_bam,
                                  args.input_control_bai, args.genome_file, args.chr_sizes, args.output_folder,
                                  final_parameterset_dict["PureCLIP"], args.threads)

    print("[NOTE] Running Piranha")
    call_ctl.peakcalling_piranha(args.input_signal_bam, args.input_control_bam,
                                 args.output_folder, final_parameterset_dict["Piranha"], args.chr_sizes)

    print("[NOTE] Running PEAKachu")
    call_ctl.peakcalling_peakachu(signal_bam_files_peakachu, control_bam_files_peakachu,
                                 args.output_folder, final_parameterset_dict["PEAKachu"], args.chr_sizes)
else:
    print("[NOTE] Running PureCLIP")
    call.peakcalling_pureclip(args.input_signal_bam, args.input_signal_bai, args.genome_file,
                              args.chr_sizes, args.output_folder, final_parameterset_dict["PureCLIP"],
                              args.threads)

    print("[NOTE] Running Piranha")
    call.peakcalling_piranha(args.input_signal_bam, args.output_folder,
                             final_parameterset_dict["Piranha"], args.chr_sizes)

    print("[NOTE] Running PEAKachu")
    call.peakcalling_peakachu(signal_bam_files_peakachu, control_bam_files_peakachu,
                              final_parameterset_dict["PEAKachu"], args.chr_sizes)


print("[NOTE] Peakcaller Comparison")
peakcomparison.peakcaller_comparison(args.output_folder)

print("[NOTE] IDR")
peakcomparison.idr(args.output_folder, args.seed, 1000, args.threads, len(args.input_signal_bam))

if (args.refinement):
   print("[NOTE] Run Refinement with adjusted p-value of <= {}".format(args.adj_pval))
   refi.refinement(args.output_folder, args.adj_pval)

print("[END]")