
import os
import subprocess as sb

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.mkdir(file_path)

def peakcalling_pureclip(signal_bam, signal_bai, ctl_bam, ctl_bai, genome_fa, outputpath, pureclip_par, threads):

    outputpath_pureclip = outputpath + "/pureclip"
    ensure_dir(outputpath_pureclip)

    for i in range(0, len(signal_bam)):
        for j in range(0, len(ctl_bam)):

            signal_name = signal_bam[i].split("/").pop()
            ctl_name = ctl_bam[i].split("/").pop()

            crosslinks_file = "{}/{}_vs_{}_crosslinkind_sites.bed".format(outputpath_pureclip, signal_name, ctl_name)
            binding_regions = "{}/{}_vs_{}_binding_regions.bed".format(outputpath_pureclip, signal_name, ctl_name)
            pureclip_parameter_file = "{}/{}_vs_{}_parameters.txt".format(outputpath_pureclip, signal_name, ctl_name)

            sb.Popen("pureclip -i {} -bai {} -ibam {} -ibai {} -g {} -o {} -or {} -p {} -nt {} {}".format(
                signal_bam[i], signal_bai[i], ctl_bam[j], ctl_bai[j], genome_fa, crosslinks_file, binding_regions,
                pureclip_parameter_file, threads, pureclip_par
            ), shell=True).wait()

    # # collect all peaks from pureclip
    # sb.Popen("cat {}/*_binding_regions.bed > {}/pureclip_peaks_raw.bed".format(outputpath_pureclip, outputpath), shell=True).wait()

    # # extend peak regions by n basepairs
    # n = 10
    # sb.Popen("bedtools slop {} -i {}/pureclip_peaks_raw.bed -g {} > {}/pureclip_peaks_extended.bed".format(n, outputpath, genome_fa, outputpath), shell=True).wait()

    # # adding peakid to peakfile
    # file_raw = open("{}/pureclip_peaks_extended.bed".format(outputpath), "r")
    # file_final = open("{}/pureclip_peaks.bed".format(outputpath), "w")
    #
    # i = 1
    # for line in file_raw:
    #     data = line.split("\t")
    #     file_final.write("{1} \t {2} \t {3} \t pureclip_{4} \t {5} \t {6} \n".format(
    #                      data[0], data[1], data[3], i, data[5], data[6]))
    # file_final.close()

def peakcalling_piranha(signal_bam, ctl_bam, outputpath, piranha_par=""):

    ##################
    #### Piranha #####
    ##################

    outputpath_piranha = outputpath + "/piranha"
    ensure_dir(outputpath_piranha)

    for i in range(0, len(signal_bam)):
        signal_name = signal_bam[i].split("/").pop()
        sb.Popen("Piranha {} {} -o {}/{}_piranha_peaks_signal.bed".format(piranha_par,
                signal_bam[i], outputpath_piranha, signal_name), shell=True).wait()

    for i in range(0, len(ctl_bam)):
        ctl_name = ctl_bam[i].split("/").pop()
        sb.Popen("Piranha {} {} -o {}/{}_piranha_peaks_control.bed".format(piranha_par,
                 ctl_bam[i], outputpath_piranha, ctl_name), shell=True).wait()

    # intersect peaks from control and peaks from signal to filter false positives
    sb.Popen("cat {}/*_piranha_peaks_signal.bed > {}/piranha_peaks_signal.bed".format(outputpath_piranha, outputpath_piranha)).wait()
    sb.Popen("cat {}/*_piranha_peaks_control.bed > {}/piranha_peaks_control.bed".format(outputpath_piranha, outputpath_piranha)).wait()
    sb.Popen("bedtools intersect -a {}/piranha_peaks_signal.bed -b ".format(outputpath_piranha) +
              "{}/piranha_peaks_control.bed -s -v > {}/piranha_peaks_raw.bed".format(outputpath_piranha, outputpath)).wait()

    # # extend peak regions by n basepairs
    # n = 10
    # sb.Popen("bedtools slop {} -i {}/piranha_peaks_raw.bed -g {} > {}/piranha_peaks_extended.bed".format(n, outputpath, genome_fa, outputpath)).wait()
    #
    # # adding peakid to peakfile
    # file_raw = open("{}/piranha_peaks_extended.bed".format(outputpath), "r")
    # file_final = open("{}/piranha_peaks.bed".format(outputpath), "w")
    #
    # i = 1
    # for line in file_raw:
    #     data = line.split("\t")
    #     file_final.write("{1} \t {2} \t {3} \t piranha_{4} \t {5} \t {6} \n".format(
    #                      data[0], data[1], data[3], i, data[5], data[6]))
    # file_final.close()


def peakcalling_peakachu(signal_bam, ctl_bam, outputpath, peakachu_par=""):

    ###################
    #### PEAKachu #####
    ###################

    outputpath_peakachu = outputpath + "/peakachu"
    ensure_dir(outputpath_peakachu)

    argumentlist_signal = " ".join(signal_bam)
    argumentlist_ctl = " ".join(ctl_bam)

    sb.Popen("peakachu adaptive --exp_libs {} --ctr_libs {} --output_folder {} {}".format(argumentlist_signal,
                                                            argumentlist_ctl, outputpath_peakachu, peakachu_par)).wait()

    # generate outputfile peaks_peakachu.bed
    sb.Popen("sed -n 1p ./tmp_output/peak_tables/*csv > {}/peakachu_peaks.tsv".format(outputpath)).wait()
    sb.Popen("tail -n +2 -q ./tmp_output/peak_tables/*.csv >> {}/peakachu_peaks.tsv".format(outputpath)).wait()
    sb.Popen('''awk '/peak/ {print $0}' ./tmp_output/peak_annotations/*.gff > {}/peakachu_peaks.gff'''.format(outputpath_peakachu)).wait()
    sb.Popen("mv ./tmp_output/plots/Initial*.png {}/ma.png".format(outputpath_peakachu)).wait()

    # # adding peakid to peakfile
    # file_raw = open("{}/peakachu_peaks.gff".format(outputpath), "r")
    # file_final = open("{}/peakachu_peaks.bed".format(outputpath), "w")
    #
    # i = 1
    # for line in file_raw:
    #     data = line.split("\t")
    #     file_final.write("{1} \t {2} \t {3} \t peakachu_{4} \t {5} \t {6} \n".format(
    #                      data[0], data[3], data[4], i, data[5], data[6]))
    # file_final.close()

