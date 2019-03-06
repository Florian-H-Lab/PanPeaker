
import os
import subprocess as sb

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.mkdir(file_path)

def peakcalling_pureclip(signal_bam, signal_bai, genome_fa, chr_sizes, outputpath, pureclip_par, threads):

    outputpath_pureclip = outputpath + "/pureclip"
    ensure_dir(outputpath_pureclip)

    # Go over each signal and each control sample and call crosslinking sites and
    # binding binding regions with PureCLIP. Use binding regions as peaks.
    print("... do PureCLIP")
    for i in range(0, len(signal_bam)):

        signal_name = signal_bam[i].split("/").pop()

        crosslinks_file = "{}/{}_crosslinkind_sites.bed".format(outputpath_pureclip, signal_name)
        binding_regions = "{}/{}_binding_regions.bed".format(outputpath_pureclip, signal_name)
        pureclip_parameter_file = "{}/{}_parameters.txt".format(outputpath_pureclip, signal_name)

        sb.Popen("pureclip -i {} -bai {} -g {} -o {} -or {} -p {} -nt {} {}".format(
            signal_bam[i], signal_bai[i], genome_fa, crosslinks_file, binding_regions,
            pureclip_parameter_file, threads, pureclip_par
        ), shell=True).wait()

    print("... edit results")

    # Collect all peaks from pureclip.
    # Do not merge replicates (intersect). I want all possible peaks called from replicates, even those
    # which might only appear in just one replicate. Because of the IDR later on, peaks that appear only in one
    # replicate might still be robust when we will look if the peak overlaps with other peakcaller's results.
    sb.Popen("cat {0}/*_binding_regions.bed > {0}/pureclip_peaks_raw.bed".format(outputpath_pureclip), shell=True).wait()

    # extend peak regions by n basepairs
    n = 10
    sb.Popen("bedtools slop -header -b {0} -i {1}/pureclip_peaks_raw.bed -g {2} > {1}/pureclip_peaks_extended.bed".format(
              n, outputpath_pureclip, chr_sizes), shell=True).wait()

    # adding peakid to peakfile
    file_raw = open("{}/pureclip_peaks_extended.bed".format(outputpath_pureclip), "r")
    file_final = open("{}/pureclip_peaks.bed".format(outputpath), "w")

    i = 1
    for line in file_raw:
        data = line.strip("\n").split("\t")
        file_final.write("{0}\t{1}\t{2}\tpureclip_{3}\t{4}\t{5}\n".format(
                         data[0], data[1], data[2], i, data[4], data[5]))
        i += 1
    file_final.close()
    file_raw.close()

def peakcalling_piranha(signal_bam, outputpath, piranha_par, chr_sizes):

    ##################
    #### Piranha #####
    ##################

    outputpath_piranha = outputpath + "/piranha"
    ensure_dir(outputpath_piranha)

    print("... do Piranha")

    # call Piranha for signal samples
    for i in range(0, len(signal_bam)):
        signal_name = signal_bam[i].split("/").pop()
        sb.Popen("Piranha {} {} -o {}/{}_piranha_peaks_signal.bed".format(piranha_par,
                signal_bam[i], outputpath_piranha, signal_name), shell=True).wait()

    print("... edit results")

    # Intersect peaks from control and peaks from signal to filter false positives.
    # Do not merge replicates of signal and control (intersect). I want all possible peaks called from replicates, even those
    # which might only appear in just one replicate. Because of the IDR later on, peaks that appear only in one
    # replicate might still be robust when we will look if the peak overlaps with other peakcaller's results.

    # first cat signals together
    sb.Popen("cat {0}/*_piranha_peaks_signal.bed > {0}/piranha_peaks_raw.bed".format(outputpath_piranha), shell=True).wait()

    # extend peak regions by n basepairs
    n = 10
    sb.Popen("bedtools slop -header -b {0} -i {1}/piranha_peaks_raw.bed -g {2} > {1}/piranha_peaks_extended.bed".format(
             n, outputpath_piranha, chr_sizes), shell=True).wait()

    # adding peakid to peakfile
    file_raw = open("{}/piranha_peaks_extended.bed".format(outputpath_piranha), "r")
    file_final = open("{}/piranha_peaks.bed".format(outputpath), "w")

    i = 1
    for line in file_raw:
        data = line.strip("\n").split("\t")
        file_final.write("{0}\t{1}\t{2}\tpiranha_{3}\t{4}\t{5}\n".format(
                         data[0], data[1], data[2], i, data[6], data[5]))
        i+=1
    file_final.close()
    file_raw.close()

def peakcalling_peakachu(signal_bam, outputpath, peakachu_par, chr_sizes):

    ###################
    #### PEAKachu #####
    ###################

    outputpath_peakachu = outputpath + "/peakachu"
    ensure_dir(outputpath_peakachu)

    argumentlist_signal = " ".join(signal_bam)

    print("... do PEAKachu")

    # call PEAKachu
    sb.Popen("peakachu adaptive --exp_libs {} --output_folder {} {}".format(argumentlist_signal,
            outputpath_peakachu, peakachu_par), shell=True).wait()

    print("... edit results")

    # generate outputfile peaks_peakachu.bed
    sb.Popen("sed -n 1p {0}/peak_tables/*csv > {0}/peakachu_peaks.tsv".format(outputpath_peakachu), shell=True).wait()
    sb.Popen("tail -n +2 -q  {0}/peak_tables/*.csv >> {0}/peakachu_peaks.tsv".format(outputpath_peakachu), shell=True).wait()
    sb.Popen("cat {0}/peak_annotations/*.gff | grep peak > {0}/peakachu_peaks.gff".format(outputpath_peakachu), shell=True).wait()

    # extend peak regions by n basepairs
    n = 10
    sb.Popen("bedtools slop -header -b {} -i {}/peakachu_peaks.gff -g {} > {}/peakachu_peaks_extended.gff".format(
             n, outputpath_peakachu, chr_sizes, outputpath_peakachu), shell=True).wait()

    # adding peakid to peakfile
    file_raw = open("{}/peakachu_peaks_extended.gff".format(outputpath_peakachu), "r")
    file_final = open("{}/peakachu_peaks.bed".format(outputpath), "w")

    i = 1
    for line in file_raw:
        data = line.strip("\n").split("\t")
        file_final.write("{0}\t{1}\t{2}\tpeakachu_{3}\t{4}\t{5}\n".format(
                         data[0], data[3], data[4], i, data[5], data[6]))
        i += 1
    file_final.close()
    file_raw.close()

