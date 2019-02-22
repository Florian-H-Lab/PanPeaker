
import os

def peakcalling(signal_bam, signal_bai, ctl_bam, ctl_bai, genome_fa,
                outputpath, outputname, peakachu_par, piranha_par, pureclip_par, threads):

    ###################
    #### PURECLIP #####
    ###################

    outputpath_pureclip = outputpath + "/pureclip"
    tmp_folder = outputpath_pureclip + "/tmp/"

    for i in range(0, len(signal_bam)):
        for j in range(0, len(ctl_bam)):

            signal_name = signal_bam[i].split("/").pop()
            ctl_name = ctl_bam[i].split("/").pop()

            crosslinks_file = "{}/{}_vs_{}_crosslinkind_sites.bed".format(outputpath_pureclip, signal_name, ctl_name)
            binding_regions = "{}/{}_vs_{}_binding_regions.bed".format(outputpath_pureclip, signal_name, ctl_name)
            pureclip_parameter_file = "{}/{}_vs_{}_parameters.txt".format(outputpath_pureclip, signal_name, ctl_name)

            os.system("pureclip -i {0} -bai {1} -ibam {2} -ibai {3} -g {4} -o {5} -tmp {6} -or {7} -p {8} -nt {9} {10}".format(
                signal_bam[i], signal_bai[i], ctl_bam[j], ctl_bai[j], genome_fa, crosslinks_file, tmp_folder, binding_regions,
                pureclip_parameter_file, threads, pureclip_par
            ))

    ##################
    #### Piranha #####
    ##################

    outputpath_piranha = outputpath + "/piranha"

    for i in range(0, len(signal_bam)):
        signal_name = signal_bam[i].split("/").pop()
        os.system("Piranha {} {} -o {}/{}_piranha_peaks.bed".format(piranha_par,  signal_bam[i], outputpath_piranha, signal_name))

    for i in range(0, len(ctl_bam)):
        ctl_name = ctl_bam[i].split("/").pop()
        os.system("Piranha {} {} -o {}/{}_piranha_peaks.bed".format(piranha_par, ctl_bam[i], outputpath_piranha, ctl_name))

    ## intersect peaks from control and peaks from signal to filter false positives


    ###################
    #### PEAKachu #####
    ###################

    outputpath_peakachu = outputpath + "/peakachu"

    argumentlist_signal = " ".join(signal_bam)
    argumentlist_ctl = " ".join(ctl_bam)

    os.system("peakachu adaptive --exp_libs {} --ctr_libs {} --output_folder {} {}".format(argumentlist_signal,
                                                            argumentlist_ctl, outputpath_peakachu, peakachu_par))

    ### generate outputfile peaks_peakachu.bed