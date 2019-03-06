
# The script will be used to:
# a) Generete a Venndiagram to compare the individual peakcallers
# b) Generate a Robust-Peak-File with all peaks validated for their
#    robustness (Np/Nt-ratio)

import subprocess as sb
import matplotlib.pyplot as plt
import os
import random
import numpy
import scipy.stats as sci
import time
import statsmodels.stats.multitest as multi
import multiprocessing

from matplotlib_venn import venn3_unweighted


def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.mkdir(file_path)

def get_line_count(file):
    count = 0
    for line in file:
        count += 1
    return count

def add_peak_to_dict(peaks, peak_dict):
    file = open(peaks, "r")
    i = 1
    if ( len(peak_dict.keys()) != 0 ):
        i += len(peak_dict.keys())
    for line in file:
        id = "peak_" + str(i)
        peak_dict[id] = line
        i += 1
    file.close()

def generate_count_file_single(flist, label1, label2, outputpath):
    options = "-s -c"
    out = "{}/{}_vs_{}.bed".format(outputpath, label1, label2)
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(flist[0], flist[1], options, out), shell=True).wait()

def generate_count_file(flist, label, outputpath):
    options = "-s -c"

    ab_out = "{}/{}_vs_b.bed".format(outputpath, label)
    ac_out = "{}/{}_vs_c.bed".format(outputpath, label)

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(flist[0], flist[1], options, ab_out), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(flist[0], flist[2], options, ac_out), shell=True).wait()

def get_intersect_counts(label, outputpath):

    ab_out = "{}/{}_vs_b.bed".format(outputpath, label)
    ac_out = "{}/{}_vs_c.bed".format(outputpath, label)

    peak_dict = dict()

    # read in files from before
    file_ab = open(ab_out, "r")
    file_ac = open(ac_out, "r")

    # adding counts together
    for line in file_ab:
        data = line.split("\t")

        # +1 to account for the peak itself
        if data[3] not in peak_dict:
            peak_dict[data[3]] = int(data[6])

    for line in file_ac:
        data = line.split("\t")

        # +1 to account for the peak itself
        if data[3] not in peak_dict:
            peak_dict[data[3]] = int(data[6])
        else:
            peak_dict[data[3]] = peak_dict[data[3]] + int(data[6])

    return(peak_dict)

# Generate pseudo_pool file from random choices of the peaks pool.
# Choices are with replacement (like int bootstrapping).
def generate_pseudo_pool(peaks, peak_keys_list, output_file):
    num_peaks = len(peak_keys_list)

    # Do some random choices, where the sample size is as big as the set size.
    random_peaks = [random.choice(peak_keys_list) for x in range(0, int(num_peaks))]

    # Generate pseudo_pool file.
    pseudo_pool_file = open(output_file, "w")

    for key in random_peaks:
        # write the peak (random choice)
        pseudo_pool_file.write(peaks[key])

    pseudo_pool_file.close()

# Adds counts of different dictionaries together
def add_pool_counts(pool_dict, sample_dict):
    for key in sample_dict:
        if key not in pool_dict:
            pool_dict[key] = sample_dict[key]
        else:
            pool_dict[key] = pool_dict[key] + sample_dict[key]

# calculate a p-value for robust peaks
def idr(outputpath, seed, num_pools, threads):
    random.seed(seed)

    # object needed to share data with processes
    manager = multiprocessing.Manager()

    # INPUT
    peakachu = "{}/peakachu_peaks.bed".format(outputpath)
    pureclip = "{}/pureclip_peaks.bed".format(outputpath)
    piranha = "{}/piranha_peaks.bed".format(outputpath)

    ## OUTPUT
    outputpath_idr = outputpath + "/idr"
    ensure_dir(outputpath_idr)

    outputpath_tmp = outputpath_idr + "/tmp"
    ensure_dir(outputpath_tmp)

    all_peaks = "{}/robust_peaks.bed".format(outputpath)

    ## Nt = number of peaks consistent between true replicates, we know were the peaks came from, specfic reproducibility
    ## High count only occur for peaks that poverlaps with the same peaks of the other peakcaller or with peaks that
    ## fall slightly in the interval

    print("[NOTE] Calculate Nt")

    generate_count_file([peakachu, pureclip, piranha], "peakachu", outputpath_idr)
    generate_count_file([piranha, peakachu, pureclip], "piranha", outputpath_idr)
    generate_count_file([pureclip, peakachu, piranha], "pureclip", outputpath_idr)

    peakachu_peak_dict = get_intersect_counts("peakachu", outputpath_idr)
    piranha_peak_dict = get_intersect_counts("piranha", outputpath_idr)
    pureclip_peak_dict = get_intersect_counts("pureclip", outputpath_idr)

    Nt_dict = dict()
    Nt_dict.update(peakachu_peak_dict)
    Nt_dict.update(piranha_peak_dict)
    Nt_dict.update(pureclip_peak_dict)

    ## Np = number of peaks conistent between pseudoreps, we are oblivious were the peak came from, generel reproducibility
    ## High number of overlaps occur when the same peaks occur quite often and not just because it overlaps just by chance with other peaks
    ## So you test for the general rbustness regardless of the peakcaller
        # - make a giant pool of peaks (all PEAKachu, Piranha, PureCLIP peaks)
        # then generate two random pools (rep1, rep2, each holding n=0.1*total number of peaks) from this giant pool
        # - if a peak is randomly selected more than once, then just merge the two instances, else I run into problems of the 
        # intersection ... If a peak occurs more than once, it wouldnt work ... And allowing to pick a peak only once, would decrease
        # the pool size with every choice which would include a bias for the other peaks with every step
        # - look of how often a peak intersect between these two pools
        # - do the random generation 1000 to get a distribution for each peak

    print("[NOTE] Calculate Np")

    # generate a big pool library
    all_peaks_raw_dict = dict()

    add_peak_to_dict(peakachu, all_peaks_raw_dict)
    add_peak_to_dict(piranha, all_peaks_raw_dict)
    add_peak_to_dict(pureclip, all_peaks_raw_dict)

    all_peaks_raw_peak_key_list = list(all_peaks_raw_dict.keys())

    Np_Nt_dict = dict()

    for key in Nt_dict:
        if key not in Np_Nt_dict:
            Np_Nt_dict[key] = [0] * num_pools

    # pools_list = [x for x in range(0, num_pools)]
    # pool = multiprocessing.Pool(int(threads))
    # for i in pools_list:
    #     pool.apply_async(np_nt_quotient_calculation, args=(i, Np_Nt_dict, Nt_dict, all_peaks_raw_dict,
    #                                                        all_peaks_raw_peak_key_list, outputpath_tmp))
    # pool.close()
    # pool.join()

    for i in range(0, num_pools):

        if ( i % 100 == 0 ):
            print("... generate pools " + str(i))

        pseudo_pool_1 = "{}/pool1.bed".format(outputpath_tmp)
        pseudo_pool_2 = "{}/pool2.bed".format(outputpath_tmp)
        pseudo_pool_3 = "{}/pool3.bed".format(outputpath_tmp)

        generate_pseudo_pool(all_peaks_raw_dict, all_peaks_raw_peak_key_list, pseudo_pool_1)
        generate_pseudo_pool(all_peaks_raw_dict, all_peaks_raw_peak_key_list, pseudo_pool_2)
        generate_pseudo_pool(all_peaks_raw_dict, all_peaks_raw_peak_key_list, pseudo_pool_3)

        labels = ["pool1", "pool2", "pool3"]
        tuple_files = [ [[pseudo_pool_1, pseudo_pool_2], [pseudo_pool_1, pseudo_pool_3]],
                        [[pseudo_pool_2, pseudo_pool_1], [pseudo_pool_2, pseudo_pool_3]],
                        [[pseudo_pool_3, pseudo_pool_1], [pseudo_pool_3, pseudo_pool_2]]]

        #start = time.time()
        pool = multiprocessing.Pool(int(threads))
        for j in range(0, len(labels)):
            pool.apply_async(generate_count_file_single, args=(tuple_files[j][0], labels[j], "b", outputpath_tmp))
            pool.apply_async(generate_count_file_single, args=(tuple_files[j][1], labels[j], "c", outputpath_tmp))
        pool.close()
        pool.join()
        #end = time.time()
        #print(end - start)

        pseudo_pool_1_dict = get_intersect_counts("pool1", outputpath_tmp)
        pseudo_pool_2_dict = get_intersect_counts("pool2", outputpath_tmp)
        pseudo_pool_3_dict = get_intersect_counts("pool3", outputpath_tmp)

        Np_dict = dict()
        add_pool_counts(Np_dict, pseudo_pool_1_dict)
        add_pool_counts(Np_dict, pseudo_pool_2_dict)
        add_pool_counts(Np_dict, pseudo_pool_3_dict)

        for key in Np_dict:
            if (Nt_dict[key] != 0):
                Np_Nt_dict[key][i] = Np_dict[key] / Nt_dict[key]
            else:
                Np_Nt_dict[key][i] = 0

    # Only peak with Nt_dict[key] != 0 are true robust peaks.
    # The feature Nt_dict[key] = 0 corresponds to peaks which only appear for the individual peakcaller.
    true_peaks = list()
    for key in Nt_dict:
        if( Nt_dict[key] != 0 ):
            true_peaks.append(key)

    ## Np/Nt < p-val
    ## A robust peak should have a high quotient
    print("[NOTE] Calculate p-value")

    # generate p-value dict and distribution of the means
    pval_dict = dict()
    mean_dict = dict()
    std_dict = dict()

    for key in true_peaks:
        mu, std = sci.norm.fit(Np_Nt_dict[key])
        mean_dict[key] = mu
        std_dict[key] = std

    casted_list_means = list(mean_dict.values())

    print("... print plot")
    # Plot the distribution of means.
    f = plt.figure()
    plt.hist(casted_list_means, bins=100, density=True, alpha=0.6, color='g')
    xmin, xmax = plt.xlim()
    x = numpy.linspace(xmin, xmax, 100)
    mu_means, std_means = sci.norm.fit(casted_list_means)
    p = sci.norm.pdf(x, mu_means, std_means)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu_means, std_means)
    plt.title(title)
    f.savefig(outputpath_idr + "/expected_mean_pdf_of_Np_Nt_quotient.pdf", bbox_inches='tight')

    print("... perform one sided ttest")
    for key in true_peaks:
        # p/2 --> to get one-sided ttest
        # I just want the significance that my quotient is bigger than the expected quotient. This
        # represents a higher and more robust peaks.
        if ( mean_dict[key] >= mu_means ):
            pval_dict[key] = (sci.ttest_ind(Np_Nt_dict[key], casted_list_means)[1])/2
        else:
            pval_dict[key] = 1.0

    # Plot distribution of p-vals
    print("... plot pval distribution")
    f = plt.figure()
    plt.hist(pval_dict.values(), bins=50, density=True, alpha=0.6, color='g')
    f.savefig(outputpath_idr + "/pval_distribution.pdf", bbox_inches='tight')

    # Benjamin hochberg correctur for multiple hypothesis testing
    print("... perform BH correction")
    pvals_corrected = multi.multipletests(list(pval_dict.values()), alpha=0.05, method='fdr_bh',
                                          is_sorted=False, returnsorted=False)

    i = 0
    pvals_corrected_dict = dict()
    for key in pval_dict:
        pvals_corrected_dict[key] = pvals_corrected[1][i]
        i += 1

    # Plot distribution of corrected p-vals
    f = plt.figure()
    plt.hist(list(pvals_corrected_dict.values()), bins=50, density=True, alpha=0.6, color='g')
    f.savefig(outputpath_idr + "/corrected_pval_distribution.pdf", bbox_inches='tight')

    # write total output file
    print("[NOTE] Generate output")
    file_robust_peaks = open(all_peaks, "w")

    for peak in all_peaks_raw_peak_key_list:
        line = all_peaks_raw_dict[peak].strip("\n")
        id = line.split("\t")[3]

        # For every peak that is not listed in the mean list, just do a line with a pval of 1.0.
        # These lines are the peaks that are not considered for the expected mean distribution and
        # p-value correction. These peaks do not overlap with any other peak from the other peakcallers.
        if ( id in mean_dict ):
            file_robust_peaks.write("{}\t{}\t{}\t{}\n".format(line, mean_dict[id], pval_dict[id], pvals_corrected_dict[id]))
        else:
            file_robust_peaks.write("{}\t{}\t{}\t{}\n".format(line, "0", "1.0", "1.0"))
    file_robust_peaks.close()

#TODO make parellel version
#TODO safe parameter list
def peakcaller_comparison(outputpath):
    print("[NOTE] Run Comparison")

    options="-wa -s -u"

    # INPUT
    peakachu="{}/peakachu_peaks.bed".format(outputpath)
    pureclip="{}/pureclip_peaks.bed".format(outputpath)
    piranha="{}/piranha_peaks.bed".format(outputpath)

    ## OUTPUT
    # for parallelisation
    # outputpath_comparison = outputpath + "/peak_comparison_{i}"
    outputpath_comparison = outputpath + "/peak_comparison"
    ensure_dir(outputpath_comparison)

    peakachu_pureclip = "{}/peakachu_pureclip_peaks.bed".format(outputpath_comparison)
    peakachu_piranha = "{}/peakachu_piranha_peaks.bed".format(outputpath_comparison)
    piranha_pureclip = "{}/piranha_pureclip_peaks.bed".format(outputpath_comparison)
    all = "{}/peakachu_piranha_pureclip_peaks.bed".format(outputpath_comparison)

    ## RUN
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(peakachu, pureclip, options, peakachu_pureclip), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(peakachu, piranha, options, peakachu_piranha), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(piranha, pureclip, options, piranha_pureclip), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(peakachu_piranha, pureclip, options, all), shell=True).wait()

    ## STATS
    f_peakachu = open(peakachu, "r")
    f_piranha = open(piranha, "r")
    f_pureclip = open(pureclip, "r")
    f_peakachu_pureclip = open(peakachu_pureclip, "r")
    f_piranha_pureclip = open(piranha_pureclip, "r")
    f_peakachu_piranha = open(peakachu_pureclip, "r")
    f_all = open(all, "r")

    np_peakachu = get_line_count(f_peakachu)
    np_piranha = get_line_count(f_piranha)
    np_pureclip = get_line_count(f_pureclip)
    np_peakachu_pureclip = get_line_count(f_peakachu_pureclip)
    np_piranha_pureclip = get_line_count(f_piranha_pureclip)
    np_peakachu_piranha = get_line_count(f_peakachu_piranha)
    np_all = get_line_count(f_all)

    f_peakachu.close()
    f_piranha.close()
    f_pureclip.close()
    f_peakachu_pureclip.close()
    f_piranha_pureclip.close()
    f_peakachu_piranha.close()
    f_all.close()

    stats_file = open("{}/stats_comparison.txt".format(outputpath_comparison), "w")
    stats_file.write("peakcaller\tpeakcount\n")
    stats_file.write("PEAKachu\t{}\n".format(np_peakachu))
    stats_file.write("Piranha\t{}\n".format(np_piranha))
    stats_file.write("PureCLIP\t{}\n".format(np_pureclip))
    stats_file.write("PEAKachu_PureCLIP\t{}\n".format(np_peakachu_pureclip))
    stats_file.write("Piranha_PureCLIP\t{}\n".format(np_piranha_pureclip))
    stats_file.write("PEAKachu_Piranha\t{}\n".format(np_peakachu_piranha))
    stats_file.write("ALL\t{}\n".format(np_all))
    stats_file.close()

    ## VENN DIAGRAM
    #(A, B, AB, C, AC, BC, ABC)
    f = plt.figure()
    venn3_unweighted(subsets=(np_peakachu, np_piranha, np_peakachu_piranha,
                              np_pureclip, np_peakachu_pureclip, np_piranha_pureclip, np_all),
                     set_labels=('PEAKachu', 'Piranha', 'PureCLIP'),
                     set_colors=('r', 'y', 'b'))
    plt.show()
    f.savefig(outputpath_comparison + "/comparison_plot.pdf", bbox_inches='tight')
    return(np_all)


