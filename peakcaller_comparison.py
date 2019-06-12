
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

from simple_venn import venn3


def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.mkdir(file_path)

def get_line_count(file):
    f = open(file, "r")
    count = 0
    for line in f:
        count += 1
    f.close()
    return count


def add_peaks_to_dict(peaks, peak_dict, peaklength_dict):
    file = open(peaks, "r")
    for line in file:
        data = line.split("\t")
        peaklength_dict[data[3]] = abs(int(data[2]) - int(data[1]))
        peak_dict[data[3]] = line
    file.close()

def generate_count_file_single(flist, label1, label2, outputpath):
    options = "-s -c"
    out = "{}/{}_vs_{}.bed".format(outputpath, label1, label2)
    # uniq --> Remove multiple lines (peaks) in each file because of the bootstrapping.
    # They need to be removed !after! the counting to calculate the quotient correctly.
    # -k 4 --> column four are the peak ids
    sb.Popen("bedtools intersect -a {} -b {} {} | sort -k 4 | uniq > {}".format(flist[0], flist[1], options, out),
             shell=True).wait()

# Count hof often a peak occurs in different peaksets of different peakcallers,
# thus account for robustness of different replciates and different peaks.
def get_variance_counts_dict(label, outputpath):

    ab_out = "{}/{}_vs_b.bed".format(outputpath, label)
    ac_out = "{}/{}_vs_c.bed".format(outputpath, label)

    peak_dict = dict()

    # read in files from before
    file_ab = open(ab_out, "r")
    file_ac = open(ac_out, "r")

    # adding counts together
    for line in file_ab:
        data = line.split("\t")

        if data[3] not in peak_dict:
            peak_dict[data[3]] = int(data[6])

    for line in file_ac:
        data = line.split("\t")

        if data[3] not in peak_dict:
            peak_dict[data[3]] = int(data[6])
        else:
            peak_dict[data[3]] = peak_dict[data[3]] + int(data[6])

    return(peak_dict)

def generate_intersect_file(flist, label, outputpath):
    options = "-s -c"

    ab_out = "{}/{}_vs_b.bed".format(outputpath, label)
    ac_out = "{}/{}_vs_c.bed".format(outputpath, label)

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(flist[0], flist[1], options, ab_out), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(flist[0], flist[2], options, ac_out), shell=True).wait()

# Just count the occurence of a peak in different peakcallers
def get_bias_counts_dict(label, outputpath):

    ab_out = "{}/{}_vs_b.bed".format(outputpath, label)
    ac_out = "{}/{}_vs_c.bed".format(outputpath, label)

    peak_dict = dict()

    # read in files from before
    file_ab = open(ab_out, "r")
    file_ac = open(ac_out, "r")

    # adding counts together
    for line in file_ab:
        data = line.split("\t")

        if data[3] not in peak_dict:
            peak_dict[data[3]] = 1

        if int(data[6]) != 0:
            peak_dict[data[3]] = peak_dict[data[3]] + 1

    for line in file_ac:
        data = line.split("\t")

        if data[3] not in peak_dict:
            peak_dict[data[3]] = 1

        if int(data[6]) != 0:
            peak_dict[data[3]] = peak_dict[data[3]] + 1

    return(peak_dict)

# Generate pseudo_pool file from random choices of the peaks pool.
# Choices are with replacement (like int bootstrapping), to tackle replicates.
def generate_pseudo_pool(peaks, peak_keys_list, output_file):

    num_peaks = len(peak_keys_list)

    # Do some random choices, where the sample size is as big as the set size.
    random_peaks = [random.choice(peak_keys_list) for x in range(0, num_peaks)]

    # Generate pseudo_pool file.
    pseudo_pool_file = open(output_file, "w")

    for key in random_peaks:
        # a) Write the random chosen peak.
        # b) To a random shifting of the peaks to tackle also peaks that intersect
        # to two peak just because is lies inside the borders of two peaks.
        data = peaks[key].strip("\n").split("\t")
        shift = random.randint(0,10)
        start = int(data[1]) + shift
        end = int(data[2]) + shift
        pseudo_pool_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(data[0], str(start), str(end),
                                                                data[3], data[4], data[5]))
    pseudo_pool_file.close()

# Adds counts of different dictionaries together
def add_pool_counts(pool_dict, sample_dict):
    for key in sample_dict:
        if key not in pool_dict:
            pool_dict[key] = sample_dict[key]
        else:
            pool_dict[key] = pool_dict[key] + sample_dict[key]

# calculate a p-value for robust peaks
def idr(outputpath, seed, num_pools, threads, num_replicates):
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

    generate_intersect_file([peakachu, pureclip, piranha], "peakachu", outputpath_idr)
    generate_intersect_file([piranha, peakachu, pureclip], "piranha", outputpath_idr)
    generate_intersect_file([pureclip, peakachu, piranha], "pureclip", outputpath_idr)

    Nt_pekachu_dict = get_bias_counts_dict("peakachu", outputpath_idr)
    Nt_piranha_dict = get_bias_counts_dict("piranha", outputpath_idr)
    Nt_pureclip_dict = get_bias_counts_dict("pureclip", outputpath_idr)

    Nt_dict = dict()
    Nt_dict.update(Nt_pekachu_dict)
    Nt_dict.update(Nt_piranha_dict)
    Nt_dict.update(Nt_pureclip_dict)

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

    peaklength_dict = dict()

    # generate a big pool library
    peakachu_peaks_dict = dict()
    add_peaks_to_dict(peakachu, peakachu_peaks_dict, peaklength_dict)
    peakachu_peak_key_list = list(peakachu_peaks_dict.keys())

    piranha_peaks_dict = dict()
    add_peaks_to_dict(piranha, piranha_peaks_dict, peaklength_dict)
    piranha_peak_key_list = list(piranha_peaks_dict.keys())

    pureclip_peaks_dict = dict()
    add_peaks_to_dict(pureclip, pureclip_peaks_dict, peaklength_dict)
    pureclip_peak_key_list = list(pureclip_peaks_dict.keys())

    Np_Nt_dict = dict()

    for key in Nt_dict:
        if key not in Np_Nt_dict:
            Np_Nt_dict[key] = [0] * num_pools

    for i in range(0, num_pools):

        if ( i % 100 == 0 ):
            print("... generate pools " + str(i))

        pseudo_pool_1 = "{}/pool_peakachu.bed".format(outputpath_tmp)
        pseudo_pool_2 = "{}/pool_piranha.bed".format(outputpath_tmp)
        pseudo_pool_3 = "{}/pool_pureclip.bed".format(outputpath_tmp)

        generate_pseudo_pool(peakachu_peaks_dict, peakachu_peak_key_list, pseudo_pool_1)
        generate_pseudo_pool(piranha_peaks_dict, piranha_peak_key_list, pseudo_pool_2)
        generate_pseudo_pool(pureclip_peaks_dict, pureclip_peak_key_list, pseudo_pool_3)

        labels = ["pool_peakachu", "pool_piranha", "pool_pureclip"]
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

        pseudo_pool_1_dict = get_variance_counts_dict("pool_peakachu", outputpath_tmp)
        pseudo_pool_2_dict = get_variance_counts_dict("pool_piranha", outputpath_tmp)
        pseudo_pool_3_dict = get_variance_counts_dict("pool_pureclip", outputpath_tmp)

        Np_dict = dict()
        add_pool_counts(Np_dict, pseudo_pool_1_dict)
        add_pool_counts(Np_dict, pseudo_pool_2_dict)
        add_pool_counts(Np_dict, pseudo_pool_3_dict)

        # Take the length of the peak into account. Broad peaks intersectin with lots of
        # other peaks from different peakcallers are not good.
        for key in Np_dict:
            Np_Nt_dict[key][i] = (Np_dict[key] * Nt_dict[key] * 1000000) / peaklength_dict[key]

    # First constraint: Only peak with Nt_dict[key] != 1 are true robust peaks.
    # The feature Nt_dict[key] = 1 corresponds to peaks which only appear for the individual peakcaller.
    # Second constraint: Takes peaks out that are shared between peakcallers, but does not intersect with
    # other peaks if you do the bootstrapping, so the peak is not good supported by other peaks (replicates,
    # number of peaks) of the other peakcallers.
    true_peaks = list()
    for key in Np_Nt_dict:
        if( Nt_dict[key] != 1 and not all(v == 0 for v in Np_Nt_dict[key]) ):
            true_peaks.append(key)

    # Create a distribution for Nt_dict[key] == 1 which is my negative model
    false_peaks = list()
    for key in Nt_dict:
        if( Nt_dict[key] == 1 and not all(v == 0 for v in Np_Nt_dict[key]) ):
            false_peaks.append(key)

    ## Np/Nt < p-val
    ## A robust peak should have a high quotient
    print("[NOTE] Calculate p-value")

    # generate p-value dict and distribution of the means
    mean_dict = dict()

    for key in true_peaks:
        mean_dict[key] = numpy.mean(Np_Nt_dict[key])

    casted_list_means = list(mean_dict.values())
    log_casted_list_means = numpy.log1p(casted_list_means)

    print("... print plot")
    # Plot the distribution of means.
    f = plt.figure()
    plt.hist(log_casted_list_means, bins=100, density=True, alpha=0.6, color='g')
    xmin, xmax = plt.xlim()
    x = numpy.linspace(xmin, xmax, 100)
    mu_means, std_means = sci.norm.fit(log_casted_list_means)
    p = sci.norm.pdf(x, mu_means, std_means)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu_means, std_means)
    plt.title(title)
    f.savefig(outputpath_idr + "/expected_mean_pdf_of_Np_Nt_quotient.pdf", bbox_inches='tight')

    mean_false_dict = dict()

    for key in false_peaks:
        mean_false_dict[key] = numpy.mean(Np_Nt_dict[key])

    casted_list_means_false = list(mean_false_dict.values())
    log_casted_list_means_false = numpy.log1p(casted_list_means_false)

    # Plot the distribution of means of the false sample set.
    f_false = plt.figure()
    plt.hist(log_casted_list_means_false, bins=100, density=True, alpha=0.6, color='g')
    xmin, xmax = plt.xlim()
    x_false = numpy.linspace(xmin, xmax, 100)
    mu_means_false, std_means_false = sci.norm.fit(log_casted_list_means_false)
    p_false = sci.norm.pdf(x_false, mu_means_false, std_means_false)
    plt.plot(x_false, p_false, 'k', linewidth=2)
    title_false = "Fit results: mu = %.2f,  std = %.2f" % (mu_means_false, std_means_false)
    plt.title(title_false)
    f_false.savefig(outputpath_idr + "/expected_mean_false_pdf_of_Np_Nt_quotient.pdf", bbox_inches='tight')

    print("... perform one sided ttest")
    #all_pvals_dict = dict()
    true_pvals_dict = dict()
    for key in true_peaks:
        # p/2 --> to get one-sided ttest
        # I just want the significance that my quotient is bigger than the expected quotient. This
        # represents a higher and more robust peaks.
        if ( mean_dict[key] >= mu_means_false ):
        #    all_pvals_dict[key] = (sci.ttest_ind(Np_Nt_dict[key], casted_list_means_false)[1])/2
            true_pvals_dict[key] = (sci.ttest_ind(numpy.log1p(Np_Nt_dict[key]), log_casted_list_means_false)[1])/2
        #else:
        #    all_pvals_dict[key] = 1.0

    # Plot distribution of p-vals
    print("... plot pval distribution")
    f = plt.figure()
    plt.hist(list(true_pvals_dict.values()), bins=50, density=True, alpha=0.6, color='g')
    f.savefig(outputpath_idr + "/pval_distribution.pdf", bbox_inches='tight')

    # Benjamin hochberg correctur for multiple hypothesis testing
    print("... perform BH correction")
    pvals_corrected = multi.multipletests(list(true_pvals_dict.values()), alpha=0.05, method='fdr_bh',
                                          is_sorted=False, returnsorted=False)

    i = 0
    pvals_corrected_dict = dict()
    for key in true_pvals_dict:
        pvals_corrected_dict[key] = pvals_corrected[1][i]
        i += 1

    # Plot distribution of corrected p-vals
    f = plt.figure()
    plt.hist(list(pvals_corrected_dict.values()), bins=50, density=True, alpha=0.6, color='g')
    f.savefig(outputpath_idr + "/corrected_pval_distribution.pdf", bbox_inches='tight')

    # write total output file
    print("[NOTE] Generate output")
    file_robust_peaks = open(all_peaks, "w")

    # whole entire set of peaks
    all_peaks_raw_dict = dict()
    all_peaks_raw_dict.update(peakachu_peaks_dict)
    all_peaks_raw_dict.update(piranha_peaks_dict)
    all_peaks_raw_dict.update(pureclip_peaks_dict)

    for peak in all_peaks_raw_dict:
        line = all_peaks_raw_dict[peak].strip("\n")
        id = line.split("\t")[3]

        # For every peak that is not listed in the mean list, just do a line with a pval of 1.0.
        # These lines are the peaks that are not considered for the expected mean distribution and
        # p-value correction. These peaks do not overlap with any other peak from the other peakcallers.
        if ( id in true_pvals_dict ):
            file_robust_peaks.write("{}\t{}\t{}\t{}\n".format(line, mean_dict[id], true_pvals_dict[id], pvals_corrected_dict[id]))
        else:
            file_robust_peaks.write("{}\t{}\t{}\t{}\n".format(line, "0", "1.0", "1.0"))
    file_robust_peaks.close()

def peakcaller_comparison(outputpath):
    print("[NOTE] Run Comparison")

    options="-s -u"

    # INPUT
    n_a = "{}/peakachu_peaks.bed".format(outputpath)
    n_b = "{}/piranha_peaks.bed".format(outputpath)
    n_c = "{}/pureclip_peaks.bed".format(outputpath)

    ## OUTPUT
    # for parallelisation
    # outputpath_comparison = outputpath + "/peak_comparison_{i}"
    outputpath_comparison = outputpath + "/peak_comparison"
    ensure_dir(outputpath_comparison)

    n_ab = "{}/peakachu_piranha_peaks.bed".format(outputpath_comparison)
    n_ac = "{}/peakachu_pureclip_peaks.bed".format(outputpath_comparison)

    n_ba = "{}/piranha_peakachu_peaks.bed".format(outputpath_comparison)
    n_bc = "{}/piranha_pureclip_peaks.bed".format(outputpath_comparison)

    n_ca = "{}/pureclip_peakachu_peaks.bed".format(outputpath_comparison)
    n_cb = "{}/pureclip_piranha_peaks.bed".format(outputpath_comparison)

    n_abc = "{}/peakachu_piranha_pureclip_peaks.bed".format(outputpath_comparison)
    n_bac = "{}/piranha_peakachu_pureclip_peaks.bed".format(outputpath_comparison)
    n_cab = "{}/pureclip_peakachu_piranha_peaks.bed".format(outputpath_comparison)

    ## RUN
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_a, n_b, options, n_ab), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_a, n_c, options, n_ac), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_b, n_a, options, n_ba), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_b, n_c, options, n_bc), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_c, n_a, options, n_ca), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_c, n_b, options, n_cb), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ab, n_c, options, n_abc), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ba, n_c, options, n_bac), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ca, n_b, options, n_cab), shell=True).wait()

    ## STATS
    abc = get_line_count(n_abc)
    bac = get_line_count(n_bac)
    cab = get_line_count(n_cab)

    ac = get_line_count(n_ac)
    ab = get_line_count(n_ab)
    ba = get_line_count(n_ba)
    bc = get_line_count(n_bc)
    ca = get_line_count(n_ca)
    cb = get_line_count(n_cb)

    a = get_line_count(n_a)
    b = get_line_count(n_b)
    c = get_line_count(n_c)

    stats_file = open("{}/stats_comparison.txt".format(outputpath_comparison), "w")
    stats_file.write("peakcaller\tpeakcount\n")
    stats_file.write("PEAKachu\t{}\n".format(a))
    stats_file.write("Piranha\t{}\n".format(b))
    stats_file.write("PureCLIP\t{}\n".format(c))
    stats_file.write("PEAKachu_Piranha\t{}\n".format(ab))
    stats_file.write("PEAKachu_PureCLIP\t{}\n".format(ac))
    stats_file.write("Piranha_PEAKachu\t{}\n".format(ba))
    stats_file.write("Piranha_PureCLIP\t{}\n".format(bc))
    stats_file.write("PureCLIP_PEAKachu\t{}\n".format(ca))
    stats_file.write("PureCLIP_PEAKachu\t{}\n".format(cb))
    stats_file.write("ALL_PEAKachu\t{}\n".format(abc))
    stats_file.write("ALL_Piranha\t{}\n".format(bac))
    stats_file.write("ALL_PureCLIP\t{}\n".format(cab))
    stats_file.close()

    # Adjust values for venn diagram
    ab -= abc
    ac -= abc

    ba -= bac
    bc -= bac

    ca -= cab
    cb -= cab

    a -= (abc + ab + ac)
    b -= (bac + ba + bc)
    c -= (cab + ca + cb)

    np_all = abc + bac + cab

    #data = [np_peakachu, np_piranha, np_peakachu_piranha, np_pureclip, np_peakachu_pureclip, np_piranha_pureclip, np_all]
    a = "A:" + str(a)
    b = "B:" + str(b)
    c = "C:" + str(c)
    ab = "A:" + str(ab) + "\nB:" + str(ba)
    ac = "A:" + str(ac) + "\nC:" + str(ca)
    bc = "B:" + str(bc) + "\nC" + str(cb)
    abc = "A:" + str(abc) + "\nB:" + str(bac) + "\nC:" + str(cab)

    data = [a, b, c, ab, ac, bc, abc]

    ## VENN DIAGRAM
    #[A, B, C, AB, AC, BC, ABC]
    f = plt.figure(figsize=(8, 8))
    plt.rcParams["font.family"] = "serif"
    venn3(subsets=data, set_labels=('PEAKachu', 'Piranha', 'PureCLIP'), set_colors=('r', 'y', 'b'))
    plt.show()
    f.savefig(outputpath_comparison + "/comparison_plot.pdf", bbox_inches='tight', dpi=300)
    return(np_all)