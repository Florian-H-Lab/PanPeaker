
import subprocess as sb
import matplotlib.pyplot as plt
import os
import random
import numpy
import scipy.stats as sci
import statsmodels.stats.multitest as multi

from matplotlib_venn import venn3

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.mkdir(file_path)

def generate_count_file(afile, bfile, cfile, label, outputpath):
    options = "-s -c"

    ab_out = "{}/{}_vs_b".format(outputpath, label)
    ac_out = "{}/{}_vs_c".format(outputpath, label)

    ## RUN
    peak_dict = dict()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(afile, bfile, options, ab_out)).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(afile, cfile, options, ac_out)).wait()

    # read in files from before
    file_ab = open(ab_out, "r")
    file_ac = open(ac_out, "r")

    # adding counts together
    for line in file_ab:
        data = line.split("\t")

        # +1 to account for the peak itself
        if data[3] not in peak_dict:
            peak_dict[data[3]] = data[6] + 1

    for line in file_ac:
        data = line.split("\t")

        # +1 to account for the peak itself
        if data[3] not in peak_dict:
            peak_dict[data[3]] = data[6] + 1
        else:
            peak_dict[data[3]] = peak_dict[data[3]] + data[6]

    return(peak_dict)

def generate_pseudo_pool(peaks, output_file):
    ## OUTPUT
    pseudo_pool_file = open(output_file, "w")

    ## RUN
    # randomly sample peaks
    num_peaks = sum(1 for line in open(peaks))

    if( num_peaks < 100 ):
        sb.Popen("cp {} {}".format(peaks, pseudo_pool_file),shell=True)
        return 0

    peaks_list = [x for x in range(0, num_peaks)]
    random_peaks = [random.choice(peaks_list) for x in range(0, int(0.1*num_peaks))]
    random_peaks.sort()

    i = 0
    for line in peaks:
        if ( i in random_peaks ):
            pseudo_pool_file.write(line)

    pseudo_pool_file.close()


def add_pool_counts(pool_dict, sample_dict):
    for key in sample_dict:
        if key not in pool_dict:
            pool_dict[key] = sample_dict[key]
        else:
            pool_dict[key] = pool_dict[key] + sample_dict[key]

# calculate a p-value for robust peaks
def idr(outputpath, seed, n):
    random.seed(seed)

    # INPUT
    peakachu = "{}/peakachu_peaks.bed".format(outputpath)
    pureclip = "{}/pureclip_peaks.bed".format(outputpath)
    piranha = "{}/piranha_peaks.bed".format(outputpath)

    ## OUTPUT
    outputpath_idr = outputpath + "/idr"
    ensure_dir(outputpath_idr)

    outputpath_tmp = outputpath_idr + "/tmp"
    ensure_dir(outputpath_idr)

    pseudo_pool_1 = "{}/p1.bed".format(outputpath_tmp)
    pseudo_pool_2 = "{}/p2.bed".format(outputpath_tmp)
    pseudo_pool_3 = "{}/p3.bed".format(outputpath_tmp)

    all_peaks = "{}/robust_peaks.bed".format(outputpath)

    ## Nt = number of peaks consistent between true replicates, we know were the peaks came from, specfic reproducibility
    ## High count only occur for peaks that poverlaps with the same peaks of the other peakcaller or with peaks that
    ## fall slightly in the interval

    peakachu_peak_dict = generate_count_file(peakachu, pureclip, piranha, "peakachu", outputpath_idr)
    piranha_peak_dict = generate_count_file(piranha, peakachu, pureclip, "piranha", outputpath_idr)
    pureclip_peak_dict = generate_count_file(pureclip, peakachu, piranha, "pureclip", outputpath_idr)

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

    # generate a big pool file
    sb.Popen("cat {} {} {} > {}".format(peakachu, piranha, pureclip, all_peaks), shell=True)

    Np_Nt_dict = dict()

    for key in Nt_dict:
        if key not in Np_Nt_dict:
            Np_Nt_dict[key] = [-1] * n

    # to the random pooling n-times
    for i in range (0, n):
        generate_pseudo_pool(all_peaks, pseudo_pool_1)
        generate_pseudo_pool(all_peaks, pseudo_pool_2)
        generate_pseudo_pool(all_peaks, pseudo_pool_3)

        pseudo_pool_1_dict = generate_count_file(pseudo_pool_1, pseudo_pool_2, pseudo_pool_3, "pool1", outputpath_tmp)
        pseudo_pool_2_dict = generate_count_file(pseudo_pool_2, pseudo_pool_1, pseudo_pool_3, "pool2", outputpath_tmp)
        pseudo_pool_3_dict = generate_count_file(pseudo_pool_3, pseudo_pool_1, pseudo_pool_2, "pool3", outputpath_tmp)

        Np_dict = dict()
        add_pool_counts(Np_dict, pseudo_pool_1_dict)
        add_pool_counts(Np_dict, pseudo_pool_2_dict)
        add_pool_counts(Np_dict, pseudo_pool_3_dict)

        for key in Np_dict:
            Np_Nt_dict[key][i] = Np_dict[key] / Nt_dict[key]

    ## Np/Nt < p-val 
    ## A robust peak should have a high quotient
    standard_data = sci.norm.rvs(3., 1., size=1000)
    expected_norm_dist = sci.norm.fit(standard_data)

    # Plot the histogram.
    plt.hist(expected_norm_dist, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = numpy.linspace(xmin, xmax, 100)
    p = sci.norm.pdf(x, 3., 1.)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (3., 1.)
    plt.title(title)

    # generate p-value dict and distribution of the means
    pval_dict = dict()
    mean_list = [-1] * len(Np_Nt_dict.keys())

    i = 0
    for key in Np_Nt_dict:
        mu, std = sci.norm.fit(Np_Nt_dict[key])
        mean_list[i] = mu
        i += 1
        pval_dict[key] = sci.ttest_ind(Np_Nt_dict[key],standard_data)[1]

    # Plot distirbution of means
    plt.hist(mean_list, bins=25, density=True, alpha=0.6, color='g')
    xmin, xmax = plt.xlim()
    x = numpy.linspace(xmin, xmax, 100)
    mu, std = sci.norm.fit(mean_list)
    p = sci.norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)

    # Benjamin hochberg correctur for multiple hypothesis testing
    pvals_corrected_dict = dict()
    pvals_corrected = multi.multipletests(pval_dict.values(), alpha=0.05, method='fdr_bh',
                                          is_sorted=False, returnsorted=False)
    i = 0
    for key in pval_dict:
        pvals_corrected_dict[key] = pvals_corrected[i]
        i += 1

    # write total output file
    file_robust_peaks = open(all_peaks, "w")
    for line in all_peaks:
        peak_id = line.split("\t")[4]
        file_robust_peaks.write("{} \t {} \t {} \t {} \n".format(line.strip("\n"), Np_Nt_dict[peak_id],
                                                                 pval_dict[key], pvals_corrected_dict[key]))
    file_robust_peaks.close()

def peakcaller_comparison(outputpath):
    options="-wa -s -u"

    # INPUT
    peakachu="{}/peakachu_peaks.gtf".format(outputpath)
    pureclip="{}/pureclip_peaks.bed".format(outputpath)
    piranha="{}/piranha_peaks.bed".format(outputpath)

    ## OUTPUT
    outputpath_comparison = outputpath + "/peak_comparison"
    ensure_dir(outputpath_comparison)

    peakachu_pureclip = "{}/peakachu_pureclip_peaks.bed".format(outputpath_comparison)
    peakachu_piranha = "{}/peakachu_piranha_peaks.bed".format(outputpath_comparison)
    piranha_pureclip = "{}/piranha_pureclip_peaks.bed".format(outputpath_comparison)
    all = "{}/peakachu_piranha_pureclip_peaks.bed".format(outputpath_comparison)

    ## RUN
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(peakachu, pureclip, options, peakachu_pureclip)).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(peakachu, piranha, options, peakachu_piranha)).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(piranha, pureclip, options, piranha_pureclip)).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(peakachu_piranha, pureclip, options, all)).wait()

    ## STATS
    np_peakachu = sum(1 for line in open(peakachu))
    np_piranha = sum(1 for line in open(piranha))
    np_pureclip = sum(1 for line in open(pureclip))
    np_peakachu_pureclip = sum(1 for line in open(peakachu_pureclip))
    np_piranha_pureclip = sum(1 for line in open(piranha_pureclip))
    np_peakachu_piranha = sum(1 for line in open(peakachu_piranha))
    np_all = sum(1 for line in open(all))

    stats_file = open("{}/stats_comparison.txt".format(outputpath_comparison), "w")
    stats_file.write("peakcaller \t peakcount\n")
    stats_file.write("PEAKachu \t {}\n".format(np_peakachu))
    stats_file.write("Piranha \t {}\n".format(np_piranha))
    stats_file.write("PureCLIP \t {}\n".format(np_pureclip))
    stats_file.write("PEAKachu_PureCLIP \t {}\n".format(np_peakachu_pureclip))
    stats_file.write("Piranha_PureCLIP \t {}\n".format(np_piranha_pureclip))
    stats_file.write("PEAKachu_Piranha \t {}\n".format(np_peakachu_piranha))
    stats_file.write("ALL \t {}\n".format(np_all))
    stats_file.close()

    ## VENN DIAGRAM
    #(A, B, AB, C, AC, BC, ABC)

    f = plt.figure()
    venn3(subsets=(np_peakachu, np_piranha, np_peakachu_piranha,
                   np_pureclip, np_peakachu_pureclip, np_piranha_pureclip, np_all),
                   set_labels=('PEAKachu', 'Piranha', 'PureCLIP'))
    plt.show()
    f.savefig(outputpath_comparison + "/comparison_plot.pdf", bbox_inches='tight')


