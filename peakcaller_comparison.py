
import subprocess as sb
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

def peakcaller_comparison(outputpath):
    options="-wa -s -u"

    # INPUT
    peakachu="{}/peakachu_peaks.gtf".format(outputpath)
    pureclip="{}/pureclip_peaks.bed".format(outputpath)
    piranha="{}/piranha_peaks.bed".format(outputpath)

    ## OUTPUT
    peakachu_pureclip = "{}/peakachu_pureclip_peaks.bed".format(outputpath)
    peakachu_piranha = "{}/peakachu_piranha_peaks.bed".format(outputpath)
    piranha_pureclip = "{}/piranha_pureclip_peaks.bed".format(outputpath)
    all = "{}/peakachu_piranha_pureclip_peaks.bed".format(outputpath)

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

    stats_file = open("{}/stats_comparison.txt".format(outputpath), "w")
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
    venn3(subsets=(np_peakachu, np_piranha, np_peakachu_piranha,
                   np_pureclip, np_peakachu_pureclip, np_piranha_pureclip, np_all),
                   set_labels=('PEAKachu', 'Piranha', 'PureCLIP'))
    plt.show()
