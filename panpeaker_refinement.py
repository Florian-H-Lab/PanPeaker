import subprocess as sb
import matplotlib.pyplot as plt
import os

from simple_venn import venn4

def get_line_count(file):
    f = open(file, "r")
    count = 0
    for line in f:
        count += 1
    f.close()
    return count

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.mkdir(file_path)

def refinement(outputpath, adj_pval):

    # Refine Panpeaker peaks
    input_file = outputpath + "/robust_peaks.bed"
    sort_file = outputpath + "/robust_peaks_sorted.bed"
    filtered_file = outputpath + "/robust_peaks_sorted_filtered.bed"
    refined_file_tmp = outputpath + "/robust_peaks_refined_tmp.bed"
    refined_file = outputpath + "/robust_peaks_refined.bed"

    sb.Popen("bedtools sort -i {} > {}".format(input_file, sort_file), shell=True).wait()
    sb.Popen("awk '$9<=" + str(adj_pval) + "{print $0}' " + sort_file + " > " + filtered_file, shell=True).wait()
    sb.Popen("bedtools merge -i {} -s -c 6,9 -o distinct,mean > {}".format(filtered_file, refined_file_tmp), shell=True).wait()
    sb.Popen("awk -v OFS='\t' '{print $1,$2,$3,\"peak_\"NR,$5,$4}' " + refined_file_tmp + " > " + refined_file, shell=True).wait()

    sb.Popen("rm {}".format(sort_file), shell=True).wait()
    sb.Popen("rm {}".format(filtered_file), shell=True).wait()
    sb.Popen("rm {}".format(refined_file_tmp), shell=True).wait()

    print("[NOTE] Run Comparison with PanPeaker")

    options="-s -u"

    # INPUT
    n_a = "{}/peakachu_peaks.bed".format(outputpath)
    n_b = "{}/piranha_peaks.bed".format(outputpath)
    n_c = "{}/pureclip_peaks.bed".format(outputpath)
    n_d = refined_file

    p = dict()
    p["a"] = "PEAKachu"
    p["b"] = "Piranha"
    p["c"] = "PureCLIP"
    p["d"] = "PanPeaker"

    ## OUTPUT
    # for parallelisation
    # outputpath_comparison = outputpath + "/peak_comparison_{i}"
    outputpath_comparison = outputpath + "/peak_comparison_with_panpeaker"
    ensure_dir(outputpath_comparison)

    n_abcd = "{}/{}_{}_{}_{}_peaks.bed".format(outputpath_comparison, p["a"], p["b"], p["c"], p["d"])
    n_bacd = "{}/{}_{}_{}_{}_peaks.bed".format(outputpath_comparison, p["b"], p["a"], p["c"], p["d"])
    n_cabd = "{}/{}_{}_{}_{}_peaks.bed".format(outputpath_comparison, p["c"], p["a"], p["b"], p["d"])
    n_dabc = "{}/{}_{}_{}_{}_peaks.bed".format(outputpath_comparison, p["d"], p["a"], p["b"], p["c"])

    n_abc = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["a"], p["b"], p["c"])
    n_abd = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["a"], p["b"], p["d"])
    n_acd = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["a"], p["c"], p["d"])

    n_bac = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["b"], p["a"], p["c"])
    n_bad = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["b"], p["a"], p["d"])
    n_bcd = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["b"], p["c"], p["d"])

    n_cab = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["c"], p["a"], p["c"])
    n_cad = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["c"], p["a"], p["d"])
    n_cbd = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["c"], p["b"], p["d"])

    n_dab = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["d"], p["a"], p["b"])
    n_dac = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["d"], p["a"], p["c"])
    n_dbc = "{}/{}_{}_{}_peaks.bed".format(outputpath_comparison, p["d"], p["b"], p["c"])

    n_ab = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["a"], p["b"])
    n_ac = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["a"], p["c"])
    n_ad = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["a"], p["d"])

    n_ba = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["b"], p["a"])
    n_bc = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["b"], p["c"])
    n_bd = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["b"], p["d"])

    n_ca = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["c"], p["a"])
    n_cb = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["c"], p["b"])
    n_cd = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["c"], p["d"])

    n_da = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["d"], p["a"])
    n_db = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["d"], p["b"])
    n_dc = "{}/{}_{}_peaks.bed".format(outputpath_comparison, p["d"], p["c"])

    ## RUN
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_a, n_b, options, n_ab), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_a, n_c, options, n_ac), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_a, n_d, options, n_ad), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_b, n_a, options, n_ba), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_b, n_c, options, n_bc), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_b, n_d, options, n_bd), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_c, n_a, options, n_ca), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_c, n_b, options, n_cb), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_c, n_d, options, n_cd), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_d, n_a, options, n_da), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_d, n_b, options, n_db), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_d, n_c, options, n_dc), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ab, n_c, options, n_abc), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ab, n_d, options, n_abd), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ac, n_d, options, n_acd), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ba, n_c, options, n_bac), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ba, n_d, options, n_bad), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_bc, n_d, options, n_bcd), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ca, n_b, options, n_cab), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_ca, n_d, options, n_cad), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_cb, n_d, options, n_cbd), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_da, n_b, options, n_dab), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_da, n_c, options, n_dac), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_db, n_c, options, n_dbc), shell=True).wait()

    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_abc, n_d, options, n_abcd), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_bac, n_d, options, n_bacd), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_cab, n_d, options, n_cabd), shell=True).wait()
    sb.Popen("bedtools intersect -a {} -b {} {} > {}".format(n_dab, n_c, options, n_dabc), shell=True).wait()

    ## STATS
    abcd = get_line_count(n_abcd)
    bacd = get_line_count(n_bacd)
    cabd = get_line_count(n_cabd)
    dabc = get_line_count(n_dabc)

    abc = get_line_count(n_abc)
    abd = get_line_count(n_abd)
    acd = get_line_count(n_acd)

    bac = get_line_count(n_bac)
    bad = get_line_count(n_bad)
    bcd = get_line_count(n_bcd)

    cab = get_line_count(n_cab)
    cad = get_line_count(n_cad)
    cbd = get_line_count(n_cbd)

    dab = get_line_count(n_dab)
    dac = get_line_count(n_dac)
    dbc = get_line_count(n_dbc)

    ab = get_line_count(n_ab)
    ac = get_line_count(n_ac)
    ad = get_line_count(n_ad)

    ba = get_line_count(n_ba)
    bc = get_line_count(n_bc)
    bd = get_line_count(n_bd)

    ca = get_line_count(n_ca)
    cb = get_line_count(n_cb)
    cd = get_line_count(n_cd)

    da = get_line_count(n_da)
    db = get_line_count(n_db)
    dc = get_line_count(n_dc)

    a = get_line_count(n_a)
    b = get_line_count(n_b)
    c = get_line_count(n_c)
    d = get_line_count(n_d)

    stats_file = open("{}/stats_comparison.txt".format(outputpath_comparison), "w")
    stats_file.write("peakcaller\tpeakcount\n")
    stats_file.write("{}\t{}\n".format(p["a"], a))
    stats_file.write("{}\t{}\n".format(p["b"], b))
    stats_file.write("{}\t{}\n".format(p["c"], c))
    stats_file.write("{}\t{}\n".format(p["d"], d))

    stats_file.write("{}_{}\t{}\n".format(p["a"], p["b"], ab))
    stats_file.write("{}_{}\t{}\n".format(p["a"], p["c"], ac))
    stats_file.write("{}_{}\t{}\n".format(p["a"], p["d"], ad))

    stats_file.write("{}_{}\t{}\n".format(p["b"], p["a"], ba))
    stats_file.write("{}_{}\t{}\n".format(p["b"], p["c"], bc))
    stats_file.write("{}_{}\t{}\n".format(p["b"], p["d"], bd))

    stats_file.write("{}_{}\t{}\n".format(p["c"], p["a"], ca))
    stats_file.write("{}_{}\t{}\n".format(p["c"], p["b"], cb))
    stats_file.write("{}_{}\t{}\n".format(p["c"], p["d"], cd))

    stats_file.write("{}_{}\t{}\n".format(p["d"], p["a"], da))
    stats_file.write("{}_{}\t{}\n".format(p["d"], p["b"], db))
    stats_file.write("{}_{}\t{}\n".format(p["d"], p["c"], dc))

    stats_file.write("{}_{}_{}\t{}\n".format(p["a"], p["b"], p["c"], abc))
    stats_file.write("{}_{}_{}\t{}\n".format(p["a"], p["b"], p["d"], abd))
    stats_file.write("{}_{}_{}\t{}\n".format(p["a"], p["c"], p["d"], acd))

    stats_file.write("{}_{}_{}\t{}\n".format(p["b"], p["a"], p["c"], bac))
    stats_file.write("{}_{}_{}\t{}\n".format(p["b"], p["a"], p["d"], bad))
    stats_file.write("{}_{}_{}\t{}\n".format(p["b"], p["c"], p["d"], bcd))

    stats_file.write("{}_{}_{}\t{}\n".format(p["c"], p["a"], p["b"], cab))
    stats_file.write("{}_{}_{}\t{}\n".format(p["c"], p["a"], p["d"], cad))
    stats_file.write("{}_{}_{}\t{}\n".format(p["c"], p["b"], p["d"], cbd))

    stats_file.write("{}_{}_{}\t{}\n".format(p["d"], p["a"], p["c"], dab))
    stats_file.write("{}_{}_{}\t{}\n".format(p["d"], p["a"], p["c"], dac))
    stats_file.write("{}_{}_{}\t{}\n".format(p["d"], p["b"], p["c"], dbc))

    stats_file.write("ALL_{}\t{}\n".format(p["a"], abcd))
    stats_file.write("ALL_{}\t{}\n".format(p["b"], bacd))
    stats_file.write("ALL_{}\t{}\n".format(p["c"], cabd))
    stats_file.write("ALL_{}\t{}\n".format(p["d"], dabc))
    stats_file.close()

    # Adjust values for venn diagram.
    abc -= abcd
    abd -= abcd
    acd -= abcd

    bac -= bacd
    bad -= bacd
    bcd -= bacd

    cab -= cabd
    cad -= cabd
    cbd -= cabd

    dab -= dabc
    dac -= dabc
    dbc -= dabc

    ab -= (abcd + abc + abd)
    ac -= (abcd + abc + acd)
    ad -= (abcd + abd + acd)

    ba -= (bacd + bac + bad)
    bc -= (bacd + bac + bcd)
    bd -= (bacd + bcd + bad)

    ca -= (cabd + cab + cad)
    cb -= (cabd + cab + cbd)
    cd -= (cabd + cbd + cad)

    da -= (dabc + dab + dac)
    db -= (dabc + dab + dbc)
    dc -= (dabc + dbc + dac)

    a -= (abcd + abc + abd + acd + ab + ac + ad)
    b -= (bacd + bac + bad + bcd + ba + bc + bd)
    c -= (cabd + cab + cad + cbd + ca + cb + cd)
    d -= (dabc + dab + dac + dbc + da + db + dc)

    a = "A:" + str(a)
    b = "B:" + str(b)
    c = "C:" + str(c)
    d = "D:" + str(d)

    ab = "A:{}\nB:{}".format(str(ab), str(ba))
    ac = "A:{}\nC:{}".format(str(ac), str(ca))
    ad = "A:{}\nD:{}".format(str(ad), str(da))

    bc = "B:{}\nC:{}".format(str(bc), str(cb))
    bd = "B:{}\nD:{}".format(str(bd), str(db))

    cd = "C:{}\nD:{}".format(str(cd), str(dc))

    abc = "A:{}\nB:{}\nC:{}".format(str(abc), str(bac), str(cab))
    abd = "A:{}\nB:{}\nD:{}".format(str(abd), str(bad), str(dab))
    acd = "A:{}\nC:{}\nD:{}".format(str(acd), str(cad), str(dac))
    bcd = "B:{}\nC:{}\nD:{}".format(str(bcd), str(cbd), str(dbc))

    abcd = "A:{}\nB:{}\nC:{}\nD:{}".format(str(abcd), str(bacd), str(cabd), str(dabc))

    data = [a, b, c, d, ab, ac, ad, bc, bd, cd, abc, abd, acd, bcd, abcd]

    ## VENN DIAGRAM
    f = plt.figure(figsize=(10, 10))
    plt.rcParams["font.family"] = "serif"
    venn4(data, set_labels=('PEAKachu', 'Piranha', 'PureCLIP', 'PanPeaker'), set_colors=('r', 'y', 'b', 'g'))
    plt.show()
    f.savefig(outputpath_comparison + "/comparison_plot.pdf", bbox_inches='tight', dpi=300)