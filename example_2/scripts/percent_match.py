#!/usr/bin/env python

import sys
import subprocess
import numpy as np

haplotypes = open("../haplotypes.polarized.vcf").readlines()
if sys.argv[-2] == '-h':
    haplotypes = open(sys.argv[-1]).readlines()
    sys.argv = sys.argv[:len(sys.argv)-2]


num_to_report = int(sys.argv[1])

header_haplotypes = haplotypes[0].split()

haps = haplotypes[0].split()[9:]
hap_inds = range(0,len(haps))

for file in sys.argv[2:]:
    observations = open(file).readlines()
    simple_similarity_scores = [0.] * len(haps)

    contig_similarity_scores = [0.] * len(haps)
    prev_mod = [0.] * len(haps)

    for i in range(1, min(len(haplotypes), len(observations))):
        observation = observations[i].split()[1]
        haplotype = haplotypes[i].split()
        for hap in hap_inds:
            if observation == haplotype[hap+9] and observation != ".":
                simple_similarity_scores[hap] += 1
                contig_similarity_scores[hap] += 1 + prev_mod[hap]
                prev_mod[hap] += 1
            else:
                prev_mod[hap] = 0

    # print(haps)
    # print("Naive:")
    # print([score/(len(haplotypes)-1) for score in simple_similarity_scores])
    # print("Contig:")
    # print([score/(len(haplotypes)-1) for score in contig_similarity_scores])
    naive_vals = np.sort(simple_similarity_scores)[::-1][0:num_to_report]
    contig_vals = np.sort(contig_similarity_scores)[::-1][0:num_to_report]

    naive = [haps[i] for i in np.argsort(simple_similarity_scores)[::-1][0:num_to_report]]
    naive_vals = np.sort(simple_similarity_scores)[::-1][0:num_to_report]
    # naive = [header_haplotypes.index(i) + 1 for i in naive_names]
    # for i in range(0,len(naive)):
    #     print(naive[i], ": ", naive_vals[i])
    contig = [haps[i] for i in np.argsort(contig_similarity_scores)[::-1][0:num_to_report]]
    contig_vals = np.sort(contig_similarity_scores)[::-1][0:num_to_report]
    # contig = [header_haplotypes.index(i) + 1 for i in contig_names]
    # for i in range(0,len(contig)):
    #     print(contig[i], ": ", contig_vals[i])

    print("Naive is Contig? ", naive == contig)
    naive_file  = "_".join(file.split(".")[0:2]) + "_top_" + str(sys.argv[1]) + "_naive_founder_matches.dat"
    contig_file = "_".join(file.split(".")[0:2]) + "_top_" + str(sys.argv[1]) + "_contig_founder_matches.dat"
    naive_write = open(naive_file, "w")
    contig_write = open(contig_file, "w")

    header = "Chrom\tPos\tRef\tAlt"
    to_cut = "1,2,4,5"
    for i in naive:
        header += "\t"+i
        # print("To add to cut" ,i)
        # print("header_haplotypes: ", header_haplotypes)
        j = header_haplotypes.index(i) + 1
        # print("This is j", j)
        to_cut += ","+str(j)
    header += "\n"
    naive_write.write(header)
    print(header)
    naive_write.close()
    to_run = "grep "+ file.split(".")[1] + " ../haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + naive_file
    print("to_run:", to_run)
    # print("File name:",naive_file)
    # # print("Chrom: ", cur_chrom)
    # # print("Line: ",line)
    # # to_run = to_run.split()
    # result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")

    header = "Chrom\tPos\tRef\tAlt"
    to_cut = "1,2,4,5"
    for i in contig:
        header += "\t"+i
        j = header_haplotypes.index(i) + 1
        to_cut += ","+str(j)
    header += "\n"
    contig_write.write(header)
    print(header)
    contig_write.close()
    to_run = "grep "+ file.split(".")[1]+ " ../haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + contig_file
    print("to_run:", to_run)
    # print("File name:", contig_file)
    # # print("Chrom: ", cur_chrom)
    # # print("Line: ",line)
    # # to_run = to_run.split()
    # result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
