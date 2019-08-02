#!/usr/bin/env python

import os
import sys
import subprocess

for file in sys.argv:
    if file.endswith("mlp"):
        with open("haplotypes.polarized.vcf") as f:
            header_haplotypes = f.readline().split()
        with open(file) as f:
            header = f.readline()
            ind_num = file[0]
            file_name = ind_num+"_2L_true_founders.dat"
            w  = open(file_name, "w")
            cur_chrom = "2L"
            founder_list = []
            for line in f.readlines():
                line = line.split()
                if line[0] != cur_chrom:
                    header = "Chrom\tPos\tRef\tAlt"
                    to_cut = "1,2,4,5"
                    for i in founder_list:
                        header += "\t"+i
                        j = header_haplotypes.index(i) + 1
                        to_cut += ","+str(j)
                    founder_list = []
                    header += "\n"
                    w.write(header)
                    to_run = "grep "+ cur_chrom+ " haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + file_name
                    print("to_run:", to_run)
                    print("File name:",file)
                    print("Chrom: ", cur_chrom)
                    print("Line: ",line)
                    # to_run = to_run.split()
                    result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
                    subprocess.run("sed -i 's/\t0/\t./g' "+file_name, shell=True)
                    subprocess.run("sed -i 's/\t1/\t0/g' "+file_name, shell=True)
                    subprocess.run("sed -i 's/\t-1/\t1/g' "+file_name, shell=True)
                    cur_chrom = line[0]
                    file_name = ind_num+ "_"+line[0] + "_true_founders.dat"
                    w.close()
                    w = open(file_name, "w")
                founder_list.append(line[1])
            header = "Chrom\tPos\tRef\tAlt"
            to_cut = "1,2,4,5"
            for i in founder_list:
                header += "\t"+i
                j = header_haplotypes.index(i) + 1
                to_cut += ","+str(j)
            header += "\n"
            w.write(header)
            to_run = "grep "+ cur_chrom+ " haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + file_name
            print("to_run:", to_run)
            print("File name:",file)
            print("Chrom: ", cur_chrom)
            print("Line: ",line)
            # to_run = to_run.split()
            result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
            cur_chrom = line[0]
            file_name = ind_num+ "_"+line[0] + "_true_founders.dat"
            w.close()

        w.close()
    elif file.endswith("counts"):
        with open(file) as f:
            header = f.readline().split()
            cur_chrom = "2L"
            file_name = file[0] + "_observations_"+cur_chrom+".dat"
            w = open(file_name, "w")
            w.write("Pos\tobs\n")
            for i in f.readlines():
                if i[0] == "@" or i[0] == "C":
                    continue
                else:
                    line = i.split()
                    if line[0] != cur_chrom:
                        w.close()
                        cur_chrom = line[0]
                        file_name = file[0] + "_observations_"+cur_chrom+".dat"
                        w = open(file_name, "w")
                        w.write("Pos\tobs\n")
                    obs = "."
                    if line[2] != 0 and line[3] == 0:
                        obs = "0"
                    elif line[3] != 0 and line[2] == 0:
                        obs = "1"
                    w.write(line[1] + '\t' + obs + '\n')

            w.close()
