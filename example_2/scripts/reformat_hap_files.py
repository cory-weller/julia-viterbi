#!/usr/bin/env python

import os
import sys
import subprocess
import random

for index in range(1, len(sys.argv)):
    file = sys.argv[index]

    if file.endswith("founder_matches.dat") or file.endswith("founders_wide.dat"): #Optimize when you get a minute
        print("1: ", sys.argv[index])
        print("2: ", sys.argv[index+1])
        keys = [line.split()[0].strip() for line in open(sys.argv[index+1]).readlines()[1:]]
        print("keys: ", keys[0:10])
        key_ind = 0
        with open(file+".stripped", "w") as dest:
            lines = open(file).readlines()
            dest.write(lines[0])
            for line in lines[1:]:
                # print(line)
                # print(keys[key_ind])
                if line.split()[1].strip() == keys[key_ind]:
                    dest.write(line)
                    key_ind += 1

    if file.endswith("key"):
        continue

    if file.endswith("mlp"): #Only get the foudners which are relevant
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
                    w.close()
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
                    w = open(file_name, "w")
                if line[3] not in founder_list:
                    founder_list.append(line[3])
                if line[4] not in founder_list:
                    founder_list.append(line[4])
            header = "Chrom\tPos\tRef\tAlt"
            to_cut = "1,2,4,5"
            for i in founder_list:
                header += "\t"+i
                j = header_haplotypes.index(i) + 1
                to_cut += ","+str(j)
            header += "\n"
            w.write(header)
            w.close()
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

    # if file.endswith("mlp"): #Only get the foudners which are relevant
    #     with open("haplotypes.polarized.vcf") as f:
    #         header_haplotypes = f.readline().split()
    #     with open(file) as f:
    #         header = f.readline()
    #         ind_num = file[0]
    #         file_name = ind_num+"_2L_true_founders.dat"
    #         w  = open(file_name, "w")
    #         cur_chrom = "2L"
    #         founder_list = []
    #         for line in f.readlines():
    #             line = line.split()
    #             if line[0] != cur_chrom:
    #                 header = "Chrom\tPos\tRef\tAlt"
    #                 to_cut = "1,2,4,5"
    #                 for i in founder_list:
    #                     header += "\t"+i
    #                     j = header_haplotypes.index(i) + 1
    #                     to_cut += ","+str(j)
    #                 founder_list = []
    #                 header += "\n"
    #                 w.write(header)
    #                 w.close()
    #                 to_run = "grep "+ cur_chrom+ " haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + file_name
    #                 print("to_run:", to_run)
    #                 print("File name:",file)
    #                 print("Chrom: ", cur_chrom)
    #                 print("Line: ",line)
    #                 # to_run = to_run.split()
    #                 result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
    #                 subprocess.run("sed -i 's/\t0/\t./g' "+file_name, shell=True)
    #                 subprocess.run("sed -i 's/\t1/\t0/g' "+file_name, shell=True)
    #                 subprocess.run("sed -i 's/\t-1/\t1/g' "+file_name, shell=True)
    #                 cur_chrom = line[0]
    #                 file_name = ind_num+ "_"+line[0] + "_true_founders.dat"
    #                 w = open(file_name, "w")
    #             founder_list.append(line[1])
    #         header = "Chrom\tPos\tRef\tAlt"
    #         to_cut = "1,2,4,5"
    #         for i in founder_list:
    #             header += "\t"+i
    #             j = header_haplotypes.index(i) + 1
    #             to_cut += ","+str(j)
    #         header += "\n"
    #         w.write(header)
    #         w.close()
    #         to_run = "grep "+ cur_chrom+ " haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + file_name
    #         print("to_run:", to_run)
    #         print("File name:",file)
    #         print("Chrom: ", cur_chrom)
    #         print("Line: ",line)
    #         # to_run = to_run.split()
    #         result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
    #         cur_chrom = line[0]
    #         file_name = ind_num+ "_"+line[0] + "_true_founders.dat"
    #
    #     w.close()

    if file.endswith("mlp"): #make an equal number of random potential founders - ANTAGONIZE
        with open("haplotypes.polarized.vcf") as f:
            header_haplotypes = f.readline().split()
        with open(file) as f:
            header = f.readline()
            ind_num = file[0]
            file_name = ind_num+"_2L_true_founders_with_random.dat"
            w  = open(file_name, "w")
            cur_chrom = "2L"
            founder_list = []
            for line in f.readlines():
                line = line.split()
                if line[0] != cur_chrom:
                    header = "Chrom\tPos\tRef\tAlt"
                    to_cut = "1,2,4,5"
                    num_real_founders = len(founder_list)
                    print("Real founders: ", num_real_founders)
                    num_rand_added = 0
                    while num_rand_added < num_real_founders:
                        # print("Sampling...")
                        sample = random.sample(header_haplotypes,1)[0]
                        # print("Sample: ", sample)
                        if sample[0:3] == "RAL":
                            if sample not in founder_list:
                                # print("Adding random...")
                                founder_list.append(sample)
                                num_rand_added += 1

                    for i in founder_list:
                        header += "\t"+i
                        j = header_haplotypes.index(i) + 1
                        to_cut += ","+str(j)
                    founder_list = []
                    header += "\n"
                    w.write(header)
                    w.close()
                    to_run = "grep "+ cur_chrom+ " haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + file_name
                    print("to_run:", to_run)
                    print("File name:",file)
                    print("Chrom: ", cur_chrom)
                    print("Line: ",line)
                    # to_run = to_run.split()
                    result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
                    subprocess.run(str("sed -i 's/\t0/\t./g' "+file_name), shell=True)
                    subprocess.run("sed -i 's/\t1/\t0/g' "+file_name, shell=True)
                    subprocess.run("sed -i 's/\t-1/\t1/g' "+file_name, shell=True)
                    cur_chrom = line[0]
                    file_name = ind_num+ "_"+line[0] + "_true_founders_with_random.dat"
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
            w.close()
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

    if file.endswith("wide.haps"): #Currently only configured for 2L
        chrom_list = ["2L"]
        with open("../haplotypes.polarized.vcf") as f:
            header_haplotypes = f.readline().split()
        with open(file) as f:
            header = f.readline()
            ind_num = file.split(".")[0]
            file_name = ind_num+"_2L_true_founders_wide.dat"
            w  = open(file_name, "w")
            cur_chrom = "2L"
            founder_list = []
            for line in f.readlines():
                line = line.split()
                if line[0] != cur_chrom:
                    if (line[0] not in chrom_list):
                        break
                    header = "Chrom\tPos\tRef\tAlt"
                    to_cut = "1,2,4,5"
                    num_real_founders = len(founder_list)
                    # print("Real founders: ", num_real_founders)
                    num_rand_added = 0
                    while num_rand_added < num_real_founders:
                        # print("Sampling...")
                        sample = random.sample(header_haplotypes,1)[0]
                        # print("Sample: ", sample)
                        if sample[0:3] == "RAL":
                            if sample not in founder_list:
                                # print("Adding random...")
                                founder_list.append(sample)
                                num_rand_added += 1

                    for i in founder_list:
                        header += "\t"+i
                        j = header_haplotypes.index(i) + 1
                        to_cut += ","+str(j)
                    founder_list = []
                    header += "\n"
                    w.write(header)
                    w.close()
                    to_run = "grep "+ cur_chrom+ " ../haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + file_name
                    print("to_run:", to_run)
                    print("File name:",file)
                    # print("Chrom: ", cur_chrom)
                    # print("Line: ",line)
                    # to_run = to_run.split()
                    result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
                    subprocess.run(str("sed -i 's/\t0/\t./g' "+file_name), shell=True)
                    subprocess.run("sed -i 's/\t1/\t0/g' "+file_name, shell=True)
                    subprocess.run("sed -i 's/\t-1/\t1/g' "+file_name, shell=True)
                    cur_chrom = line[0]
                    file_name = ind_num+ "_"+line[0] + "_true_founders_wide.dat"
                    w = open(file_name, "w")
                if (line[3] not in founder_list):
                    founder_list.append(line[3])
                if (line[4] not in founder_list):
                    founder_list.append(line[4])
            header = "Chrom\tPos\tRef\tAlt"
            to_cut = "1,2,4,5"
            for i in founder_list:
                header += "\t"+i
                j = header_haplotypes.index(i) + 1
                to_cut += ","+str(j)
            header += "\n"
            w.write(header)
            w.close()
            to_run = "grep "+ cur_chrom+ " ../haplotypes.polarized.vcf | cut -f "+ to_cut + " >> " + file_name
            print("to_run:", to_run)
            print("File name:",file)
            # print("Chrom: ", cur_chrom)
            # print("Line: ",line)
            # to_run = to_run.split()
            result = subprocess.run(to_run, stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8")
            cur_chrom = line[0]
            file_name = ind_num+ "_"+line[0] + "_true_founders_wide.dat"

        w.close()

    elif file.endswith("counts"):
        #5892 for four
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
                    if int(line[2]) != 0 and int(line[3]) == 0:
                        obs = "0"
                    elif int(line[3]) != 0 and int(line[2]) == 0:
                        obs = "1"
                    w.write(line[1] + '\t' + obs + '\n')

            w.close()

    elif file.endswith("vcf"):
        print("File: ", file + "  ")#, end='\r')
        #5892 for four
        with open(file) as f:
            header = f.readline().split()
            cur_chrom = "2L"
            file_name = file.split(".")[0] + "."+cur_chrom+".dat"
            print("File name: ", file_name)
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
                        file_name = file[0] + "."+cur_chrom+".dat"
                        w = open(file_name, "w")
                        w.write("Pos\tobs\n")
                    obs_1 = "."
                    obs_2 = "."
                    if line[2] == "1/1":
                        obs_1 = "1"
                    elif line[2] == "0/0":
                        obs_1 = "0"
                    if line[3] == "1/1":
                        obs_2 = "1"
                    elif line[3] == "0/0":
                        obs_2 = "0"
                    w.write(line[1] + '\t' + obs_1 + '\t' + obs_2 + '\n')

            w.close()
