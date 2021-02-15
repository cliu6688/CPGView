import re
import os
import csv
import sys
from Bio import SeqIO

# Part one: import data
arch    = sys.argv[1]
file1   = sys.argv[2]
file2   = sys.argv[3]
record  = SeqIO.parse(arch, "genbank")
rec     = next(record)


def get_path(path1):
    r = os.path.abspath(path1)
    return r


new_list1 = []
new_list2 = []
new_list3 = []
new_list4 = []
new_list5 = []
new_list6 = []

# Part two: select gene and subgene
y_nums = 0
for f in rec.features:
    if f.type != 'source' and f.type != 'gene':
        if "join" in str(f.location) and f.qualifiers["gene"][0].lower() != 'rps12':
            y_nums += 1
            location1 = str(f.location).replace("join{", "").replace("}", "").replace("](-)", "").replace("](+)",
                                                                                                          "").replace(
                "[", "").replace(" ", "").split(",")  # delete the direction
            sign = list(set(re.findall(r'[(](.*?)[)]', str(f.location))))[0]

            if sign == "-":
                direction = -1
                strand = "reverse"
            elif sign == "+":
                direction == 1
                strand = "forward"

            if "(+)" in str(f.location):

                start1 = int(location1[0].split(":")[0]) + 1
                end1 = int(location1[-1].split(":")[1])
                line = [str(y_nums), int(start1), int(end1), strand, direction, f.qualifiers["gene"][0]]
                new_list1.append(line)

                for jj in location1:
                    substart1 = int(jj.split(":")[0]) + 1
                    subend1 = jj.split(":")[1]
                    line3 = [str(y_nums), int(start1), int(end1), strand, direction, "Exon", int(substart1),
                             int(subend1), f.qualifiers["gene"][0]]
                    new_list2.append(line3)

            elif "(-)" in str(f.location):

                start2 = int(location1[-1].split(":")[0]) + 1
                end2 = int(location1[0].split(":")[1])
                line1 = [str(y_nums), int(start2), int(end2), strand, direction, f.qualifiers["gene"][0]]
                new_list1.append(line1)

                for jj in location1:
                    substart2 = int(jj.split(":")[0]) + 1
                    subend2 = jj.split(":")[1]

                    line2 = [str(y_nums), int(start2), int(end2), strand, direction, "Exon", int(substart2),
                             int(subend2), f.qualifiers["gene"][0]]
                    new_list2.append(line2)

        elif "join" in str(f.location) and f.qualifiers["gene"][0].lower() == 'rps12':

            location2 = str(f.location).replace("join{", "").replace("}", "").replace("]", "").replace("]", "").replace(
                "[", "").replace(" ", "").split(",")  # Remove the position of the direction
            sign2 = list(set(re.findall(r'[(](.*?)[)]', str(f.location))))
            unms = 1

            new_location = []
            for mie in location2:
                mie1 = mie + str(unms)
                new_location.append(mie1)
                unms += 1

            str_new_location = " ".join(new_location)

            if "+" in str_new_location and "-" in str_new_location:
                new_location1 = str_new_location.split(" ")
                for mxs in new_location1:
                    if "+" in mxs:
                        direction1 = 1
                        strand1 = "forward"
                        mx1 = mxs.split("(")

                        label_start1 = int(mx1[0].split(":")[0]) + 1
                        label_end1 = int(mx1[0].split(":")[1])
                        num1 = mx1[1].split(")")[1]

                        line3 = ["2 rps12", "exon" + str(num1), strand1, direction1, int(label_start1), int(label_end1)]
                        new_list3.append(line3)

                        new_line3 = ["1 Transcript 1", "exon" + str(num1), "Transcript 1", direction1,
                                     int(label_start1),
                                     int(label_end1), "exon" + str(num1)]

                        if new_line3[1] == "exon2" or new_line3[1] == "exon3":
                            new_line3[1] = new_line3[1] + " (IRa)"
                        new_list4.append(new_line3)

                    elif "-" in mxs:
                        direction2 = -1
                        strand2 = "reverse"
                        mx2 = mxs.split("(")
                        label_start2 = int(mx2[0].split(":")[0]) + 1
                        label_end2 = int(mx2[0].split(":")[1])
                        num2 = mx2[1].split(")")[1]
                        line4 = ["2 rps12", "exon" + str(num2), strand2, direction2, int(label_start2), int(label_end2)]
                        new_list3.append(line4)

                        new_line4 = ["1 Transcript 1", "exon" + str(num2), "Transcript 1", 1, int(label_start2),
                                     int(label_end2), "exon" + str(num2)]
                        new_list4.append(new_line4)

            elif "-" in str_new_location and "+" not in str_new_location:
                new_location2 = str_new_location.split(" ")
                for mxs1 in new_location2:
                    direction3 = -1
                    strand3 = "reverse"
                    mx3 = mxs1.split("(")
                    label_start3 = int(mx3[0].split(":")[0]) + 1
                    label_end3 = int(mx3[0].split(":")[1])
                    num3 = mx3[1].split(")")[1]
                    line5 = ["2 rps12", "exon" + str(num3), strand3, direction3, int(label_start3), int(label_end3)]
                    new_list3.append(line5)

                    new_line5 = ["3 Transcript 2", "exon" + str(num3), "Transcript 2", 1, int(label_start3),
                                 int(label_end3), "exon" + str(num3)]
                    if new_line5[1] == "exon2" or new_line5[1] == "exon3":
                        new_line5[1] = new_line5[1] + " (IRb)"
                    new_list5.append(new_line5)

            elif "+" in str_new_location and "-" not in str_new_location:
                new_location3 = str_new_location.split(" ")
                for mxs2 in new_location3:
                    direction4 = 1
                    strand4 = "forward"
                    mx4 = mxs2.split("(")
                    label_start4 = int(mx4[0].split(":")[0]) + 1
                    label_end4 = int(mx4[0].split(":")[1])
                    num4 = mx4[1].split(")")[1]
                    line6 = ["2 rps12", "exon" + str(num4), strand4, direction4, int(label_start4), int(label_end4)]
                    new_list3.append(line6)
                    new_line6 = ["3 Transcript 2", "exon" + str(num4), "Transcript 2", 1, int(label_start4),
                                 int(label_end4), "exon" + str(num4)]
                    if new_line6[1] == "exon2" or new_line6[1] == "exon3":
                        new_line6[1] = new_line6[1] + " (IRa)"
                    new_list6.append(new_line6)

new_list5 = (new_list6 if new_list6 != [] else new_list5)

# Part three: delete  genes and errors in genes
new_list1 = sorted(new_list1, key=lambda x: x[1])  # Sort by location of gene
new_list11 = new_list1.copy()
for i in range(0, len(new_list1)):
    for ii in range(0, len(new_list1)):
        if new_list1[ii][1] < new_list1[i][1] < new_list1[ii][2]:
            if new_list1[i] in new_list11 and new_list1[ii] in new_list11:
                new_list11.remove(new_list1[i])
                new_list11.remove(new_list1[ii])
            else:
                continue

new_list2 = sorted(new_list2, key=lambda x: x[6])
new_list22 = new_list2.copy()

for m in range(0, len(new_list2)):
    for mm in range(0, len(new_list2)):
        if new_list2[mm][1] < new_list2[m][1] < new_list2[mm][2]:
            if new_list2[m] in new_list22 and new_list2[mm] in new_list22:
                new_list22.remove(new_list2[m])
                new_list22.remove(new_list2[mm])
            else:
                continue

# Part four: Dealing with the problem of gene exon coordinate overlap
new_list23 = []
for aq in new_list22:
    total_distance = (aq[2] - aq[1]) / 10
    sub_gene_distance = aq[-2] - aq[-3]
    cha = sub_gene_distance - total_distance

    if cha >= 0:
        qs1 = [aq[0], aq[1], aq[2], aq[3], aq[4], aq[5], aq[6], aq[7], aq[6], aq[7], aq[-1]]
        new_list23.append(qs1)

    if cha < 0:
        cha1 = total_distance - sub_gene_distance
        if aq[1] == aq[6] and aq[2] != aq[7]:
            new_end = aq[7] + cha1
            qs2 = [aq[0], aq[1], aq[2], aq[3], aq[4], aq[5], aq[6], new_end, aq[6], aq[7], aq[-1]]

            new_list23.append(qs2)

        elif aq[1] != aq[6] and aq[2] == aq[7]:
            new_start = aq[6] - cha1
            qs3 = [aq[0], aq[1], aq[2], aq[3], aq[4], aq[5], new_start, aq[7], aq[6], aq[7], aq[-1]]

            new_list23.append(qs3)

        elif aq[1] != aq[6] and aq[2] != aq[7]:
            new_start = aq[6] - (cha1 / 2)
            new_end = aq[7] + (cha1 / 2)
            qs4 = [aq[0], aq[1], aq[2], aq[3], aq[4], aq[5], new_start, new_end, aq[6], aq[7], aq[-1]]
            new_list23.append(qs4)

for ow in new_list23:
    ow[6] = int(ow[6])
    ow[7] = int(ow[7])

# Part 5: add intron to subgene
new_list23 = sorted(new_list23, key=lambda x: x[8])

test_value1 = new_list23[0]
test_value2 = new_list23[1]

if test_value1[0] == test_value2[0] and test_value1[1] == test_value2[1] and test_value1[2] == test_value2[2]:
    intron = [test_value1[0], test_value1[1], test_value1[2], test_value1[3], test_value1[4], "Intron", test_value1[7],
              test_value2[6], test_value1[9], test_value2[8], test_value1[-1]]

else:
    if test_value1[1] == test_value1[-2]:
        intron = [test_value1[0], test_value1[1], test_value1[2], test_value1[3], test_value1[4], "Intron",
                  test_value1[7], test_value1[2], test_value1[9], test_value1[2], test_value1[-1]]
    elif test_value1[2] == test_value1[-1]:
        intron = [test_value1[0], test_value1[1], test_value1[2], test_value1[3], test_value1[4], "Intron",
                  test_value1[1], test_value1[6], test_value1[1], test_value1[8], test_value1[-1]]

new_list23.insert(0, intron)
new_list23 = sorted(new_list23, key=lambda x: x[8])

# The sixth part: processing the data of RPS12
new_raw_list0 = []
new_raw_list1 = []
new_raw_list2 = []
new_list33 = new_list3.copy()
list1 = list(set(list(map(tuple, new_list33))))
list2 = list(map(list, list1))
list2 = sorted(list2, key=lambda x: x[-2])
list3 = []
list4 = []
for lw in list2:
    if lw[3] == 1:
        list3.append(lw)
    elif lw[3] == -1:
        list4.append(lw)

list3 = sorted(list3, key=lambda x: x[-2])  # Positive column
list4 = sorted(list4, key=lambda x: x[-2])  # Negative column

loca11 = [[1260, 1410], [1490, 1640]]
local2 = [[600, 750], [820, 970], [1040, 1190]]
local3 = [[100, 135], [255, 290], [310, 345]]
local4 = [[2100, 2135], [2210, 2245], [2310, 2345]]

for xq, yq in zip(list3, loca11):
    xq.insert(2, yq[0])
    xq.insert(3, yq[1])
    xq.insert(8, xq[1])

    xq[4] = "Genome  (+)"
    if xq[1] == "exon2":
        xq[1] = "exon2" + " (IRa)"
    elif xq[1] == "exon3":
        xq[1] = "exon3" + " (IRa)"
    new_raw_list1.append(xq)

for xq1, yq1 in zip(list4, local2):
    xq1.insert(2, yq1[0])
    xq1.insert(3, yq1[1])
    xq1.insert(8, xq1[1])
    xq1[4] = "(-)"
    if xq1[1] == "exon2":
        xq1[1] = "exon2" + " (IRb)"
    elif xq1[1] == "exon3":
        xq1[1] = "exon3" + " (IRb)"
    new_raw_list1.append(xq1)

for xq2, yq2 in zip(new_list4, local3):
    xq2.insert(2, yq2[0])
    xq2.insert(3, yq2[1])
    new_raw_list0.append(xq2)

for xq3, yq3 in zip(new_list5, local4):
    xq3.insert(2, yq3[0])
    xq3.insert(3, yq3[1])
    new_raw_list2.append(xq3)


# The seventh part: save data to CSV file

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
        print("")
        print("---  New folder (", str(path).strip(), ") has been generated  ---")
        return 1
    else:
        print("")
        print("---  New folder (", str(path).strip(), ") already exists here. Please delete it first!  ---")
        return -1

#file1 = "data_cis_splicing_gene"
folder1 = os.path.exists(file1)
#file2 = "data_trans_splicing_gene"
folder2 = os.path.exists(file2)

if not folder1 and not folder2:

    os.makedirs(file1)
    #os.makedirs(folder1)
    print("---  New folder (data_cis_splicing_gene) has been generated  ---")
    os.makedirs(file2)
    #os.makedirs(folder2)
    print("---  New folder (data_trans_splicing_gene) has been generated  ---")

    ab_path1 = str(get_path('./' + file1))
    ab_path2 = str(get_path('./' + file2))
    name1 = ab_path1 + "/cis_splicing_gene.csv"
    name2 = ab_path1 + "/cis_splicing_subgene.csv"
    name3 = ab_path2 + "/trans_splicing_gene.csv"
    name4 = ab_path2 + "/trans_splicing_subgene.csv"
    csvfile1 = open(name1, 'w', newline='', encoding="utf-8")
    writer1 = csv.writer(csvfile1)
    writer1.writerow(["Gene", "start", "end", "strand", "direction", "Gene_label"])
    writer1.writerows(new_list11)
    csvfile1.close()

    csvfile2 = open(name2, 'w', newline='', encoding="utf-8")
    writer2 = csv.writer(csvfile2)
    writer2.writerow(["Gene", "start", "end", "strand", "direction", "Subgene", "from", "to", "label_from", "label_to",
                      "Gene_label"])
    writer2.writerows(new_list23)
    csvfile2.close()

    csvfile3 = open(name3, 'w', newline='', encoding="utf-8")
    writer3 = csv.writer(csvfile3)
    writer3.writerow(["Gene", "Exon", "start", "end", "strand", "direction", "label_from", "label_to", "label_label"])
    writer3.writerows(new_raw_list0)
    writer3.writerows(new_raw_list1)
    writer3.writerows(new_raw_list2)
    csvfile3.close()

    csvfile4 = open(name4, 'w', newline='', encoding="utf-8")
    writer4 = csv.writer(csvfile4)
    writer4.writerow(["Gene", "Exon", "start", "end", "strand", "direction", "label_from", "label_to", "label_label"])
    writer4.writerows(new_raw_list0)
    writer4.writerows(new_raw_list2)
    csvfile4.close()
    print("")
    print("---  All files have been saved to folders!  ---")

else:
    print(
        "---  New folder (data_cis_splicing_gene) or (trans_cis_splicing_gene) already exists here. Please delete it "
        "or them first!  ---")
