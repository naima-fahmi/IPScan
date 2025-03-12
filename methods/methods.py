import glob, os, csv, sys, time
from operator import attrgetter
import numpy as np
import pandas as pd
import bisect
import warnings
from bisect import bisect_left
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
csv.field_size_limit(sys.maxsize)

from methods.print_and_log import print_and_log

warnings.simplefilter(action="ignore", category=FutureWarning)

def bi_contains(lst, item):
    return bisect_left(lst, item)

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()

def count_total_read_count(chrom, start, end, bam_list, position_row):
	totalCount = 0
	length = end - start + 1
	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		if(int(bam_list[pos2][0]) != end):
			pos2 = pos2 - 1

		for t in range(pos1, pos2+1):
			read = int(bam_list[t][1])
			totalCount += read
		
	return totalCount, length


def make_full_dictionary(ann_df, chromosomes):
    ChromDict = {}
    for chrom in chromosomes:
        GeneDict = {}
        chr_rows = ann_df[ann_df["chrom"] == chrom]
        gene_list = list(set(chr_rows["name2"]))
        for gene in gene_list:
            gene_rows = chr_rows[chr_rows["name2"] == gene]
            exList = []
            for index, row in gene_rows.iterrows():
            	strand = row["strand"]
            	exonCount = row["exonCount"]
            	exonStarts = list(filter(None, row["exonStarts"].split(",")))
            	exonEnds = list(filter(None, row["exonEnds"].split(",")))
            	for i in range(exonCount):
            		st, en = int(exonStarts[i]), int(exonEnds[i])
            		if (st, en) not in exList:
            			exList.append((st, en))

            GeneDict[gene.strip().upper()] = (exList, strand)

        ChromDict[chrom] = GeneDict

    return ChromDict


def merge_intervals(inputlist):
	n = len(inputlist)
	inputlist.sort(key = lambda x: x[0], reverse = False)

	st = Stack()
	st.push(inputlist[0])

	for i in range(1,n):
		stacktop_start, stacktop_end = st.top()
		if inputlist[i][0]<= stacktop_end:
			st.pop()
			stacktop_end = max(stacktop_end, inputlist[i][1])
			st.push((stacktop_start, stacktop_end))
		else:
			st.push(inputlist[i])

	mergedList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		mergedList.append(stacktop)
		st.pop()

	mergedList.sort(key = lambda x: x[0], reverse = False)
	return mergedList

def merge_ChromDict(ChromDict, chromosomes):
    ChromDict_merged = {}
    for chrom in chromosomes:
        GeneDict_merged = {}
        GeneDict = ChromDict[chrom]
        for gene in GeneDict.keys():
            (exonList, strand) = GeneDict[gene.upper()]
            mergedExonList = merge_intervals(exonList)
            GeneDict_merged[gene.upper()] = (mergedExonList, strand)
        ChromDict_merged[chrom] = GeneDict_merged

    return ChromDict_merged

def get_intron_list(exonList):
	intronList = []
	start = exonList[0][1] + 1
	for item in exonList[1:]:
		end = item[0] - 1
		intronList.append((int(start), int(end)))
		start = item[1] + 1

	return intronList


def generate_ipscan_events(species, annotation_file, sample, output_dir, logger):
	g1_name = os.path.basename(sample)
	print_and_log("----------------------------------------------------------------", logger)
	print_and_log("| STEP 1: Detecting IPA events for sample: "+g1_name+"...                                  |", logger)
	print_and_log("----------------------------------------------------------------", logger)

	start = time.perf_counter()

	chromosomes_h = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY"
    ]

	chromosomes_m = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chrX",
        "chrY"
    ]

	# Determine the chromosomes based on the species
	if species == "hg38" or species == "hg19":
	    chromosomes = chromosomes_h
	elif species == "mm10" or species == "mm39":
	    chromosomes = chromosomes_m
	else:
	    print("Species not found. Please select among hg38, hg19, mm10 or mm39")
	    sys.exit()


	ann_df = pd.read_csv(annotation_file, delimiter = '\t')
	# Convert the whole annotation into a dictionary for faster use
	print_and_log("Creating the chromosome dictionary...", logger)
	start = time.perf_counter()
	ChromDict = make_full_dictionary(ann_df, chromosomes)
	print_and_log(f"Chromosome dictionary created in {time.perf_counter() - start:.2f} seconds", logger)


	# Merge the Exon Intervals
	print_and_log("Merging the Exon intervals...", logger)
	start = time.perf_counter()
	ChromDict_merged = merge_ChromDict(ChromDict, chromosomes)
	print_and_log(f"Exon intervals merged in {time.perf_counter() - start:.2f} seconds", logger)

	peak_df = pd.read_csv("./anno/"+species+"_polyA_sites.csv", delimiter='\t')

	writer_list = []
	output_columns = [
			'Chrom', 
			'Gene Name', 
			'Strand', 
			'Intron Start', 
			'Intron End', 
			'Peak position', 
			'Ratio_1',
			'Ratio_2',
			'Ratio_1/Ratio_2', 
			'R1',
			'R2',
			'R3',
			'L1',
			'L2',
			'L3',
			'C1',
			'C2',
			'C3'
	]

	for chrom in chromosomes:
		#print("Starting:",chrom)
		st = time.time()

		if os.path.getsize(os.path.join(sample, chrom + ".txt")) > 0:
			bam_df = pd.read_csv(os.path.join(sample, chrom + ".txt"), delimiter="\t")
			position_row = bam_df.iloc[:, 0].tolist()
			bam_list = bam_df.values.tolist()

			flag = 0
			GeneDict_merged = ChromDict_merged[chrom]
			gene_count = 0
			for geneID in GeneDict_merged.keys():
				pos_flag, pos_flag2 = [], []
				peak_rows = peak_df.loc[(peak_df['Chrom']==chrom) & (peak_df['Gene Name']==geneID)]
				#print("len(peak_rows)", len(peak_rows))

				if len(peak_rows) > 0:
					(mergedExonList, strand) = GeneDict_merged[geneID]
					range_start, range_end = min(mergedExonList)[0], max(mergedExonList)[1]
					if (range_end - range_start) <= 2500000:     ### maximum gene length of human. to stop genes occurring in two regions in the annotation, which makes this range very big
						gene_count += 1
						intronList = get_intron_list(mergedExonList)
						
						R1, L1 = 0, 0
						for (exonStart, exonEnd) in mergedExonList:
							Re, Le = count_total_read_count(chrom, exonStart, exonEnd, bam_list, position_row)
							#print(exonStart, exonEnd, Re, Le)
							R1 += Re
							L1 += Le
						
						for row in peak_rows.itertuples():
							peak_pos_list = (str(row[4]).strip('[ ]')).split(',')
							peak_pos_list = [int(p) for p in peak_pos_list]
							peak_pos_list.sort()

							for ind in range(len(intronList)):
								(intron_start, intron_end) = intronList[ind]
								for peak_pos in peak_pos_list:
									peak_pos = int(peak_pos)
									if (peak_pos >= intron_start) and (peak_pos<=intron_end+1):
										
										if strand == '+':
											R2, L2 = count_total_read_count(chrom, intron_start, peak_pos, bam_list, position_row)
											R3, L3 = count_total_read_count(chrom, peak_pos+1, intron_end, bam_list, position_row)
										else:
											R2, L2 = count_total_read_count(chrom, peak_pos, intron_end, bam_list, position_row)
											R3, L3 = count_total_read_count(chrom, intron_start, peak_pos-1, bam_list, position_row)
										
										C1, C2, C3 = 0, 0, 0
										ratio_1, ratio_2 = 0, 0
										if 0 not in [L1,L2,L3]:
											#print([R1,R2,L1,L2,L3])
											if R1>0:
												C1 = float(R1)/L1
											if R2 > 0:
												C2 = float(R2)/L2
											if R3 > 0:
												C3 = float(R3)/L3

											if C1 > 0:
												ratio_1 = float(C2)/C1
												ratio_2 = float(C3)/C1

											if ratio_2>0:
												found += 1
												writer_list.append(
													(
														chrom, 
														geneID, 
														strand, 
														intron_start, 
														intron_end, 
														peak_pos, 
														ratio_1, 
														ratio_2, 
														round(ratio_1/ratio_2, 3), 
														R1, 
														R2, 
														R3, 
														L1, 
														L2, 
														L3, 
														C1, 
														C2, 
														C3
													)
												)
												#print("Found", chrom, geneID, strand, intron_start, intron_end, peak_pos, ratio_1, ratio_2, round(ratio_1/ratio_2, 3), R1, R2, R3, L1, L2, L3, C1, C2, C3)


	df_out = pd.DataFrame(writer_list, columns=output_columns)
	df_out.to_csv(os.path.join(output_dir, "IPScan_"+sample+"_full.csv"), sep="\t", index=False)

	C2_cutoff = 10
	df_out['C2/C1'] = df['C2']/df['C1']
	newdf = df.loc[(df['C2']>C2_cutoff)&(df['C2/C1']<=0.2)&(df['C1']>df['C2'])&(df['C1']>df['C3'])&(df['C2']>df['C3'])]
	#print(len(newdf))
	newdf.to_csv(os.path.join(output_dir, "IPScan_"+sample+".csv"), sep='\t', header=True, index=False)
	print_and_log(f"IPA sites detected in {time.perf_counter() - start:.2f} seconds", logger)