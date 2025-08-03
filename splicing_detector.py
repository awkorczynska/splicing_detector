#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import time


#this function reads the given options from command line and displays help if necessary
def parse():
	parser = argparse.ArgumentParser(description = 'This program finds alternative splicing sites on a gff annotation')
    #adding arguments with the description that will be shown in help
	parser.add_argument('-a', '--annotation', help = 'name (if in the same catalog) or path to a .gff file with annotation')
	parser.add_argument('-g', '--groups', help='name (if in the same catalog) or path to a .lst file with groups of genes IDs')
	parser.add_argument('-m', '--minimal_number', type = int, help = 'the minimal number required to match the groups', default = 1)
	parser.add_argument("-t", '--program_type', type = str, help = 'Type of the program. Should be either \'identical\' or \'overlap\'', default = 'identical')
	args = parser.parse_args()	#only arguments shortcut

	if args.annotation  == None:
		print('Please provide an annotation file')
	if args.groups == None:
		print('Please provide a groups file')

	print(f"You run the program with parameters: annotation = {args.annotation}, groups = {args.groups}, minimal_number = {args.minimal_number}, type = {args.program_type}")

	annotation_path = args.annotation
	groups_path = args.groups
	minimal_number = args.minimal_number
	program_type = args.program_type

	return annotation_path, groups_path, minimal_number, program_type


def group_introns(annotation_gff, groups, minimal_number = 1, program_type = 'identical'):
	'''

	:param annotation_gff: a string path to a .gff file
	This is the file with annotations
	:param groups: a string path to an .lst format of genes grouped in rows.
	One row is one group.
	Genes in one row are overlapping each other and are potentially the same gene
	:param minimal_number: a number of identical/overlapping introns in genes
	:param program_type: The type of common introns. Default identical, could be also overlap
	:return:
	'''
	assert type(minimal_number) == int, "Minimal number of common introns must be an integer"
	assert program_type in ['identical', 'overlap'], "Program type must be either 'identical' or 'overlap'"

	import gffpandas.gffpandas as gffpd
	#import pandas as pd
	import numpy as np
	#from collections import defaultdict
	from itertools import combinations
	# Loading file with annotations to analyze
	annotation = gffpd.read_gff3(annotation_gff)

	print()
	print("---------------------Input data:")

	# Changing the format for columns instead of attributes to have an easier access by an index
	annotation_columns = annotation.attributes_to_columns()
	column_index = {column: index for index, column in enumerate(annotation_columns.columns)}

	print("Data columns: ", annotation_columns.columns)
	# Changing the format for a numpy array so that we can iterate and modify the data at the same time
	data_array = annotation_columns.values

	# Creating arrays with the selected types of records
	genes = annotation_columns[annotation_columns['type'] == 'gene'].values
	transcripts = annotation_columns[annotation_columns['type'] == 'mRNA'].values
	exons = annotation_columns[annotation_columns['type'] == 'exon'].values
	introns = annotation_columns[annotation_columns['type'] == 'intron'].values
	CDS = annotation_columns[annotation_columns['type'] == 'CDS'].values

	print("Number of genes in input file: ", len(genes))
	print("Number of transcripts in input file: ", len(transcripts))
	print("Number of exons in input file: ", len(exons))
	print("Number of intorns in input file: ", len(introns))

	# Loading the list of genes grouped by overlapping each other and potentially being the same gene
	overlap_data = open(groups)

	# creating a list of lists in witch each inner list would be a list of gene IDs from each row from
	# a groups file
	overlap_list = []
	j = 0
	for i in overlap_data:
		i_list = i.split('\t')
		i_list[len(i_list) - 1] = i_list[len(i_list) - 1].strip()
		# Adding a list from a row to the general list of lists
		overlap_list.append(i_list)
		j += 1

	# variables for statistics
	not_found_genes = set()
	divided_groups = 0#######################################################################
	how_many_transcripts_from_one_gene = {}
	how_many_transcripts_from_one_gene[0] = []
	genes_with_detected_alternative_splicing = set()
	splicing_in_how_many_genes = 0

	#reloading the input data into dictionaires and lists for easier operating

	# genes_dict dictionary
	# key - gene ID
	# value - gene record
	genes_dict = {}

	# genes_introns_dict
	# key - gene ID
	#value - list of transcript record and set of introns ([transcript, set(introns from this thanscript)])
	genes_introns_dict = {}

	# transcripts_dict
	# key - transcripts' (mRNA records') ids
	# values - list of id of the gene the transcript came from and record of the transcript
	transcripts_dict = {}

	# introns_dict dictionary
	# key - intron's id,
	# value - whole numpy.array record of an intron
	introns_dict = {}

	#reloading the data
	for gene in genes:
		gene_id = gene[column_index["ID"]]
		genes_introns_dict[gene_id] = []
		genes_dict[gene_id] = gene

	for transcript in transcripts:
		parent = transcript[column_index["Parent"]]
		transcript_id = transcript[column_index["ID"]]
		genes_introns_dict[parent].append(transcript)
		# in the first set, there wil be introns' ids stored
		genes_introns_dict[parent].append(set())
		# in the second, tuples containing the beginning and the end of an intron
		genes_introns_dict[parent].append(set())
		transcripts_dict[transcript_id] = [parent, transcript]

	for intron in introns:
		intron_id = intron[column_index["ID"]]
		parent_mRNA = intron[column_index["Parent"]]
		parent_gene = transcripts_dict[parent_mRNA][0]
		intron_start = intron[column_index["start"]]
		intron_end = intron[column_index["end"]]
		tuple = (intron_start, intron_end)
		# Adding a currently analysed intron to the set in a dictionary of its parent's key
		genes_introns_dict[parent_gene][1].add(intron_id)
		# Adding a tuple of the introns' edges
		genes_introns_dict[parent_gene][2].add(tuple)
		introns_dict[intron_id] = intron

	print('Number of groups in input file: ', len(overlap_list))
	print("Number of records in input file: ", len(data_array))
	print()

	# now all the data from data_array (loaded .gff file) is loaded into the program,
	# so we can modify data_array without loosing any information


	# analyzing each group in order to determine if all the transcripts are really from the same gene
	# and if the alternative splicing potentially occurs
	# We will use upper triangular matrices filled with numbers of common introns between the genes in the group

	alternative_splicing_sites = []
	group_number = 0
	for group in overlap_list:
		# Creation of a matrix of zeros
		n = len(group)
		#the ids of genes in which alternative splicing is found
		alternative_splicing_sites_genes = set()
		#only creating the matrix if there are at least two genes to compare
		if n >=2:
			matrix = np.zeros((n, n), dtype=int)
			# The number in matrix[i, j] would be the number of overlapping introns in i and j gene in a group

			# Creating a list of pairs of genes (from one group) which will be compared
			# by searching for common (identical or overlapping) introns
			for gene_index1, gene_index2 in list(combinations(range(n), r=2)):
				# at the end of the for loop, the matrix could be symmetrical
				# but to optimize the time and memory we are comparing the genes only
				# when the first index is lower than the second (upper triangular matrix)

				gene1 = group[gene_index1]
				gene2 = group[gene_index2]
				record1 = 0
				record2 = 0

				# checking if the genes' ids are in the loaded data
				if gene1 in genes_introns_dict:
					record1 = genes_introns_dict[gene1]
				else:
					# if not, adding its ids to the not found genes
					not_found_genes.add(gene1)
				if gene2 in genes_introns_dict:
					record2 = genes_introns_dict[gene2]
				else:
					not_found_genes.add(gene2)

				#if both records are found
				if type(record1) != int and type(record2) != int:

					# N is a number of identical introns found between two genes
					N = count_identical_introns(record1, record2)
					# M is a number of overlapping introns in a group
					M, overlapping_introns_set = overlapping_introns(record1, record2)
					# print(type(transcripts))
					if program_type == 'identical':
						matrix[gene_index1, gene_index2] = N
					if program_type == 'overlap':
						matrix[gene_index1, gene_index2] = M

					answer = set()

					#if a found number of overlapping/identical introns is at least minimal, we search for alternative splicing sites
					if matrix[gene_index1, gene_index2] >= minimal_number:
						# for further examination will be taken only the overlapping genes
						answer = find_alt_splicing(overlapping_introns_set)
						# print("answer: ", answer)
						scaffold = genes_dict[gene1][column_index["seq_id"]]
						#we add genes to the detected set if the answer is not empty
						if answer != set():
							ans = (scaffold, gene1, gene2, answer)
							alternative_splicing_sites.append(ans)
							alternative_splicing_sites_genes.add(gene1)
							alternative_splicing_sites_genes.add(gene2)
				else:
					continue

			alternative_splicing_sites.append(['group_number{}'.format(group_number), group_number, 'group'])
			group_number += 1

			#now after finding alternative splicing sites we are changing the annotation so as the transcripts in the same place in the genome
			#have the same parent gene
			genes_dict, transcripts_dict, genes_introns_dict, not_found_genes, common_number_gene_dict = changing_parents(group, column_index, matrix, genes_dict, transcripts_dict, genes_introns_dict, not_found_genes, minimal_number)

			#counting in how many genes potential alternative splicing is found
			for gene in group:
				if gene not in not_found_genes:
					if genes_dict[gene][column_index["type"]] == 'gene':
						if gene in alternative_splicing_sites_genes:
							scaffold = genes_dict[gene][column_index["seq_id"]]
							genes_with_detected_alternative_splicing.add((scaffold, gene))
							splicing_in_how_many_genes += 1


			#adding the information about the genes in the index of how many genes have common introns with them
			for number in common_number_gene_dict:
				if number not in how_many_transcripts_from_one_gene:
					how_many_transcripts_from_one_gene[number] = []
				for gene in common_number_gene_dict[number]:
					how_many_transcripts_from_one_gene[number].append(gene)

	#counting average number of transcripts per gene
	sum = 0
	for number in how_many_transcripts_from_one_gene:
		sum += len(how_many_transcripts_from_one_gene[number])
	avg_number_of_transcripts_per_gene = sum/len(genes)

	# counting in how many different groups alternative splicing has been detected
	new_return = []
	for record in alternative_splicing_sites:
		if len(record[2]) > 0:
			new_return.append(record)


	# adding every record (and every type) to the array with results
	# Total number of records: count from all sources
	num_records = len(genes_dict) + len(transcripts_dict) + len(exons) + len(introns) + len(CDS)

	# Create empty array with shape (num_records, 9) to match gff format
	new_data_array = np.empty((num_records, 9), dtype=object)

	i = 0
	for gene in genes_dict:
		new_data_array[i] = genes_dict[gene][0:9]
		# print(genes_dict[gene][0:9])
		i += 1
	for transcript in transcripts_dict:
		new_data_array[i] = transcripts_dict[transcript][1][0:9]
		i += 1
	for exon in exons:
		new_data_array[i] = exon[0:9]
		i += 1
	for intron in introns:
		new_data_array[i] = intron[0:9]
		i += 1
	for cd in CDS:
		new_data_array[i] = cd[0:9]
		i += 1

	# saving to a file
	# changing numpy.array type to dataFrame
	new_gff_data = change_to_dataframe(new_data_array)
	# data to be saved has to be in pd.DataFrame format
	save_gff_to_file(new_gff_data)

	# saving alternative splicing sites to a file
	with open('detected_splicing_sites.csv', 'w') as f:
		for record in new_return:
			f.write(f"{record}\n")

	# saving not found genes
	nfg = list(not_found_genes)
	with open('not_found_genes.csv', 'w') as f:
		for record in nfg:
			f.write(f"{record}\n")

	# saving statistics
	statistics =[]

	statistics.append(['Statistics on file {}'.format(annotation_gff) + ' and {}'.format(groups)])
	statistics.append(['Parameters of the file: minimal_number = ' +str(minimal_number) + 'and  sensitivity: {}'.format(program_type)])
	statistics.append(['Number of genes in the input file : ', len(genes)])
	genes_after = new_gff_data[new_gff_data['type'] == 'gene']
	statistics.append(['Number of genes in the output file : ', len(genes_after)])
	statistics.append(['Number of genes that were assigned to an other parent gene (duplicates or short fragmented genes): ', len(genes)-len(genes_after)])
	statistics.append(['Number of genes with only one transcript', len(how_many_transcripts_from_one_gene[0])])
	statistics.append(['Average number of transcripts per gene: ', avg_number_of_transcripts_per_gene])
	#Zapisać jeszcze:
	#genes_with_detected_alternative_splicing
	#how_many_transcripts_from_one_gene


	print("-------------------------Results:")

	print('Number of genes in output file (after compressing): ', len(genes_after))
	print('Number of deserted genes: ', len(genes)-len(genes_after))
	print('Number of genes with only one transcript: ', len(how_many_transcripts_from_one_gene[0]))
	print('Number of genes with detected potential alternative splicing: ', splicing_in_how_many_genes)
	print("Number of records in output file: ", len(new_data_array), " (should be the same as in the input)")
	# print('Liczba grup na wejściu: ', len(overlap_list))
	# print('Liczba podzielonych grup: ', divided_groups)

	# saving statistics to a file
	with open('how_many_transcripts_from_one_gene.csv', 'w') as f:
		for number in how_many_transcripts_from_one_gene:
			f.write(f"List of genes with {number} detected transcripts\n")
			f.write(f"{how_many_transcripts_from_one_gene[number]}\n")

	# saving statistics to a file
	with open('genes_with_detected_alternative_splicing.csv', 'w') as f:
		for record in genes_with_detected_alternative_splicing:
			f.write(f"{record}\n")

	#saving statistics to a file
	with open('statistics.csv', 'w') as f:
		for record in statistics:
			f.write(f"{record}\n")
	return 0


def count_identical_introns(gene1_dict_record, gene2_dict_record):
	'''
	:param gene1_dict_record: a list with a transcript's id at the [0], set of introns' ids at[1] and set of introns'
	start and end tuples at [2] for a key gene of a given record
	:param gene2_dict_record: a list with a transcript's id at the [0], set of introns' ids at[1] and set of introns'
	start and end tuples at [2] for a key gene of a given record
	:return: int: number of identical introns in the given records
	'''
	# under index 2, the set of tuples of beginnings and ends is stored
	introns_set1 = gene1_dict_record[2]
	introns_set2 = gene2_dict_record[2]
	# Number of found identical introns in the two sets
	N = len(introns_set1.intersection(introns_set2))
	return N


def find_alt_splicing(overlapping_introns_set):
	'''
	:param overlapping_introns_set: set of tuples of tuples with start and end of an intron.
	Each tuple have two overlapping introns starts and ends
	:return: set of tuples of pairs of introns' ends that differ one from another
	'''

	different_ends = set()
	pairs_list = list(overlapping_introns_set)

	# comparing introns' ends
	for pair in pairs_list:
		# saving starts and ends in variables for easier naming
		start1 = pair[0][0]
		start2 = pair[1][0]
		end1 = pair[0][1]
		end2 = pair[1][1]

		# this will be the counter of identical ends
		identical_start = 0
		identical_end = 0

		if start1 == start2:
			identical_start += 1
		if end1 == end2:
			identical_end += 1

		identical = identical_start + identical_end

		# it means both ends are the same, so there is no alternative splicing
		if identical == 2:
			continue

		# if one end is identical, we add the tuple of the second ends indices to the "different_ends" set
		if identical == 1:
			if identical_start == 1:
				different_ends.add((end1, end2))
			if identical_end == 1:
				different_ends.add((start1, start2))

		# if introns have no common ends and they overlap each other (which is checked earlier in other function)
		# there could be alternative splicing on both ends
		if identical == 0:
			different_ends.add((start1, start2))
			different_ends.add((end1, end2))
		assert identical <= 2, "There are more than 2 identical ends detected (impossible)"

	return different_ends


def overlap(intron1_tuple, intron2_tuple):
	'''
	:param gene1_tuple: tuple of integer start and of an intron
	:param gene2_tuple: tuple of integer start and of an intron
	:return: True/False value of overlap between two introns
	'''

	# taking beginnings and ends of introns from the tuples
	start1 = intron1_tuple[0]
	start2 = intron2_tuple[0]
	end1 = intron1_tuple[1]
	end2 = intron2_tuple[1]

	# checking if they overlap
	if start1 < end2 and start2 < end1:
		return True
	else:
		return False


def overlapping_introns(gene1_dict_record, gene2_dict_record):
	'''
	:param gene1_dict_record: a list with a transcript's id at the [0], set of introns' ids at[1] and set of introns'
	start and end tuples at [2] for a key gene of a given record
	:param gene2_dict_record: a list with a transcript's id at the [0], set of introns' ids at[1] and set of introns'
	start and end tuples at [2] for a key gene of a given record
	:return: int: number of overlapping introns in the given records
	'''

	from itertools import product
	# under index 2, the set of tuples of beginnings and ends is stored
	introns_set1 = gene1_dict_record[2]
	introns_set2 = gene2_dict_record[2]

	overlapping_introns_set = set()

	# we will compare every intron from one set with every from the second
	# we are using indices here to be able to link an id of an intron with its start and end positions
	introns_pairs = list(product(introns_set1, introns_set2))
	# here will be stored a number of overlapping introns
	M = 0

	for pair in introns_pairs:
		intron1_tuple = pair[0]
		intron2_tuple = pair[1]
		#if the introns overlap each other
		if overlap(intron1_tuple, intron2_tuple):
			M += 1
			overlapping_introns_set.add((intron1_tuple, intron2_tuple))

	return M, overlapping_introns_set


def changing_parents(group, column_index, matrix, genes_dict, transcripts_dict, genes_introns_dict,  not_found_genes, minimal_number):
	'''
	:param group: an iterable list of genes IDs
	:param column_index: dictionary of names of columns (from input gff file) with values of these names indices
	:param matrix: upper triangular matrix. In [i][j] (i<j) a number of common introns between i and j genes from group is stored
	:param genes_dict: a dictionary of genes IDs with values of the genes record from input file
	:param transcripts_dict: a dictionary of transcripts IDs with values of [parent, transcript record]
	:param genes_introns_dict: a dictionary of genes IDs with values of [transcript record, set(introns IDs), set((start, end of an intron))]
	:param not_found_genes: set of genes IDs not found in input annotation file
	:param minimal_number: minimal number of common introns defined during starup
	:return: returns
					changed genes_dict,
					changed transcripts_dict,
					changed genes_introns_dict,
					changed not_found_genes,
					common_number_group_dict
	'''

	# creating a dictionary with True/False values if it has a minimal_number or more common introns with
	# others for each gene in the group
	gene_index_and_length = {}
	for i in range(len(group)):
		gene = group[i]
		if gene in not_found_genes:
			gene_index_and_length[i] = [i,0]
		if gene not in not_found_genes:
			gene_record = genes_dict[gene]
			gene_start = gene_record[column_index["start"]]
			gene_end = gene_record[column_index["end"]]
			gene_length = gene_end - gene_start
			gene_index_and_length[i] = [i, gene_length]

	# in this table there will be the index and length of the potential parent of the group stored
	the_longest_parent_to_be = [0, gene_index_and_length[0][0]]
	# we create a nem nxn matrix which will store the lengths of genes which fulfill the minimal number criteria
	n = len(group)
	matrix_for_parent = np.empty((n, n), dtype=object)

	# The matrix is nxn so the length of the first row is also the number of the columns
	for i in range(len(matrix[0])):
		# we place the length of the current gene in the place of it's index
		matrix_for_parent[i][i] = gene_index_and_length[i]


		for j in range(i+1, n):
			#if minimal number criterium id fulfilled
			if matrix[i][j] >= minimal_number:
				length_i = gene_index_and_length[i][1]
				length_j = gene_index_and_length[j][1]
				#if the current gene is longer then the one we are chceckig
				if length_i > length_j:
					matrix_for_parent[i][j] = gene_index_and_length[i]
				#if the length is the same, we choose to register the greater index
				if length_i <= length_j:
					matrix_for_parent[i][j] = gene_index_and_length[j]
			#if the minimal criteria are not met
			else:
				#we make a list of checked index and 0, because this record will no longer be taken into consideration
				#(0 is lower than any length and we are searching for the longest gene)
				matrix_for_parent[i][j] = [j,0]

	# in common_numbers_group_dict on i key will be the number of genes which the gene have an i number of common intorns
	common_number_group_dict = {}

	# now we have to check for which genes the criteria is met and put their lengths in the new matrix as well
	for i in range(n):
		the_longest_parent_to_be = matrix_for_parent[i][i]
		common_number = 0
		for j in range(n):
			if matrix_for_parent[i][j] is not None:
				#if the length of the index j gene is not zero
				if matrix_for_parent[i][j][1] > 0:
					common_number += 1
				#if the length of the index j gene is greater than the current parent_to_be
				if matrix_for_parent[i][j][1] > the_longest_parent_to_be[1]:
					the_longest_parent_to_be = matrix_for_parent[i][j]
				#if the length of the genes is the same
				if matrix_for_parent[i][j][1] == the_longest_parent_to_be[1]:
					if matrix_for_parent[i][j][0] > the_longest_parent_to_be[0]:
						the_longest_parent_to_be = matrix_for_parent[i][j]
			if matrix_for_parent[j][i] is not None:
				if matrix_for_parent[j][i][1] > 0:
					common_number += 1
				#and now we check the same thing, but in i column
				if matrix_for_parent[j][i][1] > the_longest_parent_to_be[1]:
					the_longest_parent_to_be = matrix_for_parent[j][i]
				#if the length of the genes is the same
				if matrix_for_parent[j][i][1] == the_longest_parent_to_be[1]:
					if matrix_for_parent[j][i][0] > the_longest_parent_to_be[0]:
						the_longest_parent_to_be = matrix_for_parent[j][i]

		#we have to substract the current gene (added in row and column) from the number of the other genes that have sth in common
		common_number = common_number - 2
		#at the index of the common number we add the name of the gene
		if common_number not in common_number_group_dict:
			common_number_group_dict[common_number] = []
		common_number_group_dict[common_number].append(group[i])


		#now we check if the found longest parent's index is differernt then ours
		#and if the answer is yes, we change parents
		if the_longest_parent_to_be[0] != group[i] and group[i] not in not_found_genes:
			new_parent = group[the_longest_parent_to_be[0]]
			current_gene = group[i]
			if new_parent != current_gene:
				transcript_id = genes_introns_dict[current_gene][0][column_index["ID"]]
				old_parent = current_gene
				transcripts_dict[transcript_id][1][column_index["Parent"]] = new_parent

				# changing attributes of the transcript and gene
				old_attributes = transcripts_dict[transcript_id][1][column_index["attributes"]]
				old_gene_attributes = genes_dict[old_parent][column_index["attributes"]]
				transcript_attributes_list = old_attributes.split(';')
				transcript_attributes_dict = {}
				for i in transcript_attributes_list:
					attribute = i.split('=')
					transcript_attributes_dict[attribute[0]] = attribute[1]

				gene_attributes_list = old_gene_attributes.split(';')
				transcript_attributes_dict['Parent'] = new_parent
				transcript_attributes_dict['Old_Parent'] = old_parent
				gene_attributes_list.append('Taken_over_by={}'.format(new_parent))

				#making a list of changed and added attributes
				transcript_attributes_list = []
				for attribute_name in transcript_attributes_dict:
					value = attribute_name + '=' + transcript_attributes_dict[attribute_name]
					transcript_attributes_list.append(value)

				#changing the attributes of the transcript and the old_parent gene
				transcripts_dict[transcript_id][1][column_index["attributes"]] = ";".join(transcript_attributes_list)
				genes_dict[old_parent][column_index["attributes"]] = ";".join(gene_attributes_list)
				# changing type of duplicated gene
				genes_dict[old_parent][column_index["type"]] = 'deserted'

	#returning changed annotation records and other modified variables
	return genes_dict, transcripts_dict, genes_introns_dict, not_found_genes, common_number_group_dict


def change_to_dataframe(data_list):
	'''
	:param data_list: data in <class 'numpy.ndarray'> format
	:return: given array data in pandas.Dataframe format
	'''
	new_dataframe = pd.DataFrame(data=data_list, index=None,
								 columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
										  'attributes'], dtype=object, copy=None)

	return new_dataframe


def save_gff_to_file(new_dataframe):
	'''
	:param new_dataframe: pandas.DataFrame
	:return: saves a file in GFF format in a current folder
	'''
	new_dataframe.to_csv(
		'changed_annotation.gff',
		sep='\t',
		header=False,
		index=False,
		quoting=3  # Avoid quoting strings unnecessarily
	)


if __name__ == '__main__':
	# Save timestamp
	start = time.time()
	annotation_path, groups_path, minimal_number, program_type = parse()
	group_introns(annotation_path, groups_path, minimal_number, program_type)

	# Save timestamp
	end = time.time()
	print()
	print("Run time: ", end - start, " s, ", (end - start) / 60, " min")