import gffpandas.gffpandas as gffpd
import pandas as pd
import numpy as np
from collections import defaultdict

import time
# Save timestamp
start = time.time()

def group_introns(annotation_gff, groups):
    '''
    :param annotation_gff: a string path to a .gff file
    This is the file with annotations
    :param groups: a string path to an .lst format of genes grouped in rows.
    One row is one group.
    Genes in one row are overlapping each other and are potentially the same gene
    :return:
    '''

    import gffpandas.gffpandas as gffpd
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    from itertools import combinations
    # Loading file with annotations to analyze
    annotation = gffpd.read_gff3(annotation_gff)

    print()
    print("---------------------Dane na wejściu:")

    #Changing the format for colums instead of attributes to have an easier acces by an index
    annotation_columns = annotation.attributes_to_columns()
    #print(type(annotation_columns))
    print("Pola w danych: ", annotation_columns.columns)
    # Changing the format for a numpy array so that we can iterate and modify the data at the same time
    data_array = annotation_columns.values

    # Creating arrays with the selected types of records
    genes = annotation_columns[annotation_columns['type'] == 'gene'].values
    transcripts = annotation_columns[annotation_columns['type'] == 'mRNA'].values
    exons = annotation_columns[annotation_columns['type'] == 'exon'].values
    introns = annotation_columns[annotation_columns['type'] == 'intron'].values
    CDS = annotation_columns[annotation_columns['type'] == 'CDS'].values

    print("Liczba genów: ", len(genes))
    print("Liczba transkryptów: ", len(transcripts))
    print("Liczba eksonów: ", len(exons))
    print("Liczba intronów: ", len(introns))
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

    #variables for statistics
    not_found_genes = set()
    # gdzieś tu zliczanie ile miejsc w macierzy z zerami i bez
    # czyli ile genów z jednym transkryptem
    genes_with_one_transcript = set()
    how_many_transcripts_for_a_gene = []
    divided_groups = 0
    how_many_transcripts_from_one_gene =  {}
    csv_statistics = []

    #Przepisywanie danych do wygodnych struktur

    #stworzymy słownik. Kluczem jest jeden gen (z danych genes), wartością lista [transkrypt, set(intronów)]

    genes_dict = {}
    genes_introns_dict = {}
    #This dictionary is created to simplify the process of finding introns for each gene
    #Keys are transcripts' (mRNA records') ids and the key's values are ids of the gene the transcripts came from
    #and record of the transcript
    transcripts_dict = {}
    #key - intron's id, value - whole numpy.array record of an intron
    introns_dict = {}
    for gene in genes:
        gene_id = gene[9]
        genes_introns_dict[gene_id] = []
        genes_dict[gene_id] = gene
    # Every transcript in loaded annotation data will be added to the list[0] field in the value of it's parent's key
    for transcript in transcripts:
        parent = transcript[11]
        transcript_id = transcript[9]
        genes_introns_dict[parent].append(transcript)
        #in the first set, there wil be introns' ids stored
        genes_introns_dict[parent].append(set())
        #in the second, tuples containing the beginning and the end of an intron
        genes_introns_dict[parent].append(set())
        transcripts_dict[transcript_id] = [parent, transcript]
    # # Checking if gene/transcript data was correctly loaded to the dictionary. If yes, nothing should be printed
    # # All the genes should have only one transcript at this point
    # for gene in genes_introns_dict:
    #     if len(genes_introns_dict[gene]) != 1: print(gene, genes_introns_dict[gene])

    for intron in introns:
        intron_id = intron[9]
        parent_mRNA = intron[11]
        parent_gene = transcripts_dict[parent_mRNA][0]
        intron_start = intron[3]
        intron_end = intron[4]
        tuple = (intron_start, intron_end)
        # Adding a currently analysed intron to the set in a dictionary of its parent's key
        genes_introns_dict[parent_gene][1].add(intron_id)
        #genes_transcripts_disct[parent_gene].append()
        # Adding a tuple of the introns' edges
        genes_introns_dict[parent_gene][2].add(tuple)
        introns_dict[intron_id] = intron


    #print(genes_introns_dict)
    print('Liczba genów: ', len(genes))
    print('Liczba transkryptów: ', len(transcripts))
    print('Długość słownika (ile tam genów): ', len(genes_introns_dict))
    print()

    #now all the data from data_array (loaded .gff file) is loaded into the program,
    # so we can modify data_array without loosing any information

    # teraz przechodzimy po wszystkich grupach, robimy macierze i wyszukujemy splicing
    #analyzing each group in order to determine if all the transcripts are really from the same gene
    # and if the alternative splicing potentially occurs
    alternative_splicing_sites = []
    group_number = 0
    for group in overlap_list:
        # Creation of a matrix of zeros
        n = len(group)
        matrix = np.zeros((n, n), dtype=int)
        # The number in matrix[i, j] would be the number of overlapping introns in i and j gene in a group
        group_parent = group[0]
        empty_ids = []

        # print("Grupa: ", group)
        # print("parent grupy: ", group_parent)

        # Creating a list of pairs of genes (from one group) which will be compared
        # by searching for identical introns
        for gene_index1, gene_index2 in list(combinations(range(n), r=2)):
            # at the end of the for loop, the matrix could be symmetrical
            # but to optimize the time and memory we are comparing the genes only
            # when the first index is lower than the second (upper triangular matrix)

            #print('porównywane indeksy: ', gene_index1, gene_index2)

            gene1 = group[gene_index1]
            gene2 = group[gene_index2]

            #checking if the genes' ids are in the loaded data
            if gene1 in genes_introns_dict:
                record1 = genes_introns_dict[gene1]
            else:
                #if not, adding its ids to the not found genes
                not_found_genes.add(gene1)
                break
            if gene2 in genes_introns_dict:
                record2 = genes_introns_dict[gene2]
            else:
                not_found_genes.add(gene2)
                break

            N = count_identical_introns(record1, record2)
            M, overlapping_introns_set = overlapping_introns(record1, record2)
            # print(type(transcripts))
            matrix[gene_index1, gene_index2] = M

            #print("overlapping_introns_set: ", overlapping_introns_set)
            answer = set()

            #for further examination will be taken only the overlapping genes
            answer = find_alt_splicing(overlapping_introns_set)
            ans = (gene1, gene2, answer)
            #tutaj można też od razu zapisywać do pliku
            alternative_splicing_sites.append(ans)

        alternative_splicing_sites.append(['group_number{}'.format(group_number), group_number, "group" ])
        group_number += 1
        #     print("answer: ", answer)
        #
        # print(matrix)
        # print(alternative_splicing_sites)

        # if we find a whole row full of zeros, it means, that the gene group[index_of_the_row] has none introns
        # in common with the other genes in the group. We will add the id of the gene to a set of genes with this feature
        # to create a statistics at the end of the programme
        #The matrix is nxn so the lenght of the first row is also the number of the columns
        for i in range(len(matrix[0])):
            row = matrix[i]
            numbers = set(row)

            for j in range(len(matrix[0])):
                numbers.add(matrix[j][i])
            #if the row is full of zeros, the only number in the set would be 0
            if len(numbers) == 1:
                #group[i] stores the id of the gene
                genes_with_one_transcript.add(group[i])
                ###############################trzeba zobaczyć jak działa, gdy są same zera
                #if the group existed earlier, but it's ended now
                if group_parent != 0:
                    divided_groups += 1
                #reset of a group parent
                group_parent = 0
                #tutaj trzeba dodać zapisywanie, któe indeksy są puste
                empty_ids.append(i)
            #otherwise, we change the parent of the transcript
            if len(numbers) > 1:
                #if the last row and column was empty, we need to define a new parent
                if group_parent == 0:
                    # we make a current record from a group a new parent
                    group_parent = group[i]
                    break
                #if the group_parent already exists and has overlapping introns with
                #changing parents of transcripts

                transcript_id = genes_introns_dict[group[i]][0][9]
                old_record = transcripts_dict[transcript_id][1]
                old_parent = transcripts_dict[transcript_id][1][11]

                transcripts_dict[transcript_id][1][11] = group_parent

                #if we changed the parent above
                if old_parent != group_parent:
                    # changing attributes of the gene
                    old_attributes = transcripts_dict[transcript_id][1][8]
                    #print("old attributes: ", old_attributes)
                    attributes_list = old_attributes.split(';')
                    #print("attributes_list: ", attributes_list)
                    attributes_list[2] = 'Parent={}'.format(group_parent)
                    #print("attributes_list: ", attributes_list)
                    transcripts_dict[transcript_id][1][8] = ";".join(attributes_list)
                    #print(transcripts_dict)[1]


                    #changing type of duplicated gene
                    genes_dict[old_parent][2] = 'deserted'


            assert len(numbers)>= 1, "Error in the matrix {}".format(matrix)

        #ile średnio transkryptów na gen

    #counting in how many different groups alternative splicing has been detected
    splicing_in_how_many_genes = 0
    new_return = []
    for record in alternative_splicing_sites:
        if len(record[2]) > 0:
            new_return.append(record)
    for i in range(len(new_return)-1):
        record = new_return[i]
        next_record = new_return[i+1]
        if type(record[2])== str and type(next_record[2])==set:
            splicing_in_how_many_genes += 1

    print("-------------------------Wyniki:")

    print('liczba genów: ', len(genes))
    print('Liczba genów z jednym transkryptem: ', len(genes_with_one_transcript))
    print('Liczba genów analizowanych dalej: ', len(genes)-len(genes_with_one_transcript))
    print('Liczba grup na wejściu: ', len(overlap_list))
    print('Liczba podzielonych grup: ', divided_groups)
    print('Liczba genów z potencjalnym alt. splicingiem: ', splicing_in_how_many_genes)

    #print(genes_dict)

    #adding every record (and every type) to the array with results
    new_data_array = np.empty_like(data_array)
    i = 0
    for gene in genes_dict:
        new_data_array[i] = genes_dict[gene]
        i += 1
    for transcript in transcripts_dict:
        new_data_array[i] = transcripts_dict[transcript][1]
        i += 1
    for exon in exons:
        new_data_array[i] = exon
        i += 1
    for intron in introns:
        new_data_array[i] = intron
        i += 1
    for cd in CDS:
        new_data_array[i] = cd
        i += 1

    print("liczba rekordów w nowym pliku: ", len(new_data_array), " powinna być taka sama jak w pliku wejściowym")

    #ogólnie będziemy zwracać:
    # saving to a file
    #changing numpy.array type to dataFrame
    new_gff_data = change_to_dataframe(new_data_array)
    # data to be saved has to be in pd.DataFrame format
    save_gff_to_file(new_gff_data)

    #saving alternative splicing sites to a file
    with open('detected_splicing_sites.csv', 'w') as f:
        for record in new_return:
            f.write(f"{record}\n")

    #saving not found genes
    nfg = list(not_found_genes)
    with open('not_found_genes.csv', 'w') as f:
        for record in nfg:
            f.write(f"{record}\n")
    # csv ze statystykami genów
    # statystyki dla Rafała

    return 0

def count_identical_introns(gene1_dict_record, gene2_dict_record):
    '''
    :param gene1_dict_record: a list with a transcript's id at the [0], set of introns' ids at[1] and set of introns'
    start and end tuples at [2] for a key gene of a given record
    :param gene2_dict_record: a list with a transcript's id at the [0], set of introns' ids at[1] and set of introns'
    start and end tuples at [2] for a key gene of a given record
    :return: int: number of identical introns in the given records
    '''
    #under index 2, the set of tuples of begginings and ends is stored
    introns_set1 = gene1_dict_record[2]
    introns_set2 = gene2_dict_record[2]
    #Number of found identical introns in the two sets
    N = len(introns_set1.intersection(introns_set2))

    return N

def find_alt_splicing(overlapping_introns_set):
    '''

    :param overlapping_introns_set: set of tuples of tuples with start and end of an intron.
    Each tuple have two overlapping introns starts and ends
    :return:
    '''

    different_ends = set()
    pairs_list = list(overlapping_introns_set)

    #comparing introns' ends
    #if the number of identical genes is equal of the general number of introns from one or both of the genes,
    #there's no reason to compare this pair. There is no alternative splicing it such a case
    #Because of that, we compare the introns only if the number of introns from both genes is greater then N (the number of identical introns)

    for pair in pairs_list:

        #saving starts and ends in variables for easier naming
        start1 = pair [0][0]
        start2 = pair [1][0]
        end1 = pair [0][1]
        end2 = pair [1][1]


        #this will be the counter of identical ends
        identical_start = 0
        identical_end = 0

        #print("identical przed pętlą: ", identical_start, identical_end)

        if start1 == start2:
            identical_start += 1
        if end1 == end2:
            identical_end += 1

        identical = identical_start + identical_end

        #print("identical po pętli: ", identical_start, identical_end, identical)

        #it means both ends are the same, so there is no alternative splicing
        if identical == 2:
            break

        #if one end is identical, we add the tuple of the second ends indices to the "different_ends" set
        if identical == 1:
            if identical_start == 1:
                different_ends.add((end1, end2))
            if identical_end == 1:
                different_ends.add((start1,start2))

        #if introns have no common ends and they overlap each other (which is checked earlier)
        #there could be alternative splicing on both ends
        if identical == 0:
            different_ends.add((start1,start2))
            different_ends.add((end1,end2))
        assert identical <= 2, "There are more than 2 identical ends detected (impossible)"


    return different_ends

def overlap (gene1_tuple, gene2_tuple):
    '''
    :param gene1_tuple: tuple of integer start and of an intron
    :param gene2_tuple: tuple of integer start and of an intron
    :return: True/False value of overlap between two introns
    '''

    #taking beginnings and ends of introns from the tuples
    start1 = gene1_tuple[0]
    start2 = gene2_tuple[0]
    end1 = gene1_tuple[1]
    end2 = gene2_tuple[1]

    #checking if they overlap
    if start1<end2 and start2<end1:
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

    #print(gene1_dict_record)

    # range1 = range(len(gene1_dict_record[1]))
    # range2 = range(len(gene2_dict_record[1]))
    overlapping_introns_set = set()

    #we will compare every intron from one set with every from the second
    #we are using indices here to be able to link an id of an intron with its start and end positions
    introns_pairs = list(product(introns_set1, introns_set2))
    #here will be stored a number of overlapping introns
    N = 0

    for pair in introns_pairs:
        intron1_tuple = pair[0]
        intron2_tuple = pair[1]
        #print('para: ', pair[0], pair[1])
        if overlap(intron1_tuple, intron2_tuple):
            #print("nachodzą: ", pair[0], pair[1])
            N += 1
            overlapping_introns_set.add((intron1_tuple, intron2_tuple))


    return N, overlapping_introns_set

def change_to_dataframe(data_list):
    '''
    :param data_list: data in <class 'numpy.ndarray'> format
    :return: given array data in pandas.Dataframe format
    '''
    new_dataframe = pd.DataFrame(data=data_list, index=None, columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
           'attributes', 'ID', 'Name', 'Parent', 'Target', 'coverage', 'identity',
           'indels', 'matches', 'mismatches', 'unknowns'], dtype=object, copy=None)

    return new_dataframe

def save_gff_to_file(new_dataframe):
    '''
    :param new_dataframe: pandas.DataFrame
    :return: saves a file in GFF format in a current folder
    '''
    new_dataframe.to_csv(
        'changed_.gff',
        sep='\t',
        header=False,
        index=False,
        #quoting=3  # Avoid quoting strings unnecessarily
    )

group_introns('C:\\Users\\ania\\Documents\\Studia\\Licencjat\\multiexon.gff', 'C:\\Users\\ania\\Documents\\Studia\\Licencjat\\groups.lst')


# Save timestamp
end = time.time()
print()
print("Czas wykonania: ", end - start, " s, ", (end-start)/60, " min")