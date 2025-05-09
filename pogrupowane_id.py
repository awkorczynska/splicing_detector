import gffpandas.gffpandas as gffpd
import pandas as pd
import numpy as np
from collections import defaultdict

#wczytuję plik z adnotacjami
annotation = gffpd.read_gff3('C:\\Users\\ania\\Documents\\Studia\\Licencjat\\overlap_multiexon_100.gff')

#zamieniam format danych na taki z kolumnami zamiast atrybutów (żeby łatwiej dotrzeć do poszczególnych pól
annotation_columns = annotation.attributes_to_columns()
print("Pola w danych: ", annotation_columns.columns)
#zamieniłam dane na NUMPY.ARRAY, bo nie znalazłam sposobu, żeby iterować i modyfikować jednocześnie
#atrybuty w Dataframe
data_list = annotation_columns.values
#print(data_list)

#############################################nowe dane od Maćka
#wczytuję plik z adnotacjami dopasowanymi do tych samych genów
overlap_list = open('C:\\Users\\ania\\Documents\\Studia\\Licencjat\\overlap_groups.lst')
#/\ to jest lista wierszy, z których każdy zawiera po kilka rekordów 'Name' pochodzących z tego samego genu
#Nazwy 'Name' są różne, ale efeltywnie pochodzą z tego samego genu

#Tworzę słownik z tej listy. Kluczem jest pierwszy rekord z wiersza.
#Na podstawie tych rekordów-kluczy dopasuję je do genu, z którego pochodzą
multiple_dict = {}
i_list =[]
for i in overlap_list:
    i_list = i.split('\t')
    i_list[len(i_list)-1] = i_list[len(i_list)-1].strip()
    #pierwsze ID będzie kluczem, a reszta wartościami dla tego klucza
    multiple_dict[i_list[0]] = i_list[1:]


#ten kawałek zamienia atrybut "Parent" na podstawie tabelki z odpowiadaniem genów
#po wykonaniu całej pętli wszystkie transkrypty pochodzące z tego samego miejsca (jeden wiersz w pliku overlap_groups)
#mają jako rodzica wpisany gen o tym samym ID. Jest ID genu z pierwszej kolumny pliku overlap_groups
for record in data_list:
    #pole record[2] zawiera typ rekordu
    if record[2] == 'mRNA':
        #pole record[11] zawiera id rodzica
        parent = record[11]
        for key in multiple_dict.keys():
            if parent in multiple_dict[key]:
                record[11] = key
                break

#Na tym etapie w data_list mamy geny, które są w polach "Parent" transkryptów i takie geny,
#które zostały podmienione i nie mają już żadnych transkryptów jako "dzieci". Te drugie geny
#nie są nam potrzebne, więc usuniemy je z tabeli.
#Wystarczy znaleźć i usunąć geny, których ID znajduje się w tabelach wartośći w słowniku multiple_dict

#data_list TO JEST JEDNAK ARRAY TRZEBA POZMIENIAĆ WSZĘDZIE---------------------------------------------------------------------------------------
# licznik usuniętych wierszy
j = 0
for i in range(len(data_list)):
    #przesuwam indeks o tyle, ile wierszy zostało już usuniętych
    i = i-j
    for key in multiple_dict.keys():
        #pole 9 zawiera ID
        if data_list[i][9] in multiple_dict[key]:
            data_list = np.delete(data_list, (i), axis=0)
            #ponieważ usuwam wiersz z danych, na których pracuję, muszę cofnąć indeks o jeden
            j += 1

#teraz chcę pogrupować dane w taki sposób, by w jednej grupie były wszystkie transkrypty
#i ich exony. Potem w obrębie tych grup sprawdzę, czy wszystkie exony zaczynają się w tych samych
#miejscach, czy występują między nimi jakieś przesunięcia

#z danych data_list wybieram tylko te rekordy, kóre są transkryptami
transcripts = []
for i in data_list:
    if i[2] == 'mRNA':
        transcripts.append(i)
transcripts = np.array(transcripts)

# Tworzę słownik, w którym kluczami są ID genów, a wartościami rekordy transkryptów,
# dla których dany gen jest rodzicem
grouped_by_parent = defaultdict(list)
#uzupełniam słownik grupując po polu "Parent"
for record in transcripts:
    parent = record[11]  # Index 11 jest tożsamy z polem 'Parent'
    #dodaję tylko ID tych transkryptów
    grouped_by_parent[parent].append(record[9])
# Można tu jeszcze pousuwać grupy z tylko jednym ID transkryptu, czyli takie, że
# grouped_by_parent[parent]==1
# Na razie zostawiam tak, bo łatwo będzie sprawdzić czy program dobrze działa -
# nie powinien wykryć alt. spl. tam gdzie jest tylko jeden transkrypt

# zamiana list arrayów na array arrayów. Nie wiem jeszcze, który format będzie lepszy
grouped_by_parent = {k: np.array(v) for k, v in grouped_by_parent.items()}
print('grouped_by_parent', grouped_by_parent)
# print(len(transcripts))
# print(len(grouped_by_parent))


#znajduję wszystkie exony
exons = []
for i in data_list:
    if i[2] == 'exon':
        exons.append(i)
exons = np.array(exons)

# Teraz zróbmy listy exonów w obrębie grup
# Jak już będzie działało, wrzucę to wszystko w pętlę, aby wykonywało się dla każdego
# "zestawu" exonów dla każdego genu.
# ten sposób oszczędzimy trochę pamięci. Nie trzeba będzie pamiętać wyników dla
# poprzednich grup, bo będą już zapisane w odpowiednim folderze



# Tworzę słownik, w którym kluczem jest gen danej grupy, a wartościami
# lista wszystkich exonów pochodzących od wszystkich transkryptów od każdego genu
exons_by_gene = defaultdict(list)
#    uzupełniam słownik grupując po polu "Parent". Tym razem rodzicami będą transkrypty
# o zgromadzonych w liście grouped_by_parent[parent] ID

for gene in grouped_by_parent.keys():
    exons_by_gene[gene] = []
    print("gene: ", gene)
    print("długość: ", len(grouped_by_parent[gene]))
    for transcript_id in grouped_by_parent[gene]:
        for exon in exons:
            if exon[11] == transcript_id: #pole 11 jest tożsame z atrybutem "Parent"
                exons_by_gene[gene].append(exon)
                print('Parent: ', exon[11])

#print('exons_by_gene', exons_by_gene)

# #Zmiana formatu i zapisanie do pliku gff
#
# #zamieniam zmienioną tabelę z powrotem na dataframe
# new_dataframe = pd.DataFrame(data=data_list, index=None, columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
#        'attributes', 'ID', 'Name', 'Parent', 'Target', 'coverage', 'identity',
#        'indels', 'matches', 'mismatches', 'unknowns'], dtype=object, copy=None)
#
# #Działające zapisywanie wyników do pliku gff
# new_dataframe.to_csv(
#     'podmienione_i _usuniete_geny_multiexon_100.gff',
#     sep='\t',
#     header=False,
#     index=False,
#     quoting=3  # Avoid quoting strings unnecessarily
# )
#
# # wyniki = gffpd.read_gff3('C:\\Users\\ania\\Documents\\Studia\\Licencjat\\Splicing\\podmienione_geny_multiexon_100.gff')
# # print(wyniki)