from gffpandas import gffpandas as gffpd
import pandas as pd

#wczytuję plik z adnotacjami
annotation = gffpd.read_gff3('C:\\Users\\ania\\Documents\\Studia\\Licencjat\\multiexon.gff')


# Changing the format for colums instead of attributes to have an easier acces by an index
annotation_columns = annotation.attributes_to_columns()
# print(type(annotation_columns))
print("Pola w danych: ", annotation_columns.columns)
# Changing the format for a numpy array so that we can iterate and modify the data at the same time
data_array = annotation_columns.values

print(type(annotation))
print(type(annotation_columns))

data_scafold_0 = annotation_columns[annotation_columns['seq_id']=='scaffold_0']

# Select only the original 9 GFF3 columns
gff_columns = ['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
gff_df = data_scafold_0[gff_columns].copy()

# Replace NaNs in 'attributes' with empty strings
gff_df['attributes'] = gff_df['attributes'].fillna('')

# Save in GFF3-compliant format
gff_df.to_csv(
    'multiexon_0.gff',
    sep='\t',
    header=False,
    index=False,
    quoting=3  # optional: avoids adding quotes around strings
)

print(gff_df.head())
print(gff_df.dtypes)

