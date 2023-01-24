from Bio import SeqIO
import pandas as pd 

def get_highly_expressed_genes(path_to_csv) -> list:
    '''
    read in .xlsx file and return a list of genes that are categorised as highly expressed
    '''
    df = pd.read_csv(path_to_csv)
    genes_of_interest = df[df['Category'] == 'Expressed in all high']['Ensembl gene id'].tolist()
    return genes_of_interest

def subset_fasta(fasta_path:str, genes_of_interest:list, output_path:str):
    '''
    subset a fasta file to only include genes of interest
    '''
    with open(output_path, 'w') as f:
        for record in SeqIO.parse(fasta_path, 'fasta'):
            if record.id.split('|')[1].split('.')[0] in genes_of_interest:
                f.write(f'>{record.id}\n{record.seq}\n')

def main():
    genes_of_interest = get_highly_expressed_genes('/home/jack/gene_expression_list.csv')
    subset_fasta('/home/jack/Downloads/gencode.v39.transcripts.fa', genes_of_interest[:1000], 'example/highly_expressed_genes.fa')

if __name__ == '__main__':
    main()
