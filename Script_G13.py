
""" 
Bioinformatics and Computational Biology Project
2022/2023 - Group 13
Ana Araújo, nº 59457
Cláudia Afonso, nº 36273, 
João Faia - nº 47051
"""

# Importing the required modules and libraries
import gzip
from Bio import SeqIO, Entrez, Medline
from Bio.Blast import NCBIXML
import pandas as pd


# Setting the maximum width of a pandas dataframe column to diplay the full content
pd.set_option('display.max_colwidth', None)


######################### PROJECT GOAL 1 ##################################

################### DESCRIPTION OF THE GENOME FILES: ######################
###### NUMBER OF RECORDS, IDs, NAMES AND DESCRIPTION OF EACH RECORD #######

print()
print("PROJECT GOAL 1 - DESCRIPTION OF THE GENOME FILES")
print()

print('Using GenBank Files...')
print()


# Creating useful functions to avoid repeating code

def parse_gzip_files(fileName):

    """
    Opens a gzip-compressed file and parses its contents as GenBank records.

    This function takes a gzip file name as input and reads its contents 
    using SeqIO.parse() from the Biopython library. The contents of the 
    file are expected to be in GenBank format. The function returns a 
    list of parsed GenBank records.

    Args:
        fileName (str): The name of the gzip-compressed file to open and 
        parse.

    Returns:
        list: A list of parsed GenBank records from the file.
    """

    records = []

    with gzip.open(fileName, "rt") as handle:
        # Use SeqIO.parse() to parse the sequence from the compressed file
        records = list(SeqIO.parse(handle, "genbank"))

    return records


def create_records_dataframe(records):
    """
    Creates a pandas dataframe from a list of records

    This function takes a list of records as SeqRecord objects and extracts their 
    respective ID, name, and description. It then creates a pandas DataFrame with 
    these attributes as columns and returns the DataFrame.

    Args:
        records (list): A list of SeqRecord objects from GenBank records

    Returns:
        pandas.dataframe: A dataframe containing the ID, name and description of each 
        SeqRecord objects
    """

    ids = [record.id for record in records]
    names = [record.name for record in records]
    descriptions = [record.description for record in records]

    return pd.DataFrame({"ID": ids, "Name": names, "Description": descriptions})


# For Escherichia coli
print("Species 1: Escherichia coli")
print()

records_Coli = parse_gzip_files("GCF_000005845.2_ASM584v2_genomic.gbff.gz")
print()

# Print the number of records in the reference genome of Escherichia coli
print("There is %i record for the reference genome of Escherichia coli" % len(records_Coli))
print()

# Printing a dataframe with the ID, name and description of the single record in the reference genome of Escherichia coli
print(create_records_dataframe(records_Coli))
print()

print(26*"-")
print()

# For Vibrio cholerae
print("Species 2: Vibrio cholerae")
print()

records_Cholerae = parse_gzip_files("GCF_008369605.1_ASM836960v1_genomic.gbff.gz")

# Print the number of records in the reference genome of Vibrio cholerae
print("There are %i records for the reference genome of Vibrio cholerae" % len(records_Cholerae))
print()

# Printing a dataframe with the ID, name and description of the single record in the reference genome of Vibrio cholerae
print(create_records_dataframe(records_Cholerae))
print()


######################### PROJECT GOAL 2 ##################################

################## SEQUENCE ALIGNMENT OF THE GENOMES: #####################
###### ALIGNMENT OF THE SEQUENCES USING THE ALIGNMENT TOOLS AND ###########
######### IDENTIFICATION OF SIMILARITIES BETWEEN THE SEQUENCES ############

print()
print("PROJECT GOAL 2 - SEQUENCE ALIGNMENT OF THE GENOMES")
print()


# Reading the XML file resulting from the BLAST alignment of the two genomes

result_handle = open("6FFWMUDU114-Alignment.xml")
blast_record = NCBIXML.read(result_handle)

scores = []

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        score = 100*(hsp.identities/hsp.align_length)
        scores.append(score)

total_score = 0

for score in scores:
    total_score += score

mean_score = total_score / len(scores)

print(f'The similarity score between Escherichia coli and Vibrio cholerae is {mean_score:.2f}', '%')
print()

e_value_threshold = 0.000000000000000000000000000001
e_value_threshold_string = f"{e_value_threshold:.2e}"

total_hsp = 0

print('Sequence title, length and e-value for the alignments below a threshold of', e_value_threshold_string)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < e_value_threshold:    
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            total_hsp += 1

print()
print("The total number of alignments below an e-value threshold of", e_value_threshold_string, "is", total_hsp)
print()


######################### PROJECT GOAL 3 ##################################

######### RETRIEVAL OF RESEARCH PAPERS RELATED WITH EACH SPECIES ##########

print()
print("PROJECT GOAL 3A - RETRIEVAL OF RESEARCH PAPERS RELATED WITH EACH SPECIES")
print()


# Creating useful functions to avoid repeating code

def search_pubmed_records(query_term, retmax=10, sort="pub_date"):
    """
    Accesses Entrez to perform a search query to retrieve PMIDs.

    This function takes a query term and performs a search in the PubMed database
    using the Entrez.esearch() function from the Biopython library. It retrieves the
    PubMed IDs (PMIDs) of the search results and returns them as a sorted list in
    descending order. By default, it retrieves a maximum of 10 PMIDs, but you can 
    specify a different value using the `retmax` parameter.
    
    Args:
        query_term (str): The search term to query in PubMed.
        retmax (int, optional): The maximum number of PMIDs to retrieve. Defaults to 10.
        sort (string, optional): The method used to sort PMIDs in the ESearch output. Defaults to publication date.

    Returns:
        list: A sorted list of PMIDsof the search results  in descending order.
    """

    # Perform a search query using Entrez.esearch()
    handle = Entrez.esearch(db="pubmed", term=query_term, retmax=retmax, sort=sort)
    record = Entrez.read(handle)

    id_list = record["IdList"]

    return id_list


def retrieve_medline_records(id_list):

    """
    Retrieves and prints the details of PubMed records as Medline in plain-text format.

    This function takes a list of PubMed IDs (PMIDs) and uses the Bio.Entrez.efetch() 
    function from the Biopython library to retrieve the corresponding Medline records in
    plain-text format. It then iterates over the retrieved records and prints the article 
    title, authors, journal, and publication date for each record.

    Args:
        id_list (list): A list of PMIDs to retrieve and print the details of.

    Returns:
        None
    """

    # Using Bio.Entrez.efetch to retrieve the records
    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")

    # Reading Medline records one by one from the handle
    records = Medline.parse(handle)

    # Iterate over the PubMed records and print the article details
    i = 1
    for record in records:
        print("Article Number", i)
        print('Article Title:', record["TI"])
        print('Authors:', ", ".join(record["AU"]))
        print('Journal:', record["JT"])
        print('Publication Date:', record["DP"])
        print("---")
        i += 1  


Entrez.email = "fc36273@alunos.fc.ul.pt"  # Telling NCBI who we are


# Searching in PubMed for publications that include Escherichia coli in their title
print('Searching for the 10 most recent publications that include Escherechia coli in their title...')
print()

retrieve_medline_records(search_pubmed_records("Escherichia coli[title]"))
print()


# Searching in PubMed for publications that include Vibrio cholerae in their title
print('Searching for the 10 most recent publications that include Vibrio cholerae in their title...')
print()

retrieve_medline_records(search_pubmed_records("Vibrio cholerae[title]"))


######### SELECTION OF TWO GENES FROM EACH SPECIES ###############
########### AND RETRIEVAL OF THEIR 10 MOST RECENT PAPERS #########


print()
print("PROJECT GOAL 3B - RETRIEVAL OF RESEARCH PAPERS RELATED WITH TWO GENES OF EACH SPECIES")
print()


# Creating useful functions to avoid repeating code

def retrieve_gene_names(records):
    """
    Retrieves gene names from a list of SeqRecord objects obtained from GenBank records.

    This function takes a list of SeqRecord objects obtained from GenBank records, and 
    iterates over them to extract the names of the genes. It checks each feature in the records 
    and appends the gene name to a list if the feature type is 'gene' and if a gene name is
    available.

    Args:
        records (list): A list of SeqRecord objects from GenBank records.

    Returns:
        list: A list of all gene names contained within the GenBank records.
    """

    gene_names = []
    
    # Iterate over the records to see their genes
    for record in records:
        # Iterate over the features
        for feature in record.features:
            # Check if the feature is a gene
            if feature.type == "gene":
                # Append the name of the gene to the list
                gene_name = feature.qualifiers.get("gene", [""])[0]
                if gene_name:
                    gene_names.append(gene_name)
    
    return gene_names


def select_genes(possible_genes, gene_names):
    """
    Selects genes from a list of possible genes based on their presence in a list of gene names.

    This function takes a list of possible gene names and checks if each name exists in the provided
    list of gene names. The function returns a list of selected genes that are present in the provided
    list of gene names.

    Args:
        possible_genes (list): A list of gene names to consider for selection.
        gene_names (list): A list of gene names against which to check the possible genes of a given species.

    Returns:
        list: A list of selected genes that are present in the provided list of gene names.
    """

    selected_genes = []
    
    for possible_gene in possible_genes:
        if possible_gene in gene_names:
            if possible_gene not in selected_genes:
                selected_genes.append(possible_gene)
    
    return selected_genes


# For Escherichia coli

gene_names_Coli = retrieve_gene_names(records_Coli)

#print('List of gene names in Escherichia coli: ', gene_names_Coli)


# Selecting two genes from the Escherichia coli genome

possible_gene_names_Coli = ['recA', 'lacZ']

selected_genes_Coli = select_genes(possible_gene_names_Coli, gene_names_Coli)

print('The two selected genes of Escherichia coli are:', selected_genes_Coli[0], 'and', selected_genes_Coli[1])
print()


# Searching in PubMed for publications that include the selected genes in their title

print("Searching in PubMed for publications that include the gene 'recA' in their title...")
print()

# Searching in PubMed for publications that include recA in their title
retrieve_medline_records(search_pubmed_records("recA[title]"))
print()


print("Searching in PubMed for publications that include the gene 'lacZ' in their title...")
print()

# Searching in PubMed for publications that include lacZ in their title
retrieve_medline_records(search_pubmed_records("lacZ[title]"))
print()

# For Vibrio cholerae

gene_names_Cholerae = retrieve_gene_names(records_Cholerae)

#print('List of gene names in Vibrio cholerae: ', gene_names_Cholerae)

possible_gene_names_Cholerae = ['hlyA','rtxA']

selected_genes_Cholerae = select_genes(possible_gene_names_Cholerae, gene_names_Cholerae)

# Selecting two genes from the Vibrio Cholerae genome

print('The two selected genes of Vibrio cholerae are:', selected_genes_Cholerae[0], 'and', selected_genes_Cholerae[1])
print()


# Searching in PubMed for publications that include the selected genes in their title

print("Searching in PubMed for publications that include the gene 'hlyA' in their title...")
print()

# Searching in PubMed for publications that include hlyA in their title
retrieve_medline_records(search_pubmed_records("hlyA[title]"))
print()

print("Searching in PubMed for publications that include the gene 'rtxA' in their title...")
print()

# Searching in PubMed for publications that include rtxA in their title
retrieve_medline_records(search_pubmed_records("rtxA[title]"))

