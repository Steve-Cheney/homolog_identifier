import argparse
import re
import time
import requests
from Bio.Seq import Seq
from Bio import SeqIO


def getEnsemblId(gene: str ="MC1R") -> str:
    """
    Given a gene id, return the ensembl id within the human species.

    Params
    ------
    gene
        The gene id to lookup.
   
    Returns
    -------
    The ensembl id.
    """   

    server = "http://mygene.info/v3"
    params = {
        "q": gene,
        "species": "human"
    }
    endpoint = "/query"
    r = requests.get(server+endpoint, params=params)
    if r.status_code == 200:
        r = r.json()
        e_gene_id = r["hits"][0]["_id"]
        r = requests.get(server+"/gene/"+e_gene_id)
        if r.status_code == 200:
            r = r.json()
            return r["ensembl"]["gene"]
        else:
            print("No results found for given entrez id - Code:", r.status_code)
            return None
    else:
        print("No results found for given gene - Code:", r.status_code)
        return None


def getDNASeq(ensembl_gene_id: str) -> str:
    """
    Given a an ensembl gene id, return DNA sequence in FASTA format as a string.

    Params
    ------
    ensembl_gene_id
        The ensembl gene id to lookup.
   
    Returns
    -------
    The DNA sequence in FASTA format.
    """   

    json_header = {'Content-Type': 'application/json'}
    ensembl_server = 'http://rest.ensembl.org'
    endpoint = f'/sequence/id/{ensembl_gene_id}?content-type=text/plain'
    
    r = requests.get(ensembl_server+endpoint, headers = json_header)
    
    if r.status_code == 200:
        r = r.json()
        return f">{ensembl_gene_id}\n{r['seq']}"
    else:
        print("No results found for given gene - Code:", r.status_code)
        return None


def getLongestORF(fasta):
    """
    Given a fasta, return DNA sequence in FASTA format as a string.
    Adapted from https://stackoverflow.com/a/31758161

    Params
    ------
    ensembl_gene_id
        The ensembl gene id to lookup.
   
    Returns
    -------
    The DNA sequence in FASTA format.
    """   
    lines = fasta.split('\n')
    fasta_no_header = lines[1]
    return max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)', fasta_no_header), key = len)


def findHomologs(ensembl_gene_id: str, gene: str = "M1CR") -> None:
    """
    Given an ensembl gene id, create a list of homologous species.

    Params
    ------
    ensembl_gene_id
        The ensembl gene id to lookup.
    gene
        The gene name given as the initial argument
        
    Returns
    -------
    A text file called <gene>_homology_list.txt with a list of unique homolog species names.
    """   

    ensembl_server = f"https://rest.ensembl.org/homology/id/{ensembl_gene_id}"

    r = requests.get(ensembl_server, headers={"Content-Type": "application/json"})

    if r.status_code == 200:
        homologs = r.json()
        species_list = []

        for homology in homologs["data"][0]["homologies"]:
            #print(homology)
            
            target= homology["target"]
            species = target["species"]
            if species not in species_list:
                species_list.append(species)
            

        # Write species names to a text file
        with open(f"{gene}_homology_list.txt", "w") as f:
            for species_name in species_list:
                f.write(species_name + "\n")
    else:
        print("Failed to retrieve homologs - Code:", r.status_code)

def main():
    start_time = time.perf_counter()
    # Define arguments
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-g", "--gene", required=True, help="Gene id")
    
    args = parser.parse_args()

    print("Getting ensembl ID...")
    e_id = getEnsemblId(args.gene)
    print("Getting FASTA...")
    fasta = getDNASeq(e_id)
    print("Getting ORF...")
    seq = getLongestORF(fasta)
    # Create a Seq object
    seq = Seq(seq)
    # Translate the DNA sequence to amino acids
    aa_seq = seq.translate()

    fasta_info = f">{e_id}\n{aa_seq}"
    fasta_file_path = f"{e_id}.fasta"

    # Write the FASTA information to the file
    with open(fasta_file_path, "w") as fasta_file:
        fasta_file.write(fasta_info)

    print("Finding homologs...")
    findHomologs(e_id)

    end_time = time.perf_counter()
    print(f"\n\nProcess completed in {round(end_time-start_time, 3)} seconds.")

if __name__ == "__main__":
    main()
