# Homologous Species Identifier

This script is designed to retrieve a gene, return a FASTA file of its protein sequence and identify homologous species.

## Installation

Ensure you have Python installed on your system. This script requires the following dependencies:

- `argparse`
- `requests`
- `BioPython`

You can install these dependencies using pip:

```pip install argparse requests biopython```

## Usage

You can run the script by executing the following command:

```python pipeline.py -g <gene_id>```


Replace `<gene_id>` with the desired gene identifier.

## Functions

### `getEnsemblId(gene: str) -> str`

Given a gene ID, this function returns the Ensembl ID within the human species.

### `getDNASeq(ensembl_gene_id: str) -> str`

Given an Ensembl gene ID, this function returns the DNA sequence in FASTA format as a string.

### `getLongestORF(fasta: str) -> str`

Given a FASTA sequence, this function identifies the longest open reading frame (ORF) and returns it in FASTA format.

### `findHomologs(ensembl_gene_id: str, gene: str = "M1CR") -> None`

Given an Ensembl gene ID, this function creates a list of homologous species and writes it to a text file named `<gene>_homology_list.txt`.

### `main()`

The main function orchestrates the entire process, including parsing command-line arguments, executing the necessary functions, and printing the process completion time.

## Output

`<ensembl_gene_id>.fasta`

- A FASTA file of the amino acid sequence

`<gene>_homology_list.txt`

- A list of unique homologous species

## Example

```python pipeline.py -g MC1R```


This command retrieves information about the MC1R gene, including Ensembl ID, protein sequence, and homologous species.


Steve Cheney
RBIF100 A4