# PanSearch

# Installation using conda

Clone this repo and run the following commands

`conda env create -f requirements.yaml`

`conda activate pansearch`

# Usage

Run the following command to view the help menu.

`python pansearch.py -h` 

OR 

`python pansearch.py -h`

The following arguments are mandatory:

1. --pangenome : the path to the folder containing your genome files in FASTA format. For example, if your files are in a folder named 'genomes' and it is in your current directory, you can provide the relative path (genomes) or absolute path (usually something like - /home/user/current_directory/genomes/).
2. --query : the path to the file containing the query sequences in FASTA format. For example, if the file named 'query.fa' is in your current directory, you can just provide the name of the file (query.fa) or its absolute path (usually something like - /home/user/current_directory/query.fa).
3. If you're running the script using the --nhmmer option, make sure you provide the path to an alignment file of your query sequences in stockholm (.sto) format or use the -a option to allow pansearch to generate one for you using your query file or a file specified with -i. You can convert .aln files to .sto files using the CLUSTAL alignment tool.

Default usage of the script containing only the required arguments will look like this:

`python pansearch.py -g genome_dir/ -q query.fa `

`python pansearch.py -g genome_dir/ -q query.fa -n True -st query.sto`

# Output

The script generates FASTA files containing the top hits for each query as a separate file for each genome. You can import them to geneious or any other alignment tool and perform multiple sequence alignments with them. There are also corresponding BED files and summary files for each FASTA file if you'd like more information on their identity and positional coordinates.

