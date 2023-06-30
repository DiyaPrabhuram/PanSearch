#import statements
import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#arguments
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--pangenome", help="enter the path of the directory your genome files are stored", required=True)
parser.add_argument("-q", "--query", help="enter the path to the query file containing sequences you'd like to search against the genomes", required=True)
parser.add_argument("-a", "--align", help="add True to this argument if you'd like a clustal alignment of your final hits", default="False")
parser.add_argument("-i", "--alignment_input", help="enter the path to the file that you would like to use to generate a stockholm alignment", default="0")
parser.add_argument("-e", "--evalue", help="enter the expect value you want blast to use (in decimal, for example: 0.00001)", default=0.0000001)
parser.add_argument("-o", "--output_path", help="please specify the path of your desired output folder")
parser.add_argument("-n", "--nhmmer", help="add True to this argument if you'd like to run nhmmer", default="False")
parser.add_argument("-ex", "--extension", help="specify the number of base pairs you'd like to extend your hits", default="750")
parser.add_argument("-qc", "--query_coverage", help="specify the query coverage cut-off value for BLAST. (as an integer, for example: 80)", default=80, type=int)
parser.add_argument("-pi", "--percent_identity", help="specify the percent identity cut-off value for BLAST. (as an integer, for example: 70)", default=70, type=int)
parser.add_argument("-ne", "--nhmmer_evalue", help="enter the expect value you want nhmmer to use (in decimal, for example: 0.00001)", default="0.000001")
parser.add_argument("-st", "--stockholm", help="please provide an alignment of your query sequences in stockholm format")
args = parser.parse_args()

if (args.nhmmer) and not (args.stockholm or args.align):
	sys.exit("ERROR: please provide an alignment file in stockholm format if you're running nHMMER.\n")

'''
Notes for README:
1. arguments- genome.fa path, 3-letter identifier, query.fa path, evalue, ncbiblast bin path
2. install bedtools using sudo apt-get install bedtools
3. install biopython
4. install ncbiblast
5. to run nhmmer, make sure the query contains CDS sequences

Issues to sort out:
5. separate directories for each output type maybe
'''

def generate_table(blast_output, idt, pident, qcover, hmm=False):

	file = open(blast_output, 'r')

	query_name = []
	top_hit_name = []
	percent_identity = []
	query_coverage = []
	start = []
	stop = []
	strand = []
	for line in file.readlines():
		if not line.startswith("#"):
			attributes = line.split(',')
			if not hmm:
				if attributes[0] not in query_name:
					query_name.append(attributes[0])
					top_hit_name.append(attributes[1])
					percent_identity.append(attributes[5])
					query_coverage.append(attributes[4])
					if attributes[10].strip() == "plus":
						start.append(int(attributes[8]))
						stop.append(attributes[9])
						strand.append("+")
					else:
						start.append(int(attributes[9]))
						stop.append(attributes[8])
						strand.append("-")
			else:
				if float(attributes[5]) >= pident and int(attributes[4]) > qcover:
					query_name.append(attributes[0])
					top_hit_name.append(attributes[1])
					percent_identity.append(attributes[5])
					query_coverage.append(attributes[4])
					if attributes[10].strip() == "plus":
						start.append(int(attributes[8]))
						stop.append(attributes[9])
						strand.append("+")
					else:
						start.append(int(attributes[9]))
						stop.append(attributes[8])
						strand.append("-")

	if not hmm:
	
		df1 = pd.DataFrame({"QUERY_NAME": query_name, "TOP_HIT_NAME": top_hit_name, "START": start, "STOP": stop, "STRAND": strand, "PERCENT_IDENTITY": percent_identity, "QUERY_COVERAGE": query_coverage})
		df2 = pd.DataFrame({"TOP_HIT_NAME": top_hit_name, "START": start, "STOP": stop, "STRAND": strand})
		df1.to_csv(f'{idt}_summary.txt', index=None, header=True, sep='\t')
		df2_minus = df2[df2['STRAND'] == '-']
		df2_plus = df2[df2['STRAND'] == '+']
		df2_minus.to_csv(f'{idt}_reverse.bed', index=None, header=False, sep='\t')
		df2_plus.to_csv(f'{idt}_forward.bed', index=None, header=False, sep='\t')

		file.close()

		return "A summary of your hits and forward/reverse BED files have been generated."

	else:

		df2 = pd.DataFrame({"TOP_HIT_NAME": top_hit_name, "START": start, "STOP": stop, "STRAND": strand})
		df2.to_csv(f'{idt}.bed', index=None, header=False, sep='\t')

		file.close()

		return "A BED file has been generated."



def reverse_complement(fasta_file, idt):

	output_handle = open(f'{idt}_reverse_complemented.fa', "w")
	for record in SeqIO.parse(fasta_file, 'fasta'):
		new_record = SeqRecord(seq = record.seq.reverse_complement(), \
				 id = record.id, \
				 description = "reverse complement")
		SeqIO.write(new_record, output_handle, "fasta")
	output_handle.close()

	return "The reverse hits have been reverse complemented"



def identifier_generator(filename):

	return str(filename)[:3]



def nhmmer_bed(filename, idt):
	
	file = open(filename, 'r')
	bed_lines = []
	for line in file.readlines():
		if not line.startswith("#"):
			attributes = line.split()
			if int(attributes[6]) < int(attributes[7]):
				bed_lines.append([attributes[0], attributes[6], attributes[7], '+'])
			elif int(attributes[6]) > int(attributes[7]):
				bed_lines.append([attributes[0], attributes[7], attributes[6], '-'])
	df = pd.DataFrame(bed_lines, columns=["CHROM", "START", "STOP", "STRAND"])
	df_minus = df[df['STRAND'] == '-']
	df_plus = df[df['STRAND'] == '+']
	df_minus.to_csv(f'{idt}_reverse_hmmer.bed', index=None, header=False, sep='\t')
	df_plus.to_csv(f'{idt}_forward_hmmer.bed', index=None, header=False, sep='\t')

	return "nhmmer output converted to BED file"

	

def extend_positions(df, idt, bp):

	df["START"] = df["START"] - int(bp)
	df["STOP"] = df["STOP"] + int(bp)
	df.to_csv(f"{idt}_extended.bed", index=False, sep="\t", header=False)

	return "positions have been extended"



def filter_tuples(lst):
	dump_lst = []
	new_lst = []
	for indx, tup in enumerate(lst):
		for indx1, tup1 in enumerate(lst):
			if indx != indx1:
				if tup[0] in range(tup1[0], tup1[1]) and tup[1] in range(tup1[0], tup1[1]):
					dump_lst.append(tup)
				elif tup[0] < tup1[0] and tup[1] in range(tup1[0], tup1[1]):
					if abs(tup[0] - tup[1]) < abs(tup1[0] - tup1[1]):
						dump_lst.append(tup)
					else:
						dump_lst.append(tup1)
				elif tup[0] < tup1[0] and tup[1] == tup1[1]:
					dump_lst.append(tup1)
				elif tup[0] == tup1[0] and tup[1] in range(tup1[0], tup1[1]):
					dump_lst.append(tup)
	for item in lst:
		if item not in dump_lst:
			new_lst.append(item)

	return new_lst
	


def filter_df(filename, idt):
	
	df = pd.read_csv(filename, sep='\t', names=['CHROM', 'START', 'STOP', 'STRAND'])

	dataframes = []

	for chrom in df['CHROM'].drop_duplicates():

		for direc in df['STRAND'].drop_duplicates():

			df_temp = df[(df['CHROM'] == chrom) & (df['STRAND'] == direc)]
			df_temp['COORDINATES'] = df_temp[['START', 'STOP']].apply(tuple, axis=1)
			coord_list = df_temp['COORDINATES'].tolist()
			updated_coord_list = filter_tuples(coord_list)
			no_dup_list = list(set([i for i in updated_coord_list]))
			df_temp_new = df_temp.loc[df_temp['COORDINATES'].isin(no_dup_list)]
			no_dup_df = df_temp_new.drop_duplicates()
			new_df = no_dup_df[['CHROM', 'START', 'STOP', 'STRAND']]
			dataframes.append(new_df)

	overall_df = pd.concat(dataframes)
	return overall_df

	
####### PIPELINE IMPLEMENTATION #######	

if __name__ == "__main__":

	if args.align == "True":

		if args.alignment_input == "0":
			os.system(f'clustalo -i {args.query} -o seqs.sto --outfmt=st --output-order=tree-order --force')
		else:
			os.system(f'clustalo -i {args.alignment_input} -o seqs.sto --outfmt=st --output-order=tree-order --force')

	for genome in os.listdir(args.pangenome):

		identifier = identifier_generator(genome)

		print("making BLAST database")
		os.system(f'makeblastdb -in {args.pangenome}{genome} -out {identifier}_db -parse_seqids -dbtype nucl')

		print("running nucleotide BLAST")
		os.system(f'blastn -db {identifier}_db -query {args.query} -task blastn -evalue {args.evalue} -outfmt "7 delim=, qacc sacc evalue bitscore qcovus pident qstart qend sstart send sstrand" > {identifier}_out')

		os.system(f'rm *_db.*')

			####### BLAST ONLY #######

		if args.nhmmer == "False":

			summary_message = generate_table(f'{identifier}_out', identifier, pident=args.percent_identity, qcover=args.query_coverage)
			print(summary_message)

			print("getting sequences of top hit set")
			os.system(f'bedtools getfasta -fi {args.pangenome}{genome} -bed {identifier}_forward.bed -s -fo {identifier}_forward_hits.fa')
			os.system(f'bedtools getfasta -fi {args.pangenome}{genome} -bed {identifier}_reverse.bed -s -fo {identifier}_reverse_hits.fa')

			os.system(f'rm {args.pangenome}*.fai')

			reverse_complement_message = reverse_complement(f'{identifier}_reverse_hits.fa', identifier)
			print(reverse_complement_message)

			os.system(f'mv {identifier}* {args.output_path}')

			####### HMMER ONLY #######

		else:

			summary_message = generate_table(f'{identifier}_out', identifier, pident=args.percent_identity, qcover=args.query_coverage, hmm=True)
			print(summary_message)

			filtered_df = filter_df(f'{identifier}.bed', identifier)

			extension_message = extend_positions(filtered_df, identifier, args.extension)
			print(extension_message)

			print("fetching fasta sequences for further analysis")
			os.system(f'bedtools getfasta -fi {args.pangenome}{genome} -bed {identifier}_extended.bed -s -fo {identifier}_filtered.fa')

			os.system(f'rm {args.pangenome}*.fai')

			print("running hmmer")
			if args.align == "False":
				os.system(f'nhmmer -E {args.nhmmer_evalue} --acc --tblout {identifier}_hmmer.out {args.stockholm} {identifier}_filtered.fa')
			else:
				os.system(f'nhmmer -E {args.nhmmer_evalue} --acc --tblout {identifier}_hmmer.out seqs.sto {identifier}_filtered.fa')

			print("parsing hmmer output")
			nhmmer_output_message = nhmmer_bed(f'{identifier}_hmmer.out', identifier)

			print("getting sequences from hmmer output")
			os.system(f'bedtools getfasta -fi {identifier}_filtered.fa -bed {identifier}_forward_hmmer.bed -s -fo {identifier}_forward_hmmer.fa')
			os.system(f'bedtools getfasta -fi {identifier}_filtered.fa -bed {identifier}_reverse_hmmer.bed -s -fo {identifier}_reverse_hmmer.fa')

			os.system(f'rm *.fai')

			reverse_complement_message = reverse_complement(f'{identifier}_reverse_hmmer.fa', identifier)
			print(reverse_complement_message)
				
			os.system(f'mv {identifier}* {args.output_path}')