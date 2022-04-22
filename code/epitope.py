#!/usr/bin/python3
from epitope_helper import *
import pandas
import requests
import re
import argparse
import sys
import csv

def query_to_alphafold(filename, species):

    #Check if output folder exists.
    does_folder_exist("outputs")

    #read in input file.
    tf = pandas.read_csv(filename)
    
    # fetch uniprot sequence from AlphaFold website
    uniprot_pages = [fetch_uniprot(p, species) for p in tf

    #fetch transcripts from ensembl API.
    ts = list(map(fetch_transcripts, tf['ensemblID']))
    
    #select longest sequence.
    max_ts = list(map(get_max_str, ts)) 

    #convert string to SeqIO fasta object.
    fastas = list(map(string_to_fasta, max_ts))

    #add gene names to fasta object.
    blast_queries0 = name_fastas(fastas, tf["geneName"])

    #filter out queries where no transcript was returned.
    blast_queries = [ele for ele in blast_queries0 if str(ele.seq)!=""]

    #Write blast queries to file.
    SeqIO.write(blast_queries, "blastqueries.fasta", "fasta")

    #blast against folder of databases.
    multi_blast("../data/blastdbs", "blastqueries.fasta", evalue = 1e-3, outfmt = 5)


    #parse blast xml output and fetch full transcripts for hits.
    output_list = []
    for file in list_files("outputs", ".xml"):
        print("parsing " + file)
        seqs = list(map(fetch_transcripts, xml_parse(file)))
        max_seqs = list(map(get_max_str, as_list(seqs)))
        fastas = list(map(string_to_fasta, max_seqs))
        output_list.append(fastas)

    #combine matching transcripts.
    matches = tuple(zip(blast_queries, *output_list))

    descriptions = "Danio_renio", "Xenopus_tropicalis", "Mus_musculus", "Gallus_gallus", "Homo_sapiens", "Takifugu_rubripes"]
    desc_matches = desc_fastas(matches, descriptions)
    #write matching transcripts to individual fastas.
    #using protein sequence
    to_fasta(desc_matches, protein = True)

    #MSA with muscle.
    for fasta in list_files("outputs", ".fasta"):
        muscle_cline = MuscleCommandline(
        input=fasta, 
        clw = True, 
        out=fasta.split(".fasta")[0] + ".clw")
        stdout, stderr = muscle_cline()

   	#parse muscle txt output and calculate shannon entropy for each column.
   	#write result to csv file.
    with open('outputs/entropies.csv', 'w', encoding='UTF8', newline='') as f:
    	writer = csv.writer(f)
    	for file in list_files("outputs", ".clw"):
        	#print("calculating entropy of" + file)
        	res = shannon_entropy_list_msa(AlignIO.read(file, "clustal"))
        	res.insert(0, file.split("/")[1].split(".")[0])
        	writer.writerow(res)
    f.close()

    print("No results returned for:")
    for i in neg_indices:
        print(tf.geneName[i])
        	


###
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= 
    'Takes an input file of genes and associated Ensembl IDs, retrieves \
    associated sequences from the Ensembl REST API, BLASTs the sequences against \
    local blasteable databases, and then performs a MSA with the homologous sequences.')

    parser.add_argument("-i", "--input", dest="filename", required=True,
                    help="single column text file with gene name queries", metavar="txt file",
                    type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-c", "--config", dest="config", required = True, 
                    help = "path to config.cfg file with parameter values", metavar="cfg file",
                    default = "code/config.yml",
                    type=lambda x: is_valid_file(parser, x))

    args = parser.parse_args() 

    with open(args.config, 'r') as file:
        params = yaml.safe_load(file)

    query_to_msa(args.filename, species)

