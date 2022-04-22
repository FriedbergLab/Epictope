#!/usr/bin/python3
from epitope_helper import *
import pandas
import requests
import re
import argparse
import sys
import csv

def query_to_alphafold(filename, species):

    ## read in input file
    tf = pandas.read_csv(filename)

    # search alphafold for protein entry from species
    print("Searching Alphafold...")
    alphafold_pages = [search_alphafold(p, params["species"]) for p in tf["geneName"]]

    # search for links to uniprot entry
    uniprot_links = list(map(lambda x: search_link(x, link = "http://www.uniprot.org/uniprot/.[^\s\"]*")+".fasta", alphafold_pages))

    # if no results found, remove from queries
    neg_indices = [idx for idx, s in enumerate(uniprot_links) if 'No result' in s]
    delete_multiple_element(alphafold_pages, neg_indices)
    delete_multiple_element(uniprot_links, neg_indices)

    # search for links for alphafold entry
    entry_links = list(map(lambda x: re.sub("\.", "https://alphafold.ebi.ac.uk/", search_link(x, link = "./entry/.[^\s\"]*")), alphafold_pages))

    # search alphafold page for pdb links
    entry_pages = [render(x) for x in entry_links]
    pdb_links = list(map(lambda x: search_link(x, link = "https.[^\s]*.pdb"), entry_pages))


    # download sequence from uniprot fasta
    print("Downloading fasta files...")
    fastas = list(map(fetch_content, uniprot_links))
    named_fastas = dict(zip(tf['geneName'], fastas))

    # write sequences to file
    seq_list = []
    for key, value in named_fastas.items():
        value = value.decode("utf-8")
        seq = re.search("[A-Z]{10,}", value.replace("\n", "")).group()
        rec = SeqRecord(Seq(seq), 
                id = key, 
                name = key, 
                description = "")
        seq_list.append(rec)

    print("Writing fasta sequences to file...")
    SeqIO.write(seq_list, params["data_file"]+"blast_queries.fasta", "fasta")

    # download pdb files for queried proteins
    does_folder_exist(params["data_file"]+"pdb", create = True)
    print("Downloading pdb entries...")
    pdb = list(map(fetch_content, pdb_links))
    named_pdb = dict(zip(tf['geneName'], pdb))

    for key, value in named_pdb.items():
        f = open(params["data_file"]+"pdb/"+key+".pdb", 'wb')
        f.write(value)
        f.close()

    # blast sequences
    print("BLASTing for similar proteins...")
    multi_blast(params["data_file"]+"blastdb", params["data_file"]+"blast_queries.fasta", evalue = 1e-3, outfmt = 5)

    #parse blast xml output and fetch full transcripts for hits.
    output_list = []
    for file in list_files("../outputs", ".xml"):
        print("parsing " + file)
        seqs = list(map(fetch_transcripts, xml_parse(file)))
        max_seqs = list(map(get_max_str, as_list(seqs)))
        fastas = list(map(string_to_fasta, max_seqs))
        output_list.append(fastas)

    #combine matching transcripts.
    matches = tuple(zip(seq_list, *output_list))

    # add animal descriptions
    descriptions = list(map(string_strip, list_files("../outputs", ".xml")))
    descriptions.insert(0, params["species"])
    desc_matches = desc_fastas(matches, descriptions)

    # write protein sequences for multiple species to file
    to_fasta(desc_matches, protein = True)

    #MSA protein sequences with muscle.
    print("Performing MSA with muscle...")
    for fasta in list_files("../outputs", ".fasta"):
        muscle_cline = MuscleCommandline(
        input=fasta, 
        clw = True, 
        out=fasta.split(".fasta")[0] + ".clw")
        stdout, stderr = muscle_cline()

    # print gene names where no results.
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
    parser.add_argument("-c", "--config", dest="config", required = False, 
                    help = "path to config.cfg file with parameter values", metavar="cfg file",
                    default = "code/config.yaml",
                    type=lambda x: is_valid_file(parser, x))

    args = parser.parse_args() 

    params = yaml.safe_load(args.config)

    query_to_alphafold(args.filename, params["species"])
    print("Operation complete.")



