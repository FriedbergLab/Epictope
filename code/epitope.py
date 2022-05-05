#!/usr/bin/python3
from epitope_helper import *
import pandas
import requests
import re
import argparse
import sys
import csv

def query_to_alphafold(filename, params):

    tf = pandas.read_csv(filename)

    if(len(tf) < 1):
        print("No gene names found.")
        sys.exit()

    # convert gene names to uniprot ID
    swissprot_list = [uniprot_mapping(p, taxa_id = params["species"]) for p in tf["geneName"]]

    swissprot_df = uniprot_to_df(swissprot_list, queries = tf["geneName"])

    pos = pandas.DataFrame(columns = ['id_no', "geneName"])
    pos = pandas.concat([pos, swissprot_df[swissprot_df['id_no'].notnull()]])

    neg = pandas.DataFrame(columns = ['id_no', "geneName"])
    neg = pandas.concat([neg, swissprot_df[swissprot_df['id_no'].isnull()]])

    if(len(neg) > 0):
        alt_list = [uniprot_mapping(p, totype = "AC", taxa_id = params["species"]) for p in neg["geneName"]]
        alt_df = uniprot_to_df(alt_list, queries = neg["geneName"])
        pos = pandas.concat([pos, alt_df[alt_df["id_no"].notnull()]])

    fastas = [fetch_content("http://www.uniprot.org/uniprot/"+i+".fasta") for i in pos["id_no"]]
    named_fastas = dict(zip(tf['geneName'], fastas))

    # write sequences to file
    seq_list = []
    for key, value in named_fastas.items():
        value = value.decode("utf-8")
        start = re.search("GN", value).start()
        result = value[start:].split("=")[1].strip()
        seq = re.search("[A-Z]{10,}", value.replace("\n", "")).group()
        rec = SeqRecord(Seq(seq), 
                id = key, 
                name = result, 
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
    for file in list_files(params["output_file"], ".xml"):
        print("parsing " + file)
        seqs = list(map(fetch_transcripts, xml_parse(file)))
        max_seqs = list(map(get_max_str, as_list(seqs)))
        fastas = list(map(string_to_fasta, max_seqs))
        output_list.append(fastas)

    #combine matching transcripts.
    matches = tuple(zip(seq_list, *output_list))

    # add animal descriptions
    descriptions = list(map(string_strip, list_files(params["output_file"], ".xml")))
    descriptions.insert(0, params["species"])
    desc_matches = desc_fastas(matches, descriptions)

    # write protein sequences for multiple species to file
    to_fasta(desc_matches, protein = True, params = params)

    #MSA protein sequences with muscle.
    print("Performing MSA with muscle...")
    for fasta in list_files(params["output_file"], ".fasta"):
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

    parameters = yaml.safe_load(args.config)

    query_to_alphafold(args.filename, params = parameters)
    print("Operation complete.")



