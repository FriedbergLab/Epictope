from epitope_helper import *
import pandas
import requests
import re

def query_to_msa():
    tf = pandas.read_csv("../data/zebra_transcription_factors.csv")
    
    #fetch transcripts from ensembl API
    ts = list(map(fetch_transcripts, tf['ensemblID']))
    
    #select longest sequence
    max_ts = list(map(get_max_str, ts)) 

    #convert string to SeqIO fasta object
    fastas = list(map(string_to_fasta, max_ts))

    #add gene names to fasta object
    blast_queries0 = name_fastas(fastas, tf["geneName"])

    #filter out queries where no transcript was returned
    blast_queries = [ele for ele in blast_queries0 if str(ele.seq)!=""]

    #Write blast queries to file
    SeqIO.write(blast_queries, "blastqueries.fasta", "fasta")

    #blast against folder of databases
    multi_blast("../data/blastdbs", "blastqueries.fasta", evalue = 1e-3, outfmt = 5)

    #parse blast xml output and fetch full transcripts for hits
    output_list = []
    for file in list_files("../outputs", ".xml"):
        print("parsing" + file)
        seqs = list(map(fetch_transcripts, xml_parse(file)))
        max_seqs = list(map(get_max_str, as_list(seqs)))
        fastas = list(map(string_to_fasta, max_seqs))
        output_list.append(fastas)

    #combine matching transcripts
    matches = tuple(zip(blast_queries, *output_list))

    descriptions = ["Denio_renio", "Xenopus_tropicalis", "Mus_msculus", "Gallus_gallus", "Homo_sapiens", "Takifugu_rubripes"]
    desc_matches = desc_fastas(matches, descriptions)
    #write matching transcripts to individual fastas
    #using protein sequence
    to_fasta(desc_matches, protein = True)

    #MSA with muscle.
    for fasta in list_files("../outputs", ".fasta"):
        muscle_cline = MuscleCommandline(
        input=fasta, 
        out=fasta.split(".fasta")[0] + ".txt")
        stdout, stderr = muscle_cline()

###
if __name__ == '__main__':
    query_to_msa()

