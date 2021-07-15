import pandas
import Bio
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ensemblrest import EnsemblRest
import re
import requests
import os
import os.path


#Fetch sequences using ensembl API.  
#Inputs:
#server - server url.
#ex - resource extension.
#seq_type - type of sequence.  
#content_type - return data format.  
#multi_seq - boolean for return multiple sequences.

#Outputs: list object of fetched sequences as string.
def fetch_transcripts(query, 
    server = "https://rest.ensembl.org", 
    ext = "/sequence/id/", 
    seq_type = "cds",
    content_type = "fasta",
    multi_seq = "1"
):
    url = server + ext + query + \
        "?multiple_sequences=" + multi_seq + \
        ";type=" + seq_type + \
        ";content_type=" + content_type
    r = requests.get(url, allow_redirects = True, headers = {'Content-type': 'text/x-fasta'})
    
    if re.search("error", r.text):
        print("No sequence found for " + query)
        transcripts = ""
    else:
        transcripts = str(r.text).split(">")[1:]
    return(transcripts)


#function to return max element length from a list.  
#Input: list.
#Output: longest element from input list.
def get_max_str(lst):
    return max(lst, key=len) if lst else ""

#function to add name attribute to list of SeqIO objects.  
#Input: 
#fastas - list of SeqIO objects.  
#names - list of names (string).  
#Output: list of SeqIO objects with edited name attribute.
def name_fastas(fastas, names):
    i = 0
    named_fastas = []
    for fasta in fastas:
        fasta.name = names[i]
        i += 1
        named_fastas.append(fasta)
    return(named_fastas)    

def desc_fastas(fastas, desc):
    desc_fastas = []
    for fasta in fastas:
        i = 0
        for entry in fasta:
                entry.description = desc[i]
                i += 1
        desc_fastas.append(fasta)
    return(desc_fastas)    



#function that identifies all files of specific type in a folder
#Inputs:
#folder - folder to search through.  
#filetype - type of file to search for.  
#Outputs: list of filenames of specified type in folder.
def list_files(folder, filetype):
    folder_list = []
    for dirpath, dirnames, filenames in os.walk(folder):
        for filename in [f for f in filenames if f.endswith(filetype)]:
            folder_list.append(os.path.join(dirpath, filename))
    return folder_list


#function to perform blast against local databases
#Inputs:
#db_path - path to local reference databases.
#fasta_file - fasta file of blast queries.
#evalue - tuning parameter for stringency of blast search.
#outfmt - blast output file format (5 = xml)
#Outputs - xml files with blast results
def multi_blast(db_path, fasta_file, evalue = 1e-3, outfmt = 5):
    db_list = list_files(db_path, ".fa")
    # Loop to blast query fasta against any databbases in db folder.
    for db in db_list:
        blastn_cline = NcbiblastnCommandline(
            query = fasta_file, 
            db = db, 
            evalue = evalue, 
            outfmt = outfmt, 
            word_size = 11,
            gapopen = 2,
            gapextend = 2,
            reward = 1, 
            penalty = -3, 
            out = "outputs/" + db.split("/")[3] + "_blast.xml",)
        stdout, stderr  = blastn_cline()

#function to return inputs as objects as single element lists.
#Input: object
#Output: single element list
def as_list(arg):
    return arg if isinstance(arg, list) else [arg]

#function to parse blast xml output files
#Input: path to xml file.
#Output: list of first hits from blast query.
def xml_parse(xmlpath):
    # Initate ensembleRest object and lists
    id_list = []
    # Parse records from xml output.
    for record in NCBIXML.parse(open(xmlpath)): 
        try:
            id_temp = record.alignments[0].title.split(" ")[1].split(".")[0]
        # If not hit from blast query, will return IndexError.  
        # Append "" if no hit to maintain index order of queries.
        except IndexError:
            id_temp = ""
        id_list.append(id_temp)
    return(id_list)

#function to write list of SeqIO objects to fasta file with
#option to convert nucleotide sequence to protein sequence. 
#Input
#Output
def to_fasta(seq_list, protein = False, name_dict = ""):
   # Write matching sequences to fastas for msa.
    for match in seq_list:
        filtered_match = list(filter(None, match))
        if len(filtered_match) > 1:
            if protein == False:
                SeqIO.write(filtered_match, 
                    "outputs/" + match[0].name + ".fasta", "fasta")
            else:
                for newseq in filtered_match: 
                    temp = newseq.translate()
                    newseq.seq = temp.seq
                SeqIO.write(filtered_match,
                    "outputs/" + match[0].name + ".fasta", "fasta")

#function to convert a character string consisting of fasta information into a 
#Biopython SeqIO object.
#Input: character string consisting of fasta ID and sequence.
#Output Biopython SeqIO object.
def string_to_fasta(string, desc = "", name = ""):
    if string == "":
        rec = SeqRecord(Seq(""), id = "", name = name, description = desc)
    else:  
        temp_string = string.replace("\n", "")
        seq_id = re.search("ENS[A-Z]+[0-9]+", temp_string).group()
        seq = re.search("[A,C,T,G]{3,}", temp_string).group()
        rec = SeqRecord(Seq(seq), 
                id = seq_id, 
                name = name, 
                description = desc)
    return(rec)

#function to determine if an input file exists.
#Input:  a filename
#Output: opens file if exists, reports if file does not exist. 
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

#function to determine if output folder exists, if it does not exist
#create the output folder.
#Input: folder name
#Output: No output if folder exists, creates folder if does not exist.
def does_folder_exist(folder):
    if os.path.exists(folder) == False:
        os.mkdir(folder)
        print("Created " + folder + " folder.")
###

