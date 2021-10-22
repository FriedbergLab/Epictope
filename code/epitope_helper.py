#!/usr/bin/python3
import pandas
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from ensemblrest import EnsemblRest
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
def to_fasta(seq_list, name_dict = ""):
   # Write matching sequences to fastas for msa.
    for match in seq_list:
        filtered_match = list(filter(None, match))
        if len(filtered_match) > 1:
            if match[0].name == "<unknown id>":
                SeqIO.write(filtered_match, 
                    "../outputs/fastas/" + "unknown" + ".fasta", "fasta")
            else:
                SeqIO.write(filtered_match,
                    "../outputs/fastas/" + match[0].name + ".fasta", "fasta")

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


# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i 
# Gaps and N's are included in the calculation
def shannon_entropy(list_input):
    import math
    unique_base = set(list_input)                           # Get only the unique bases in a column
    #unique_base = unique_base.discard("-")
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base)                        # Number of residues of type i                   
        P_i = n_i/float(M)                                  # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    return sh_entropy

#function to calculate shannon entropy for all columns in a MSA.
def shannon_entropy_list_msa(alignment_file):
    shannon_entropy_list = []
    for col_no in range(len(list(alignment_file[0]))):
        list_input = list(alignment_file[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))
    return shannon_entropy_list


def read_alpha(file):
    data = pandas.read_csv(file, sep = "\t")
    data.columns = ["chain_id", "residue_no", "residue", "secondary_struct", "surface_area"]
    return(data)


def alpha_seq(named_data):
    seq = "".join(list(named_data[0]["residue"]))
    rec = SeqRecord(Seq(seq), id = "AlphaFold Sequence", name = named_data[1], description = "")
    return(rec)


def seq_translate(SeqRecord):
    SeqRecord.seq = SeqRecord.translate().seq
    return(SeqRecord)


def add_alphas(matches, alpha_seqs):
    res = list()
    for match in matches:
        for alpha in alpha_seqs:
            if match[0].name + ".txt" == alpha.name :
                temp = list(map(seq_translate, list(match)))
                temp.append(alpha)
                res.append(temp)
    return(res)


def fetch_alphas(alphafold_folder):
    alphafold_files = list_files(alphafold_folder, "txt")
    alphafold_data = list(map(read_alpha, alphafold_files))
    names = os.listdir(alphafold_folder)
    named_data = list(zip(alphafold_data, names))
    alphafold_seqs = list(map(alpha_seq, named_data))
    return(alphafold_seqs)
