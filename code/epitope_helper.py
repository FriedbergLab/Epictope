#!/usr/bin/python3
import pandas
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import configparser
import re
import requests
import os
import os.path
import selenium
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import requests
import time
import yaml

def uniprot_mapping(identifier, fromtype = "GENENAME", totype = "SWISSPROT", taxa_id = 9606):
    """Takes an identifier, and types of identifier 
    (to and from), and calls the UniProt mapping service"""
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from':fromtype,
                'to':totype,
                'format':'tab',
                'query':identifier,
                'taxon': taxa_id,
    }
    #urllib turns the dictionary params into an encoded url suffix
    data = urllib.parse.urlencode(params)
    #construct the UniProt URL
    url = base+'/'+tool+'?'+data
    #and grab the mapping
    response = urllib.request.urlopen(url)
    #response.read() provides tab-delimited output of the mapping
    string = str(response.read())
    #clean string
    clean_string = re.sub(r'\\t|\\n|b\'|\'','_', string).split("_")
    #remove empty elements at the beginning
    full_string = list(filter(None, clean_string))
    #lengthen string to 4 elements
    while(len(full_string) < 4):
        full_string.append(None)
    return(full_string)

def uniprot_to_df(sp_list, queries):
    '''Converts string output from swissprot id mapping to 
    Pandas Dataframe'''
    sp_df = pandas.DataFrame(columns = ['id_no', "geneName"])
    sp_df = pandas.DataFrame(sp_list).iloc[:, [3]]
    sp_df.columns = ["id_no"]
    sp_df["geneName"] = list(queries)
    return(sp_df)

def fetch_content(url):
    """retrives content
    from URL
    """
    res = requests.get(url, allow_redirects = True)
    return(res.content)

def fetch_content(url):
    """retrives content
    from URL
    """
    res = requests.get(url, allow_redirects = True)
    return(res.content)

def search_link(soup, link):
    """Searches a beautifulsoup object
    for a link and returns matches.
    """
    all_links = ""
    for a in soup.find_all('a', href = True):
        all_links += a['href']+" "
    search = re.search(link, all_links)
    if search is None: return("No result found")
    return(search.group())

def delete_multiple_element(list_object, indices):
    """Deletes elements of a list given a list and
    a list of elements to delete.
    """
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)

def does_folder_exist(folder, create = False):
    """ Checks if folder exists,
    optional argument to create file if 
    desired.
    """
    if os.path.exists(folder) == False:
        print(folder+" does not exist")
        if create == True:
            os.mkdir(folder)
            print("Created " + folder + " folder.")
    else:
        print(folder+" exists.")


def is_valid_file(parser, arg):
    """checks if file path exists.
    """
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

def list_files(folder, filetype):
    """lists all of the files
    of a specific filetype in a 
    given folder.
    """
    folder_list = []
    for dirpath, dirnames, filenames in os.walk(folder):
        for filename in [f for f in filenames if f.endswith(filetype)]:
            folder_list.append(os.path.join(dirpath, filename))
    return folder_list


def multi_blast(db_path, fasta_file, evalue = 1e-3, outfmt = 5):
    """performs a blast search for each
    sequence listed in a fasta file against a
    database of organisms.
    """
    db_list = list_files(db_path, ".fa")
    # Loop to blast query fasta against any databases in db folder.
    for db in db_list:
        blastp_cline = NcbiblastpCommandline(
            query = fasta_file, 
            db = db, 
            evalue = evalue, 
            outfmt = outfmt, 
            out = "outputs/" + db.split("/")[3] + "_blast.xml",)
        stdout, stderr  = blastp_cline()

def xml_parse(xmlpath):
    """Parses the xml output from
    NCBI blast output.
    """
    # Initate ensembleRest object and lists
    id_list = []
    # Parse records from xml output.
    for record in NCBIXML.parse(open(xmlpath)): 
        try:
            id_temp = record.alignments[0].title.split(" ")[0].split(".")[0]
        # If not hit from blast query, will return IndexError.  
        # Append "" if no hit to maintain index order of queries.
        except IndexError:
            id_temp = ""
        id_list.append(id_temp)
    return(id_list)


def fetch_transcripts(query, 
    server = "https://rest.ensembl.org", 
    ext = "/sequence/id/", 
    seq_type = "cds",
    content_type = "fasta",
    multi_seq = "1"
):
    """fetch the coding sequence
    for a protein from ensembl
    """
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


def get_max_str(lst):
    """determine which string 
    is longest from a list of strings
    """
    return max(lst, key=len) if lst else ""

def as_list(arg):
    """convert an object 
    to a list.
    """
    return arg if isinstance(arg, list) else [arg]

def string_to_fasta(string, desc = "", name = ""):
    """converts a string object
    to a Bio.SeqRecord object
    """
    if string == "":
        rec = SeqRecord(Seq(""), id = "", name = name, description = desc)
    else:  
        temp_string = string.replace("\n", "")
        seq_id = re.search("ENS[A-Z]+[0-9]+", temp_string).group()
        seq = re.search("[A-Z]{10,}", temp_string).group()
        rec = SeqRecord(Seq(seq), 
                id = seq_id, 
                name = name, 
                description = desc)
    return(rec)

def string_strip(string):
    """strip the folder argument from a string.
    """
    res = string.split("/")[1].split("_blast")[0]
    return(res)

def desc_fastas(fastas, desc):
    """adds a description attribute to 
    fasta SeqRecord objects
    """
    desc_fastas = []
    for fasta in fastas:
        i = 0
        for entry in fasta:
                entry.description = desc[i]
                i += 1
        desc_fastas.append(fasta)
    return(desc_fastas)   


def to_fasta(seq_list, params, protein = False, name_dict = ""):
   # Write matching sequences to fastas for msa.
    for match in seq_list:
        filtered_match = list(filter(None, match))
        if protein == True:
            SeqIO.write(filtered_match, 
                params["output_file"] + match[0].id + ".fasta", "fasta")
        else:
            for newseq in filtered_match: 
                temp = newseq.translate()
                newseq.seq = temp.seq
            SeqIO.write(filtered_match,
                params["output_file"] + match[0].id + ".fasta", "fasta")


def shannon_entropy(list_input):
    """Function to calcuate the Shannon's entropy per alignment column
    """
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


def shannon_entropy_list_msa(alignment_file):
    """Function to calculate shannon entropy for all columns in a MSA.
    """
    shannon_entropy_list = []
    for col_no in range(len(list(alignment_file[0]))):
        list_input = list(alignment_file[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))
    return shannon_entropy_list
