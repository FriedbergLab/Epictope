import pandas
import re
import os
import selenium
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastpCommandline

# read in transcription factors from file
filename = "../data/tfs.txt"
tf = pandas.read_csv(filename)


# pull uniprot ID from alphafold search

def render(url, driver_loc = '/usr/bin/chromedriver'):
	s = Service(driver_loc)
	opts = Options()
	opts.add_argument("--disable-extensions")
	opts.add_argument("--disable-gpu")
	opts.add_argument("--no-sandbox") # linux only
	opts.add_argument("--headless")
	browser = webdriver.Chrome(service=s, options = opts)
	browser.get(url)
	soup = BeautifulSoup(browser.page_source, features = "html.parser")
	browser.quit()
	return(soup)

def fetch_uniprot(query, organism, server = "https://alphafold.ebi.ac.uk/search/text/"):
	if organism != "" :
		organism_ext = "?organismScientificName=" + re.sub(" ", "%20", organism)
	else:
		organism_ext = ""
	alphafold_url = server + query + organism_ext
	res = render(alphafold_url)
	return(res)

uniprot_pages = [fetch_uniprot(p, "Danio rerio") for p in tf['geneName']]

def search_link(soup, link = "http://www.uniprot.org/uniprot/.*"):
	text = re.search(link, soup.prettify()).group()
	url = re.sub(" target=_blank>", "", text.replace("\"", "")) + ".fasta"
	return(url)

def fetch_text(url):
	res = requests.get(url).text
	return(res)

uniprot_links = list(map(search_link, uniprot_pages))

fastas = list(map(fetch_text, uniprot_links))


named_fastas = dict(zip(tf['geneName'], fastas))

seq_list = []
for key, value in named_fastas.items():
	seq = re.search("[A-Z]{10,}", value.replace("\n", "")).group()
	rec = SeqRecord(Seq(seq), 
                id = key, 
                name = key, 
                description = "")
	seq_list.append(rec)

SeqIO.write(seq_list, "blast_queries.fasta", "fasta")

def list_files(folder, filetype):
    folder_list = []
    for dirpath, dirnames, filenames in os.walk(folder):
        for filename in [f for f in filenames if f.endswith(filetype)]:
            folder_list.append(os.path.join(dirpath, filename))
    return folder_list

def multi_blast(db_path, fasta_file, evalue = 1e-3, outfmt = 5):
    db_list = list_files(db_path, ".fa")
    # Loop to blast query fasta against any databbases in db folder.
    for db in db_list:
        blastp_cline = NcbiblastpCommandline(
            query = fasta_file, 
            db = db, 
            evalue = evalue, 
            outfmt = outfmt, 
            out = "../outputs/" + db.split("/")[3] + "_blast.xml",)
        stdout, stderr  = blastp_cline()


multi_blast("../data/blastdb", "blast_queries.fasta", evalue = 1e-3, outfmt = 5)


def xml_parse(xmlpath):
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

def get_max_str(lst):
    return max(lst, key=len) if lst else ""

def as_list(arg):
    return arg if isinstance(arg, list) else [arg]

def string_to_fasta(string, desc = "", name = ""):
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


def string_strip(string):
	res = string.split("/")[2].split("_blast")[0]
	return(res)

descriptions = list(map(string_strip, list_files("../outputs", ".xml")))
descriptions.insert(0, "Danio_rerio")

def desc_fastas(fastas, desc):
    desc_fastas = []
    for fasta in fastas:
        i = 0
        for entry in fasta:
                entry.description = desc[i]
                i += 1
        desc_fastas.append(fasta)
    return(desc_fastas)   

desc_matches = desc_fastas(matches, descriptions)
#write matching transcripts to individual fastas.
#using protein sequence
def to_fasta(seq_list, name_dict = ""):
   # Write matching sequences to fastas for msa.
    for match in seq_list:
        filtered_match = list(filter(None, match))
        if len(filtered_match) > 1:
                SeqIO.write(filtered_match, "../outputs/" + match[0].name + ".fasta", "fasta")
to_fasta(desc_matches)

#MSA with muscle.
for fasta in list_files("../outputs", ".fasta"):
    muscle_cline = MuscleCommandline(
    input=fasta, 
    clw = True, 
    out=fasta.split(".fasta")[0] + ".clw")
    stdout, stderr = muscle_cline()


