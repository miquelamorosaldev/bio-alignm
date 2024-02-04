from Bio import Entrez, SeqIO
import os

Entrez.email = "david@xtec.dev"

def print_matrix(name):
    ""
    

def sequence(id):
    
    db = "nuccore"
    if id.startswith("NP"):
        db = "protein"

    dir = "data/"
    if not os.path.isdir(dir):
        os.makedirs(dir)

    file = "{}/{}.fasta".format(dir, id)
    if not os.path.exists(file):
        print(f"Downloading file {id}...")
        with Entrez.efetch(db=db, id=id, rettype="fasta") as response:
            with open(file, "w") as out:
                out.write(response.read())

    return SeqIO.read(file, "fasta")


