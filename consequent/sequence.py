"""
Reading Uniprot sequences by CID code.

"""
from requests import post
from Bio import SeqIO
from io import StringIO
from sys import exit
from os import path

def getUniprotSeq(cID):

    if path.exists(cID+".fasta") :
        print("Reading file '{}' from disk".format(cID+".fasta"))
        pSeq = list(SeqIO.parse(cID+".fasta", "fasta"))
    else :
        print("Fetching file '{}' from Uniprot".format(cID+".fasta"))
        baseUrl="http://www.uniprot.org/uniprot/"
        currentUrl=baseUrl+cID+".fasta"
        response = post(currentUrl)
        cData=''.join(response.text)
        Seq=StringIO(cData)
        pSeq=list(SeqIO.parse(Seq,'fasta'))
        if len(pSeq) == 0 :
            print("Error getting sequence '{}'.".format(cID))
        
        # Write the file on local directory to save further URL access
        SeqIO.write(pSeq[0],cID+".fasta","fasta")
    return pSeq[0]

