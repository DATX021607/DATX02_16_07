#!/usr/bin/env python

from sys import *
from Bio.PDB import *
	
#Send dpb file as input together with the program
#
#
#Readig the input pdb code
pdbfile=str(argv[1]).lower();
#Creates path to download for later scan
path_after_dl='Unparsed-pdb/pdb'+pdbfile+'.ent';
#Empty sequence array
sequence = [];
#Create a downloadlist and a parser
pdbl    = PDBList();
parser  = PDBParser();
#Downloads the pdb file
pdbl.retrieve_pdb_file(pdbfile, obsolete=False, pdir='Unparsed-pdb'); 

#Parse the data from the downloaded pdb
structure = parser.get_structure(pdbfile, path_after_dl)
#Creates a PollyPeptideBuilder
ppb=PPBuilder()
#Build the peptide
ppb.build_peptides(structure) 
#Loop through the object and create a sequence list
for pp in ppb.build_peptides(structure): 
	pp.get_sequence();

#Store the sequence lsit
seq = pp.get_sequence()
#Open a file, and write the sequence to it, then close
f = open("Parsed-pdb/"+pdbfile, 'w')
f.write(str(seq))
f.close

