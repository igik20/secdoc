from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import subprocess
from itertools import groupby
import sys

class SecStruc():
    def __init__(self, pdb_path, *, dssp_path="dssp", chain="A"):
        """
        Create a DSSP object from an input PDB and extract some useful variables.
        """
        # get version of the DSSP executable
        dssp_ver = subprocess.run(["dssp", "--version"], capture_output=True, text=True).stdout.strip().split()[2]
        
        # create a Biopython DSSP dict object
        self.dssp = dssp_dict_from_pdb_file(pdb_path, DSSP=dssp_path, dssp_version=dssp_ver)[0]
        
        # filter the entries in the dict and convert to a nicer format
        self.dssp = dict([(k[1][1], v[1]) for (k,v) in self.dssp.items() if k[0] == chain])
        
        # extract the secondary structure into a string
        self.secseq = ''.join(self.dssp.values())
    
    def threestate(self, *, helix="GHI", sheet="BE"):
        """
        Convert a DSSP secondary structure string into three-state format (helix/sheet/coil).
        DSSP states considered helices and sheets can be changed with the helix and sheet arguments.
        Pure - returns the threestate string.
        """
        out = ""
        for c in self.secseq:
            if c in helix:
                out += 'H'
            elif c in sheet:
                out += 'E'
            else:
                out += 'C'
        return out

    def to_threestate(self, *, helix="GHI", sheet="BE"):
        """
        Convert the structure string to three-state in-place. All arguments as in the pure function.
        """
        self.secseq = self.threestate(helix=helix, sheet=sheet)

    def runlength(self):
        """
        Convert the structure string to run-length encoding. Pure.
        """
        out = ""
        counts = ((x, sum(1 for _ in y)) for (x,y) in groupby(self.secseq))
        for (x,y) in counts:
            out += x
            out += str(y)
        return out

    def featurelist(self):
        """
        Create a list of pairs with secondary structures and their lengths. Pure.
        """
        return list((x, sum(1 for _ in y)) for (x,y) in groupby(self.secseq))

    def featurereport(self, path=1):
        """
        Pretty-print a report of secondary structure to stdout or a file.
        """
        features = self.featurelist()
        with open(path, 'w') as f:
            for feature in features:
                if feature[0] == 'H':
                    f.write("HELIX " + str(feature[1]) + '\n')
                elif feature[0] == 'E':
                    f.write("SHEET " + str(feature[1]) + '\n')
                elif feature[0] == 'C':
                    f.write("COIL  " + str(feature[1]) + '\n')



