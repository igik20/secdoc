from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import subprocess

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

    def threestate(seq, helix="GHI", sheet="BE"):
        """
        Convert a DSSP secondary structure string into three-state format (helix/sheet/coil).
        """
        out = ""
        for c in seq:
            if c in helix:
                out.append('H')
            elif c in sheet:
                out.append('E')
            else:
                out.append('C')
        return out
