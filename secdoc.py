from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import subprocess
from itertools import groupby

class SecFeat():
    def __init__(self, _type, _start, _end):
        self.type = _type
        self.start = _start
        self.end = _end
        self.length = self.end - self.start + 1

    def __len__(self):
        return self.length
    
    def __str__(self):
        labels = {'H': 'HELIX', 'E': 'SHEET', 'C': 'COIL'}
        return f"{labels[self.type]}\t{self.start}\t{self.end}\t{len(self)}"

class SecStruc():
    def __init__(self, pdb_path, *, dssp_path="dssp", chain="A"):
        """
        Create a DSSP object from an input PDB and extract some useful variables.
        """
        dssp_ver = subprocess.run(["dssp", "--version"], capture_output=True, text=True).stdout.strip().split()[2]
        
        self.dssp = dssp_dict_from_pdb_file(pdb_path, DSSP=dssp_path, dssp_version=dssp_ver)[0]
        
        # filter the entries in the dict and convert to a nicer format
        self.dssp = dict([(k[1][1], v[1]) for (k,v) in self.dssp.items() if k[0] == chain])
        
        self.secseq = ''.join(self.dssp.values())
        self.to_threestate()

        self.features = self.make_features()

    def make_features(self):
        pos = 1
        feats = []
        for (c, group) in groupby(self.secseq):
            L = sum([1 for _ in group])
            feats.append(SecFeat(c, pos, pos + L - 1))
            pos += L
        return feats
    
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

    def shorthand(self):
        """
        Return a shorthand of secondary structure features, without lengths. Pure.
        """
        return ''.join(f.type for f in self.features)

    def runlength(self):
        """
        Convert the structure string to run-length encoding. Pure.
        """
        return ''.join(f.type + str(len(f)) for f in self.features)

    def print_report(self):
        """
        Pretty-print a report of secondary structure to stdout.
        """
        labels = {'H': 'HELIX', 'E': 'SHEET', 'C': 'COIL'}
        print("TYPE\tLENGTH\tSTART\tEND")
        for feature in self.features:
            print(f"{labels[feature.type]}\t{str(len(feature))}\t{str(feature.start)}\t{str(feature.end)}")

    def write_report(self, path):
        """
        Pretty-print a report of secondary structure to a file.
        """
        labels = {'H': 'HELIX', 'E': 'SHEET', 'C': 'COIL'}
        with open(path, 'w') as f:
            f.write("TYPE\tLENGTH\tSTART\tEND\n")
            for feature in self.features:
                f.write(f"{labels[feature.type]}\t{str(len(feature))}\t{str(feature.start)}\t{str(feature.end)}\n")

    def find_features(self, feature_type, feature_start = lambda x: True, feature_end = lambda x: True, feature_len = lambda x: True):
        """
        Find all features in a secondary structure that are of a certain type. Optionally, requirements on start position, end position or length can be added.
        """
        if not all([callable(feature_start), callable(feature_end), callable(feature_len)]):
            raise TypeError("feature_start, feature_end and feature_len must be callable")
        out_features = []
        for feature in self.features:
            if feature.type == feature_type and feature_start(feature.start) and feature_end(feature.end) and feature_len(len(feature)):
                out_features.append(feature)
        return out_features


