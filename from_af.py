import requests
import tempfile
from secdoc import SecStruc

def fetch_af_struc(crossref):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{crossref}-F1-model_v4.pdb"
    req = requests.get(url, allow_redirects = True)
    temp = tempfile.NamedTemporaryFile()
    temp.write(req.content)
    s = SecStruc(temp.name)
    return s