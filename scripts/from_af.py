import requests
import tempfile
import subprocess
from secdoc import SecStruc

def fetch_af_struc(crossref):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{crossref}-F1-model_v4.pdb"
    req = requests.get(url, allow_redirects = True)
    temp = tempfile.NamedTemporaryFile()
    temp.write(req.content)
    res = subprocess.run(["grep", "-v", "DBREF", temp.name], capture_output = True)
    temp2 = tempfile.NamedTemporaryFile()
    temp2.write(res.stdout)
    s = SecStruc(temp2.name)
    return s