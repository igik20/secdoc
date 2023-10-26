from from_af import fetch_af_struc

def main():
    with open("testdata/xrefs.txt") as f:
        for line in f:
            xref = line.strip()
            try:
                s = fetch_af_struc(xref)
                s.write_report(f"tables/{xref}.tsv")
            except Exception as e:
                print(f"Error while processing {xref}: \n {e}")

if __name__ == "__main__":
    main()