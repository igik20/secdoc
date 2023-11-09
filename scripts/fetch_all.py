from from_af import fetch_af_struc

def main():
    succ = 0
    fail = 0
    with open("testdata/xrefs.txt") as f:
        for line in f:
            xref = line.strip()
            try:
                s = fetch_af_struc(xref)
                s.write_report(f"tables/{xref}.tsv")
                succ += 1
            except Exception as e:
                print(f"Error while processing {xref}")
                fail += 1
    print(f"Success: {succ}")
    print(f"Fail: {fail}")

if __name__ == "__main__":
    main()

"""
Success: 83
Fail: 70
"""