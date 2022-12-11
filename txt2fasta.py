import string
import sys


if __name__ == "__main__":
    file_path = sys.argv[1]
    id_prefix = sys.argv[2]
    out_path = sys.argv[3]


    with open(file_path, mode="r") as rfp:
        lines = rfp.readlines()
    
    seqs = [x.strip() for x in lines if len(x.strip()) > 0]
    ids_ = [f"{id_prefix}{i}" for i,_ in enumerate(seqs)]


    out_text = ""
    for i in range (len(seqs)):
        new_text = f">{ids_[i]}\n{seqs[i]}\n"
        out_text += new_text

    with open (out_path,'w') as wfp:
        wfp.write(out_text)
