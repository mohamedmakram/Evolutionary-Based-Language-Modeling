import plotly.express as px

fileInput = open("protein_seqs.fasta", "r")
fileOutput = open("clean_protein_seqs.fasta", "w")

count = 1
l=[]

for line in fileInput:
    line = line.rstrip("\n")
    #fileOutput.write(line + "\n")
    l.append(len(line))
    if len(line)>352:
      continue
    if len(line) < 288:
       continue
    fileOutput.write(line + "\n")
    count = count + 1
fig = px.box(l)

fig.show()
fileInput.close()
fileOutput.close()