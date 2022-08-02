import cpplib
import matplotlib.pyplot as plt
import os
import json

Bio = cpplib.Bio()

userInp = input("do you want to process an entire genome or a sequence(1:genone, 2:sequence): ")
fileInp = input("Please enter the path to your genomic file: ")
if userInp == "genome" or "1":
    Bio.genomeProcessor(fileInp)
else:
    gene = ""
    with open(fileInp, "r") as geneFile:
        gene = geneFile.read().splitlines(True)
    lines = []
    for line in gene:
        if not(line.startswith(">")):
            lines.append(line)
        gene = "".join(lines)
        gene = gene.replace("\n", "")

    minNuc = int(input("please enter the amount of minimun nucleotides: "))

    print("Processing....................")

    for i in range(6):
        Bio.ORFfinder(gene, minNuc, frame=i);

print("finding longest ORF...............")
longestORF = max(Bio.totalORFs, key=len)

print(f"Longest ORF: {longestORF} \n Total ORF's found {len(Bio.totalORFs)}")

if not os.path.exists('Data'):
    os.makedirs('Data')
print("writting data to a file.............")
print(f"Data/genome{str(len(os.listdir('Data')))}.json")
data = {
    "Id":len(os.listdir("Data")),
    "FileName": fileInp,
    "LongestORF":longestORF,
    "protien":Bio.translate(longestORF)
}

currentFile = str(len(os.listdir('Data')))

with open(f"Data/genome{currentFile}.json", "w") as jsonFile:
    json.dump(data, jsonFile)

nucCount = Bio.nucCounter(longestORF);
nucleotideLabels = []
nucleotideY = []
for key, value in nucCount.items():
    nucleotideLabels.append(key)
    nucleotideY.append(value)

aminoAcids = Bio.translate(longestORF)
aminoAcidsCount = Bio.aminoAcidsCounter(aminoAcids)
aminoAcidsPieChartLables = []
aminoAcidsY = []
for key, value in aminoAcidsCount.items():
    aminoAcidsPieChartLables.append(key)
    aminoAcidsY.append(value)

codons = Bio.toCodon(longestORF)
codonCount = Bio.codonCounter(codons);
codonPieChartLabels = []
codonY = []
for key, value in codonCount.items():
    codonPieChartLabels.append(key)
    codonY.append(value)

plt.subplot(1, 3, 1)
plt.pie(nucleotideY, labels = nucleotideLabels, textprops = {"fontsize":7}, autopct="%.2f%%")
plt.title("Nucleotide Frequency")
plt.legend(nucleotideY,title="nucleotide frequency", loc="lower right", prop={"size":7}, bbox_to_anchor=(1.3, -0.1), ncol=2)

plt.subplot(1, 3, 2)
plt.pie(codonY, labels = codonPieChartLabels, textprops = {"fontsize":6}, autopct="%.2f%%")
plt.title("Codon Frequency")
plt.legend(codonY,title="codons", loc="lower right", prop={"size":6}, bbox_to_anchor=(1.2, -0.5), ncol=4)

plt.subplot(1, 3, 3)
plt.pie(aminoAcidsY, labels = aminoAcidsPieChartLables, textprops = {"fontsize":8}, autopct="%.2f%%")
plt.title("Amino Acids")
plt.legend(aminoAcidsY,title="Amino Acid Frequency", loc="lower right", prop={"size":5}, bbox_to_anchor=(1.3, -0.4), ncol=2)


plt.savefig(f"Data/genome{currentFile}.png", dpi=1200)