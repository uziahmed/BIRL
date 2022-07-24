import cpplib
import matplotlib.pyplot as plt

Bio = cpplib.Bio()

fileInp = input("Please enter the path to your genomic file: ")
gene = ""
file = open(fileInp, "r")
lines = []
gene = file.read().splitlines(True)
for line in gene:
    if (line.startswith("A" or "T" or "C" or "G")):
        lines.append(line)
    gene = "".join(lines)
    gene = gene.replace("\n", "")

minNuc = int(input("please enter the amount of minimun nucleotides: "))

print("Processing....................")

ORFs = Bio.ORFfinder(gene, minNuc);

longestORF = max(ORFs, key=len)

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

plt.savefig("fig.png")
plt.show()