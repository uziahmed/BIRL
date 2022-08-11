#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm> 
#include <fstream>
#include <iostream>

#include "include/json.hpp"

using json = nlohmann::json;

std::ifstream aminoAcidsFile("amino_acids.json");

class Bio{
    public:
        std::vector<std::string> TotalORF;
        json aminoAcidsData = json::parse(aminoAcidsFile);

        Bio(){}

        std::vector<std::string> toCodons(std::string gene, int start = 0){
            std::vector<std::string> output;
            for (int i = start; i < gene.length(); i += 3) {
                if (gene.substr(i, 3).length() == 3)
                {
                    output.push_back(gene.substr(i, 3));
                }
            }
            return output;
        }

        std::string reverse_comp(std::string gene){
            std::unordered_map<char, char> rev_comp;
            rev_comp['T'] = 'A';
            rev_comp['A'] = 'T';
            rev_comp['C'] = 'G';
            rev_comp['G'] = 'C';
            std::vector<char> rev;
            for (int i = 0; i < gene.length(); i++)
            {
                rev.push_back(rev_comp[gene[i]]);
            }
            reverse(rev.begin(), rev.end());
            std::string output(rev.begin(), rev.end());
            return output;
        }

        std::string compliment(std::string gene){
            std::unordered_map<char, char> rev_comp;
            rev_comp['T'] = 'A';
            rev_comp['A'] = 'T';
            rev_comp['C'] = 'G';
            rev_comp['G'] = 'C';
            std::vector<char> rev;
            for (int i = 0; i < gene.length(); i++)
            {
                rev.push_back(rev_comp[gene[i]]);
            }
            std::string output(rev.begin(), rev.end());
            return output;
        }

        std::string transcribe(std::string gene){
            gene = compliment(gene);
            for (int i = 0; i < gene.length(); i++)
            {
                if (gene[i] == 'T')
                {
                    gene[i] = 'U';
                }
            }
            return gene;
        }

        std::vector<std::string> ORF_finder(std::string gene, int minNuc = 70, int frame = 0, bool isGenome = false){
            std::vector<std::string> singleCodon;
            if (isGenome == false)
            {
                if(frame <= 2) singleCodon = toCodons(gene, frame);
                else singleCodon = toCodons(reverse_comp(gene), frame-3);
            }else singleCodon = toCodons(gene);
            
            std::vector<std::string> currentORF;
            bool startCodonFound = false;
            for (int i = 0; i < singleCodon.size(); i++)
            {
                if (singleCodon[i] == "ATG")
                {
                    currentORF.push_back(singleCodon[i]);
                    startCodonFound = true;
                }
                else if(singleCodon[i] == "TAA" || singleCodon[i] == "TAG" || singleCodon[i] == "TGA"){
                    if (startCodonFound)
                    {
                        currentORF.push_back(singleCodon[i]);
                        std::string out = std::accumulate(currentORF.begin(), currentORF.end(), std::string(""));
                        TotalORF.push_back(out);
                        currentORF.clear();
                        startCodonFound = false;
                    }
                }else if(singleCodon[i] == "NNN"){
                    currentORF.clear();
                    startCodonFound = false;
                }
                else{
                    if (startCodonFound)
                    {
                        currentORF.push_back(singleCodon[i]);
                    }
                }
            }
            return TotalORF;
        }

        std::unordered_map<std::string, int> codonCounter(std::vector<std::string> gene){
            std::unordered_map<std::string, int> output;
            for (std::string i : gene)
            {
                output[i]++;
            }
            return output;
        }

        std::unordered_map<char, int> aminoAcidsCounter(std::string protien){
            std::unordered_map<char, int> output;
            for (char i : protien)
            {
                output[i]++;
            }
            return output;
        }

        std::string translate(std::string Seq){
            std::vector<std::string> gene = toCodons(Seq);
            std::vector<std::string> output;
            for (int i = 0; i < gene.size(); i++)
            {
                for (int j = 0; j < aminoAcidsData.size(); j++)
                {
                    for (int k = 0; k < aminoAcidsData[j]["codons"].size(); k++)
                    {
                        if (gene[i] == aminoAcidsData[j]["codons"][k])
                        {
                            output.push_back(aminoAcidsData[j]["SLA"]);
                        }
                    }
                }
            }
            std::string out = std::accumulate(output.begin(), output.end(), std::string(""));
            return out;
            
        }

        std::unordered_map<char, int> nucleotideFrequency(std::string gene){
            std::unordered_map<char, int> nucTable;
            nucTable['T'] = 0;
            nucTable['G'] = 0;
            nucTable['A'] = 0;
            nucTable['C'] = 0;
            for (int i = 0; i < gene.size(); i++)
            {
                if (gene[i] == 'A') nucTable['A']++;
                if (gene[i] == 'C') nucTable['C']++;
                if (gene[i] == 'T') nucTable['T']++;
                if (gene[i] == 'G') nucTable['G']++;
            }
            return nucTable;
        }

        std::vector<std::string> genomeProcessor(std::string fileInp){
            std::ifstream bigFile(fileInp);
            constexpr size_t size = 10;
            constexpr size_t bufferSize = 1024 * (1024 * size);
            std::unique_ptr<char[]> buffer(new char[bufferSize]);
            int counter = 0;
            std::cout<<"started processing.......\n";
            while (bigFile)
            {
                bigFile.read(buffer.get(), bufferSize);
                std::string Seq = std::string(buffer.get(), bigFile.gcount());
                Seq.erase(std::remove(Seq.begin(), Seq.end(), '\n'), Seq.end());
                std::transform(Seq.begin(), Seq.end(), Seq.begin(), ::toupper);
                counter += size;
                ORF_finder(Seq, 70, 0, true);
                std::cout<<counter<< " MB's processed"<<std::endl;
            }
            return TotalORF;
        }
};
