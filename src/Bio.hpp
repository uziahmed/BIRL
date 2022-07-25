#include<string>
#include<vector>
#include<unordered_map>
#include <algorithm> 
#include<fstream>

#include "include/json.hpp"

using json = nlohmann::json;

std::ifstream aminoAcidsFile("amino_acids.json");

class Bio{
    public:
        std::vector<std::string> codons;
        std::string gene;
        std::vector<std::string> TotalORF;
        json aminoAcidsData = json::parse(aminoAcidsFile);

        Bio(){}

        std::vector<std::string> toCodons(std::string gene){
            std::vector<std::string> output;
            for (int i = 0; i < gene.length(); i += 3) {
                if (gene.substr(i, 3).length() == 3)
                {
                    output.insert(output.end(),gene.substr(i, 3));
                    codons = output;
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

        std::vector<std::string> ORF_finder(std::string gene, int minNuc = 70){
            std::vector<std::string> singleCodon = toCodons(gene);
            for (int i = 0; i < singleCodon.size(); i++)
            {
                int currentIteration = 0;
                if (i < currentIteration){continue;}
                if (singleCodon[i] == "ATG")
                {
                    std::vector<std::string> currentORF;
                    currentORF.insert(currentORF.end(),singleCodon[i]);
                    for(int j = i; j < singleCodon.size(); j++)
                    {
                        if (singleCodon[j] == "TGA" || singleCodon[j] == "TAA" || singleCodon[j] == "TAG")
                        {
                            currentORF.insert(currentORF.end(),singleCodon[j]);
                            currentIteration = j;
                            break;
                        }else{
                            currentORF.insert(currentORF.end(),singleCodon[j]);
                            continue;
                        }
                    }
                    currentORF.erase(currentORF.begin());
                    std::string out = std::accumulate(currentORF.begin(), currentORF.end(), std::string(""));
                    if (out.size() >= minNuc) TotalORF.insert(TotalORF.end(),out);
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

};
