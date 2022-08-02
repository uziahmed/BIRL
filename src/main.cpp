#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <fstream>

#include "Bio.hpp"
#include "include/json.hpp"

using json = nlohmann::json;

int main(int argc, char const *argv[])
{
    std::unique_ptr<Bio> bio(new Bio());

    std::string fileInp;
    std::cout<<"Please enter the path to your file: ";
    std::cin>>fileInp;

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
        bio->ORF_finder(Seq, 70, 0, true);
        std::cout<<counter<< " MB's processed"<<std::endl;
    }

    // for(auto orf : bio->TotalORF){
    //     std::cout<<orf<<std::endl;
    // }
    std::cout<<"finding largest ORF........"<<std::endl;
    auto longestORF = std::max_element(bio->TotalORF.begin(), bio->TotalORF.end(),
                            [](const auto& a, const auto& b) {
                                return a.size() < b.size();
                            });

    std::filesystem::create_directory("Data");

    int fileCount = 0;
    for(const auto & file : std::filesystem::directory_iterator("Data")){
        fileCount++;
    }

    std::string outFile = "Data/seq"+std::to_string(fileCount)+".json";
    std::ofstream MyFile(outFile);

    json data = {
        {"Id", fileCount},
        {"GenomeName", fileInp},
        {"FileName", outFile},
        {"LongestORF", *longestORF}
    };

    std::cout<<"writting data to file"<<std::endl;
    MyFile << data;

    MyFile.close();
    std::cout<<outFile<<" saved."<<std::endl;

    return 0;
}
