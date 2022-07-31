#include <iostream>
#include <vector>
#include <memory>
#include <string>

#include "Bio.hpp"

int main(int argc, char const *argv[])
{
    std::unique_ptr<Bio> bio(new Bio());

    std::string fileInp;
    std::cout<<"Please enter the path to your file: ";
    std::cin >>fileInp;

    std::ifstream bigFile(fileInp);
    constexpr size_t bufferSize = 1024 * 4096;
    std::unique_ptr<char[]> buffer(new char[bufferSize]);
    int counter = 0;
    while (bigFile)
    {
        bigFile.read(buffer.get(), bufferSize);
        counter += 4;
        bio->ORF_finder(buffer.get(), 70, 0, false, false);
        std::cout<<counter<< " MB's processed"<<std::endl;
    }

    for (int i = 0; i < bio->TotalORF.size(); i++)
    {
        std::cout<< bio->TotalORF[i]<<std::endl;
    }
    
    return 0;
}
