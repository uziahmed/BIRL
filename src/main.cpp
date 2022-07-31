#include <iostream>
#include <vector>
#include <memory>
#include <string>

int main(int argc, char const *argv[])
{
    std::unique_ptr<Bio> bio(new Bio());

    std::string fileInp;
    std::cout<<"Please enter the path to your file: ";
    std::cin >>fileInp;

    std::vector<std::string> allORF;
    std::ifstream bigFile(fileInp);
    constexpr size_t bufferSize = 1024 * 4096;
    std::unique_ptr<char[]> buffer(new char[bufferSize]);
    while (bigFile)
    {
        bigFile.read(buffer.get(), bufferSize);
        bio->ORF_finder(buffer.get());
    }

    for (int i = 0; i < bio->TotalORF.size(); i++)
    {
        std::cout<< bio->TotalORF[i]<<std::endl;
    }
    
    return 0;
}
