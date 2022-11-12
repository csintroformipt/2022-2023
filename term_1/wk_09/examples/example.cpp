#include <fstream>
#include <iostream>
#include <limits>

#include "json.hpp"

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    if (argc != 2){
        std::cout << "Error: invalid number of CMD arguments: " << argc <<" != 2. Aborting." << std::endl;
        exit(1);
    }
    else{
        std::cout << "Using input file " << argv[1] << std::endl;
    }
    
    std::ifstream f(argv[1]);
    if (!f.good()){
        std::cout << "invalid input file: " << argv[1] << " Aborting." << std::endl;
        exit(2);
    }

    json data = json::parse(f);

    std::cout << "json content:" << std::endl;

    std::cout << data << std::endl;

    int some_int_var = data["key1_int"];
    std::string some_string = data["key2_string"];
    double some_double_var = data["dict_key"]["dict_value"];

    std::cout << "here's some values from variables based on the json:" << std::endl;

    std::cout.precision(std::numeric_limits< double >::max_digits10 + 2); // PLEASE TRY COMMENTING ME

    std::cout << "\tint var:" << some_int_var << std::endl << 
                 "\tdouble var: " << some_double_var << std::endl <<
                 "\tstring var: " << some_string << std::endl;
    
    return 0;
}
