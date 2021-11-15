#include <istream>
#include <fstream>      // std::ifstream
#include <string>
#include <vector>

std::vector<std::string> readCSVRow(const std::string &row);

std::vector<std::vector<std::string>> readCSV(std::istream &in);

std::vector<std::vector<std::string>> readControlFile(std::string control_file);
  
  