#include <istream>
#include <fstream>      // std::ifstream
#include <string>
#include <vector>

#include "../include/read_control.h"


enum class CSVState {
  UnquotedField,
  QuotedField,
  QuotedQuote
};

std::vector<std::string> readCSVRow(const std::string &row) {
  CSVState state = CSVState::UnquotedField;
  std::vector<std::string> fields {""};
  size_t i = 0; // index of the current field
  for (char c : row) {
    switch (state) {
    case CSVState::UnquotedField:
      switch (c) {
      case ',': // end of field
        fields.push_back(""); i++;
        break;
      case '"': state = CSVState::QuotedField;
        break;
      default:  fields[i].push_back(c);
      break; }
      break;
    case CSVState::QuotedField:
      switch (c) {
      case '"': state = CSVState::QuotedQuote;
        break;
      default:  fields[i].push_back(c);
      break; }
      break;
    case CSVState::QuotedQuote:
      switch (c) {
      case ',': // , after closing quote
        fields.push_back(""); i++;
        state = CSVState::UnquotedField;
        break;
      case '"': // "" -> "
        fields[i].push_back('"');
        state = CSVState::QuotedField;
        break;
      default:  // end of quote
        state = CSVState::UnquotedField;
      break; }
      break;
    }
  }
  return fields;
}

/// Read CSV file, Excel dialect. Accept "quoted fields ""with quotes"""
std::vector<std::vector<std::string>> readCSV(std::istream &in) {
  std::vector<std::vector<std::string>> table;
  std::string row;
  while (!in.eof()) {
    std::getline(in, row);
    if (in.bad() || in.fail()) {
      break;
    }
    auto fields = readCSVRow(row);
    table.push_back(fields);
  }
  return table;
}





//' Read control file into armadillo matrix
//' 
//' @param control_file path to the control file
//' @return armadillo matrix with contents of control file
//' @details control file is a csv, with column headers, with one row per chromosome. 
//'     The csv control file has ten columns: chromosome, small effect summary file, large effect summary file, 
//'     bfile for reference data, sample size of the summary data, mafMax, number of SNPs, block information, heritability, number of threads,  
//'     and output path for SNP effect estimates.
//' @reference https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
//'     
std::vector<std::vector<std::string>> readControlFile(std::string control_file){
  // make istream object
  std::ifstream is;
  is.open(control_file);
  //read csv file from stream
  std::vector<std::vector<std::string>> out = readCSV(is);
  return (out);
} 
