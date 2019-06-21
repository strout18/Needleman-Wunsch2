#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "nw_sim_mat.hh"
#include <cassert>

Needleman_sm::Needleman_sm(std::string seq_a, std::string seq_b, int indel, const char* filemat)
: Needleman(seq_a, seq_b, indel, 0, 0) {
    std::ifstream infile(filemat);

    assert(infile.is_open());
    std::getline(infile, letters_); //store order of letters in first row
    letters_.erase(std::remove(letters_.begin(), letters_.end(), '\t'), letters_.end()); //trim tab characters
    std::string currline;
    std::stringstream ss;
    while (std::getline(infile, currline)) {
        std::vector<int> row;
        ss.str(currline);
        std::string currwd;
        ss >> currwd;    //skip first column entry
        while (ss >> currwd) {
            row.push_back(stoi(currwd));
        }
        mat_.push_back(row);
        ss.clear(); //clear any flags
    }
    
}

int Needleman_sm::get_match_score(char a, char b) {
    int a_pos = std::find(letters_.begin(), letters_.end(), a) - letters_.begin(); //index of a in letters_
    int b_pos = std::find(letters_.begin(), letters_.end(), b) - letters_.begin(); //index of b in letters_
    return mat_[a_pos][b_pos];
}
