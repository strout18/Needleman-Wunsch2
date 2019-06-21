#pragma once
#include "nw.hh"

class Needleman_sm: public Needleman{
    public:
      Needleman_sm(std::string seq_a, std::string seq_b, int indel, const char* filemat);
      virtual ~Needleman_sm() = default;
    private:
        std::vector< std::vector<int> > mat_;
        std::string letters_;
        int get_match_score(char a, char b) override;
};
