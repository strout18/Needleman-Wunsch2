#pragma once
#include <vector>
#include <string>
#include <utility>

class Alt_Needleman{
    public:
      Alt_Needleman(std::string seq_a, std::string seq_b, int indel, int match, int mismatch);
      virtual ~Alt_Needleman() = default;
      void run();

    protected:
      void make_grid();
      virtual int get_match_score(char a, char b);
      int score_grid();
      std::pair<int, int> pos_from_align(std::string aligna, std::string alignb);
      void trace(std::string subaligna, std::string subalignb);     
      void print_alignments();
      std::vector<char> get_arrows(int row, int col);
      std::vector< std::vector<int> > grid_;   //used to store scores
      std::vector< std::pair<std::string, std::string> > align_;    //stores alignments
      std::string seq_a_;
      std::string seq_b_;
      int len_a_;
      int len_b_;
      int indel_;
      int match_;
      int mismatch_;
};
