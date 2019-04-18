#pragma once
#include <vector>
#include <string>
#include <utility>

class Needleman{
    public:
      Needleman(std::string seq_a, std::string seq_b, int match, int indel, int mismatch);
      ~Needleman() = default;
      void make_grid();
      int score_grid();
      void get_dirs(std::string dirlst);
      void trace();     
      void print_alignments();
      std::pair<int, int> pos_from_path(std::string path);
      void run();

    private:
      std::vector< std::vector<int> > grid_;   //used to store scores
      std::vector< std::vector<std::string> > trace_;  //used to store directions
      std::vector< std::pair<std::string, std::string> > align_;    //stores alignments
      std::vector<std::string> paths_;
      std::string seq_a_;
      std::string seq_b_;
      int len_a_;
      int len_b_;
      int match_;
      int indel_;
      int mismatch_;
};
