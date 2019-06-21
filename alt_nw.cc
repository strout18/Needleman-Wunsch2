#include "alt_nw.hh"
#include <algorithm>
#include <iostream>
#include <cassert>

Alt_Needleman::Alt_Needleman(std::string seq_a, std::string seq_b, int indel, int match, int mismatch)
: seq_a_(seq_a), seq_b_(seq_b), len_a_(seq_a.length()), len_b_(seq_b.length()), match_(match), indel_(indel), mismatch_(mismatch) {
}

void Alt_Needleman::make_grid() {
    //makes first row of grid
    std::vector<int> line1;
    for (auto i = 0; i < len_a_ + 1; i++) { //first row in form <0, indel, 2 * indel, ..., len_a_ * indel>
        line1.push_back(indel_ * i);
    }
    grid_.push_back(line1);
}

int Alt_Needleman::get_match_score(char a, char b) {
    if (a == b) {
        return match_;
    }
    else {return mismatch_;}
}

int Alt_Needleman::score_grid() {
    for (auto row = 1; row < len_b_ + 1; row++) { //by row

        std::vector<int> new_vect = {indel_ * row};
        grid_.push_back(new_vect);

        for (auto col = 1; col < len_a_ + 1; col++) {
            assert (row < grid_.size() && col == grid_[row].size());          //make sure we're in vector bounds and j is about to be newest elem in vector
            int top = grid_[row - 1][col] + indel_;       //cell above
            int left = grid_[row][col - 1] + indel_;      //cell to the left
            int diag = grid_[row-1][col-1];
            diag += get_match_score(seq_b_[row-1], seq_a_[col-1]); //everything is shifted top left by 2 place bc extra column + row
            int maximum = std::max({top, left, diag});
            grid_[row].push_back(maximum);
        }
    }
    return grid_[len_a_][len_b_];
}

std::vector<char> Alt_Needleman::get_arrows(int row, int col) { //get which cell (t, l, d) was responsible for current cell score
    std::vector<char> ret;
    if (row && grid_[row - 1][col] + indel_ == grid_[row][col]) {
        ret.push_back('t');
    };
    if (col && grid_[row][col - 1] + indel_ == grid_[row][col]) {
        ret.push_back('l');
    };
    if (col && row && grid_[row-1][col-1] + get_match_score(seq_b_[row-1], seq_a_[col-1]) == grid_[row][col]) {
        ret.push_back('d');
    };
    return ret;
};

std::pair<int, int> Alt_Needleman::pos_from_align(std::string aligna, std::string alignb) {
    int row = len_b_;
    int col = len_a_;
    assert (aligna.length() == alignb.length());
    for (auto i = 0; i < aligna.length(); i++) {
        row -= 1;
        col -= 1;
        if (aligna[i] == '-') {
            col += 1;
        };
        if (alignb[i] == '-') {     //should it be an else? is it even possible for both to align with a gap?
            row += 1;
        };
    };
    return std::make_pair(row, col);
}

void Alt_Needleman::trace(std::string subaligna, std::string subalignb) {
    std::pair position = pos_from_align(subaligna, subalignb);
    int& row = position.first;  //DOUBLE CHECK
    int& col = position.second;
    //while ((row || col) && (row >=0 && col >= 0)) {
    while (row || col) {
        std::vector<char> arrows = get_arrows(row, col);
        char arrow;
        auto orig = arrows.size();
        for (auto i = 0; i < orig; i++) {
            arrow = arrows.back();
            if (arrow == 't') {
                if (arrows.size() > 1) {
                    arrows.pop_back();
                    trace("-" + subaligna, seq_b_[row-1] + subalignb);
                }
                else {
                    subaligna = "-" + subaligna;
                    subalignb = seq_b_[row-1] + subalignb;
                }
            }
            else if (arrow == 'l') {
                if (arrows.size() > 1) {
                    arrows.pop_back();
                    trace(seq_a_[col-1] + subaligna, "-" + subalignb);
                }
                else {
                    subaligna = seq_a_[col-1] + subaligna;
                    subalignb = "-" + subalignb;
                }
            }
            else { 
                if (arrows.size() > 1) {
                    arrows.pop_back();
                    trace(seq_a_[col-1] + subaligna, seq_b_[row-1] + subalignb);
                }
                else {
                    subaligna = seq_a_[col-1] + subaligna;
                    subalignb = seq_b_[row-1] + subalignb;
                }
            }
             //continue in while loop
            arrow = arrows.back();
        };
        if (arrow != 'l') {         //arrow == t or d
            row -= 1;
        }
        if (arrow != 't') {
            col -= 1;             //arrow == l or d
        }
    };
    align_.push_back(std::make_pair(subaligna, subalignb));
    return;
}

void Alt_Needleman::print_alignments() {
    for (auto alignments: align_) {
        std::cout << alignments.first << "\n" << alignments.second << "\n";
        std::cout << "---" << "\n";
    }
    return;
}

void Alt_Needleman::run() {
    make_grid();
    //std::cout << score_grid() << "\n";
    score_grid();
    trace("", "");
    print_alignments();
    return;
}