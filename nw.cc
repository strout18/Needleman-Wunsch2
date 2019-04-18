#include "nw.hh"
#include <algorithm>
#include <iostream>
#include <cassert>

Needleman::Needleman(std::string seq_a, std::string seq_b, int match, int indel, int mismatch)
: seq_a_(seq_a), seq_b_(seq_b), len_a_(seq_a.length()), len_b_(seq_b.length()), match_(match), indel_(indel), mismatch_(mismatch) {
}

void Needleman::make_grid() {
    //makes first row of grid
    std::vector<int> line1;
    std::vector<std::string> line1dir(len_a_ + 1, "l");         //all going left from top row
    for (auto i = 0; i < len_a_ + 1; i++) { //first row in form <0, indel, 2 * indel, ..., len_a_ * indel>
        line1.push_back(indel_ * i);
    }
    grid_.push_back(line1);
    trace_.push_back(line1dir);

    //makes rest of grid, which will initially be len_b_ + 1 rows where each row is populated with len_a_ of the row's first number
    for (auto i = 1; i < len_b_ + 1; i++) {        //starts from 1 bc we've already done first row
        std::vector<int> new_vect = {indel_ * i};
        grid_.push_back(new_vect);
    }
    assert (grid_.size() == len_b_ + 1);
}

//seq_a on top, seq_b on side
int Needleman::score_grid() {                           //returns alignment score
    for (auto row = 1; row < len_b_ + 1; row++) { //by row
        //can only go top from first in row
        std::vector<std::string> row_vect = {"t"};                 //e.g. < <'t', 'td', ..., 'd'>, <'d', ..., 't'>, ... , <'l', 'l', ... 'l'> > (stored by rows)
        for (auto col = 1; col < len_a_ + 1; col++) {
            assert (row < grid_.size() && col == grid_[row].size());          //make sure we're in vector bounds and j is about to be newest elem in vector
            int top = grid_[row - 1][col] + indel_;       //cell above
            int left = grid_[row][col - 1] + indel_;      //cell to the left
            int diag = grid_[row-1][col-1];
            if (seq_b_[row-1] == seq_a_[col-1]) {           //everything is shifted top left by 2 place bc extra column + row
                diag += match_;
            }
            else {
                diag += mismatch_;
            }
            int maximum = std::max({top, left, diag});
            grid_[row].push_back(maximum);

            std::string dirstring;
            //assigning directions for traceback
            if (maximum == top) {
                dirstring += 't';
            };
            if (maximum == left) {
                dirstring += 'l';
            }
            if (maximum == diag) {
                dirstring += 'd';
            }

            row_vect.push_back(dirstring);
        }
        trace_.push_back(row_vect);
    }
    return grid_[len_a_][len_b_];
}

std::pair<int, int> Needleman::pos_from_path(std::string path) { //returns the grid position the path leads to
    int row = len_b_;         //start at bottom left corner
    int col = len_a_;
    for (char c: path) {
        if (c != 'l') {         //c == t or d
            row -= 1;
        }
        if (c != 't') {
            col -= 1;             //c == l or d
        }
    }
    return std::make_pair(row, col);
}

void Needleman::get_dirs(std::string dirs) {        //helper function for trace
    auto pos = pos_from_path(dirs);
    if (pos == std::pair(0, 0)) {             //exit condition
        paths_.push_back(dirs);
        return;
    }
    else {
        std::string traceval = trace_[pos.first][pos.second];
        for (char c: traceval) {
            get_dirs(dirs + c);
        }
    }
}

void Needleman::trace() {
    get_dirs("");
    for (auto str: paths_) {
        std::string aligna;
        std::string alignb;

        int col = len_a_;         //start at bottom left corner
        int row = len_b_;
        for (char c: str) {
            if (c == 'l') {     //indel
                aligna = seq_a_[col-1] + aligna;
                alignb = "-" + alignb;
                col -= 1;
            }
            else if (c == 't') {    //indel
                aligna = "-" + aligna;
                alignb = seq_b_[row-1] + alignb;
                row -= 1;
            }
            else {                  //match or mismatch
                aligna = seq_a_[col-1] + aligna;
                alignb = seq_b_[row-1] + alignb;
                col -= 1;
                row -= 1;
            }
        }
        align_.push_back(std::make_pair(aligna, alignb));
    }
}

void Needleman::print_alignments() {
    for (auto alignments: align_) {
        std::cout << alignments.first << "\n" << alignments.second << "\n";
        std::cout << "---" << "\n";
    }
    return;
}

void Needleman::run() {
    make_grid();
    //std::cout << score_grid() << "\n";
    score_grid();
    trace();
    print_alignments();
    return;
}
