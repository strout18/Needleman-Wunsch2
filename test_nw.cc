#include "nw_sim_mat.hh"
#include <iostream>

int main(int argc, char** argv) {
    if (argc == 5) {
        Needleman_sm test(argv[1], argv[2], atoi(argv[3]), argv[4]);
        test.run();
    }
    else if (argc == 6) {
        Needleman test(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        test.run();
    }
    else {
        std::cerr << "Wrong number of arguments! \n";
        return -3;
    }
   // Needleman test1("GCATGCU", "GATTACA", -1, 1, -1);
    return 0;
}
