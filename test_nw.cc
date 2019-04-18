#include "nw.hh"

int main() {
    Needleman test1("GCATGCU", "GATTACA", 1, -1, -1);
    test1.run();
    return 0;
}
