#include "src/solve.h"
// #include "ortoo"


int main()
{
    SolveILP ex = SolveILP("spec_files/5_block.ilp", 10,  true);
    vector<float> x_i, y_i, w_i, h_i;
    vector<int> z_i;
    float Y;

    // ex.create_constraints();
    tie(Y, x_i, y_i, z_i, w_i, h_i) = ex.solve(2);
    vector<float> W, H;
    tie(W, H) = ex.visualize(Y, x_i, y_i, z_i, w_i, h_i);

    
    // cout << "Y = " << Y << endl;

    return 0;
}

