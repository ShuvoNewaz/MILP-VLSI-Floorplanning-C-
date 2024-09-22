#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <math.h>

#include "fusion.h"
using namespace std;
using namespace mosek::fusion;
using namespace monty;

int factorial(int n)
{
    if(n == 0)
    {
        return 1;
    }
    else
    {
        return n * factorial(n - 1);
    }
}

int NcR(int n, int r)
{
    return factorial(n) / (factorial(r) * factorial(n - r));
}

int position(int n, int i, int j)
{
    int pos = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
    return pos;
}


shared_ptr<monty::ndarray<double, 1> > nonOverlapConstraint(int total_variables,
                        int ind_1, int ind_2, int ind_3, int ind_4, int ind_5,
                        float val_3, float val_4, float val_5)
{
    vector<double> A(total_variables);
    A[ind_1] = 1; 
    A[ind_2] = -1;
    A[ind_3] = val_3; 
    A[ind_4] = val_4;
    A[ind_5] = val_5;

    return new_array_ptr<double>(A);
}

void greaterThanZeroConstraint(int total_variables, int ind_1, Model::t M, Variable::t X)
{
    vector<double> A(total_variables);
    A[ind_1] = 1;
    auto Coefficients = new_array_ptr<double>(A);
    M->constraint(Expr::dot(Coefficients, X), Domain::greaterThan(0)); // x, y >= 0
}

void softWidthConstraint(int total_variables, int ind_1, float w_min, float w_max, Model::t M, Variable::t X)
{
    vector<double> A(total_variables);
    A[ind_1] = 1;
    auto Coefficients = new_array_ptr<double>(A);
    M->constraint(Expr::dot(Coefficients, X), Domain::inRange(w_min, w_max));
}

void chipDimensionConstraint(int total_variables, int ind_1, int ind_2, float dim_1, float dim_2, Model::t M, Variable::t X)
{
    vector<double> A(total_variables);
    A[ind_1] = 1;
    A[ind_2] = dim_1 - dim_2;
    A[total_variables - 1] = -1;  // Y = -1
    auto Coefficients = new_array_ptr<double>(A);
    M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(-dim_2));
}

void binaryConstraint(int total_variables, int ind_1, Model::t M, Variable::t X)
{
    vector<double> A(total_variables);
    A[ind_1] = 1;
    auto Coefficients = new_array_ptr<double>(A);
    M->constraint(Expr::dot(Coefficients, X), Domain::inRange(0, 1)); // 0 <= x_ij <= 1
}

void writeFile(ofstream &output_file, vector<float> variable, unsigned short int count)
{
    for(unsigned short int i=0; i<count-1; i++)
    {
        output_file << to_string(variable[i]) + ",";
    }
    output_file << to_string(variable[count-1]) + "\n";
}