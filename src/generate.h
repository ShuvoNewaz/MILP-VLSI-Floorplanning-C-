#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <math.h>
#include "fusion.h"
using namespace mosek::fusion;
using namespace monty;

using namespace std;

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

class GenerateProblem
{
    public:
        bool hard_exists;
        bool soft_exists;
        unsigned short int num_hard_modules, num_soft_modules, num_total_modules;
        vector<float> hard_module_width, hard_module_height;
        vector<float> area, min_aspect, max_aspect;
        vector<vector<float>> soft_module_width_range, soft_module_height_range;
        vector<float> gradient, intercept;
        float bound;

        string file;
        int num_blocks;
        bool underestimation;
        vector<string> lines;
        string output {"constraints.lp"};
        ofstream constraint_file;

        GenerateProblem(string, int, bool); // Constructor is defined outside the class

        tuple<int, int> total_modules()
        /*
            Returns a tuple (number of hard modules, number of soft modules)
            found in the input file.
        */
        {

            for(string line : lines)
            {
                if(line.substr(0, 4) == "hard")
                {
                    hard_exists = true;
                    num_hard_modules = stoi(line.substr(7, line.size() - 7));
                }
                else if(line.substr(0, 4) == "soft")
                {
                    soft_exists = true;
                    num_soft_modules = stoi(line.substr(7, line.size() - 7));
                }
                else if(hard_exists && !soft_exists)
                {
                    num_soft_modules = 0;
                }
                else if(!hard_exists && soft_exists)
                {
                    num_hard_modules = 0;
                }
            }

            return make_tuple(num_hard_modules, num_soft_modules);
        }

        tuple<vector<float>, vector<float>> hard_module_dimension()
        /*
            Returns tuple containing the dimensions (width, height)
            of the hard modules.
        */
        {
            tuple<vector<float>, vector<float>> hard_dimension;
            if (hard_exists)
            {
                unsigned short int i = 0;

                for(string line : lines)
                {
                    if(line.substr(0, 4) == "hard")
                    {
                        continue;
                    }
                    
                    unsigned short int comma_index = line.find(",");
                    float width, height;
                    width = stof(line.substr(0, comma_index+1));
                    height = stof(line.substr(comma_index+1, line.size()-comma_index));
                    get<0>(hard_dimension).push_back(width);
                    get<1>(hard_dimension).push_back(height);
                    i += 1;
                    if(i >= num_hard_modules)
                    {
                        break;
                    }
                }
            }
            else
            {
                get<0>(hard_dimension).push_back(0);
                get<1>(hard_dimension).push_back(0);;
            }

            return hard_dimension;
        }

        tuple<vector<float>, vector<float>, vector<float>> soft_module_properties()
        /*
            Returns a 2-d array (vector) containing the properties
            (area, minimum aspect ration, maximum aspect ration)
            of the soft modules.
        */
        {
            bool soft_encountered {false};
            tuple<vector<float>, vector<float>, vector<float>> soft_properties;
            if (num_soft_modules != 0)
            {
                unsigned short int i = 0;
                for(string line : lines)
                {
                    if(line.substr(0, 4) == "soft")
                    {
                        soft_encountered = true;
                        continue;
                    }
                    if(soft_encountered)
                    {
                        unsigned short int comma_index_1 = line.find(",");
                        unsigned short int comma_index_2 = line.rfind(",");
                        float area, min_aspect, max_aspect;
                        area = stof(line.substr(0, comma_index_1+1));
                        min_aspect = stof(line.substr(comma_index_1+1, comma_index_2-comma_index_1));
                        max_aspect = stof(line.substr(comma_index_2+1, line.size()-comma_index_2));
                        get<0>(soft_properties).push_back(area);
                        get<1>(soft_properties).push_back(min_aspect);
                        get<2>(soft_properties).push_back(max_aspect);
                        i += 1;
                        if(i >= num_soft_modules)
                        {
                            break;
                        }
                    }
                }
            }
            else
            {
                get<0>(soft_properties).push_back(0);
                get<1>(soft_properties).push_back(0);
                get<2>(soft_properties).push_back(0);
            }

            return soft_properties;
        }
        
        tuple<vector<vector<float>>, vector<vector<float>>> soft_module_dimension_range()
        /*
            Returns a tuple containing the minimum and maximum allowable widths and
            heights for the soft modules.
        */
        {
            vector<vector<float>> width_range, height_range;
            if(soft_exists)
            {
                float a, min_w, max_w, min_h, max_h;
                for(unsigned short int i = 0; i < num_soft_modules; i++)
                {
                    a = area[i];
                    min_w = sqrtf(a * min_aspect[i]);
                    max_w = sqrtf(a * max_aspect[i]);
                    width_range.push_back({min_w, max_w});
                    if(underestimation)
                    {
                        min_h = a / max_w;
                        max_h = a / max_w + (max_w - min_w) * a / (max_w * max_w);
                    }
                    else
                    {
                        min_h = a / max_w;
                        max_h = a / min_w;
                    }
                    height_range.push_back(vector<float> {min_h, max_h});
                }                
            }
            else
            {
                width_range = {{0}, {0}};
                height_range = {{0}, {0}};
            }

            return make_tuple(width_range, height_range);
        }
        
        tuple<vector<float>, vector<float>> linear_approximation()
        {
            float a, g, c;
            float min_w, max_w, min_h, max_h;
            vector<float> gradients, intercepts;

            if(soft_exists)
            {
                for(short unsigned int i = 0; i < num_soft_modules; i++)
                {
                    a = area[i];
                    min_w = soft_module_width_range[i][0];
                    max_w = soft_module_width_range[i][1];
                    min_h = soft_module_height_range[i][0];
                    max_h = soft_module_height_range[i][1];
                    if(underestimation)
                    {
                        g = - a / (max_w * max_w);
                        c = 2 * a / max_w;
                    }
                    else
                    {
                        g = (max_h - min_h) / (min_w - max_w);
                        c = max_h - g * min_w;
                    }
                    gradients.push_back(g);
                    intercepts.push_back(c);
                }
            }
            else
            {
                gradients = {0};
                intercepts = {0};
            }

            return make_tuple(gradients, intercepts);
        }

        float upper_bound()
        /*
                Assuming a square chip block, this functions returns the maximum required
                length (and width) of a block that can be formed using the given modules.
        */
        {
            float W_hard{0}, H_hard{0}, W_soft{0}, H_soft{0};
            float W_block, H_block;

            for(short unsigned int i = 0; i < num_hard_modules; i++)
            {
                W_hard += max(hard_module_width[i], hard_module_height[i]);
            }
            H_hard = W_hard;

            for(short unsigned int i = 0; i < num_soft_modules; i++)
            {
                W_soft += soft_module_width_range[i][1];
                H_soft += soft_module_height_range[i][1];
            }

            W_block = W_hard + W_soft;
            H_block = H_hard + H_soft;

            return max(W_block, H_block);
        }

        /*
            Start putting together the *.lp file containing the constraints for
            the problem. This file will allow solving using an LP solver.
        */

        void objective()
        {
            constraint_file.open(output, ios::out);
            if(constraint_file.is_open())
            {
                constraint_file << "/* Objective Function */\nmin: Y;\n\n\n";
            }
            constraint_file.close();
        }

        void hard_hard_nonoverlap()
        {
            if(hard_exists)
            {
                float width_i, height_i, width_j, height_j;
                constraint_file.open(output, ios::app);
                constraint_file << "/* Non-overlap constraints hard-hard */\n";

                for(short unsigned int i = 1; i <= num_hard_modules; i++)
                {
                    for(short unsigned int j = 1; j <= num_hard_modules; j++)
                    {
                        if(j > i)
                        {
                            width_i = hard_module_width[i-1];
                            height_i = hard_module_height[i-1];
                            width_j = hard_module_width[j-1];
                            height_j = hard_module_height[j-1];
                            
                            constraint_file << "x" << to_string(i) << " + " << to_string(height_i) << " z" << to_string(i) << " + " << to_string(width_i) << " - " << to_string(width_i) << " z" << to_string(i) << " <= x" << to_string(j) << " + " << to_string(roundf(bound)) + " x" << to_string(i) << to_string(j) << " + " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                            constraint_file << "x" << to_string(i) << " - " << to_string(height_j) << " z" << to_string(j) << " - " << to_string(width_j) << " + " << to_string(width_j) << " z" << to_string(j) << " >= x" << to_string(j) << " - " << to_string(roundf(bound)*1) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " - " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                            constraint_file << "y" << to_string(i) << " + " << to_string(width_i) << " z" << to_string(i) << " + " << to_string(height_i) << " - " << to_string(height_i) << " z" << to_string(i) << " <= y" << to_string(j) << " + " << to_string(roundf(bound)*1) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " - " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                            constraint_file << "y" << to_string(i) << " - " << to_string(width_j) << " z" << to_string(j) << " - " << to_string(height_j) << " + " << to_string(height_j) << " z" << to_string(j) << " >= y" << to_string(j) << " - " << to_string(roundf(bound)*2) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " + " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n\n\n";
                        }
                    }
                }
                constraint_file.close();
            }
        }

        void hard_soft_nonoverlap()
        {
            if(hard_exists && soft_exists)
            {
                float width_h, height_h, g, c;
                constraint_file.open(output, ios::app);
                constraint_file << "/* Non-overlap constraints hard-soft */\n";

                for(short unsigned int i = 1; i <= num_hard_modules; i++)
                {
                    for(short unsigned int j = num_hard_modules+1; j <= num_total_modules; j++)
                    {
                        width_h = hard_module_width[i-1];
                        height_h = hard_module_height[i-1];

                        g = gradient[j-num_hard_modules-1];
                        c = intercept[j-num_hard_modules-1];

                        constraint_file << "x" << to_string(i) << " + " << to_string(height_h) << " z" << to_string(i) << " + " << to_string(width_h) << " - " << to_string(width_h) << " z" << to_string(i) << " <= x" << to_string(j) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " + " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                        constraint_file << "x" << to_string(i) << " - w" << to_string(j) << " >= x" << to_string(j) << " - " << to_string(roundf(bound)*1) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " - " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                        constraint_file << "y" << to_string(i) << " + " << to_string(width_h) << " z" << to_string(i) << " + " << to_string(height_h) << " - " << to_string(height_h) << " z" << to_string(i) << " <= y" << to_string(j) << " + " << to_string(roundf(bound)*1) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " - " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                        constraint_file << "y" << to_string(i) << " + " << to_string(-1*g) << " w" << to_string(j) << " - " << to_string(c) << " >= y" << to_string(j) << " - " << to_string(roundf(bound)*2) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " + " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n\n\n";
                    }
                }
                constraint_file.close();
            }
        }

        void soft_soft_nonoverlap()
        {
            if(soft_exists)
            {
                float g_i, c_i, g_j, c_j;
                
                constraint_file.open(output, ios::app);
                constraint_file << "/* Non-overlap constraints soft-soft */\n";

                for(short unsigned int i = num_hard_modules+1; i <= num_total_modules; i++)
                {
                    for(short unsigned int j = num_hard_modules+1; j <= num_total_modules; j++)
                    {
                        if(j <= i)
                        {
                            continue;
                        }
                        else
                        {
                            g_i = gradient[i-num_hard_modules-1];
                            c_i = intercept[i-num_hard_modules-1];
                            g_j = gradient[j-num_hard_modules-1];
                            c_i = intercept[j-num_hard_modules-1];

                            constraint_file << "x" << to_string(i) << " + w" << to_string(i) << " <= x" << to_string(j) << " + " << to_string(roundf(bound)*1) << " x" << to_string(i) << to_string(j) << " + " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                            constraint_file << "x" << to_string(i) << " - w" << to_string(j) << " >= x" << to_string(j) << " - " << to_string(roundf(bound)*1) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " - " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                            constraint_file << "y" << to_string(i) << " - " << to_string(-1*g_i) << " w" << to_string(i) << " + " << to_string(c_i) << " <= y" << to_string(j) << " + " << to_string(roundf(bound)*1) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " - " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n";
                            constraint_file << "y" << to_string(i) << " + " << to_string(-1*g_j) << " w" << to_string(j) << " - " << to_string(c_j) << " >= y" << to_string(j) << " - " << to_string(roundf(bound)*2) << " + " << to_string(roundf(bound)) << " x" << to_string(i) << to_string(j) << " + " << to_string(roundf(bound)) << " y" << to_string(i) << to_string(j) << ";\n\n\n";
                        }
                    }
                }
                constraint_file.close();
            }
        }

        void variable_type_constraints()
        {
            float soft_width_min, soft_width_max;
            constraint_file.open(output, ios::app);
            constraint_file << "/* variable type constraints */\n";

            auto soft_range = soft_module_dimension_range();

            for(short unsigned int i = 1; i <= num_total_modules; i++)
            {
                constraint_file << "x" << to_string(i) << " >= 0;\n";
                constraint_file << "y" << to_string(i) << " >= 0;\n";
            }
            for(short unsigned int i = num_hard_modules+1; i <= num_total_modules; i++)
            {
                soft_width_min = soft_module_width_range[i-num_hard_modules-1][0];
                soft_width_max = soft_module_width_range[i-num_hard_modules-1][1];

                constraint_file << "w" << to_string(i) << " >= " << to_string(soft_width_min) << ";\n";
                constraint_file << "w" << to_string(i) << " <= " << to_string(soft_width_max) << ";\n";
            }
            constraint_file << "\n\n";
            constraint_file.close();
        }

        void chip_width_constraints()
        {
                float width_h, height_h;
                constraint_file.open(output, ios::app);
                constraint_file << "/* chip width constraints */\n";

                for(short unsigned int i = 1; i <= num_hard_modules; i++)
                {
                    width_h = hard_module_width[i-1];
                    height_h = hard_module_height[i-1];

                    constraint_file << "x" << to_string(i) << " + " << to_string(width_h) << " - " << to_string(width_h) << " z" << to_string(i) << " + " << to_string(height_h) << " z" << to_string(i) << " <= Y;\n";
                }
                for(short unsigned int i = num_hard_modules+1; i <= num_total_modules; i++)
                {
                    constraint_file << "x" << to_string(i) << " + w" << to_string(i) << " <= Y;\n"; 
                }
                constraint_file << "\n\n";
                constraint_file.close();
        }

        void chip_height_constraints()
        {
                auto hard_dimensions = hard_module_dimension();
                auto linear_parameters = linear_approximation();
                float width_h, height_h, g, c;
                constraint_file.open(output, ios::app);
                constraint_file << "/* chip height constraints */\n";

                for(short unsigned int i = 1; i <= num_hard_modules; i++)
                {
                    width_h = hard_module_width[i-1];
                    height_h = hard_module_height[i-1];


                    constraint_file << "y" << to_string(i) << " + " << to_string(height_h) << " - " << to_string(height_h) << " z" << to_string(i) << " + " << to_string(width_h) << " z" << to_string(i) << " <= Y;\n";
                }

                for(short unsigned int i = num_hard_modules+1; i <= num_total_modules; i++)
                {
                    g = gradient[i-num_hard_modules-1];
                    c = intercept[i-num_hard_modules-1];

                    constraint_file << "y" << to_string(i) << " - " << to_string(-1*g) << " w" << to_string(i) << " + " << to_string(c) << " <= Y;\n";
                }
                constraint_file << "\n\n";
                constraint_file.close();
        }

        void binary_constraints()
        {
            constraint_file.open(output, ios::app);
            constraint_file << "/* variable type constraints */\n";
            constraint_file << "bin ";
            for(short int i = 1; i <= num_total_modules; i++)
            {
                for(short int j = 1; j <= num_total_modules; j++)
                {
                    if(j > i)
                    {
                        if(i == num_total_modules-1 && j == num_total_modules)
                        {
                            constraint_file << "x" << to_string(i) << to_string(j) << ";\n";
                        }
                        else
                        {
                            constraint_file << "x" << to_string(i) << to_string(j) << ", ";
                        }
                    }
                }
            }

            constraint_file << "bin ";
            for(short int i = 1; i <= num_total_modules; i++)
            {
                for(short int j = 1; j <= num_total_modules; j++)
                {
                    if(j > i)
                    {
                        if(i == num_total_modules-1 && j == num_total_modules)
                        {
                            constraint_file << "y" << to_string(i) << to_string(j) << ";\n";
                        }
                        else
                        {
                            constraint_file << "y" << to_string(i) << to_string(j) << ", ";
                        }
                    }
                }
            }

            constraint_file << "bin ";
            for(short int i = 1; i <= num_hard_modules; i++)
            {
                if(i == num_hard_modules)
                {
                    constraint_file << "z" << to_string(i) << ";\n";
                }
                else
                {
                    constraint_file << "z" << to_string(i) << ", ";
                }
            }
            constraint_file.close();
        }

        void create_ilp_file()
        {
            objective();
            hard_hard_nonoverlap();
            hard_soft_nonoverlap();
            soft_soft_nonoverlap();
            variable_type_constraints();
            chip_width_constraints();
            chip_height_constraints();
            binary_constraints();
        }
};

// Create constructor for GenerateProblem

GenerateProblem::GenerateProblem(string fname, int n, bool u=true)
    {
        file = fname;
        num_blocks = n;
        underestimation = u;
        ifstream f; // Read file contents into f
        f.open(file);
        string line;
        if(f.is_open())
        {
            while(getline(f, line)) // store contents of f in line
            {
                lines.push_back(line);
            }
            f.close();
        }

        tie(num_hard_modules, num_soft_modules) = total_modules();
        num_total_modules = num_hard_modules + num_soft_modules;
        tie(hard_module_width, hard_module_height) = hard_module_dimension();
        tie(area, min_aspect, max_aspect) = soft_module_properties();
        tie(soft_module_width_range, soft_module_height_range) = soft_module_dimension_range();
        tie(gradient, intercept) = linear_approximation();
        bound = upper_bound();
    }