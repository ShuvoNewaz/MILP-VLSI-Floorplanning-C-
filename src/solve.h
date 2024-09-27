#include"generate.h"
#include "fusion.h"
#include"utils.h"

using namespace mosek::fusion;
using namespace monty;


class SolveILP
{
    public:
        GenerateProblem problem;
        short unsigned int num_hard_modules, num_soft_modules, num_total_modules;

        // Hard module properties
        
        vector<float> hard_module_width, hard_module_height;

        // Soft module properties
        vector<float> soft_area, min_aspect, max_aspect;
        vector<vector<float>> soft_module_width_range, soft_module_height_range;
        vector<float> gradient, intercept;

        float bound;
        
        SolveILP(string, int, bool); // Class constructor is defined outside of class

        // Set up the MOSEK environment and variables

        Model::t M;
        vector<float> h;
        Variable::t X; // Contains N x's, N y's, N_hard z's, N_soft w's, Nc2 x_ij, Nc2 y_ij, 1 Y

        void create_constraints()
        {

            // Offset indices

            unsigned short int x_offset = num_total_modules;
            unsigned short int y_offset = x_offset + num_total_modules;
            unsigned short int z_offset = y_offset + num_hard_modules;
            unsigned short int w_offset = z_offset + num_soft_modules;
            unsigned short int x_ij_offset = w_offset + NcR(num_total_modules, 2);
            unsigned short int total_variables = 3 * num_total_modules + 2 * NcR(num_total_modules, 2) + 1;
            unsigned short int tri_flat;
            unsigned short int i, j;
            
            // Hard-hard non-overlap
            
            if(problem.hard_exists)
            {
                for(i=0; i<num_hard_modules; i++)
                {
                    for(j=0; j<num_hard_modules; j++)
                    {
                        if(j > i)
                        {
                            tri_flat = position(num_total_modules, i, j);
                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    i, j, y_offset + i, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    hard_module_height[i] - hard_module_width[i], -bound, -bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(-hard_module_width[i]));
                            }
                            
                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    i, j, y_offset + j, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    hard_module_width[j] - hard_module_height[j], -bound, bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::greaterThan(hard_module_width[j] - bound));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    x_offset + i, x_offset + j, y_offset + i, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    hard_module_width[i] - hard_module_height[i], -bound, bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(bound - hard_module_height[i]));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    x_offset + i, x_offset + j, y_offset + j, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    hard_module_height[j] - hard_module_width[j], -bound, -bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::greaterThan(hard_module_height[j] - 2 * bound));
                            }                           
                        }
                    }
                }
            }

            // Hard-soft non-overlap

            if(problem.hard_exists && problem.soft_exists)
            {
                for(i=0; i<num_hard_modules; i++)
                {
                    for(j=num_hard_modules; j<num_total_modules; j++)
                    {
                        if(j > i)
                        {
                            tri_flat = position(num_total_modules, i, j);

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    i, j, y_offset + i, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    hard_module_height[i] - hard_module_width[i], -bound, -bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(-hard_module_width[i]));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    i, j, y_offset + j, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    -1, -bound, bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::greaterThan(-bound));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    x_offset + i, x_offset + j, y_offset + i, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    hard_module_width[i] - hard_module_height[i], -bound, bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(bound - hard_module_height[i]));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    x_offset + i, x_offset + j, y_offset + j, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    -gradient[j-num_hard_modules], -bound, -bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::greaterThan(intercept[j-num_hard_modules] - 2 * bound));
                            }
                        }
                    }
                }
            }

            // Soft-soft non-overlap

            if(problem.soft_exists)
            {
                for(i=num_hard_modules; i<num_total_modules; i++)
                {
                    for(j=num_hard_modules; j<num_total_modules; j++)
                    {
                        if(j > i)
                        {
                            tri_flat = position(num_total_modules, i, j);

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    i, j, y_offset + i, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    1, -bound, -bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(0));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    i, j, y_offset + j, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    -1, -bound, bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::greaterThan(-bound));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    x_offset + i, x_offset + j, y_offset + i, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    gradient[i-num_hard_modules], -bound, bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(bound - intercept[i-num_hard_modules]));
                            }

                            {
                                auto Coefficients = nonOverlapConstraint(total_variables,
                                                                    x_offset + i, x_offset + j, y_offset + j, w_offset + tri_flat, x_ij_offset + tri_flat,
                                                                    -gradient[j-num_hard_modules], -bound, -bound);
                                M->constraint(Expr::dot(Coefficients, X), Domain::greaterThan(intercept[j-num_hard_modules] - 2 * bound));
                            }
                        }
                    }
                }
            }

            // Variable type constraints

            for(i=0; i<num_total_modules; i++)
            {
                greaterThanZeroConstraint(total_variables, i, M, X);
                greaterThanZeroConstraint(total_variables, x_offset + i, M, X);
            }
            
            for(i=0; i<num_soft_modules; i++)
            {
                double w_min = soft_module_width_range[i][0];
                double w_max = soft_module_width_range[i][1];
                softWidthConstraint(total_variables, z_offset+i, w_min, w_max, M, X);
            }

            for(i=0; i<NcR(num_total_modules, 2); i++)
            {
                binaryConstraint(total_variables, w_offset+i, M, X);
                binaryConstraint(total_variables, x_ij_offset+i, M, X);
            }

            if(problem.hard_exists)
            {
                for(i=0; i < num_hard_modules; i++)
                {
                    binaryConstraint(total_variables, y_offset+i, M, X);
                }

                for(i=0; i<num_hard_modules; i++)
                {
                    chipDimensionConstraint(total_variables, i, y_offset+i, hard_module_height[i], hard_module_width[i], M, X); // Chip width constraint
                    chipDimensionConstraint(total_variables, x_offset+i, y_offset+i, hard_module_width[i], hard_module_height[i], M, X); // Chip height constraint
                }
            }
            if(problem.soft_exists)
            {
                for(i=num_hard_modules; i<num_total_modules; i++)
                {
                    {
                        vector<double> A(total_variables);
                        A[i] = 1; // x_i = 1
                        A[y_offset + i] = 1; // w_i = 1
                        A[A.size() - 1] = -1; // Y = -1
                        auto Coefficients = new_array_ptr<double>(A);
                        M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(0));
                    }

                    {
                        vector<double> A(total_variables);
                        A[x_offset + i] = 1; // y_i = 1
                        A[y_offset + i] = gradient[i-num_hard_modules]; // w_i = m_i
                        A[A.size() - 1] = -1; // Y = -1
                        auto Coefficients = new_array_ptr<double>(A);
                        M->constraint(Expr::dot(Coefficients, X), Domain::lessThan(-intercept[i-num_hard_modules]));
                    }
                }
                

                
            }
            
        }

        tuple<float, vector<float>, vector<float>, vector<float>, vector<float>, vector<float>> solve(float run_time, bool sub_module)
        {
            vector<float> x_i, y_i, w_i, h_i, z_i;
            float Y;           
            unsigned short int total_variables = 3 * num_total_modules + 2 * NcR(num_total_modules, 2) + 1;

            // Populate the optimization model and variables

            M = new Model("Floorplan_Optimization");

            auto _M = finally([&]() { M->dispose(); });

            X = M->variable(total_variables, Domain::inRange(0, bound));
            X->slice(2 * num_total_modules, 2 * num_total_modules + num_hard_modules)->makeInteger(); // z_i are binary
            X->slice(3 * num_total_modules, total_variables - 1)->makeInteger(); // x_ij and y_ij are binary
            
            create_constraints();

            // Set max solution time

            M->setSolverParam("optimizerMaxTime", run_time);

            // Set max relative gap (to its default value)

            M->setSolverParam("mioTolRelGap", 1e-4);

            // Set max absolute gap (to its default value)

            M->setSolverParam("mioTolAbsGap", 0.1);

            vector<double> A(total_variables);
            A[A.size() - 1] = 1; // Set up the objective
            auto Coefficients = new_array_ptr<double>(A);
            M->objective("Objective", ObjectiveSense::Minimize, Expr::dot(Coefficients, X));
            M->writeTask("VLSI.ptf");
            // M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } ); // Uncomment to see details
            if(sub_module)
            {
                M->acceptedSolutionStatus(AccSolutionStatus::Feasible); // Otherwise the time limit will yield an error
            }
            
            M->solve();
            
            cout << M->getPrimalSolutionStatus();


            // Get all the parameter values
            
            auto sol = X->level();
            for(int i=0; i<num_total_modules; i++)
            {
                x_i.push_back((*sol)[i]);
                y_i.push_back((*sol)[i + num_total_modules]);
            }
            for(int i=0; i<num_hard_modules; i++)
            {
                z_i.push_back((*sol)[i + 2*num_total_modules]);
            }
            for(int i=0; i<num_soft_modules; i++)
            {
                w_i.push_back((*sol)[i + 2*num_total_modules + num_hard_modules]);
                h_i.push_back(gradient[i] * w_i[i] + intercept[i]);
            }
            
            Y = (*sol)[total_variables - 1];

            return make_tuple(Y, x_i, y_i, z_i, w_i, h_i);
        }

        float export_results(float Y, vector<float> x_i, vector<float>y_i, vector<float> z_i, vector<float> w_i, vector<float> h_i, vector<float> utilizations, string output_file_name)
        /*
            Exports the dimensions and coordinates of the optimized blocks to a text file. The results can
            later be read and plotted using Python.
        */
        {
            unsigned short int i;
            vector<float> W, H;
            float utilization = 0;

            for(i=0; i<num_hard_modules; i++)
            {
                W.push_back(hard_module_width[i]);
                H.push_back(hard_module_height[i]);
                utilization += W[i] * H[i];
            }
            for(i=0; i<num_soft_modules; i++)
            {
                W.push_back(w_i[i]);
                H.push_back(h_i[i]);
                utilization += W[i+num_hard_modules] * H[i+num_hard_modules];
            }
            for(i=0; i<utilizations.size(); i++)
            {
                utilization *= utilizations[i];
            }
            
            float chip_area = pow(Y, 2);
            utilization = (utilization / chip_area);

            ofstream output_file;
            output_file.open(output_file_name, ios::out);
            if(output_file.is_open())
            {   
                output_file << to_string(num_hard_modules) + "\n";
                output_file << to_string(num_soft_modules) + "\n";
                writeFile(output_file, x_i, num_total_modules);
                writeFile(output_file, y_i, num_total_modules);
                writeFile(output_file, W, num_total_modules);
                writeFile(output_file, H, num_total_modules);
                if(problem.hard_exists)
                {
                    writeFile(output_file, z_i, num_hard_modules);
                }
                else
                {
                    output_file << "0\n";
                }
                output_file << to_string(utilization) + "\n";
                output_file << to_string(Y) + "\n";
            }
            
            output_file.close();

            return utilization;
        }
};

// Create constructor for SolveILP

SolveILP::SolveILP(string file, int num_blocks, bool underestimation) : problem(file, num_blocks, underestimation)
    {
        num_hard_modules = problem.num_hard_modules;
        num_soft_modules = problem.num_soft_modules;
        num_total_modules = problem.num_total_modules;

        tie(hard_module_width, hard_module_height) = problem.hard_module_dimension();
        tie(soft_module_width_range, soft_module_height_range) = problem.soft_module_dimension_range();
        tie(soft_area, min_aspect, max_aspect) = problem.soft_module_properties();
        tie(gradient, intercept) = problem.linear_approximation();
        
        bound = problem.upper_bound();

        for(unsigned short int i=0; i<max<unsigned short int>(1, num_soft_modules); i++)
        {
            h.push_back(0);
        }
    }