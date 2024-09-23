#include "src/solve.h"
#include "src/augment.h"

int main()
{
    int num_blocks{30}, sub_block_size{4}, num_augmentations, i;
    bool underestimation{true}, successive_augmentation{true}, visualize_superblock{true}, lp_solve{true};
    float runtime{2};
    string src_file_path;
    string spec_files_dir = "spec_files/";
    string result_dir = "results/" + to_string(num_blocks) + "/";
    string sa_files_dir = spec_files_dir + "successive_augmentation/" + to_string(num_blocks) + "/";
    string file = to_string(num_blocks) + "_block.ilp";
    vector<float> utilizations;
    if (successive_augmentation)
    {
        system(("rm -Rf " + sa_files_dir).c_str());
        system(("mkdir -p " + sa_files_dir).c_str());
        system(("rm -Rf " + result_dir).c_str());
        system(("mkdir -p " + result_dir).c_str());
    }

    Augment aug = Augment(file);
    aug.break_problem(sub_block_size);
    num_augmentations = int(ceil(float(num_blocks) / float(sub_block_size)));
    for(i=1; i<num_augmentations+1; i++)
    {
        cout << "\nOptimizing sub block " << to_string(i) << endl;
        src_file_path = sa_files_dir + to_string(num_blocks) + "_" + to_string(i) + ".ilp";
        SolveILP problem = SolveILP(src_file_path, sub_block_size, underestimation);
        vector<float>x_i, y_i, z_i, w_i, h_i;
        float Y;
        tie(Y, x_i, y_i, z_i, w_i, h_i) = problem.solve(runtime);
        string output_file_name = result_dir + to_string(num_blocks) + "_" + to_string(i) + ".txt";
        problem.export_results(Y, x_i, y_i, z_i, w_i, h_i, output_file_name);
        // system(("python src/visualize.py -f " + output_file_name + " --glob False --sa True -show True").c_str());
    }
    return 0;
}

