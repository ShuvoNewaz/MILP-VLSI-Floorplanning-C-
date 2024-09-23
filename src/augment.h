#include<iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

int minimum(int a, int b)
{
    return a * (a < b) + b * (b <= a);
}

class Augment
{
    public:
        
        bool hard_exists, soft_exists;
        unsigned short int num_hard_modules, num_soft_modules;
        string file, spec_file;
        string spec_files_dir, sa_files_dir, sa_file_prefix;
        int num_blocks;
        vector<string> lines;

        Augment(string);

        tuple<int, int> total_modules()
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

        void break_problem(int sub_block_size)
        {
            unsigned short int num_total_modules;
            vector<float> hard_module_width, hard_module_height, area, min_aspect, max_aspect;

            tie(num_hard_modules, num_soft_modules) = total_modules();
            num_total_modules = num_hard_modules + num_soft_modules;
            tie(hard_module_width, hard_module_height) = hard_module_dimension();
            tie(area, min_aspect, max_aspect) = soft_module_properties();

            if(num_total_modules > sub_block_size)
            {
                unsigned short int soft_count {0}, i, j, k, num_subblocks, modules_in_subblock, soft_left;
                short int hard_left;
                num_subblocks = int(ceil(float(num_blocks) / float(sub_block_size)));
                soft_count = 0;
                vector<float> area_clipped, min_aspect_clipped, max_aspect_clipped;
                ofstream g;

                for(i=0; i<num_subblocks; i++)
                {
                    modules_in_subblock = minimum(sub_block_size, num_total_modules - i * sub_block_size);
                    hard_left = num_hard_modules - sub_block_size*i;
                    if(hard_left >= modules_in_subblock)
                    {
                        g.open(sa_file_prefix + "/" + to_string(num_blocks) + "_" + to_string(i+1) + ".ilp", ios::out);
                        g << "hard - " + to_string(modules_in_subblock) + "\n";
                        g.close();
                        g.open(sa_file_prefix + "/" + to_string(num_blocks) + "_" + to_string(i+1) + ".ilp", ios::app);
                        
                        for(j=0; j<modules_in_subblock; j++)
                        {
                            g << to_string(hard_module_width[i*modules_in_subblock+j]) + "," + to_string(hard_module_height[i*modules_in_subblock+j]) + "\n";
                        }
                        g.close();
                    }
                    else if(modules_in_subblock > hard_left && hard_left > 0)
                    {
                        g.open(sa_file_prefix + "/" + to_string(num_blocks) + "_" + to_string(i+1) + ".ilp", ios::out);
                        g << "hard - " + to_string(hard_left) + "\n";
                        g.close();
                        g.open(sa_file_prefix + "/" + to_string(num_blocks) + "_" + to_string(i+1) + ".ilp", ios::app);
                        for(j=0; j<hard_left; j++)
                        {
                            g << to_string(hard_module_width[i*modules_in_subblock+j]) + "," + to_string(hard_module_height[i*modules_in_subblock+j]) + "\n";
                        }
                        if(num_soft_modules > 0)
                        {
                            g << "\nsoft - " + to_string(modules_in_subblock - hard_left) + "\n";
                            for(j=0; j<modules_in_subblock - hard_left; j++)
                            {
                                g << to_string(area[j]) + "," + to_string(min_aspect[j]) + ","  + to_string(max_aspect[j])  + "\n";
                            }
                            g.close();
                            for(k=j; k<num_soft_modules; k++)
                            {
                                area_clipped.push_back(area[k]);
                                min_aspect_clipped.push_back(min_aspect[k]);
                                max_aspect_clipped.push_back(max_aspect[k]);
                            }
                            area = area_clipped;
                            min_aspect = min_aspect_clipped;
                            max_aspect = max_aspect_clipped;
                        }
                        soft_left = num_soft_modules - (modules_in_subblock - hard_left);
                    }
                    else if(hard_left <= 0)
                    {
                        soft_left = num_soft_modules - modules_in_subblock * soft_count;
                        g.open(sa_file_prefix + "/" + to_string(num_blocks) + "_" + to_string(i+1) + ".ilp", ios::out);
                        g << "soft - " + to_string(modules_in_subblock) + "\n";
                        g.close();
                        g.open(sa_file_prefix + "/" + to_string(num_blocks) + "_" + to_string(i+1) + ".ilp", ios::app);
                        if(soft_left > 0)
                        {
                            for(j=0; j<modules_in_subblock; j++)
                            {
                                g << to_string(area[soft_count*modules_in_subblock+j]) + "," + to_string(min_aspect[soft_count*modules_in_subblock+j]) + "," + to_string(max_aspect[soft_count*modules_in_subblock+j]) + "\n";
                            }
                            g.close();
                            soft_count ++;
                        }
                    }
                }
            }

        }
};

Augment::Augment(string fname)
{
    file = fname;
    num_blocks = stoi(file.substr(0, file.find("_")));
    spec_files_dir = "spec_files";
    sa_files_dir = spec_files_dir + "/successive_augmentation";
    sa_file_prefix = sa_files_dir + "/" + to_string(num_blocks);
    // system("mkdir " + sa_files_dir + " -p");
    spec_file = spec_files_dir + "/" + file;

    ifstream f; // Read file contents into f
    f.open(spec_file);
    string line;
    if(f.is_open())
    {
        while(getline(f, line)) // store contents of f in line
        {
            lines.push_back(line);
        }
        f.close();
    }
}