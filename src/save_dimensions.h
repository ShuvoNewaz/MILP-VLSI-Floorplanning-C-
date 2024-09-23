#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>

using namespace std;

void writeHard(string file_name, int num_augmentations)
{
    ofstream f(file_name);
    f << "hard - " << num_augmentations << "\n";
    f.close();
}

// void save_augmented_dimensions(string file_name)
// {
//     ofstream f(file_name);
//     f << "hard - " ;
//     f.close();
// }

void save_augmented_dimensions(string file_name, vector<float> bounds)
{
    ofstream f;
    f.open(file_name, ios_base::app);
    for(int i=0; i<bounds.size(); i++)
    {
        f << bounds[i] << "," << bounds[i] << "\n";
    }
    f.close();
}