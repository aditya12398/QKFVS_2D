#include<fstream>
#include<iostream>

using namespace std;

int main()
{
    int imax, jmax;
    int i, j;
    double x, y, z;
    ifstream infile("./gridpro.dat");
    ofstream outfile("./gridout.dat");
    infile >> imax >> jmax;
    --imax;
    for (i = 1; i <= imax; i++)
    {
        for (j = 1; j <= jmax; j++)
        {
            /*if (j == jmax)
                continue;*/
            infile >> x >> y >> z;
            outfile << x << "\t" << y << "\t" << i << "\t" << j << "\n";
        }
    }
    infile.close();
    outfile.close();
}