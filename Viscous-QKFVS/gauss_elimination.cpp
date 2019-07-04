#include <iostream>
#include <cmath>

using namespace std;
const int dim = 3;

double gauss_elimination(double matrix[dim][(dim + 1)], double solution[dim])
{
    int i, j, k;
    //double solution[dim];
    double swap;
    //Pivot the matrix for diagonal dominance
    for (i = 0; i < dim; i++)
    {
        solution[i] = 0;
        double max = fabs(matrix[i][i]);
        int pos = i;
        for (j = i; j < dim; j++)
        {
            if (fabs(matrix[j][i]) > max)
            {
                pos = j;
                max = fabs(matrix[j][i]);
            }
        }
        //Swap the elements
        for (k = 0; k <= dim; k++)
        {
            swap = matrix[i][k];
            matrix[i][k] = matrix[pos][k];
            matrix[pos][k] = swap;
        }
    }
    //Convert the matrix to a lower triangle matrix
    for (i = (dim - 1); i > 0; i--)
        for (j = i - 1; j >= 0; j--)
            for (k = dim; k >= 0; k--)
                matrix[j][k] -= (matrix[j][i] / matrix[i][i]) * matrix[i][k];

    solution[0] = matrix[0][dim] / matrix[0][0];
    for(i = 1; i < dim; i++)
    {
        solution[i] = matrix[i][dim];
        for (j = i - 1; j >= 0; j--)
        {
            solution[i] -= (solution[j] * matrix[i][j]);
        }
        solution[i] = solution[i] / matrix[i][i];
    }
    //return solution;
}

int main()
{
    double matrix[dim][(dim + 1)] = {{4, 1, 1, 10}, {1, 4, 1, 10}, {1, 1, 4, 10}};
    // = {{1, 2, 3, 4, 5, 6},{2, 3, 4, 5, 6, 7},{3, 4, 5, 6, 7, 8},{4, 5, 6, 7, 8, 9},{5, 6, 7, 8, 9, 10}};
    double sol[dim] = {0, 0, 0};
    double fx, fy, fxx, fyy, fxy;
    gauss_elimination(matrix, sol);
    /*fx = sol[0];
    fy = sol[1];
    fxx = sol[2];
    fyy = sol[3];
    fxy = sol[4];*/
    for (int i = 0; i < dim; i++)
        cout << sol[i] << ", ";
    cout <<endl;
    return 0;
}