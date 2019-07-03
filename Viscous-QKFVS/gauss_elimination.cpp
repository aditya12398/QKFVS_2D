#include <iostream>
#include <cmath>

using namespace std;

double * gauss_elimination(double matrix[5][6])
{
    int i, j, k;
    double solution[5];
    double swap;
    //Pivot the matrix for diagonal dominance
    for (i = 0; i < 5; i++)
    {
        solution[i] = 0;
        double max = fabs(matrix[i][i]);
        int pos = i;
        for (j = i; j < 5; j++)
        {
            if (fabs(matrix[j][i]) > max)
            {
                pos = j;
                max = fabs(matrix[j][i]);
            }
        }
        //Swap the elements
        for (k = 0; k <= 5; k++)
        {
            swap = matrix[i][k];
            matrix[i][k] = matrix[pos][k];
            matrix[pos][k] = swap;
        }
    }
    //Convert the matrix to a lower triangle matrix
    for (i = 4; i > 0; i--)
        for (j = i - 1; j >= 0; j--)
            for (k = 5; k >= 0; k--)
                matrix[j][k] -= (matrix[j][i] / matrix[i][i]) * matrix[i][k];

    solution[0] = matrix[0][5] / matrix[0][0];
    for(i = 1; i < 5; i++)
    {
        solution[i] = matrix[i][5];
        for (j = i - 1; j >= 0; j--)
        {
            solution[i] -= (solution[j] * matrix[i][j]);
        }
        solution[i] = solution[i] / matrix[i][i];
    }
    return solution;
}

int main()
{
    double matrix[5][6] = {{1, 2, 3, 4, 5, 6},{2, 3, 4, 5, 6, 7},{3, 4, 5, 6, 7, 8},{4, 5, 6, 7, 8, 9},{5, 6, 7, 8, 9, 10}};
    double * sol;
    double fx, fy, fxx, fyy, fxy;
    sol = gauss_elimination(matrix);
    fx = sol[0];
    fy = sol[1];
    fxx = sol[2];
    fyy = sol[3];
    fxy = sol[4];
    for (int i = 0; i < 5; i++)
        cout << sol[i] << ", ";
    cout <<endl;
    return 0;
}