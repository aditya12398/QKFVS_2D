#include <iostream>
#include <cmath>

using namespace std;

double *gauss_elimination(double matrix[5][(6)])
{
    int i, j, k;
    double *solution = new double[5];
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
    for (i = 1; i < 5; i++)
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

/*void construct_equation(int k, char var, double matrix[5][6])
{
    double sig_dxdx, sig_dydy, sig_dxdy;
    double sig_dxdxdx, sig_dydydy, sig_dxdxdy, sig_dxdydy;
    double sig_dxdxdxdx, sig_dydydydy, sig_dxdxdxdy, sig_dxdxdydy, sig_dxdydydy;
    double sig_dfdx, sig_dfdy, sig_dfdxdx, sig_dfdydy, sig_dfdxdy;
    double dx, dy, df;
    sig_dxdx = sig_dydy = sig_dxdy = sig_dxdxdx = sig_dydydy = sig_dxdxdy = sig_dxdydy = sig_dxdxdxdx = sig_dydydydy = sig_dxdxdxdy = sig_dxdxdydy = sig_dxdydydy = sig_dfdx = sig_dfdy = sig_dfdxdx = sig_dfdydy = sig_dfdxdy = 0;

    for (int i = 0; i < cell[k].nbhs; i++)
    {
        int p = cell[k].conn[i];
        dx = cell[k].cx - cell[p].cx;
        dy = cell[k].cy - cell[p].cy;
        if (var == 'u')
            df = cell[k].u1 - cell[p].u1;
        else if (var == 'v')
            df = cell[k].u2 - cell[p].u2;
        else if (var == 't')
            df = cell[k].tp - cell[p].tp;

        sig_dfdx += df * dx;
        sig_dfdy += df * dy;
        sig_dfdxdx += df * dx * dx;
        sig_dfdydy += df * dy * dy;
        sig_dfdxdy += df * dx * dy;

        sig_dxdx += dx * dx;
        sig_dxdy += dx * dy;
        sig_dydy += dy * dy;

        sig_dxdxdx += dx * dx * dx;
        sig_dxdxdy += dx * dx * dy;
        sig_dxdydy += dx * dy * dy;
        sig_dydydy += dy * dy * dy;

        sig_dxdxdxdx += dx * dx * dx * dx;
        sig_dxdxdxdy += dx * dx * dx * dy;
        sig_dxdxdydy += dx * dx * dy * dy;
        sig_dxdydydy += dx * dy * dy * dy;
        sig_dydydydy += dy * dy * dy * dy;
    }
    matrix[0][5] = sig_dfdx;
    matrix[1][5] = sig_dfdy;
    matrix[2][5] = sig_dfdxdx / 2;
    matrix[3][5] = sig_dfdydy / 2;
    matrix[4][5] = sig_dfdxdy;

    matrix[0][0] = sig_dxdx;
    matrix[1][1] = sig_dydy;
    matrix[2][2] = sig_dxdxdxdx / 4;
    matrix[3][3] = sig_dydydydy / 4;
    matrix[4][4] = sig_dxdxdydy;

    matrix[0][1] = matrix[1][0] = sig_dxdy;
    matrix[0][2] = matrix[2][0] = sig_dxdxdx / 2;
    matrix[0][3] = matrix[3][0] = sig_dxdydy / 2;
    matrix[0][4] = matrix[4][0] = sig_dxdxdy;

    matrix[1][2] = matrix[2][1] = sig_dxdxdy / 2;
    matrix[1][3] = matrix[3][1] = sig_dydydy / 2;
    matrix[1][4] = matrix[4][1] = sig_dxdydy;

    matrix[2][3] = matrix[3][2] = sig_dxdxdydy / 4;
    matrix[2][4] = matrix[4][2] = sig_dxdxdxdy / 2;

    matrix[3][4] = matrix[4][3] = sig_dxdydydy / 2;
}*/

int main()
{
    double matrix[5][(6)] = {{1, 2, 3, 4, 5, 6},{2, 3, 4, 5, 6, 7},{3, 4, 5, 6, 7, 8},{4, 5, 6, 7, 8, 9},{5, 6, 7, 8, 9, 10}};
    double *sol;
    sol = gauss_elimination(matrix);
    for (int i = 0; i < 5; i++)
        cout << sol[i] << ", ";
    cout << endl;
    delete[] sol;
    return 0;
}