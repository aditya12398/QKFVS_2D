#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

int max_edges, max_cells, max_iters;
double Mach, aoa, cfl, limiter_const;
//double L = 1, rho_ref = 1.225, pr_ref, u_ref, tr_ref, tr_inf = 288.15; //Reference Quantities
double gma = 1.4, R = 287, Prandtl_number = 0.72; //Flow and material parameters
double residue, max_res;                          //RMS Residue and maximum residue in the fluid domain
int max_res_cell;                                 //Cell number with maximum residue
double pi = 4.0 * atan(1.0);

//Structure to store edge data
struct Edge
{
    int lcell, rcell;      //Holds cell numbers on either side of the edge. Prefix l stands for left and r for right.
    double nx, ny, mx, my; //Hold point data. Prefix m stands for midpoint and n for edge normal.
    double length;         //Holds edge length data
    char status;           //Specifies whether the edge is Internal(f), Wall(w) or Outer Boundary(o).
};
struct Cell
{
    int *edge, noe, nbhs, *conn; //Holds enclosing edges and neighbour cell data.
    double area, cx, cy;         //Area of the cell and cell centre coordinates
    double rho, u1, u2, pr, tp;  //Values of density, x - velocity, y - velocity, pressure and temperature of the cell
    double u1x, u2x, u1y, u2y;   //Velocity derivatives
    double tpx, tpy;             //Temperature derivatives
    double f;
    double flux[5];               //Kinetic Fluxes. Reference: See function `void KFVS_pos_flux(...)`
    double Upold[5], Upnew[5];    //Corresponds to void forward_sweep()
    double Unew[5], Uold[5];      //Corresponds to void backward_sweep()
    int alias, cell_with_alias_k; //Aliasing of cells for implicitisation
    double q[5], qx[5], qy[5];    //Entropy Variables. Reference: See Boltzmann Equations
};
struct Cell *cell;
struct Edge *edge;
void input_data()
{
    ifstream infile("ogrid_viscous_dim");
    ifstream infile2("fp_viscous");
    //Input Flow Parameters
    infile2 >> Mach >> aoa >> cfl >> max_iters >> limiter_const;
    //Input Edge data
    infile >> max_edges;
    edge = new Edge[max_edges + 1];
    for (int k = 1; k <= max_edges; k++)
    {
        infile >> edge[k].mx >> edge[k].my >> edge[k].lcell >> edge[k].rcell >> edge[k].status >> edge[k].nx >> edge[k].ny >> edge[k].length;
    }
    //Input Cell data
    infile >> max_cells;
    cell = new Cell[max_cells + 1];
    for (int k = 1; k <= max_cells; k++)
    {
        infile >> cell[k].cx >> cell[k].cy >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr >> cell[k].area >> cell[k].noe;
        cell[k].tp = cell[k].pr / (R * cell[k].rho);
        //Set enclosing edges
        cell[k].edge = new int[cell[k].noe + 1];
        for (int r = 1; r <= cell[k].noe; r++)
            infile >> cell[k].edge[r];
        //Set neighbours
        infile >> cell[k].nbhs;
        cell[k].conn = new int[cell[k].nbhs];
        for (int r = 0; r < cell[k].nbhs; r++)
            infile >> cell[k].conn[r];
        cell[k].f = 2 * cell[k].cx * cell[k].cx + 3 * cell[k].cy * cell[k].cy + 4 * cell[k].cx * cell[k].cy;
    }
    infile.close();
    infile2.close();
} //End of the function

double *gauss_elimination(double matrix[5][6])
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
    {
        for (j = i - 1; j >= 0; j--)
        {
            double d = (matrix[j][i] / matrix[i][i]);
            for (k = 5; k >= 0; k--)
            {
                matrix[j][k] -= d * matrix[i][k];
            }
        }
    }
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

void construct_equation(int k, char var, double matrix[5][6])
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
        dx = cell[p].cx - cell[k].cx;
        dy = cell[p].cy - cell[k].cy;
        double dist = sqrt(dx * dx + dy * dy); //Distance between cell and edge midpoint
        double w = 1 / pow(dist, 2.0);                 //Least Squares weight

        /*if (var == 'u')
            df = cell[p].u1 - cell[k].u1;
        else if (var == 'v')
            df = cell[p].u2 - cell[k].u2;
        else if (var == 't')
            df = cell[p].tp - cell[k].tp;*/
        df = cell[p].f - cell[k].f;
        sig_dfdx += w * df * dx;
        sig_dfdy += w * df * dy;
        sig_dfdxdx += w * df * dx * dx;
        sig_dfdydy += w * df * dy * dy;
        sig_dfdxdy += w * df * dx * dy;

        sig_dxdx += w * dx * dx;
        sig_dxdy += w * dx * dy;
        sig_dydy += w * dy * dy;

        sig_dxdxdx += w * dx * dx * dx;
        sig_dxdxdy += w * dx * dx * dy;
        sig_dxdydy += w * dx * dy * dy;
        sig_dydydy += w * dy * dy * dy;

        sig_dxdxdxdx += w * dx * dx * dx * dx;
        sig_dxdxdxdy += w * dx * dx * dx * dy;
        sig_dxdxdydy += w * dx * dx * dy * dy;
        sig_dxdydydy += w * dx * dy * dy * dy;
        sig_dydydydy += w * dy * dy * dy * dy;
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
}

int main()
{
    double matrix[5][6]; // = {{4, 1, 1, 1, 1, 8}, {1, 4, 1, 1, 1, 8}, {1, 1, 4, 1, 1, 8}, {1, 1, 1, 4, 1, 8}, {1, 1, 1, 1, 4, 8}};
    input_data();
    construct_equation(1, 'w', matrix);
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j <= 5; j++)
        {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
    double *sol;
    sol = gauss_elimination(matrix);
    for (int i = 0; i < 5; i++)
        cout << sol[i] << ", ";
    cout << endl;
    delete[] sol;
    return 0;
}