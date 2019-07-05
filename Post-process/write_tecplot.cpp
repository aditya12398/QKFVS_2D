#include <fstream>
#include <iostream>
#include <cmath>
#include<string>

using namespace std;

int max_edges, max_cells, max_iters;
double Mach, aoa, cfl, limiter_const;
double residue, max_res; //RMS Residue and residue in the fluid domain
int max_res_cell;        //Cell number with maximum residue
double pi = 4.0 * atan(1.0);

//Structure to store edge data
struct Edge
{
    int lcell, rcell;      //Holds cell numbers on either side of the edge. Prefix l stands for left and r for right.
    double nx, ny, mx, my; //Hold point data. Prefix m stands for midpoint and n for edge normal.
    double length;         //Holds edge length data
    char status;           //Specifies whether the edge is Internal(f), Wall(w) or Outer Boundary(o).
};
//Structure to store cell data
struct Cell
{
    int *edge, noe, nbhs, *conn;  //Holds enclosing edges and neighbour cell data.
    double area, cx, cy;          //Area of the cell and cell centre coordinates
    double rho, u1, u2, pr;       //Values of density, x - velocity, y - velocity and pressure of the cell
    double flux[5];               //Kinetic Fluxes. Reference: See function `void KFVS_pos_flux(...)`
    double Upold[5], Upnew[5];    //Corresponds to void forward_sweep()
    double Unew[5], Uold[5];      //Corresponds to void backward_sweep()
    int alias, cell_with_alias_k; //Aliasing of cells for implicitisation
    double q[5], qx[5], qy[5];    //Entropy Variables. Reference: See Boltzmann Equations
};

struct Edge *edge;
struct Cell *cell;

void input_data()
{
    ifstream infile("2order-input-data");
    ifstream infile2("primitive-vector.dat");
    //Input Edge data
    int temp;
    infile >> max_edges;
    edge = new Edge[max_edges + 1];
    for (int k = 1; k <= max_edges; k++)
        infile >> edge[k].mx >> edge[k].my >> edge[k].lcell >> edge[k].rcell >> edge[k].status >> edge[k].nx >> edge[k].ny >> edge[k].length;
    //Input Cell data
    infile >> max_cells;
    cell = new Cell[max_cells + 1];
    for (int k = 1; k <= max_cells; k++)
    {
        infile2 >> temp >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr;
        infile >> cell[k].cx >> cell[k].cy >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr >> cell[k].area >> cell[k].noe;
        cell[k].edge = new int[cell[k].noe + 1];
        for (int r = 1; r <= cell[k].noe; r++)
            infile >> cell[k].edge[r];
        infile >> cell[k].nbhs;
        cell[k].conn = new int[cell[k].nbhs];
        for (int r = 0; r < cell[k].nbhs; r++)
            infile >> cell[k].conn[r];
        infile2 >> temp >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr;
    }
    infile.close();
    infile2.close();
} //End of the function

void write_tecplot()
{
    string title = "Solution_Tecplot.dat";
    ofstream tecplotfile("./" + title);
    tecplotfile << "TITLE: \"QKFVS Viscous Code - NAL\"\n";
    tecplotfile << "VARIABLES= \"X\", \"Y\", \"Density\", \"Pressure\", \"x-velocity\", \"y-velocity\", \"Velocity Magnitude\"\n";
    tecplotfile << "ZONE I= 160 J= 59 , DATAPACKING=BLOCK\n";
    for (int j = 1; j < 60; j++)
    {
        for (int i = 1; i <= 160; i++)
        {
            int k = (j - 1) * 160 + i;
            tecplotfile << cell[k].cx << "\n";
        }
    }
    tecplotfile << "\n";
    for (int j = 1; j < 60; j++)
    {
        for (int i = 1; i <= 160; i++)
        {
            int k = (j - 1) * 160 + i;
            tecplotfile << cell[k].cy << "\n";
        }
    }
    tecplotfile << "\n";
    for (int j = 1; j < 60; j++)
    {
        for (int i = 1; i <= 160; i++)
        {
            int k = (j - 1) * 160 + i;
            tecplotfile << cell[k].rho << "\n";
        }
    }
    tecplotfile << "\n";
    for (int j = 1; j < 60; j++)
    {
        for (int i = 1; i <= 160; i++)
        {
            int k = (j - 1) * 160 + i;
            tecplotfile << cell[k].pr << "\n";
        }
    }
    tecplotfile << "\n";
    for (int j = 1; j < 60; j++)
    {
        for (int i = 1; i <= 160; i++)
        {
            int k = (j - 1) * 160 + i;
            tecplotfile << cell[k].u1 << "\n";
        }
    }
    tecplotfile << "\n";
    for (int j = 1; j < 60; j++)
    {
        for (int i = 1; i <= 160; i++)
        {
            int k = (j - 1) * 160 + i;
            tecplotfile << cell[k].u2 << "\n";
        }
    }
    tecplotfile << "\n";
    for (int j = 1; j < 60; j++)
    {
        for (int i = 1; i <= 160; i++)
        {
            int k = (j - 1) * 160 + i;
            tecplotfile << sqrt(pow(cell[k].u1, 2) + pow(cell[k].u2, 2)) << "\n";
        }
    }
    tecplotfile << "\n";
    tecplotfile.close();
}

int main()
{
    input_data();
    write_tecplot();
}