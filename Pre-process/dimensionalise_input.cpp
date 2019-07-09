#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int max_cells, max_edges;

struct Edge
{
	int lcell, rcell;	  //Holds cell numbers on either side of the edge. Prefix l stands for left and r for right.
	double nx, ny, mx, my; //Hold point data. Prefix m stands for midpoint and n for edge normal.
	double length;		   //Holds edge length data
	char status;		   //Specifies whether the edge is Internal(f), Wall(w) or Outer Boundary(o).
};

struct Cell
{
    int *edge, noe, nbhs, *conn; //Holds enclosing edges and neighbour cell data.
    double area, cx, cy;         //Area of the cell and cell centre coordinates
    double rho, u1, u2, pr;      //Values of density, x - velocity, y - velocity, pressure and temperature of the cell
};

struct Cell *cell;
struct Edge * edge;
int inputdata(string fname)
{
    int temp;
    double rho, u1, u2, pr;
    ifstream infile(fname);
    //ifstream infile2("primitive-vector_implicit.dat");
    //Input Cell data
    if (!infile)
        return 1;

    infile >> max_edges;
	edge = new Edge[max_edges + 1];
	for (int k = 1; k <= max_edges; k++)
	{
		infile >> edge[k].mx >> edge[k].my >> edge[k].lcell >> edge[k].rcell >> edge[k].status >> edge[k].nx >> edge[k].ny >> edge[k].length;
	}
    infile >> max_cells;
    //Input Edge data
    cell = new Cell[max_cells + 1];
    for (int k = 1; k <= max_cells; k++)
    {
        infile >> cell[k].cx >> cell[k].cy >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr >> cell[k].area >> cell[k].noe;
        //infile2 >> temp >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr;
        //Set enclosing edges
        cell[k].edge = new int[cell[k].noe + 1];
        for (int r = 1; r <= cell[k].noe; r++)
            infile >> cell[k].edge[r];
        //Set neighbours
        infile >> cell[k].nbhs;
        cell[k].conn = new int[cell[k].nbhs];
        for (int r = 0; r < cell[k].nbhs; r++)
            infile >> cell[k].conn[r];
    }
    infile.close();
    return 0;
}

int fixdata(string fname)
{
    double rho_ref = 1.225, pr_ref, u_ref, tr_inf = 288.20; //Reference Quantities
    u_ref = sqrt(1.4 * 287 * tr_inf);
    pr_ref = rho_ref * u_ref * u_ref;

    ofstream outfile(fname);
    if (!outfile)
        return 1;

    outfile << max_edges;
	for (int k = 1; k <= max_edges; k++)
	{
		outfile <<endl << edge[k].mx << "\t" << edge[k].my <<"\t" << edge[k].lcell << "\t" << edge[k].rcell <<"\t" << edge[k].status << "\t" << edge[k].nx << "\t" << edge[k].ny << "\t" << edge[k].length;
	}
    outfile << endl << max_cells;
    //Input Edge data
    for (int k = 1; k <= max_cells; k++)
    {
        outfile << endl << cell[k].cx << "\t" << cell[k].cy << "\t" << (cell[k].rho * rho_ref) << "\t" << (cell[k].u1 * u_ref) << "\t" << (cell[k].u2 * u_ref) << "\t" << (cell[k].pr * pr_ref) << "\t" << cell[k].area << "\t" << cell[k].noe;
        //Set enclosing edges
        for (int r = 1; r <= cell[k].noe; r++)
            outfile << "\t" << cell[k].edge[r];
        //Set neighbours
        outfile << "\t" << cell[k].nbhs;
        for (int r = 0; r < cell[k].nbhs; r++)
            outfile << "\t" << cell[k].conn[r];
    }
    outfile.close();
    delete[] cell;
    return 0;
}

int main()
{
    int error = 0;
    error = inputdata("./Output/2order-input-data");
    if (error)
    {
        cout << "Error in reading file, please check the file or file title\nExiting...\n";
        exit(1);
    }
    error = 0;
    error = fixdata("./Output/naca0012_viscous_ogrid_4");
    if (error)
    {
        cout << "Error in writing file, please check the file or file title\nExiting...\n";
        exit(1);
    }
    return 0;
}