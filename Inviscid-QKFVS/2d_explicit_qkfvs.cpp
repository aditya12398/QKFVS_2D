/*
This code is a solver for 2D Euler equations for compressible domain.
It uses `Kinetic Flux Vector Splitting(KFVS)` to calculate fluxes across the edges of a cell and conserve them. The related equations can be found in the Reference: "J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations" at location ../References/Kinetic_Flux_Vector_Splitting

Strongly Recommended: Use Visual Studio Code(text editor) while understanding the code. The comments are coupled appropriately with the functions and variables for easier understanding.
*/

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>

using namespace std;
int max_edges, max_cells, max_iters;
double Mach, aoa, cfl, limiter_const;
double residue, max_res; //RMS Residue and maximum residue in the fluid domain
int max_res_cell; //Cell number with maximum residue
double pi = 4.0 * atan(1.0);

//Structure to store edge data
struct Edge
{
	int lcell, rcell; //Holds cell numbers on either side of the edge. Prefix l stands for left and r for right.
	double nx, ny, mx, my; //Hold point data. Prefix m stands for midpoint and n for edge normal.
	double length; //Holds edge length data
	char status; //Specifies whether the edge is Internal(f), Wall(w) or Outer Boundary(o).
};
//Structure to store cell data
struct Cell
{
	int *edge, noe, nbhs, *conn; //Holds enclosing edges and neighbour cell data.
	double area, cx, cy; //Area of the cell and cell centre coordinates
	double rho, u1, u2, pr; //Values of density, x - velocity, y - velocity and pressure of the cell
	double flux[5]; //Kinetic Fluxes. Reference: See function `void KFVS_pos_flux(...)`
	double Unew[5], Uold[5]; //Corresponds to void backward_sweep()
	int alias, cell_with_alias_k; //Aliasing of cells for implicitisation
	double q[5], qx[5], qy[5]; //Entropy Variables. Reference: See Boltzmann Equations
};

struct Edge *edge;
struct Cell *cell;

//-------------------Main function starts--------------------------
int main(int arg, char *argv[])
{
	void input_data();
	void initial_conditions();
	void q_variables();
	void q_derivatives();
	void evaluate_flux();
	void state_update();
	void print_output();
	void write_tecplot();
	ofstream outfile("residue_explicit");

	double res_old;

	input_data();

	int T = 250;
	for (int t = 1; t <= max_iters; t++)
	{
		initial_conditions();
		q_variables();
		q_derivatives();

		evaluate_flux();
		state_update();

		if (t == 1)
			res_old = residue;
		residue = log10(residue / res_old);
		if ((residue) < -12) {
			cout << "Convergence criteria reached.\n";
			break;
		}
		cout << t << "\t" << residue << "\t" << max_res << "\t" << max_res_cell << endl;
		outfile << t << "\t" << residue << "\t" << max_res << "\t" << max_res_cell << endl;

		if (t == T)
		{
			T = T + 250;
			print_output();
			write_tecplot();
		}
	}
	print_output();
} //End of the main function

/*
This function reads flow and grid data from the files:
Grid data: 2order-input-data
Flow Parameters: flow-parameters-qkfvs
*/
void input_data()
{
	ifstream infile("2order-input-data");
	ifstream infile2("flow-parameters-qkfvs");

	infile2 >> Mach >> aoa >> cfl >> max_iters >> limiter_const;
	//Input Edge data
	infile >> max_edges;
	edge = new Edge[max_edges + 1];
	for (int k = 1; k <= max_edges; k++)
		infile >> edge[k].mx >> edge[k].my >> edge[k].lcell >> edge[k].rcell >> edge[k].status >> edge[k].nx >> edge[k].ny >> edge[k].length;
	//Input Cell data
	infile >> max_cells;
	cell = new Cell[max_cells + 1];
	for (int k = 1; k <= max_cells; k++)
	{
		infile >> cell[k].cx >> cell[k].cy >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr >> cell[k].area >> cell[k].noe;
		cell[k].edge = new int[cell[k].noe + 1];
		for (int r = 1; r <= cell[k].noe; r++)
			infile >> cell[k].edge[r];
		infile >> cell[k].nbhs;
		cell[k].conn = new int[cell[k].nbhs];
		for (int r = 0; r < cell[k].nbhs; r++)
			infile >> cell[k].conn[r];
	}
} //End of the function

//Fluxes are evaluated in this function
void evaluate_flux()
{
	//Function Prototypes
	void linear_reconstruction(double *, int, int);
	void KFVS_pos_flux(double *, double, double, double, double, double, double);
	void KFVS_neg_flux(double *, double, double, double, double, double, double);
	void KFVS_wall_flux(double *, double, double, double, double, double, double);
	void KFVS_outer_flux(double *, double, double, double, double, double, double);

	int lcell, rcell;
	double nx, ny, s;
	double Gp[5], Gn[5];
	double Gwall[5], Gout[5];
	double lprim[5], rprim[5];
	char status;

	double rhol, u1l, u2l, prl; //Left Cell Data
	double rhor, u1r, u2r, prr; //Right Cell Data

	for (int k = 1; k <= max_edges; k++)
	{
		lcell = edge[k].lcell;
		rcell = edge[k].rcell;
		nx = edge[k].nx;
		ny = edge[k].ny;
		s = edge[k].length;
		status = edge[k].status;

		if (status == 'f') //Flux across interior edges
		{
			linear_reconstruction(lprim, lcell, k);
			linear_reconstruction(rprim, rcell, k);

			rhol = lprim[1];
			u1l = lprim[2];
			u2l = lprim[3];
			prl = lprim[4];

			rhor = rprim[1];
			u1r = rprim[2];
			u2r = rprim[3];
			prr = rprim[4];

			KFVS_pos_flux(Gp, nx, ny, rhol, u1l, u2l, prl);
			KFVS_neg_flux(Gn, nx, ny, rhor, u1r, u2r, prr);

			for (int r = 1; r <= 4; r++)
			{
				cell[lcell].flux[r] = cell[lcell].flux[r] + s * (*(Gp + r) + *(Gn + r));
				cell[rcell].flux[r] = cell[rcell].flux[r] - s * (*(Gp + r) + *(Gn + r));
			}
		}

		if (status == 'w') //Flux across wall edge
		{
			if (lcell == 0)
			{
				nx = -nx;
				ny = -ny;
				lcell = rcell;
			}

			linear_reconstruction(lprim, lcell, k);

			rhol = lprim[1];
			u1l = lprim[2];
			u2l = lprim[3];
			prl = lprim[4];

			KFVS_wall_flux(Gwall, nx, ny, rhol, u1l, u2l, prl);

			for (int r = 1; r <= 4; r++)
				cell[lcell].flux[r] = cell[lcell].flux[r] + s * (*(Gwall + r));
		}

		if (status == 'o') //Flux across outer boundary edge
		{
			if (lcell == 0)
			{
				nx = -nx;
				ny = -ny;
				lcell = rcell;
			}

			linear_reconstruction(lprim, lcell, k);

			rhol = lprim[1];
			u1l = lprim[2];
			u2l = lprim[3];
			prl = lprim[4];

			KFVS_outer_flux(Gout, nx, ny, rhol, u1l, u2l, prl);

			for (int r = 1; r <= 4; r++)
				cell[lcell].flux[r] = cell[lcell].flux[r] + s * (*(Gout + r));
		}

	} //End of k loop
} //End of the flux function

//Expressions for the kfvs - fluxes

/*
Function to find the kfvs - positive flux
The following function uses the G Flux equations.
Reference: J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations; Page 7 of 32; Book Page Number 453
Reference Location: ../References/Kinetic_Flux_Vector_Splitting/
*/
void KFVS_pos_flux(double *G, double nx, double ny, double rho, double u1, double u2, double pr)
{
	double tx, ty;
	tx = ny;
	ty = -nx;
	double un = u1 * nx + u2 * ny;
	double ut = u1 * tx + u2 * ty;

	double beta = 0.5 * rho / pr;
	double S = un * sqrt(beta);
	double A = 0.5 * (1 + erf(S));
	double B = 0.5 * exp(-S * S) / sqrt(pi * beta);

	*(G + 1) = rho * (un * A + B);
	double temp1 = (pr + rho * un * un) * A + rho * un * B;
	double temp2 = rho * (un * ut * A + ut * B);
	*(G + 2) = -ty * temp1 + ny * temp2;
	*(G + 3) = tx * temp1 - nx * temp2;
	temp1 = 0.5 * rho * (un * un + ut * ut);
	*(G + 4) = (3.5 * pr + temp1) * un * A + (3.0 * pr + temp1) * B;
} //End of G-positive flux

/*
Function to find the kfvs - negative flux
The following function uses the G Flux equations.
Reference: J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations; Page 7 of 32; Book Page Number 453
Reference Location: ../References/Kinetic_Flux_Vector_Splitting/
*/
void KFVS_neg_flux(double *G, double nx, double ny, double rho, double u1, double u2, double pr)
{
	double tx, ty;
	tx = ny;
	ty = -nx;
	double un = u1 * nx + u2 * ny;
	double ut = u1 * tx + u2 * ty;

	double beta = 0.5 * rho / pr;
	double S = un * sqrt(beta);
	double A = 0.5 * (1 - erf(S));
	double B = 0.5 * exp(-S * S) / sqrt(pi * beta);

	*(G + 1) = rho * (un * A - B);
	double temp1 = (pr + rho * un * un) * A - rho * un * B;
	double temp2 = rho * (un * ut * A - ut * B);
	*(G + 2) = -ty * temp1 + ny * temp2;
	*(G + 3) = tx * temp1 - nx * temp2;
	temp1 = 0.5 * rho * (un * un + ut * ut);
	*(G + 4) = (3.5 * pr + temp1) * un * A - (3.0 * pr + temp1) * B;
} //End of G-negative flux

/*
Function to find the kfvs - wall boundary flux
The following function uses the G Flux equations.
Reference: J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations; Page 7 of 32; Book Page Number 453
Reference Location: ../References/Kinetic_Flux_Vector_Splitting/
*/
void KFVS_wall_flux(double *G, double nx, double ny, double rho, double u1, double u2, double pr)
{
	double tx, ty, un, ut;
	tx = ny;
	ty = -nx;
	ut = u1 * tx + u2 * ty;
	un = u1 * nx + u2 * ny;

	double beta = 0.5 * rho / pr;
	double Sw = un * sqrt(beta);
	double Aw = 0.5 * (1 + erf(Sw));
	double Bw = 0.5 * exp(-Sw * Sw) / sqrt(pi * beta);

	double temp = (pr + rho * un * un) * Aw + rho * un * Bw;
	temp = 2.0 * temp;
	*(G + 1) = 0.0;
	*(G + 2) = nx * temp;
	*(G + 3) = ny * temp;
	*(G + 4) = 0.0;

} //End of the wall flux function

/*
Function to find the kfvs - outer boundary flux
The following function uses the G Flux equations.
Reference: J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations; Page 7 of 32; Book Page Number 453
Reference Location: ../References/Kinetic_Flux_Vector_Splitting/
*/
void KFVS_outer_flux(double *Gout, double nx, double ny, double rho, double u1, double u2, double pr)
{
	void KFVS_pos_flux(double *, double, double, double, double, double, double);
	void KFVS_neg_flux(double *, double, double, double, double, double, double);

	double Gp[5];
	double Gn_inf[5];
	double rho_inf, u1_inf, u2_inf, pr_inf;

	double theta = aoa * pi / 180;

	rho_inf = 1.0;
	u1_inf = Mach * cos(theta);
	u2_inf = Mach * sin(theta);
	pr_inf = 1.0 / 1.4;

	KFVS_pos_flux(Gp, nx, ny, rho, u1, u2, pr);
	KFVS_neg_flux(Gn_inf, nx, ny, rho_inf, u1_inf, u2_inf, pr_inf);

	for (int r = 1; r <= 4; r++)
		*(Gout + r) = *(Gp + r) + *(Gn_inf + r);
} //End of the function

//Function to find the conserved vector from the given primitive vector
void prim_to_conserved(double *U, int k)
{
	double rho, u1, u2, pr, e;

	rho = cell[k].rho;
	u1 = cell[k].u1;
	u2 = cell[k].u2;
	pr = cell[k].pr;
	e = 2.5 * pr / rho + 0.5 * (u1 * u1 + u2 * u2);

	*(U + 1) = rho;
	*(U + 2) = rho * u1;
	*(U + 3) = rho * u2;
	*(U + 4) = rho * e;

} //End of the function

//Function to get the primitive vector from the given conserved vector
void conserved_to_primitive(double *U, int k)
{
	double U1 = *(U + 1);
	double U2 = *(U + 2);
	double U3 = *(U + 3);
	double U4 = *(U + 4);

	cell[k].rho = U1;
	double temp = 1 / U1;
	cell[k].u1 = U2 * temp;
	cell[k].u2 = U3 * temp;
	temp = U4 - (0.5 * temp) * (U2 * U2 + U3 * U3);
	cell[k].pr = 0.4 * temp;
} //End of the function

//Function to find the delt (delta t) for each cell
double func_delt(int k)
{
	double rho, u1, u2, pr;
	double area, temp1;

	int edges, e;
	double nx, ny, delt = 0.0;

	area = cell[k].area;
	u1 = cell[k].u1;
	u2 = cell[k].u2;
	rho = cell[k].rho;
	pr = cell[k].pr;

	edges = cell[k].noe;

	temp1 = 3.0 * sqrt(pr / rho);

	for (int r = 1; r <= edges; r++)
	{
		e = cell[k].edge[r];
		double length = edge[e].length;
		nx = edge[e].nx;
		ny = edge[e].ny;
		if (edge[e].rcell == k)
		{
			nx = -nx;
			ny = -ny;
		}

		double un = u1 * nx + u2 * ny;
		double temp2 = (fabs(un) + temp1);
		temp2 = length * temp2;
		delt = delt + temp2;
	}
	delt = 2.0 * area / delt;
	return (cfl * delt);
} //End of the function

//Function to give the initial conditions
void initial_conditions()
{
	double rho, u1, u2, pr, e, U[5];

	for (int k = 1; k <= max_cells; k++)
	{
		rho = cell[k].rho;
		u1 = cell[k].u1;
		u2 = cell[k].u2;
		pr = cell[k].pr;
		e = 2.5 * pr / rho + 0.5 * (u1 * u1 + u2 * u2); //(2.5) corresponds to (1 / (gamma - 1))

		U[1] = rho;
		U[2] = rho * u1;
		U[3] = rho * u2;
		U[4] = rho * e;

		for (int r = 1; r <= 4; r++)
		{
			cell[k].flux[r] = 0.0;
			cell[k].Uold[r] = U[r];
			cell[k].Unew[r] = U[r];
		}
	}
} //End of the function

void state_update()
{
	void prim_to_conserved(double *, int);
	void conserved_to_primitive(double *, int);
	residue = 0.0;
	double delt, area;
	max_res = 0.0;

	double U[5];

	for (int k = 1; k <= max_cells; k++)
	{
		prim_to_conserved(U, k);

		double temp = U[1];

		for (int r = 1; r <= 4; r++)
		{
			cell[k].Unew[r] = cell[k].Uold[r] - func_delt(k) * cell[k].flux[r] / cell[k].area;
			U[r] = cell[k].Unew[r];
			cell[k].Uold[r] = cell[k].Unew[r];
		}

		residue = residue + (temp - U[1]) * (temp - U[1]);
		temp = fabs(temp - U[1]);

		if (max_res < temp)
		{
			max_res = temp;
			max_res_cell = k;
		}
		conserved_to_primitive(U, k);
	}
	residue = residue / max_cells;
	residue = sqrt(residue);
} //End of the function
void write_tecplot()
{
    string title = "Solution_Tecplot_Explicit.dat";
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
//Function which prints the final primitive vector into the file "primitive-vector.dat"
void print_output()
{
	ofstream outfile("primitive-vector_explicit.dat");
	for (int k = 1; k <= max_cells; k++)
		outfile << k << "\t" << cell[k].rho << "\t" << cell[k].u1 << "\t" << cell[k].u2 << "\t" << cell[k].pr << endl;
} //End of the function

/*
Second order part of the code starts from here
Second order accuracy is achieved via linear reconstruction with q-variables
*/

//Function to calculate the q-variables for all the nodes
void q_variables()
{
	int i;
	double rho, u1, u2, pr, beta;

	for (i = 1; i <= max_cells; i++)
	{
		rho = cell[i].rho;
		u1 = cell[i].u1;
		u2 = cell[i].u2;
		pr = cell[i].pr;

		beta = 0.5 * rho / pr;

		cell[i].q[1] = log(rho) + (log(beta) * 2.5) - beta * (u1 * u1 + u2 * u2);
		cell[i].q[2] = 2 * beta * u1;
		cell[i].q[3] = 2 * beta * u2;
		cell[i].q[4] = -2 * beta;
	}
} //End of the function

//Function to convert q-varibales to primitive variables
void func_qtilde_to_prim_variables(double *prim, double *qtilde)
{
	double q1, q2, q3, q4, beta, u1, u2;

	q1 = *(qtilde + 1);
	q2 = *(qtilde + 2);
	q3 = *(qtilde + 3);
	q4 = *(qtilde + 4);

	beta = -q4 * 0.5;

	double temp = 0.5 / beta;

	*(prim + 2) = q2 * temp;
	*(prim + 3) = q3 * temp;

	u1 = *(prim + 2);
	u2 = *(prim + 3);

	double temp1 = q1 + beta * (u1 * u1 + u2 * u2);
	//double temp2 = temp1 - (log(beta)/(gamma-1));

	temp1 = temp1 - (log(beta) * 2.5);

	*(prim + 1) = exp(temp1);
	*(prim + 4) = *(prim + 1) * temp;
} //End of the function

//Function to find q-derivatives using the method of least squares
void q_derivatives()
{
	//double power = 2.0;
	for (int k = 1; k <= max_cells; k++)
	{
		//Initialise the sigma values
		double sig_dxdx = 0.0;
		double sig_dydy = 0.0;
		double sig_dxdy = 0.0;

		double sig_dxdq[5], sig_dydq[5];

		for (int r = 1; r <= 4; r++)
		{
			sig_dxdq[r] = 0.0;
			sig_dydq[r] = 0.0;
		}

		for (int j = 0; j < cell[k].nbhs; j++)
		{
			int p = cell[k].conn[j];
			double delx = cell[p].cx - cell[k].cx;
			double dely = cell[p].cy - cell[k].cy;
			double dist = sqrt(delx * delx + dely * dely);
			double w = 1 / pow(dist, 2.0);

			sig_dxdx = sig_dxdx + w * delx * delx;
			sig_dydy = sig_dydy + w * dely * dely;
			sig_dxdy = sig_dxdy + w * delx * dely;

			for (int r = 1; r <= 4; r++)
			{
				double delq = cell[p].q[r] - cell[k].q[r];

				sig_dxdq[r] = sig_dxdq[r] + w * delx * delq;
				sig_dydq[r] = sig_dydq[r] + w * dely * delq;
			}
		}

		double det = sig_dxdx * sig_dydy - sig_dxdy * sig_dxdy;
		for (int r = 1; r <= 4; r++)
		{
			double detx = sig_dydy * sig_dxdq[r] - sig_dxdy * sig_dydq[r];
			double dety = sig_dxdx * sig_dydq[r] - sig_dxdy * sig_dxdq[r];

			cell[k].qx[r] = detx / det;
			cell[k].qy[r] = dety / det;
		}
	} //End of k loop
} //End of the function

//Function to find the linear reconstruction values of primitive variables at the mid point of a given edge
void linear_reconstruction(double *prim, int CELL, int edg)
{
	void limiter(double *, double *, int);
	void func_qtilde_to_prim_variables(double *, double *);
	double delx, dely, phi[5], qtilde[5];

	delx = edge[edg].mx - cell[CELL].cx;
	dely = edge[edg].my - cell[CELL].cy;

	for (int r = 1; r <= 4; r++)
		qtilde[r] = cell[CELL].q[r] + delx * cell[CELL].qx[r] + dely * cell[CELL].qy[r];

	/*limiter(qtilde,phi,CELL);

	for(int r=1;r<=4;r++)
	qtilde[r] = cell[CELL].q[r] + phi[r]*(delx*cell[CELL].qx[r] + dely*cell[CELL].qy[r]);*/

	func_qtilde_to_prim_variables(prim, qtilde);
} //End of the function

//Function to find the limiter
void limiter(double *qtilde, double *phi, int k)
{
	double maximum(int, int);
	double minimum(int, int);
	double smallest_dist(int);
	double max_q, min_q, temp;
	double q, del_neg, del_pos;

	for (int r = 1; r <= 4; r++)
	{
		q = cell[k].q[r];

		del_neg = qtilde[r] - q;

		if (fabs(del_neg) <= 10e-6)
			phi[r] = 1.0;
		else
		{
			if (del_neg > 0.0)
			{
				max_q = maximum(k, r);
				del_pos = max_q - q;
			}
			else if (del_neg < 0.0)
			{
				min_q = minimum(k, r);
				del_pos = min_q - q;
			}
			double num, den, ds;

			ds = smallest_dist(k);
			double epsi = limiter_const * ds;
			epsi = pow(epsi, 3.0);

			num = (del_pos * del_pos) + (epsi * epsi);
			num = num * del_neg + 2.0 * del_neg * del_neg * del_pos;

			den = del_pos * del_pos + 2.0 * del_neg * del_neg;
			den = den + del_neg * del_pos + epsi * epsi;
			den = den * del_neg;

			temp = num / den;

			if (temp < 1.0)
				phi[r] = temp;
			else
				phi[r] = 1.0;
		}
	}
} //End of the function

double maximum(int k, int r)
{
	double max = cell[k].q[r];

	for (int j = 0; j < cell[k].nbhs; j++)
	{
		int p = cell[k].conn[j];

		if (cell[p].q[r] > max)
			max = cell[p].q[r];
	}
	return (max);
} //End of the function

double minimum(int k, int r)
{
	double min = cell[k].q[r];

	for (int j = 0; j < cell[k].nbhs; j++)
	{
		int p = cell[k].conn[j];

		if (cell[p].q[r] < min)
			min = cell[p].q[r];
	}
	return (min);
} //End of the function

//Finding the dels; the smallest distance measured over the neighbours of a given cell
double smallest_dist(int k)
{
	int p;
	double dx, dy, ds, min;

	for (int j = 0; j < cell[k].nbhs; j++)
	{
		p = cell[k].conn[j];

		dx = cell[p].cx - cell[k].cx;
		dy = cell[p].cy - cell[k].cy;

		ds = sqrt(dx * dx + dy * dy);

		if (j == 0)
			min = ds;
		if (ds < min)
			min = ds;
	}
	return (ds);
} //End of the function