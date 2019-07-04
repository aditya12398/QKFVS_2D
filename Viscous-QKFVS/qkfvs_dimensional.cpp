/*
This code is a solver for 2D Euler equations for compressible domain.
It uses `Kinetic Flux Vector Splitting(KFVS)` to calculate fluxes across the edges of a cell and conserve them. The related equations can be found in the Reference: "J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations" at location ../References/Kinetic_Flux_Vector_Splitting

Strongly Recommended: Use Visual Studio Code(text editor) while understanding the code. The comments are coupled appropriately with the functions and variables for easier understanding.
*/
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

int max_edges, max_cells, max_iters;
double Mach, aoa, cfl, limiter_const;
//double L = 1, rho_ref = 1.225, pr_ref, u_ref, tr_ref, tr_inf = 288.15; //Reference Quantities
double gma = 1.4, R = 287, Prandtl_number = 0.72; //Flow and material parameters
double residue, max_res;						  //RMS Residue and maximum residue in the fluid domain
int max_res_cell;								  //Cell number with maximum residue
double pi = 4.0 * atan(1.0);

//Structure to store edge data
struct Edge
{
	int lcell, rcell;	  //Holds cell numbers on either side of the edge. Prefix l stands for left and r for right.
	double nx, ny, mx, my; //Hold point data. Prefix m stands for midpoint and n for edge normal.
	double length;		   //Holds edge length data
	char status;		   //Specifies whether the edge is Internal(f), Wall(w) or Outer Boundary(o).
};
//Structure to store cell data
struct Cell
{
	int *edge, noe, nbhs, *conn;  //Holds enclosing edges and neighbour cell data.
	double area, cx, cy;		  //Area of the cell and cell centre coordinates
	double rho, u1, u2, pr, tp;   //Values of density, x - velocity, y - velocity, pressure and temperature of the cell
	double u1x, u2x, u1y, u2y;	//Velocity derivatives
	double tpx, tpy;			  //Temperature derivatives
	double flux[5];				  //Kinetic Fluxes. Reference: See function `void KFVS_pos_flux(...)`
	double Upold[5], Upnew[5];	//Corresponds to void forward_sweep()
	double Unew[5], Uold[5];	  //Corresponds to void backward_sweep()
	int alias, cell_with_alias_k; //Aliasing of cells for implicitisation
	double q[5], qx[5], qy[5];	//Entropy Variables. Reference: See Boltzmann Equations
};

struct Edge *edge;
struct Cell *cell;

//--------------------------Main function starts--------------------------
int main(int arg, char *argv[])
{
	void input_data();
	void aliasing();
	void initial_conditions();
	void q_variables();
	void q_derivatives();
	void f_derivatives();
	void evaluate_flux();
	void forward_sweep();
	void backward_sweep();
	void state_update();
	void print_output();

	ofstream outfile("./Output/Residue");

	double res_old;

	input_data();
	aliasing();
	int T = 250;
	for (int t = 1; t <= max_iters; t++)
	{
		initial_conditions();
		q_variables();
		q_derivatives();
		f_derivatives();
		evaluate_flux();
		forward_sweep();
		backward_sweep();
		state_update();

		if (t == 1)
			res_old = residue;
		residue = log10(residue / res_old);
		/*if (log10(max_res) < -12 && log10(residue) < -12)
		{
			cout << "Convergence criteria reached.\n";
			break;
		}*/
		cout << t << "\t" << residue << "\t" << max_res << "\t" << max_res_cell << endl;
		outfile << t << "\t" << residue << "\t" << max_res << "\t" << max_res_cell << endl;

		if (t == T)
		{
			T = T + 250;
			print_output();
		}
	}
	print_output();
} //End of the main function

/*
This function reads flow and grid data from the files:
Grid data: ogrid_viscous_dim
Flow Parameters: fp_viscous
*/
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
	}
	infile.close();
	infile2.close();
} //End of the function

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
}

void f_derivatives()
{
	for (int k = 1; k <= max_cells; k++)
	{
		double matrix[5][6], *solution;
		//Derivative of velocity
		construct_equation(k, 'u', matrix);
		solution = gauss_elimination(matrix);
		cell[k].u1x = solution[0];
		cell[k].u1y = solution[1];
		delete[] solution;
		construct_equation(k, 'v', matrix);
		solution = gauss_elimination(matrix);
		cell[k].u2x = solution[0];
		cell[k].u2y = solution[1];
		delete[] solution;
		construct_equation(k, 't', matrix);
		solution = gauss_elimination(matrix);
		cell[k].tpx = solution[0];
		cell[k].tpy = solution[1];
		delete[] solution;
	}
}

void viscous_flux(double *G, double nx, double ny, int e, char status)
{
	int lcell, rcell;
	double uref[3], u_dash[5], tp_dash[3];
	double tauxx, tauyy, tauxy, t_ref, mu, kappa;
	mu = 1.458E-6 * pow(t_ref, 1.5) / (t_ref + 110.4);
	kappa = mu * (gma * R) / (Prandtl_number * (gma - 1));
	lcell = edge[e].lcell;
	rcell = edge[e].rcell;
	double u_inf = sqrt(gma * R * 288.20);
	double theta = aoa * pi / 180;

	double w1 = sqrt(pow((cell[lcell].cx - edge[e].mx), 2) + pow((cell[lcell].cy - edge[e].my), 2)); //Distance between left cell center and edge midpoint
	double w2 = sqrt(pow((cell[rcell].cx - edge[e].mx), 2) + pow((cell[rcell].cy - edge[e].my), 2)); //Distance between right cell center and edge midpoint
	t_ref = (w1 * cell[lcell].tp + w2 * cell[rcell].tp) / (w1 + w2);
	uref[1] = (w1 * cell[lcell].u1 + w2 * cell[rcell].u1) / (w1 + w2);
	uref[2] = (w1 * cell[lcell].u2 + w2 * cell[rcell].u2) / (w1 + w2);

	u_dash[1] = (w1 * cell[lcell].u1x + w2 * cell[rcell].u1x) / (w1 + w2);
	u_dash[2] = (w1 * cell[lcell].u1y + w2 * cell[rcell].u1y) / (w1 + w2);
	u_dash[3] = (w1 * cell[lcell].u2x + w2 * cell[rcell].u2x) / (w1 + w2);
	u_dash[4] = (w1 * cell[lcell].u2y + w2 * cell[rcell].u2y) / (w1 + w2);

	tauxx = mu * (4 / 3 * u_dash[1] - 2 / 3 * u_dash[4]);
	tauyy = mu * (4 / 3 * u_dash[4] - 2 / 3 * u_dash[1]);
	tauxy = mu * (u_dash[2] + u_dash[3]);

	//Storing viscous fluxes
	if (status == 'w')
	{
		uref[1] = uref[2] = 0;
		tp_dash[1] = tp_dash[2] = 0;
	}
	else if (status == 'o')
	{
		uref[1] = u_inf * Mach * cos(theta);
		uref[2] = u_inf * Mach * sin(theta);
	}
	G[1] = 0;
	G[2] = -(tauxx * nx + tauxy * ny);
	G[3] = -(tauxy * nx + tauyy * ny);
	G[4] = -((uref[1] * tauxx + uref[2] * tauxy + kappa * tp_dash[1]) * nx) + ((uref[1] * tauxy + uref[2] * tauyy + kappa * tp_dash[2]) * ny);
	//delete[] edge_nbs;
}

//Fluxes are evaluated in this function
void evaluate_flux()
{
	//Function Prototypes
	void linear_reconstruction(double *, int, int);
	void KFVS_pos_flux(double *, double, double, double, double, double, double);
	void KFVS_neg_flux(double *, double, double, double, double, double, double);
	void KFVS_wall_flux(double *, double, double, double, double, double, double);
	void KFVS_outer_flux(double *, double, double, double, double, double, double);

	double u_dash[5];
	int lcell, rcell;
	double nx, ny, s;
	double Gp[5], Gn[5], Gd[5];
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
			viscous_flux(Gd, nx, ny, k, status);
			for (int r = 1; r <= 4; r++)
			{
				cell[lcell].flux[r] = cell[lcell].flux[r] + s * (*(Gp + r) + *(Gn + r) + *(Gd + r));
				cell[rcell].flux[r] = cell[rcell].flux[r] - s * (*(Gp + r) + *(Gn + r) + *(Gd + r));
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
			viscous_flux(Gd, nx, ny, k, status);
			for (int r = 1; r <= 4; r++)
				cell[lcell].flux[r] = cell[lcell].flux[r] + s * (*(Gwall + r) + *(Gd + r));
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
			viscous_flux(Gd, nx, ny, k, status);
			for (int r = 1; r <= 4; r++)
				cell[lcell].flux[r] = cell[lcell].flux[r] + s * (*(Gout + r) + *(Gd + r));
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
	*(G + 4) = ((gma / (gma - 1)) * pr + temp1) * un * A + (((gma + 1) / (2 * (gma - 1))) * pr + temp1) * B;
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
	*(G + 4) = ((gma / (gma - 1)) * pr + temp1) * un * A - (((gma + 1) / (2 * (gma - 1))) * pr + temp1) * B;
} //End of G-negative flux

/*
Function to find the kfvs - wall boundary flux
The following function uses the G Flux equations.
Reference: J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations; Page 7 of 32; Book Page Number 453
Reference Location: ../References/Kinetic_Flux_Vector_Splitting/
*/
void KFVS_wall_flux(double *G, double nx, double ny, double rho, double u1, double u2, double pr)
{
	double un, ut;
	ut = 0.0; //No slip boundary condition
	un = 0.0;

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
	double rho_inf, u1_inf, u2_inf, pr_inf, u_ref;
	u_ref = sqrt(gma * R * 288.20);

	double theta = aoa * pi / 180;
	rho_inf = 1.225;
	u1_inf = u_ref * Mach * cos(theta);
	u2_inf = u_ref * Mach * sin(theta);
	pr_inf = 101325.0;

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
	e = (1 / (gma - 1)) * pr / rho + 0.5 * (u1 * u1 + u2 * u2);

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
	cell[k].tp = cell[k].pr / (R * cell[k].rho);
} //End of the function

//Function to find the delt (delta t) for each cell
double func_delt(int k)
{
	double smallest_dist(int);
	double rho, u1, u2, pr;
	double area, temp1, smin;

	int edges, e;
	double nx, ny, delt_inv = 0.0, delt_visc = 0.0;

	area = cell[k].area;
	u1 = cell[k].u1;
	u2 = cell[k].u2;
	rho = cell[k].rho;
	pr = cell[k].pr;
	smin = smallest_dist(k);
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
		delt_inv = delt_inv + temp2;
	}
	delt_inv = 2.0 * area / delt_inv;
	delt_visc = rho * ((gma * R) / (gma - 1)) * smin * smin / (0.096 * gma);
	double delt = 1 / ((1 / delt_inv) + (1 / delt_visc));
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
		e = (1 / (gma - 1)) * pr / rho + 0.5 * (u1 * u1 + u2 * u2);

		U[1] = rho;
		U[2] = rho * u1;
		U[3] = rho * u2;
		U[4] = rho * e;

		for (int r = 1; r <= 4; r++)
		{
			cell[k].flux[r] = 0.0;
			cell[k].Uold[r] = U[r];
			cell[k].Unew[r] = U[r];
			cell[k].Upold[r] = U[r];
			cell[k].Upnew[r] = U[r];
		}
	}
} //End of the function

void state_update()
{
	void prim_to_conserved(double *, int);
	void conserved_to_primitive(double *, int);
	residue = 0.0;
	max_res = 0.0;
	double U[5];

	for (int k = 1; k <= max_cells; k++)
	{
		prim_to_conserved(U, k);

		double temp = U[1];
		for (int r = 1; r <= 4; r++)
		{
			U[r] = cell[k].Unew[r];
			cell[k].Upold[r] = cell[k].Upnew[r];
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

//Function writes the data in supported tecplot format
void write_tecplot()
{
	static int count = 1;
	string title = to_string(count) + "_Solution_Tecplot.dat";
	ofstream tecplotfile("./Output/" + title);
	tecplotfile << "TITLE: \"QKFVS Viscous Code - NAL\"\n";
	tecplotfile << "VARIABLES= \"X\", \"Y\", \"Density\", \"Pressure\", \"x-velocity\", \"y-velocity\", \"Velocity Magnitude\", \"Temperature\"\n";
	tecplotfile << "ZONE I= 298 J= 179, DATAPACKING=BLOCK\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << cell[k].cx << "\n";
		}
	}
	tecplotfile << "\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << cell[k].cy << "\n";
		}
	}
	tecplotfile << "\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << cell[k].rho << "\n";
		}
	}
	tecplotfile << "\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << cell[k].pr << "\n";
		}
	}
	tecplotfile << "\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << cell[k].u1 << "\n";
		}
	}
	tecplotfile << "\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << cell[k].u2 << "\n";
		}
	}
	tecplotfile << "\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << sqrt(pow(cell[k].u1, 2) + pow(cell[k].u2, 2)) << "\n";
		}
	}
	tecplotfile << "\n";
	for (int j = 1; j < 180; j++)
	{
		for (int i = 1; i <= 298; i++)
		{
			int k = (j - 1) * 298 + i;
			tecplotfile << cell[k].tp << "\n";
		}
	}
	tecplotfile << "\n";
	tecplotfile.close();
	count++;
}

//Function which prints the final primitive vector into the file "primitive-vector.dat"
void print_output()
{
	ofstream outfile("./Output/Primitive-Vector.dat");

	for (int k = 1; k <= max_cells; k++)
	{
		outfile << k << "\t" << cell[k].rho << "\t" << cell[k].u1 << "\t" << cell[k].u2 << "\t" << cell[k].pr << endl;
	}
	write_tecplot();
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

		cell[i].q[1] = log(rho) + (log(beta) * (1 / (gma - 1))) - beta * (u1 * u1 + u2 * u2);
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
	//double temp2 = temp1 - (log(beta)/(gma-1));
	temp1 = temp1 - (log(beta) * (1 / (gma - 1)));

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

	//limiter(qtilde, phi, CELL);

	/*for (int r = 1; r <= 4; r++)
		qtilde[r] = cell[CELL].q[r] + phi[r] * (delx * cell[CELL].qx[r] + dely * cell[CELL].qy[r]);*/

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

/*Implicit part of the code starts from here*/
//Function to find the forward sweep
void forward_sweep()
{
	double get_D_inv(int);
	void get_LdU(int, double *);
	double D_inv, LdU[5], temp;

	for (int k = 1; k <= max_cells; k++)
	{
		int CELL = cell[k].cell_with_alias_k;

		D_inv = get_D_inv(CELL);

		get_LdU(CELL, LdU);

		for (int r = 1; r <= 4; r++)
		{
			temp = cell[CELL].flux[r] + LdU[r];
			temp = -D_inv * temp;

			cell[CELL].Upnew[r] = cell[CELL].Upold[r] + temp;
		}
	}
} //End of the function

//Function to find the backward sweep
void backward_sweep()
{
	double get_D_inv(int);
	void get_UdU(int, double *);
	double D_inv, UdU[5], temp;

	for (int k = max_cells; k >= 1; k--)
	{
		int CELL = cell[k].cell_with_alias_k;

		D_inv = get_D_inv(CELL);

		get_UdU(CELL, UdU);

		for (int r = 1; r <= 4; r++)
		{
			temp = cell[CELL].Upnew[r] - cell[CELL].Upold[r];
			temp = temp - D_inv * UdU[r];

			cell[CELL].Unew[r] = cell[CELL].Uold[r] + temp;
		}
	}
} //End of the function
//Lower Upper Symmetric Gauss Seidel splitting in lower triangular matrix
void get_LdU(int k, double *L)
{
	void get_delG(int, double, double, double *, char);

	double length, nx, ny, rho, u1, u2, pr, a;
	double dG[5];

	for (int r = 1; r <= 4; r++)
		L[r] = 0.0;

	for (int r = 1; r <= cell[k].noe; r++)
	{
		int e = cell[k].edge[r];
		length = edge[e].length;
		int nbh = edge[e].rcell;

		nx = edge[e].nx;
		ny = edge[e].ny;

		if (nbh == k)
		{
			nx = -nx;
			ny = -ny;
			nbh = edge[e].lcell;
		}

		int alias_of_k = cell[k].alias;
		int alias_of_nbh = cell[nbh].alias;

		if (nbh != 0 && alias_of_nbh < alias_of_k)
		{
			rho = cell[nbh].rho;
			u1 = cell[nbh].u1;
			u2 = cell[nbh].u2;
			pr = cell[nbh].pr;
			a = sqrt(gma * pr / rho);

			double sr = fabs(u1 * nx + u2 * ny) + a;

			get_delG(nbh, nx, ny, dG, 'f');

			for (int j = 1; j <= 4; j++)
			{
				double temp = cell[nbh].Upnew[j] - cell[nbh].Upold[j];
				temp = sr * temp;
				temp = dG[j] - temp;
				temp = temp * length;
				L[j] = L[j] + 0.5 * temp;
			}
		}
	}
} //End of the function
//Lower Upper Symmetric Gauss Seidel splitting in upper triangular matrix
void get_UdU(int k, double *U)
{
	void get_delG(int, double, double, double *, char);

	double rho, u1, u2, pr, a, length;
	double nx, ny, dG[5];

	for (int r = 1; r <= 4; r++)
		U[r] = 0.0;

	for (int r = 1; r <= cell[k].noe; r++)
	{
		int e = cell[k].edge[r];
		length = edge[e].length;
		int nbh = edge[e].rcell;

		nx = edge[e].nx;
		ny = edge[e].ny;

		if (nbh == k)
		{
			nx = -nx;
			ny = -ny;
			nbh = edge[e].lcell;
		}

		int alias_of_k = cell[k].alias;
		int alias_of_nbh = cell[nbh].alias;

		if (nbh != 0 && alias_of_nbh > alias_of_k)
		{
			rho = cell[nbh].rho;
			u1 = cell[nbh].u1;
			u2 = cell[nbh].u2;
			pr = cell[nbh].pr;
			a = sqrt(gma * pr / rho);

			double sr = fabs(u1 * nx + u2 * ny) + a;

			get_delG(nbh, nx, ny, dG, 'b');

			for (int j = 1; j <= 4; j++)
			{
				double temp = cell[nbh].Unew[j] - cell[nbh].Uold[j];
				temp = sr * temp;
				temp = dG[j] - temp;
				temp = temp * length;
				U[j] = U[j] + 0.5 * temp;
			}
		}
	}
} //End of the function
//Lower Upper Symmetric Gauss Seidel splitting. D stands for the diagonal matrix
double get_D_inv(int k)
{
	double func_delt(int);
	double area, delt, temp = 0.0;
	int e, nbh;
	double length, rho, u1, u2, pr, a;
	double nx, ny, h, delx, dely;
	area = cell[k].area;
	delt = func_delt(k);

	rho = cell[k].rho;
	u1 = cell[k].u1;
	u2 = cell[k].u2;
	pr = cell[k].pr;
	a = sqrt(gma * pr / rho); //a is the nondimensionalised local speed of sound;

	for (int r = 1; r <= cell[k].noe; r++)
	{
		e = cell[k].edge[r];
		length = edge[e].length;
		nbh = edge[e].rcell;
		delx = (nbh == k ? cell[edge[e].lcell].cx : cell[nbh].cx) - cell[k].cx;
		dely = (nbh == k ? cell[edge[e].lcell].cy : cell[nbh].cy) - cell[k].cy;
		h = sqrt(pow(delx * nx, 2) + pow(dely * ny, 2));

		nx = edge[e].nx;
		ny = edge[e].ny;

		if (nbh == k)
		{
			nx = -nx;
			ny = -ny;
		}
		double t_ref = cell[k].tp;
		double mu = 1.458E-6 * pow(t_ref, 1.5) / (t_ref + 110.4);

		temp = temp + (fabs(u1 * nx + u2 * ny) + a + (mu / (rho * h))) * length;
	}
	temp = (area / delt) + 0.5 * temp;
	return (1.0 / temp);
} //End of the function

//Function to get delG
void get_delG(int k, double nx, double ny, double *dG, char status)
{
	void get_G(double, double, double *, double *);
	double U[5], Gold[5], Gnew[5];

	if (status == 'f')
	{
		for (int r = 1; r <= 4; r++)
			U[r] = cell[k].Upold[r];
		get_G(nx, ny, U, Gold);

		for (int r = 1; r <= 4; r++)
			U[r] = cell[k].Upnew[r];
		get_G(nx, ny, U, Gnew);
	}
	else if (status == 'b')
	{
		for (int r = 1; r <= 4; r++)
			U[r] = cell[k].Uold[r];
		get_G(nx, ny, U, Gold);

		for (int r = 1; r <= 4; r++)
			U[r] = cell[k].Unew[r];
		get_G(nx, ny, U, Gnew);
	}

	for (int r = 1; r <= 4; r++)
		dG[r] = Gnew[r] - Gold[r];
} //End of the function

//Function to get the G_normal vector
void get_G(double nx, double ny, double *U, double *G)
{
	double rho, u1, u2, pr;
	double U1 = *(U + 1);
	double U2 = *(U + 2);
	double U3 = *(U + 3);
	double U4 = *(U + 4);
	//Get primitive forms of variables
	rho = U1;
	u1 = U2 / rho;
	u2 = U3 / rho;
	pr = 0.4 * (U4 - (0.5 / rho) * (U2 * U2 + U3 * U3));
	//Calculate convective fluxes
	G[1] = rho * (u1 * nx + u2 * ny);
	G[2] = ((pr + rho * u1 * u1) * nx) + (rho * u1 * u2 * ny);
	G[3] = ((pr + rho * u2 * u2) * ny) + (rho * u1 * u2 * nx);
	double e = (1 / (gma - 1)) * pr / rho + (0.5 * (u1 * u1 + u2 * u2));
	G[4] = (pr + rho * e) * (u1 * nx + u2 * ny);
} //End of the function

//Aliasing the cell numbers
void aliasing()
{
	int count = 1;
	int e, rcell, lcell;

	for (int k = 2 /*1*/; k <= max_cells; k++)
	{
		cell[k].alias = 0;			   //Alias number of cell k
		cell[k].cell_with_alias_k = 0; //Cell number with alias = k
	}

	cell[1].alias = 1;
	cell[1].cell_with_alias_k = 1;

	for (int k = 1; k <= max_cells; k++)
	{
		int CELL = cell[k].cell_with_alias_k;
		for (int r = 1; r <= cell[CELL].noe; r++)
		{
			e = cell[CELL].edge[r];

			rcell = edge[e].rcell;
			lcell = edge[e].lcell;

			if (rcell == CELL)
				rcell = lcell;
			if (rcell != 0)
			{
				if (cell[rcell].alias == 0)
				{
					//count = count + 1;
					cell[rcell].alias = ++count;
					cell[count].cell_with_alias_k = rcell;
				}
			}
		}
	}
} //End of the function
