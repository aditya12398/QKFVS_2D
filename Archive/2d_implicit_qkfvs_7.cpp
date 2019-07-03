/*
This code is a solver for 2D Euler equations for compressible domain.
It uses `Kinetic Flux Vector Splitting(KFVS)` to calculate fluxes across the edges of a cell and conserve them. The related equations can be found in the Reference: "J. C. Mandal, S. M. Deshpande - Kinetic Flux Vector Splitting For Euler Equations" at location ../References/Kinetic_Flux_Vector_Splitting

Strongly Recommended: Use Visual Studio Code(text editor) while understanding the code. The comments are coupled appropriately with the functions and variables for easier understanding.
*/
#include <fstream>
#include <iostream>
#include <cmath>

int max_edges, max_cells, max_iters;
double Mach, aoa, cfl, limiter_const;
double L = 4.45136E-4, rho_ref = 1.225, pr_ref, u_ref, tr_ref, tr_inf = 288.15; //Reference Quantities
double gma = 1.4, mu = 1.802E-5, R = 287, Re_inv;								//Flow and material parameters
double residue, max_res;														//RMS Residue and maximum residue in the fluid domain
int max_res_cell;																//Cell number with maximum residue
double pi = 4.0 * atan(1.0);
int imax = 298;

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
	void evaluate_flux();
	void forward_sweep();
	void backward_sweep();
	void state_update();
	void print_output();

	std::ofstream outfile("residue");

	double res_old;

	input_data();
	aliasing();
	int T = 250;
	for (int t = 1; t <= max_iters; t++)
	{
		initial_conditions();
		q_variables();
		q_derivatives();
		evaluate_flux();
		forward_sweep();
		backward_sweep();
		state_update();

		if (t == 1)
			res_old = residue;
		residue = log10(residue / res_old);
		/*if (log10(max_res) < -12 && log10(residue) < -12)
		{
			std::cout << "Convergence criteria reached.\n";
			break;
		}*/
		std::cout << t << "\t" << residue << "\t" << max_res << "\t" << max_res_cell << std::endl;
		outfile << t << "\t" << residue << "\t" << max_res << "\t" << max_res_cell << std::endl;

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
Grid data: 2order-input-data
Flow Parameters: flow-parameters-qkfvs
*/
void input_data()
{
	std::ifstream infile("ogrid_viscous");
	std::ifstream infile2("fp_viscous");
	//Input Flow Parameters
	infile2 >> Mach >> aoa >> cfl >> max_iters >> limiter_const;
	//Set Reference Values
	u_ref = sqrt(R * gma * tr_inf);
	tr_ref = gma * tr_inf;
	pr_ref = rho_ref * u_ref * u_ref;
	Re_inv = mu / (rho_ref * u_ref * L);
	//Input Edge data
	infile >> max_edges;
	edge = new Edge[max_edges + 1];
	for (int k = 1; k <= max_edges; k++)
	{
		infile >> edge[k].mx >> edge[k].my >> edge[k].lcell >> edge[k].rcell >> edge[k].status >> edge[k].nx >> edge[k].ny >> edge[k].length;
		//Non-dimensionalise x and y wrt L
		edge[k].mx = edge[k].mx / L;
		edge[k].my = edge[k].my / L;
		edge[k].length = edge[k].length / L;
	}
	//Input Cell data
	infile >> max_cells;
	cell = new Cell[max_cells + 1];
	for (int k = 1; k <= max_cells; k++)
	{
		infile >> cell[k].cx >> cell[k].cy >> cell[k].rho >> cell[k].u1 >> cell[k].u2 >> cell[k].pr >> cell[k].area >> cell[k].noe;
		cell[k].tp = cell[k].pr / cell[k].rho;
		//Set enclosing edges
		cell[k].edge = new int[cell[k].noe + 1];
		for (int r = 1; r <= cell[k].noe; r++)
			infile >> cell[k].edge[r];
		//Set neighbours
		infile >> cell[k].nbhs;
		cell[k].conn = new int[cell[k].nbhs];
		for (int r = 0; r < cell[k].nbhs; r++)
			infile >> cell[k].conn[r];
		//Non-dimensionalise x and y wrt L
		cell[k].cx = cell[k].cx / L;
		cell[k].cy = cell[k].cy / L;
		cell[k].area = cell[k].area / (L * L);
	}
	infile.close();
	infile2.close();
} //End of the function

/*
Calculate derivatives of velocity at cell edges necessary to compute the wall shear stresses
This function uses the neighbours of edge 'e' formed by `void construct_neighbours()`
*/
void u_derivatives(int *edge_nbs, int nbs, double *uref, double *u_dash, int lcell, int rcell, int e)
{
	double sig_dxdx, sig_dydy, sig_dxdy;
	double sig_dxdu1, sig_dxdu2, sig_dydu1, sig_dydu2; //Summation of (delta x (or y)) * (delta u1 (or u2))
	sig_dxdx = sig_dxdy = sig_dydy = sig_dxdu1 = sig_dxdu2 = sig_dydu1 = sig_dydu2 = 0;
	if (edge[e].status == 'w')
	{
		uref[1] = (cell[lcell].u1); //Linear Interpolation
		uref[2] = (cell[lcell].u2);
	}
	else if (edge[e].status == 'o')
	{
		uref[1] = (cell[rcell].u1); //Linear Interpolation
		uref[2] = (cell[rcell].u2);
	}
	else
	{
		double w1 = sqrt(pow((cell[lcell].cx - edge[e].mx), 2) + pow((cell[lcell].cy - edge[e].my), 2)); //Distance between left cell center and edge midpoint
		double w2 = sqrt(pow((cell[rcell].cx - edge[e].mx), 2) + pow((cell[rcell].cy - edge[e].my), 2)); //Distance between right cell center and edge midpoint
		uref[1] = (w1 * cell[lcell].u1 + w2 * cell[rcell].u1) / (w1 + w2);								 //Linear Interpolation
		uref[2] = (w1 * cell[lcell].u2 + w2 * cell[rcell].u2) / (w1 + w2);								 //Linear Interpolation
	}

	for (int j = 1; j <= nbs; j++)
	{
		int p = edge_nbs[j];
		double delx = cell[p].cx - edge[e].mx;
		double dely = cell[p].cy - edge[e].my;
		double dist = sqrt(delx * delx + dely * dely); //Distance between cell and edge midpoint
		double w = 1 / pow(dist, 2.0);				   //Least Squares weight
		double delu1 = cell[p].u1 - (uref[1]);
		double delu2 = cell[p].u2 - (uref[2]);

		sig_dxdx = sig_dxdx + w * delx * delx;
		sig_dydy = sig_dydy + w * dely * dely;
		sig_dxdy = sig_dxdy + w * delx * dely;

		sig_dxdu1 = sig_dxdu1 + w * delx * delu1;
		sig_dydu1 = sig_dydu1 + w * dely * delu1;
		sig_dxdu2 = sig_dxdu2 + w * delx * delu2;
		sig_dydu2 = sig_dydu2 + w * dely * delu2;
	}
	double det = sig_dxdx * sig_dydy - sig_dxdy * sig_dxdy;
	u_dash[1] = (sig_dydy * sig_dxdu1 - sig_dxdy * sig_dydu1) / det;
	u_dash[2] = (sig_dxdx * sig_dydu1 - sig_dxdy * sig_dxdu1) / det;
	u_dash[3] = (sig_dydy * sig_dxdu2 - sig_dxdy * sig_dydu2) / det;
	u_dash[4] = (sig_dxdx * sig_dydu2 - sig_dxdy * sig_dxdu2) / det;
} //End of the function

/*
Calculate derivatives of temperature at the wall face.
This function uses the neighbours of edge 'e' formed by `void construct_neighbours()`
*/
void tp_derivatives(int *edge_nbs, int nbs, double *tp_dash, int lcell, int rcell, int e)
{
	double sig_dxdx, sig_dydy, sig_dxdy;
	double sig_dxdtp, sig_dydtp; //Summation of (delta x (or y)) * (delta tp)
	double tpref;
	sig_dxdx = sig_dxdy = sig_dydy = sig_dxdtp = sig_dydtp = 0;
	if (lcell == 0)
	{
		tpref = cell[rcell].tp;
	}
	else if (rcell == 0)
	{
		tpref = cell[lcell].tp;
	}
	else
	{
		double w1 = sqrt(pow((cell[lcell].cx - edge[e].mx), 2) + pow((cell[lcell].cy - edge[e].my), 2)); //Distance between left cell center and edge midpoint
		double w2 = sqrt(pow((cell[rcell].cx - edge[e].mx), 2) + pow((cell[rcell].cy - edge[e].my), 2)); //Distance between right cell center and edge midpoint
		tpref = (w1 * cell[lcell].tp + w2 * cell[rcell].tp) / (w1 + w2);								 //Linear Interpolation
	}
	if (edge[e].status == 'w') //Boundary condition, adiabetic wall
		tp_dash[1] = tp_dash[2] = 0;
	else
	{
		for (int j = 1; j <= nbs; j++)
		{
			int p = edge_nbs[j];
			double delx = cell[p].cx - edge[e].mx;
			double dely = cell[p].cy - edge[e].my;
			double dist = sqrt(delx * delx + dely * dely); //Distance between cell and edge midpoint
			double w = 1 / pow(dist, 2.0);				   //Least Squares weight
			double deltp = cell[p].tp - (tpref);

			sig_dxdx = sig_dxdx + w * delx * delx;
			sig_dydy = sig_dydy + w * dely * dely;
			sig_dxdy = sig_dxdy + w * delx * dely;

			sig_dxdtp = sig_dxdtp + w * delx * deltp;
			sig_dydtp = sig_dydtp + w * dely * deltp;
		}
		double det = sig_dxdx * sig_dydy - sig_dxdy * sig_dxdy;
		tp_dash[1] = (sig_dydy * sig_dxdtp - sig_dxdy * sig_dydtp) / det;
		tp_dash[2] = (sig_dxdx * sig_dydtp - sig_dxdy * sig_dxdtp) / det;
	}
	/*std::fstream tpder;
	tpder.open("tpder.dat", std::fstream::app);
	if (e == 1)
		std::remove("tpder.dat");
	tpder << "Edge Number: " << e << " dT/dx = " << tp_dash[1] << " dT/dy = " << tp_dash[2] << "\n";
	tpder.close();*/
} //End of the function

int construct_neighbours(int *edge_nbs, int e, int lcell, int rcell, char status, int nbs)
{
	int flag, k = 0;
	if (nbs == 6)
	{
		if (status == 'f')
		{
			edge_nbs[++k] = rcell;
			for (int i = 0; i < cell[lcell].nbhs; i++)
			{
				for (int j = 0; j < cell[rcell].nbhs; j++)
				{
					if (cell[lcell].conn[i] == cell[rcell].conn[j])
					{
						edge_nbs[++k] = cell[lcell].conn[i];
					}
				}
			}
			edge_nbs[++k] = lcell;
		}
		else if (status == 'w')
		{
			for (int i = 0; i < cell[lcell].nbhs; i++)
			{
				edge_nbs[++k] = cell[lcell].conn[i];
			}
			edge_nbs[++k] = lcell;
		}
		else if (status == 'o')
		{
			for (int i = 0; i < cell[rcell].nbhs; i++)
			{
				edge_nbs[++k] = cell[rcell].conn[i];
			}
			edge_nbs[++k] = rcell;
		}
	}
	else if (nbs == 8)
	{
		for (int i = 0; i < cell[lcell].nbhs; i++)
		{
			flag = 0;
			for (int j = 0; j < cell[rcell].nbhs; j++)
			{
				if ((cell[lcell].conn[i] == cell[rcell].conn[j]))
				{
					flag = 1;
				}
			}
			if (flag == 1)
				continue;
			else
			{
				edge_nbs[++k] = cell[lcell].conn[i];
			}
		}
		for (int i = 0; i < cell[rcell].nbhs; i++)
		{
			edge_nbs[++k] = cell[rcell].conn[i];
		}
	}
	return k;
}

void viscous_flux(double *G, double nx, double ny, int e, char status)
{
	int *edge_nbs;
	int lcell, rcell;
	lcell = edge[e].lcell;
	rcell = edge[e].rcell;
	int nbs = 6;
	if (status == 'f')
	{
		if ((cell[lcell].nbhs == 5) && (cell[rcell].nbhs == 5))
			nbs = 8;
	}
	edge_nbs = new int[nbs + 1];
	int k = construct_neighbours(edge_nbs, e, lcell, rcell, status, nbs);
	if (k != 6 && k != 8)
	{
		std::cout << "Error in construct_neighbours(...). k = " << k << ", e = " << e << std::endl;
		exit(0);
	}
	//Calculating flux
	double uref[3], u_dash[5], tp_dash[3] = {0, 0, 0};
	double tauxx, tauyy, tauxy;
	double kappa = 0.024 / (L * rho_ref * R * sqrt(gma * R * tr_inf));
	//u_derivatives(edge_nbs, nbs, uref, u_dash, lcell, rcell, e);
	//shear_stress(tauxx, tauyy, tauxy, u_dash);
	tauxx = Re_inv * (4 / 3 * u_dash[1] - 2 / 3 * u_dash[4]);
	tauyy = Re_inv * (4 / 3 * u_dash[4] - 2 / 3 * u_dash[1]);
	tauxy = Re_inv * (u_dash[2] + u_dash[3]);
	tp_derivatives(edge_nbs, nbs, tp_dash, lcell, rcell, e);
	//Storing viscous fluxes
	tp_dash[1] = tp_dash[2] = 0;
	if (status == 'w')
	{
		uref[1] = uref[2] = 0;
	}
	G[1] = 0;
	G[2] = (tauxx * nx + tauxy * ny);
	G[3] = (tauxy * nx + tauyy * ny);
	G[4] = ((uref[1] * tauxx + uref[2] * tauxy + kappa * tp_dash[1]) * nx) + ((uref[1] * tauxy + uref[2] * tauyy + kappa * tp_dash[2]) * ny);
	delete[] edge_nbs;
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
	double rho_inf, u1_inf, u2_inf, pr_inf;

	double theta = aoa * pi / 180;

	rho_inf = 1.0;
	u1_inf = Mach * cos(theta);
	u2_inf = Mach * sin(theta);
	pr_inf = 1.0 / gma;

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
	cell[k].tp = cell[k].pr / cell[k].rho;
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
	//double delt, area;
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

//Function which prints the final primitive vector into the file "primitive-vector.dat"

void write_tecplot()
{
	std::ofstream tecplotfile("flow.dat");
	tecplotfile << "TITLE: \"QKFVS Viscous Code - NAL\"\n";
	tecplotfile << "VARIABLES= \"X\", \"Y\", \"Density\", \"Pressure\", \"x-velocity\", \"y-velocity\", \"Mach Number\", \"Temperature\"\n";
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
}
void print_output()
{
	std::ofstream outfile("primitive-vector.dat");

	for (int k = 1; k <= max_cells; k++)
	{
		outfile << k << "\t" << cell[k].rho << "\t" << cell[k].u1 << "\t" << cell[k].u2 << "\t" << cell[k].pr << std::endl;
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
		temp = temp + (fabs(u1 * nx + u2 * ny) + a + (Re_inv / (rho * h))) * length;
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
