#include <fstream>
#include <iostream>
#include <cmath>

int max_edges, max_cells, max_iters;
double Mach, aoa, cfl, limiter_const;
double L = 4.45136e-4, rho_ref = 1.225, pr_ref, u_ref, tr_ref, tr_inf = 288.15; //Reference Quantities
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
/*
Calculate derivatives of velocity at cell edges necessary to compute the wall shear stresses
This function uses the neighbours of edge 'e' formed by `void construct_neighbours()`
*/
void u_derivatives(int *edge_nbs, int nbs, double &u1ref, double &u2ref, double *u_dash, int lcell, int rcell, int e)
{
	double sig_dxdx, sig_dydy, sig_dxdy;
	double sig_dxdu1, sig_dxdu2, sig_dydu1, sig_dydu2; //Summation of (delta x (or y)) * (delta u1 (or u2))
	for (int j = 1; j <= nbs; j++)
	{
		int p = edge_nbs[j];
		double delx = cell[p].cx - edge[e].mx;
		double dely = cell[p].cy - edge[e].my;
		double dist = sqrt(delx * delx + dely * dely); //Distance between cell and edge midpoint
		double w = 1 / pow(dist, 2.0);				   //Least Squares weight
		if (lcell == 0)
		{
			u1ref = cell[rcell].u1;
			u2ref = cell[rcell].u2;
		}
		else if (rcell == 0)
		{
			u1ref = cell[lcell].u1;
			u2ref = cell[lcell].u2;
		}
		else
		{
			double w1 = sqrt(pow((cell[lcell].cx - edge[e].mx), 2) + pow((cell[lcell].cy - edge[e].my), 2)); //Distance between left cell center and edge midpoint
			double w2 = sqrt(pow((cell[rcell].cx - edge[e].mx), 2) + pow((cell[rcell].cy - edge[e].my), 2)); //Distance between right cell center and edge midpoint
			u1ref = (w1 * cell[lcell].u1 + w2 * cell[rcell].u1) / (w1 + w2);								 //Linear Interpolation
			u2ref = (w1 * cell[lcell].u2 + w2 * cell[rcell].u2) / (w1 + w2);								 //Linear Interpolation
		}
		double delu1 = cell[p].u1 - (u1ref);
		double delu2 = cell[p].u2 - (u2ref);

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
	for (int j = 1; j <= nbs; j++)
	{
		int p = edge_nbs[j];
		double delx = cell[p].cx - edge[e].mx;
		double dely = cell[p].cy - edge[e].my;
		double dist = sqrt(delx * delx + dely * dely); //Distance between cell and edge midpoint
		double w = 1 / pow(dist, 2.0);				   //Least Squares weight
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
	if (edge[e].status == 'w') //Boundary condition, adiabetic wall
		tp_dash[1] = tp_dash[2] = 0;
} //End of the function

//Calculate the wall shear stresses required to evaluate viscous wall fluxes
void shear_stress(double &tauxx, double &tauyy, double &tauxy, double *u_dash)
{
	tauxx = Re_inv * (4 / 3 * u_dash[1] - 2 / 3 * u_dash[4]);
	tauyy = Re_inv * (4 / 3 * u_dash[4] - 2 / 3 * u_dash[1]);
	tauxy = Re_inv * (u_dash[2] + u_dash[3]);
}

/*
Constructs neighbour cells for the given edge 'e'
Function returns the number of cells that are closest neighbours to the edge 'e'.
*/
int construct_neighbours(int edge_nbs[9], int e, int lcell, int rcell, char status)
{
	int k = 0;
	int nbs;
	if (status == 'f')
	{
		if ((cell[lcell].nbhs == 5) && (cell[rcell].nbhs == 5))
			nbs = 4;
		else //if (((cell[lcell].nbhs == 8) && (cell[rcell].nbhs == 5)) || ((cell[lcell].nbhs == 8) && (cell[rcell].nbhs == 8)) || ((cell[lcell].nbhs == 5) && (cell[rcell].nbhs == 8)))
			nbs = 6;
	}
	else if (status == 'w' || status == 'o')
	{
		nbs = 3;
	}
	//edge_nbs = new int[nbs + 1];
	if (status == 'f')
	{
		edge_nbs[++k] = rcell;
		for (int i = 1; i <= cell[lcell].nbhs; i++)
		{
			for (int j = 1; j <= cell[rcell].nbhs; j++)
			{
				if ((cell[lcell].conn[i] == cell[rcell].conn[j]))
				{
					edge_nbs[++k] = cell[lcell].conn[i];
				}
			}
		}
		edge_nbs[++k] = lcell;
	}
	else if (status == 'w')
	{
		edge_nbs[++k] = lcell;
		edge_nbs[++k] = edge[(e % imax) + 1].lcell;
		edge_nbs[++k] = (e == 1) ? (edge[imax].lcell) : (edge[e - 1].lcell);
	}
	else if (status == 'o')
	{
		edge_nbs[++k] = rcell;
		edge_nbs[++k] = (e == max_edges) ? (edge[e - (imax - 1)].lcell) : (edge[e + 1].lcell);
		edge_nbs[++k] = ((e % imax) == 1) ? (edge[max_edges].lcell) : (edge[e - 1].lcell);
	}
	return nbs;
}

void viscous_flux(double *G, double nx, double ny, int e, char status)
{
	//int construct_neighbours(int *, int, int, int, char);
	//void u_derivatives(int *, int, double, double, double *, int, int, int);
	//void shear_stress(double, double, double, double *);
	//void tp_derivatives(int *, int, double *, int, int, int);
	//Neighbouring cell identification
	int edge_nbs[9];
	int nbs;
	int lcell, rcell;
	lcell = edge[e].lcell;
	rcell = edge[e].rcell;
	nbs = construct_neighbours(edge_nbs, e, lcell, rcell, status);
	//Calculating flux
	double u1ref, u2ref, u_dash[5], tp_dash[3];
	double tauxx, tauyy, tauxy;
	double kappa = 0.024 / (L * rho_ref * R * sqrt(gma * R * tr_inf));
	u_derivatives(edge_nbs, nbs, u1ref, u2ref, u_dash, lcell, rcell, e);
	shear_stress(tauxx, tauyy, tauxy, u_dash);
	tp_derivatives(edge_nbs, nbs, tp_dash, lcell, rcell, e);
	//Storing viscous fluxes
	G[1] = 0;
	G[2] = (tauxx * nx + tauxy * ny);
	G[3] = (tauxy * nx + tauyy * ny);
	G[4] = ((u1ref * tauxx + u2ref * tauxy + kappa * tp_dash[1]) * nx) + ((u1ref * tauxy + u2ref * tauyy + kappa * tp_dash[2]) * ny);
}
