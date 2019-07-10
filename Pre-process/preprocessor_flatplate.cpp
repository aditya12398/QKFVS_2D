/*
This Pre-processor is made to handle rectangular, flatplate grids specially. Do not try to run it on an O-Grid or a C-H Grid.
*/

#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

int imax = 201;
int jmax = 51;
int max_points;
int max_cells;
int max_edges;

double M = 0.1;
double alpha = 0.0;

struct Point
{

	double x, y;
	int noc, c[15];
};

struct Cell
{
	int v1, v2, v4, v3;
	double area, cx, cy;
	double rho, u1, u2, pr;
	int noe, e[5], nbhs, conn[110];
};

struct Edge
{
	int v1, v2;
	int lcell, rcell;
	char status;
	double nx, ny, length;
	double mx, my;
};

struct Point point[154000];
struct Cell cell[154000];
struct Edge edge[400101];

int main(int arg, char *argv[])
{
	string fname = "./Flatplate/fixed_Flatplate_2.dat";
	void point_data(string);
	void print();
	void cell_data();
	void edge_data();
	void output1();
	void output2();
	void point_data_print();
	void vigie_plot();
	max_points = imax * jmax;
	//cout << "Enter the input file name of the grid file : ";
	//cin >> fname;
	point_data(fname);
	cell_data();
	edge_data();
	print();
	output1();
	output2();
	vigie_plot();
	point_data_print();
} //Main end of the program

void point_data(string title)
{
	ifstream infile(title);
	double x, y;
	int i, j, num;
//	infile >> imax >> jmax;
	for (int k = 1; k <= max_points; k++)
	{
		infile >> x >> y >> i >> j;
		num = i + (j - 1) * imax;
		point[num].x = x;
		point[num].y = y;
		if (j == 1)
		{
			point[k].noc = 2;
			if (i == 1)
			{
				point[k].c[1] = i;
				point[k].c[2] = imax;
			}
			else if (i > 1 && i <= imax)
			{
				point[k].c[1] = i;
				point[k].c[2] = i - 1;
			}
		} //End of j = 1
		else if (j > 1 && j < jmax)
		{
			point[k].noc = 4;
			if (i == 1)
			{
				point[k].c[1] = i + (j - 1) * imax;
				point[k].c[2] = j * imax;
				point[k].c[3] = point[k].c[1] - imax;
				point[k].c[4] = point[k].c[2] - imax;
			}
			else if (i > 1 && i <= imax)
			{
				point[k].c[1] = i + (j - 1) * imax;
				point[k].c[2] = point[k].c[1] - 1;
				point[k].c[3] = point[k].c[1] - imax;
				point[k].c[4] = point[k].c[3] - 1;
			}
		} //End of j>1 && j<jmax
		else if (j == jmax)
		{
			point[k].noc = 2;
			if (i == 1)
			{
				point[k].c[1] = i + (j - 2) * imax;
				point[k].c[2] = (j - 1) * imax;
			}
			else if (i > 1 && i <= imax)
			{
				point[k].c[1] = i + (j - 2) * imax;
				point[k].c[2] = point[k].c[1] - 1;
			}
		}
	}
} //End of the function

void cell_data()
{
	int get_cell(int, int);
	max_cells = (imax - 1) * (jmax - 1);
	int i, j;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;
	for (int k = 1; k <= max_cells; k++)
	{
		int quo = int(k / (imax - 1));
		int rem = int(k % (imax - 1));
		if (rem == 0)
		{
			j = quo;
			i = (imax - 1);
		}
		else if (rem != 0)
		{
			j = quo + 1;
			i = rem;
		}

		cell[k].v1 = i + (j - 1) * imax;
		cell[k].v2 = cell[k].v1 + 1;
		cell[k].v3 = cell[k].v2 + imax;
		cell[k].v4 = cell[k].v3 - 1;

		cell[k].e[1] = i + (2 * imax - 1) * (j - 1);
		cell[k].e[2] = cell[k].e[1] + (imax - 1);
		cell[k].e[3] = cell[k].e[2] + imax;
		cell[k].e[4] = cell[k].e[2] + 1;

		x1 = point[cell[k].v1].x;
		y1 = point[cell[k].v1].y;
		x2 = point[cell[k].v2].x;
		y2 = point[cell[k].v2].y;
		x3 = point[cell[k].v3].x;
		y3 = point[cell[k].v3].y;
		x4 = point[cell[k].v4].x;
		y4 = point[cell[k].v4].y;
		double Area = (x1 * y2 - x2 * y1) + (x2 * y3 - x3 * y2) + (x3 * y4 - x4 * y3) + (x4 * y1 - x1 * y4);
		Area = Area * 0.5;
		cell[k].area = fabs(Area);
		// finding the centroid of the quadrilateral
		// finding the centroid
		double cx = (x1 + x2) * (x1 * y2 - x2 * y1) + (x2 + x3) * (x2 * y3 - x3 * y2) +
					(x3 + x4) * (x3 * y4 - x4 * y3) + (x4 + x1) * (x4 * y1 - x1 * y4);

		cell[k].cx = cx / (6.0 * Area);
		double cy = (y1 + y2) * (x1 * y2 - x2 * y1) + (y2 + y3) * (x2 * y3 - x3 * y2) +
					(y3 + y4) * (x3 * y4 - x4 * y3) + (y4 + y1) * (x4 * y1 - x1 * y4);
		cell[k].cy = cy / (6.0 * Area);
		// initial conditions for each cell
		double pi = 4.0 * atan(1.0);
		double theta = alpha * pi / 180;
		double u_ref = sqrt(1.4 * 287 * 288.20);
		double Re = 5000;
		double mu = 1.461E-6 * (pow(288.20, 1.5) / (288.20 + 110.5));
		double rho = mu * Re / (u_ref * M);
		double u1 = u_ref * M * cos(theta);
		double u2 = u_ref * M * sin(theta);
		double pr = rho * 287 * 288.20;
		cell[k].rho = rho;
		cell[k].u1 =  u1;
		cell[k].u2 =  u2;
		cell[k].pr =  pr;
		// finding the connectivity for each cell
		int nbhs = 0;
		if (j == 1)
		{
			if (i == 1)
			{
				cell[k].conn[nbhs++] = get_cell(i, j + 1);
				cell[k].conn[nbhs++] = get_cell(i, j + 2);
				cell[k].conn[nbhs++] = get_cell(i + 1, j);
				cell[k].conn[nbhs++] = get_cell(i + 1, j + 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j + 2);
				cell[k].conn[nbhs++] = get_cell(i + 2, j);
				cell[k].conn[nbhs++] = get_cell(i + 2, j + 1);
				cell[k].conn[nbhs++] = get_cell(i + 2, j + 2);
			}
			else if (i > 1 && i < (imax - 1))
			{
				cell[k].conn[nbhs++] = get_cell(i - 1, j);
				cell[k].conn[nbhs++] = get_cell(i - 1, j + 1);
				cell[k].conn[nbhs++] = get_cell(i - 1, j + 2);
				cell[k].conn[nbhs++] = get_cell(i, j + 1);
				cell[k].conn[nbhs++] = get_cell(i, j + 2);
				cell[k].conn[nbhs++] = get_cell(i + 1, j);
				cell[k].conn[nbhs++] = get_cell(i + 1, j + 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j + 2);
			}
			if (i == (imax - 1))
			{
				cell[k].conn[nbhs++] = get_cell(i - 2, j);
				cell[k].conn[nbhs++] = get_cell(i - 2, j + 1);
				cell[k].conn[nbhs++] = get_cell(i - 2, j + 2);
				cell[k].conn[nbhs++] = get_cell(i - 1, j);
				cell[k].conn[nbhs++] = get_cell(i - 1, j + 1);
				cell[k].conn[nbhs++] = get_cell(i - 1, j + 2);
				cell[k].conn[nbhs++] = get_cell(i, j + 1);
				cell[k].conn[nbhs++] = get_cell(i, j + 2);
			}
		}
		else if (j > 1 && j < jmax - 1)
		{
			if (i == 1)
			{
				cell[k].conn[nbhs++] = get_cell(i, j - 1);
				cell[k].conn[nbhs++] = get_cell(i, j + 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j);
				cell[k].conn[nbhs++] = get_cell(i + 1, j + 1);
				cell[k].conn[nbhs++] = get_cell(i + 2, j - 1);
				cell[k].conn[nbhs++] = get_cell(i + 2, j);
				cell[k].conn[nbhs++] = get_cell(i + 2, j + 1);
			}
			else if (i > 1 && i < (imax - 1))
			{
				cell[k].conn[nbhs++] = get_cell(i - 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i, j - 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i - 1, j);
				cell[k].conn[nbhs++] = get_cell(i + 1, j);
				cell[k].conn[nbhs++] = get_cell(i - 1, j + 1);
				cell[k].conn[nbhs++] = get_cell(i, j + 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j + 1);
			}
			else if (i == (imax - 1))
			{
				cell[k].conn[nbhs++] = get_cell(i, j - 1);
				cell[k].conn[nbhs++] = get_cell(i, j + 1);
				cell[k].conn[nbhs++] = get_cell(i - 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i - 1, j);
				cell[k].conn[nbhs++] = get_cell(i - 1, j + 1);
				cell[k].conn[nbhs++] = get_cell(i - 2, j - 1);
				cell[k].conn[nbhs++] = get_cell(i - 2, j);
				cell[k].conn[nbhs++] = get_cell(i - 2, j + 1);
			}
		}
		else if (j == jmax - 1)
		{
			if (i == 1)
			{
				cell[k].conn[nbhs++] = get_cell(i, j - 1);
				cell[k].conn[nbhs++] = get_cell(i, j - 2);
				cell[k].conn[nbhs++] = get_cell(i + 1, j);
				cell[k].conn[nbhs++] = get_cell(i + 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j - 2);
				cell[k].conn[nbhs++] = get_cell(i + 2, j);
				cell[k].conn[nbhs++] = get_cell(i + 2, j - 1);
				cell[k].conn[nbhs++] = get_cell(i + 2, j - 2);
			}
			if (i > 1 && i < (imax - 1))
			{
				cell[k].conn[nbhs++] = get_cell(i - 1, j - 2);
				cell[k].conn[nbhs++] = get_cell(i - 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i - 1, j);
				cell[k].conn[nbhs++] = get_cell(i, j - 2);
				cell[k].conn[nbhs++] = get_cell(i, j - 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j - 2);
				cell[k].conn[nbhs++] = get_cell(i + 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i + 1, j);
			}
			if (i == (imax - 1))
			{
				cell[k].conn[nbhs++] = get_cell(i, j - 1);
				cell[k].conn[nbhs++] = get_cell(i, j - 2);
				cell[k].conn[nbhs++] = get_cell(i - 1, j);
				cell[k].conn[nbhs++] = get_cell(i - 1, j - 1);
				cell[k].conn[nbhs++] = get_cell(i - 1, j - 2);
				cell[k].conn[nbhs++] = get_cell(i - 2, j);
				cell[k].conn[nbhs++] = get_cell(i - 2, j - 1);
				cell[k].conn[nbhs++] = get_cell(i - 2, j - 2);
			}
		}
		cell[k].nbhs = nbhs;
	}
} //End of the function

// giving edge data
void edge_data()
{
	void normals(int);
	void length(int);
	max_edges = (imax - 1) * jmax + imax * (jmax - 1);
	int N = 1;
	for (int j = 1; j <= jmax; j++)
	{
		if (j == 1)
		{
			//The following loop takes care of the horizontal edges
			for (int i = 1; i < imax; i++)
			{
				edge[N].v1 = i;
				edge[N].v2 = i + 1;
				normals(N);
				length(N);
				edge[N].lcell = i;
				edge[N].rcell = 0;
				edge[N].status = 'w';
				N++;
			}
			//The following loop takes care of the vertical edges
			for (int i = 1; i <= imax; i++)
			{
				edge[N].v1 = i;
				edge[N].v2 = i + imax;
				normals(N);
				length(N);
				if (i == 1)
				{
					edge[N].lcell = 0; //imax;
					edge[N].rcell = i;
					edge[N].status = 'i';
				}
				else if (i > 1 && i < imax)
				{
					edge[N].lcell = i - 1;
					edge[N].rcell = i;
					edge[N].status = 'f';
				}
				else if (i == imax)
				{
					edge[N].lcell = i - 1;
					edge[N].rcell = 0;
					edge[N].status = 'e';
				}

				N++;
			}
		} //End of j == 1
		if (j > 1 && j < jmax)
		{
			//The following loop takes care of the horizontal edges
			for (int i = 1; i < imax; i++)
			{
				edge[N].v1 = i + (j - 1) * imax;
				edge[N].v2 = edge[N].v1 + 1;
				normals(N);
				length(N);
				edge[N].lcell = i + (j - 1) * (imax - 1);
				edge[N].rcell = i + (j - 2) * (imax - 1);
				edge[N].status = 'f';
				N++;
			}
			//The following loop takes care of the vertical edges
			for (int i = 1; i <= imax; i++)
			{
				edge[N].v1 = i + (j - 1) * imax;
				edge[N].v2 = edge[N].v1 + imax;
				normals(N);
				length(N);
				if (i == 1)
				{
					edge[N].lcell = 0; //j * imax;
					edge[N].rcell = i + (j - 1) * (imax - 1);
					edge[N].status = 'i';
				}
				else if (i > 1 && i < imax)
				{
					edge[N].lcell = (i - 1) + (j - 1) * (imax - 1);
					edge[N].rcell = edge[N].lcell + 1;
					edge[N].status = 'f';
				}
				else if (i == imax)
				{
					edge[N].lcell = (i - 1) + (j - 1) * (imax - 1);
					edge[N].rcell = 0; //edge[N].lcell + 1;
					edge[N].status = 'e';
				}
				N++;
			}
		}
		if (j == jmax)
		{
			//The following loop takes care of horizontal edges
			for (int i = 1; i < imax; i++)
			{
				edge[N].v1 = i + (j - 1) * imax;
				edge[N].v2 = edge[N].v1 + 1;
				normals(N);
				length(N);
				edge[N].lcell = 0;
				edge[N].rcell = i + (j - 2) * (imax - 1);
				edge[N].status = 'o';
				N++;
			}
		}
	} //End of j loop
} //End of the edge data function

void normals(int N)
{
	int l, r;
	double lx, ly, rx, ry;
	double mx, my;
	double ny1, ny2, nx1, nx2;
	l = edge[N].v1;
	r = edge[N].v2;
	lx = point[l].x;
	rx = point[r].x;
	mx = (lx + rx) * 0.5;
	edge[N].mx = mx;
	ly = point[l].y;
	ry = point[r].y;
	my = (ly + ry) * 0.5;
	edge[N].my = my;
	ny1 = mx - lx;
	ny2 = rx - mx;
	nx1 = my - ly;
	nx2 = ry - my;
	double nx = 0.5 * (nx1 + nx2);
	double ny = 0.5 * (ny1 + ny2);
	double det = sqrt(nx * nx + ny * ny);
	edge[N].nx = nx / det;
	edge[N].ny = -ny / det;
} //End of the function

void length(int N)
{
	int v1, v2;
	double v1x, v2x;
	double v1y, v2y;
	v1 = edge[N].v1;
	v2 = edge[N].v2;
	v1x = point[v1].x;
	v1y = point[v1].y;
	v2x = point[v2].x;
	v2y = point[v2].y;
	double length = (v2x - v1x) * (v2x - v1x) + (v2y - v1y) * (v2y - v1y);
	length = sqrt(length);
	edge[N].length = length;
} //End of the function

void print()
{
	ofstream outfile("Flatplate/print");
	for (int k = 1; k <= max_edges; k++)
		outfile << edge[k].v1 << "\t" << edge[k].v2 << "\t" << edge[k].lcell << "\t" << edge[k].rcell << endl;
} //End of the function

void output1()
{
	ofstream outfile("Flatplate/1order-input-data");
	outfile << max_edges << endl;
	for (int k = 1; k <= max_edges; k++)
		outfile << edge[k].lcell << "\t" << edge[k].rcell << "\t" << edge[k].status << "\t" << edge[k].nx << "\t" << edge[k].ny << "\t" << edge[k].length << endl;

	outfile << max_cells << endl;
	for (int k = 1; k <= max_cells; k++)
	{
		cell[k].noe = 4;
		outfile << cell[k].rho << "\t" << cell[k].u1 << "\t" << cell[k].u2 << "\t" << cell[k].pr << "\t" << cell[k].area << "\t" << cell[k].noe;
		for (int r = 1; r <= cell[k].noe; r++)
			outfile << "\t" << cell[k].e[r];
		outfile << endl;
	}
} //End of the function

void output2()
{
	ofstream outfile("Flatplate/2order-input-data");
	outfile << max_edges << endl;
	for (int k = 1; k <= max_edges; k++)
		outfile << edge[k].mx << "\t" << edge[k].my << "\t" << edge[k].lcell << "\t" << edge[k].rcell << "\t" << edge[k].status << "\t" << edge[k].nx << "\t" << edge[k].ny << "\t" << edge[k].length << endl;
	outfile << max_cells << endl;

	for (int k = 1; k <= max_cells; k++)
	{
		cell[k].noe = 4;
		outfile << cell[k].cx << "\t" << cell[k].cy << "\t" << cell[k].rho << "\t" << cell[k].u1 << "\t" << cell[k].u2 << "\t" << cell[k].pr << "\t" << cell[k].area << "\t" << cell[k].noe;
		for (int r = 1; r <= cell[k].noe; r++)
			outfile << "\t" << cell[k].e[r];
		outfile << "\t" << cell[k].nbhs;
		for (int r = 0; r < cell[k].nbhs; r++)
			outfile << "\t" << cell[k].conn[r];
		outfile << endl;
	}
	//e1<<"\t"<<cell[k].e2<<"\t"<<cell[k].e3<<"\t"<<cell[k].e4<<"\t"<<cell[k].area<<endl;
} //End of the function

//input file for vigie plot
void vigie_plot()
{
	ofstream outfile("Flatplate/vigie-file");
	outfile << max_points << endl;
	outfile << endl;
	for (int k = 1; k <= max_points; k++)
		outfile << point[k].x << "\t" << point[k].y << endl;
	outfile << endl;

	outfile << max_cells << endl;
	outfile << endl;
	for (int k = 1; k <= max_cells; k++)
		outfile << cell[k].v1 << "\t" << cell[k].v2 << "\t" << cell[k].v3 << "\t" << cell[k].v4 << "\t" << cell[k].area << endl;
}
// printing the point-data information
void point_data_print()
{
	ofstream outfile("Flatplate/point-data");
	for (int k = 1; k <= max_points; k++)
	{
		outfile << k << "\t" << point[k].noc;
		for (int r = 1; r <= point[k].noc; r++)
			outfile << "\t" << point[k].c[r];

		outfile << endl;
	}
} //End of the function
// from i and j values getting cell number
int get_cell(int i, int j)
{
	return (i + (j - 1) * (imax - 1));
}