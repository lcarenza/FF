/* Binary Fluid - coupled to hydrodynamics and P - */

//Per far ripartire il file trovare nel main #define READ_FROM_FILE 0
//Scrivere 1 al posto di 0.
//Usare il programmino Fortran starting.f per creare il file starting.dat da cui 
//rileggere i valori di Q, della phi e delle f.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* useful constants */
#define Pi 3.141592653589793
#define Piover2 1.570796326794897
#define threePiover2 4.71238898
#define twoPi 6.283185307
#define fivePiover2 7.853981634

/*
 *	Lattice Boltzmann
 */

/* program parameters */
#define Lx 1
#define Ly 128
#define Lz 128				/*32*/ /*240*/ /*100*/

#define Nmax 10			/*2000000*/ /* total number of iterations */
#define stepskip 1                /*2500*/ /* output graphics every stepskip steps */
#define tau1 2.5			/* 2.5 proportional to viscosity */
#define dt 1.0
#define densityinit 2.0			/* initial value for the density everywhere */

int iupa[Lx],idwna[Lx];
int jupa[Ly],jdwna[Ly];
int kupa[Lz],kdwna[Lz];

/* physical parameters */
#define temperature 0.5
double freeEQbulk;
double freeEphi;
double freeEel;
double freeEcol;
double freeEanc;
double freeEce;
double freeEtot;
double bodyforce = 0.0;			/* body force to set Pouiseuille flow */

double friction;
double amplitude;			/* amplitude for perturbation */

double vwtp = 0.0;
double vwbt = 0.0;
double xcm;
double ycm;
double zcm;
double sumphi;
double sumx;
double sumy;
double sumz;

/* function declarations */
void observables(double [Lx][Ly][Lz][15]);
void equilibriumdistQ(void);

void initialize(void);

void collision  (double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);
void collisionpr(double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);

void update0(double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);
void update (double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);


void plotQ(void);
void multipleplot(int n);

void plotbifurcation(void);
void plotV(int n);

void correction(int);
double globalP[3];

#define n 3
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot);
#undef n

/* lattice Boltzmann variables */
double C[Lx][Ly][Lz][15];
double Cpr[Lx][Ly][Lz][15];
double feq[Lx][Ly][Lz][15];

double f[Lx][Ly][Lz][15];
double fold[Lx][Ly][Lz][15];
double fpr[Lx][Ly][Lz][15];
double fnew[Lx][Ly][Lz][15];

int e[15][3];

/* observables */
double density[Lx][Ly][Lz];
double ux[Lx][Ly][Lz];
double uy[Lx][Ly][Lz];
double uz[Lx][Ly][Lz];
double umax;

void plotumax(int n);

/* switch */
#define CHOL 1
#define BC 0

#define READ_FROM_FILE 0

/*
 *	Binary Fluid
 */

/* binary fluid parameters */
#define alpha2phi 0.07
#define Kphi 0.14
#define K2phi 0.0

double M = 0.05;				/*0.05*/ /* need to tweak */

int radius = 32;

/* switch */
#define DROPLET 1

/* variables */
double phi[Lx][Ly][Lz];
double phiold[Lx][Ly][Lz];

double mu[Lx][Ly][Lz];		/* chemical potential */

double sxx[Lx][Ly][Lz];		/* dxphidxphi */
double sxy[Lx][Ly][Lz];
double sxz[Lx][Ly][Lz];
double syy[Lx][Ly][Lz];
double syz[Lx][Ly][Lz];
double szz[Lx][Ly][Lz];

double Ja[Lx][Ly][Lz][3];	/* advective current */
double Jd[Lx][Ly][Lz][3];	/* diffusive current */

double h_phi[Lx][Ly][Lz];
double h_phiold[Lx][Ly][Lz];

double Vx = 0.0;
double Vy = 0.0;
double Vz = 0.0;

int R1;

double phi_average;
double phi_average0;

/* functions */
void phi_function(void);

void initializephi(void);

void updatephi0(double phipr[Lx][Ly][Lz], double phiold[Lx][Ly][Lz]);
void updatephi(double phinew[Lx][Ly][Lz], double phiold[Lx][Ly][Lz]);

/* mathematical functions */
double max(double, double);
double min(double, double);
double dx_(double phi[Lx][Ly][Lz], int i, int j, int k);
double dy_(double phi[Lx][Ly][Lz], int i, int j, int k);
double dz_(double phi[Lx][Ly][Lz], int i, int j, int k);
double laplacian_(double phi[Lx][Ly][Lz], int i , int j , int k);

/* elastic parameters */
#define K 0.04
#define W 0.04

double zeta = 0.0;			/*-0.005*/ /* activity */
double xi = 0.7; //1.1 			/*1.1*/ /* |xi|>1: flow-aligning, |xi|<1: flow-tumbling; xi>0: rod-like, xi<0: disc-like */
double G = 1.0;				/* rotational diffusion constant */

/* function declarations */
void Q_function(void);

void initializeQ(void);

void updateQ0(double Qxxpr[Lx][Ly][Lz],  double Qxypr[Lx][Ly][Lz],  double Qxzpr[Lx][Ly][Lz], double Qyypr[Lx][Ly][Lz],  double Qyzpr[Lx][Ly][Lz],
	      double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz]);
void updateQ(double Qxxnew[Lx][Ly][Lz], double Qxynew[Lx][Ly][Lz], double Qxznew[Lx][Ly][Lz], double Qyynew[Lx][Ly][Lz], double Qyznew[Lx][Ly][Lz],
	     double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz]);

void add_Qperturbation(void);

/*Q tensor variables*/
double Qxx[Lx][Ly][Lz],Qxy[Lx][Ly][Lz],Qxz[Lx][Ly][Lz],Qyy[Lx][Ly][Lz],Qyz[Lx][Ly][Lz];
double Qxxold[Lx][Ly][Lz],Qxyold[Lx][Ly][Lz],Qxzold[Lx][Ly][Lz],Qyyold[Lx][Ly][Lz],Qyzold[Lx][Ly][Lz];
double DEHxx[Lx][Ly][Lz],DEHxy[Lx][Ly][Lz],DEHxz[Lx][Ly][Lz],DEHyy[Lx][Ly][Lz],DEHyz[Lx][Ly][Lz];
double DEH2xx[Lx][Ly][Lz],DEH2xy[Lx][Ly][Lz],DEH2xz[Lx][Ly][Lz],DEH2yy[Lx][Ly][Lz],DEH2yz[Lx][Ly][Lz];
double DEHxxold[Lx][Ly][Lz],DEHxyold[Lx][Ly][Lz],DEHxzold[Lx][Ly][Lz],DEHyyold[Lx][Ly][Lz],DEHyzold[Lx][Ly][Lz];
double Fh[Lx][Ly][Lz][15];
double tauxy[Lx][Ly][Lz], tauxz[Lx][Ly][Lz];
double tauyz[Lx][Ly][Lz];

double DG2xx[Lx][Ly][Lz],DG2yy[Lx][Ly][Lz],DG2xy[Lx][Ly][Lz];
double DG2xz[Lx][Ly][Lz],DG2yz[Lx][Ly][Lz],DG2zz[Lx][Ly][Lz],pG[Lx][Ly][Lz];
double Abulk = 1.0;//0.08;//0.001
double L1 = 0.03;//0.04
double phivr = 1.0;
double gammac;
double R;

/*cholesteric pitch*/

double q0=Pi/16.0;//Pi/16.0;//0.5*sin(Pi/50.);
double aa;

#define gamma 2.5 

double energy;


/*electric field variables*/

#define Inv12Pi 0.026525823848649223

double delVx=0.0;
double delVy=0.0;
double delVz=0.0;

double epsa=41.4;
double efb=0.0;//0.125;
double theta2;
double omeg=0.0;
double theta3;


double Ex[Lx][Ly][Lz],Ey[Lx][Ly][Lz],Ez[Lx][Ly][Lz];

/* Switch */
#define PERTURBED 0

/*
 *	coupling parameters
 */

double beta = 0.0;		/* 0.05, not 0.1 */

FILE *output1;
char filename1[20];
FILE *output2;
char filename2[20];
FILE *output3;
char filename3[20];


/* MPI function and definition */
#define MASTER 0
int  numtasks, taskid; //numero di task , numero del task 
int  row; //numero di cluster per lato  - number of cluster per side
int up , down , left , right , upright , upleft , downleft , downright; //variabili di vicinato - neighborood variable 
int tag_up_down = 1 , tag_down_up = 3 , tag_left_right = 2 , tag_right_left = 4; //these go as tags in passing_ function

void passing_field(double field[Lx][Ly][Lz]); // this function allows processors to communicate 
void passing_15(double field[Lx][Ly][Lz][15]);
void passing_3(double field[Lx][Ly][Lz][3]);
void passing_fields(void); 
void passing_fields_P_function(void); 

void passing_to_MASTER (double ***field_M ,double field[Lx][Ly][Lz]);

MPI_Datatype columntype;
MPI_Datatype rowtype;
MPI_Datatype columntype_M;
MPI_Datatype rowtype_M;



int n;
int i, j, k, l;



int main(int argc, char** argv)
{

//Initializing MPI environment	
MPI_Init(&argc, &argv);

	
	double starttime , endtime;			// this goes for MPI time of processing measurements
	
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);		// obtaining the rank of the process
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);	// obtain the number of processes in the communicator

	starttime = MPI_Wtime();			// starting time measurement



/****************** HERE GOES VIRTUAL TOPOLOGY **********************************/


	zeta = 0.0;
	bodyforce = 0.0;
	friction = 0.0;


    if (q0 > 0.00001) 
       aa=(2.0*cos(2.0*Pi/Ly)-2.0+4.0*Pi/Ly*sin(2.0*Pi/Ly))/(Pi/Ly*Pi/Ly);
    else
      aa=4.0;



	double tol = 1.0E-12;
	double tolerance, Vybefore, Vzbefore;

	initialize();
	initializephi();
	initializeQ();

MPI_Finalize();	
}



/* MPI function */

void passing_field(double field[Lx][Ly][Lz])
{

	double auxrow1[Lz-4] , auxrow2[Lz-4] ,auxrow3[Lz-4] ,auxrow4[Lz-4] , auxcolumn1[Ly-4] , auxcolumn2[Ly-4] ,auxcolumn3[Ly-4] ,auxcolumn4[Ly-4];

	MPI_Status Stat1; // required variable for receive routines	
	MPI_Status Stat2;
	MPI_Status Stat3;
	MPI_Status Stat4;
	MPI_Status Stat5;
	MPI_Status Stat6;
	MPI_Status Stat7;
	MPI_Status Stat8;
	MPI_Status Stat9;
	MPI_Status Stat10;
	MPI_Status Stat11;
	MPI_Status Stat12;

	/*****************************************invio dall'alto e ricevo dal basso*******************************************/

	// invio riga 2 alla penultima riga (Ly-2), 
	MPI_Send(&field[0][2][2], 1 , rowtype , up , 0 , MPI_COMM_WORLD);
	MPI_Recv(&auxrow1 , Lz-4 , MPI_DOUBLE , down , 0 , MPI_COMM_WORLD , &Stat1);
	for (int k = 0 ; k < Lz-4 ; k++) {
		field[0][Ly-2][k+2] = auxrow1[k];
	}

	// invio riga 3 all'ultima riga (Ly-1)
	MPI_Send(&field[0][3][2], 1 , rowtype , up , 1 ,MPI_COMM_WORLD);
	MPI_Recv(&auxrow2 , Lz-4 , MPI_DOUBLE , down , 1 , MPI_COMM_WORLD , &Stat2);
	for (int k = 0 ; k < Lz-4 ; k++) {
		field[0][Ly-1][k+2] = auxrow2[k];
	}		
	
	/*************************************************************************************************************************/

	
	/**************************************************invio dal basso ricevo dall'alto***************************************/
	
	//invio riga Ly-3 alla riga 1
	MPI_Send(&field[0][Ly-3][2] , 1 , rowtype , down , 2 , MPI_COMM_WORLD);
	MPI_Recv(&auxrow3 , Lz-4 , MPI_DOUBLE , up , 2 , MPI_COMM_WORLD , &Stat3);
	for (int k = 0 ; k < Lz-4 ; k++) {
		field[0][1][k+2] = auxrow3[k];
	} 

	//invio riga Ly-4 alla riga 0 
	MPI_Send(&field[0][Ly-4][2] , 1 , rowtype , down , 3 ,MPI_COMM_WORLD);
	MPI_Recv(&auxrow4 , Lz-4 , MPI_DOUBLE , up , 3 , MPI_COMM_WORLD , &Stat4);
	for (int k = 0 ; k < Lz-4 ; k++) {
		field[0][0][k+2] = auxrow4[k];
	} 
	/*************************************************************************************************************************/

	
	/***************************************************invio a sinista ricevo da destra**************************************/

	//invio la colonna 2 alla colonna Lz-2 
	MPI_Send(&field[0][2][2] , 1 , columntype , left , 4 , MPI_COMM_WORLD);	
	MPI_Recv(&auxcolumn1 , Ly-4 , MPI_DOUBLE , right , 4 , MPI_COMM_WORLD , &Stat5);
	for (int j = 0 ; j < Ly-4 ; j++) {
		field[0][j+2][Lz-2] = auxcolumn1[j];
	} 
	
	//invio la colonna 3 alla colonna Lz-1 
	MPI_Send(&field[0][2][3] , 1 , columntype , left , 5 , MPI_COMM_WORLD);
	MPI_Recv(&auxcolumn2 , Ly-4 , MPI_DOUBLE , right , 5 , MPI_COMM_WORLD , &Stat6);
	for (int j = 0 ; j < Ly-4 ; j++) {
		field[0][j+2][Lz-1] = auxcolumn2[j];
	} 
	/*************************************************************************************************************************/


	/***************************************************invio a destra ricevo da sinistra**************************************/
	
	//invio la colonna Lz-3 nella colonna 1
	MPI_Send(&field[0][2][Lz-3] , 1 , columntype , right , 6 , MPI_COMM_WORLD);
	MPI_Recv(&auxcolumn3 , Ly-4 , MPI_DOUBLE , left , 6 ,MPI_COMM_WORLD , &Stat7); 
	for (int j = 0 ; j < Ly-4 ; j++) {
		field[0][j+2][1] = auxcolumn3[j];
	}

	//invio la colonna Lz-4 nella colonna 0 
	MPI_Send(&field[0][2][Lz-4] , 1 , columntype , right , 7 , MPI_COMM_WORLD);	
	MPI_Recv(&auxcolumn4 , Ly-4 , MPI_DOUBLE , left , 7 , MPI_COMM_WORLD, &Stat8); 
	for (int j = 0 ; j < Ly-4 ; j++) {
		field[0][j+2][0] = auxcolumn4[j];
	}
	/*************************************************************************************************************************/

	//invio angolo in alto a destra in angolo in basso a sinistra
	MPI_Send(&field[0][2][Lz-3], 1 , MPI_DOUBLE , upright , 8  , MPI_COMM_WORLD);	
	MPI_Recv(&field[0][Ly-2][1] , 1 , MPI_DOUBLE , downleft ,  8 , MPI_COMM_WORLD , &Stat9 );

	//invio angolo in alto a sx ricevo dall'angolo in basso a dx
	MPI_Send(&field[0][2][2], 1 , MPI_DOUBLE , upleft , 9 , MPI_COMM_WORLD);		
	MPI_Recv(&field[0][Ly-2][Lz-2] , 1 , MPI_DOUBLE , downright , 9  , MPI_COMM_WORLD , &Stat10 );

	//invio angolo in basso a destra ricevo dall'angolo in alto a sinistra
	MPI_Send(&field[0][Ly-3][Lz-3], 1 , MPI_DOUBLE , downright , 10  ,  MPI_COMM_WORLD);
	MPI_Recv(&field[0][1][1] , 1 , MPI_DOUBLE , upleft , 10  , MPI_COMM_WORLD , &Stat11 );

	//invio angolo in basso a sx ricevo dall'angolo in alto a destra
	MPI_Send(&field[0][Ly-3][2], 1 , MPI_DOUBLE , downleft , 11  ,MPI_COMM_WORLD);	
	MPI_Recv(&field[0][1][Lz-2] , 1 , MPI_DOUBLE , upright ,  11 , MPI_COMM_WORLD , &Stat12 );

MPI_Barrier;
}

void passing_15(double field[Lx][Ly][Lz][15])
{
	double auxrow1[Lz-4] , auxrow2[Lz-4] ,auxrow3[Lz-4] ,auxrow4[Lz-4] , auxcolumn1[Ly-4] , auxcolumn2[Ly-4] ,auxcolumn3[Ly-4] ,auxcolumn4[Ly-4];

	MPI_Status Stat1; // required variable for receive routines	
	MPI_Status Stat2;
	MPI_Status Stat3;
	MPI_Status Stat4;
	MPI_Status Stat5;
	MPI_Status Stat6;
	MPI_Status Stat7;
	MPI_Status Stat8;
	MPI_Status Stat9;
	MPI_Status Stat10;
	MPI_Status Stat11;
	MPI_Status Stat12;

	for(int l = 0 ; l < 15 ; l++ ) {

		/*****************************************invio dall'alto e ricevo dal basso*******************************************/

		// invio riga 2 alla penultima riga (Ly-2), 
		MPI_Send(&field[0][2][2][l], 1 , rowtype_M , up , 0+l , MPI_COMM_WORLD);
		MPI_Recv(&auxrow1 , Lz-4 , MPI_DOUBLE , down , 0+l , MPI_COMM_WORLD , &Stat1);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][Ly-2][k+2][l] = auxrow1[k];
		}

		// invio riga 3 all'ultima riga (Ly-1)
		MPI_Send(&field[0][3][2][l], 1 , rowtype_M , up , 1+l ,MPI_COMM_WORLD);
		MPI_Recv(&auxrow2 , Lz-4 , MPI_DOUBLE , down , 1+l , MPI_COMM_WORLD , &Stat2);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][Ly-1][k+2][l] = auxrow2[k];
		}		
	
		/*************************************************************************************************************************/

	
		/**************************************************invio dal basso ricevo dall'alto***************************************/
	
		//invio riga Ly-3 alla riga 1
		MPI_Send(&field[0][Ly-3][2][l] , 1 , rowtype_M , down , 2+l , MPI_COMM_WORLD);
		MPI_Recv(&auxrow3 , Lz-4 , MPI_DOUBLE , up , 2+l , MPI_COMM_WORLD , &Stat3);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][1][k+2][l] = auxrow3[k];
		} 

		//invio riga Ly-4 alla riga 0 
		MPI_Send(&field[0][Ly-4][2][l] , 1 , rowtype_M , down , 3+l , MPI_COMM_WORLD);
		MPI_Recv(&auxrow4 , Lz-4 , MPI_DOUBLE , up , 3+l , MPI_COMM_WORLD, &Stat4);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][0][k+2][l] = auxrow4[k];
		} 
		/*************************************************************************************************************************/

	
		/***************************************************invio a sinista ricevo da destra**************************************/

		//invio la colonna 2 alla colonna Lz-2 
		MPI_Send(&field[0][2][2][l] , 1 , columntype_M , left , 4+l ,  MPI_COMM_WORLD);	
		MPI_Recv(&auxcolumn1 , Ly-4 , MPI_DOUBLE , right , 4+l , MPI_COMM_WORLD , &Stat5);
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][Lz-2][l] = auxcolumn1[j];
		} 
	
		//invio la colonna 3 alla colonna Lz-1 
		MPI_Send(&field[0][2][3][l] , 1 , columntype_M , left , 5 , MPI_COMM_WORLD);
		MPI_Recv(&auxcolumn2 , Ly-4 , MPI_DOUBLE , right , 5 , MPI_COMM_WORLD, &Stat6);
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][Lz-1][l] = auxcolumn2[j];
		} 
		/*************************************************************************************************************************/


		/***************************************************invio a destra ricevo da sinistra**************************************/
	
		//invio la colonna Lz-3 nella colonna 1
		MPI_Send(&field[0][2][Lz-3][l] , 1 , columntype_M , right , 6+l , MPI_COMM_WORLD);
		MPI_Recv(&auxcolumn3 , Ly-4 , MPI_DOUBLE , left , 6+l , MPI_COMM_WORLD , &Stat7); 
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][1][l] = auxcolumn3[j];
		}

		//invio la colonna Lz-4 nella colonna 0 
		MPI_Send(&field[0][2][Lz-4][l] , 1 , columntype_M , right , 7+l ,  MPI_COMM_WORLD);	
		MPI_Recv(&auxcolumn4 , Ly-4 , MPI_DOUBLE , left , 7+l , MPI_COMM_WORLD , &Stat8); 
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][0][l] = auxcolumn4[j];
		}
		/*************************************************************************************************************************/

		//invio angolo in alto a destra in angolo in basso a sinistra
		MPI_Send(&field[0][2][Lz-3][l], 1 , MPI_DOUBLE , upright , 8+l  , MPI_COMM_WORLD);	
		MPI_Recv(&field[0][Ly-2][1][l] , 1 , MPI_DOUBLE , downleft ,  8+l , MPI_COMM_WORLD , &Stat9 );

		//invio angolo in alto a sx ricevo dall'angolo in basso a dx
		MPI_Send(&field[0][2][2][l], 1 , MPI_DOUBLE , upleft , 9+l , MPI_COMM_WORLD);		
		MPI_Recv(&field[0][Ly-2][Lz-2][l] , 1 , MPI_DOUBLE , downright , 9 +l , MPI_COMM_WORLD , &Stat10 );

		//invio angolo in basso a destra ricevo dall'angolo in alto a sinistra
		MPI_Send(&field[0][Ly-3][Lz-3][l], 1 , MPI_DOUBLE , downright , 10+l  ,	MPI_COMM_WORLD);
		MPI_Recv(&field[0][1][1][l] , 1 , MPI_DOUBLE , upleft , 10+l  , MPI_COMM_WORLD , &Stat11 );

		//invio angolo in basso a sx ricevo dall'angolo in alto a destra
		MPI_Send(&field[0][Ly-3][2][l], 1 , MPI_DOUBLE , downleft , 11+l  ,MPI_COMM_WORLD);	
		MPI_Recv(&field[0][1][Lz-2][l] , 1 , MPI_DOUBLE , upright ,  11+l , MPI_COMM_WORLD, &Stat12 );

		MPI_Barrier;
	}
}	


void passing_3(double field[Lx][Ly][Lz][3])
{
	double auxrow1[Lz-4] , auxrow2[Lz-4] ,auxrow3[Lz-4] ,auxrow4[Lz-4] , auxcolumn1[Ly-4] , auxcolumn2[Ly-4] ,auxcolumn3[Ly-4] ,auxcolumn4[Ly-4];

	MPI_Status Stat1; // required variable for receive routines	
	MPI_Status Stat2;
	MPI_Status Stat3;
	MPI_Status Stat4;
	MPI_Status Stat5;
	MPI_Status Stat6;
	MPI_Status Stat7;
	MPI_Status Stat8;
	MPI_Status Stat9;
	MPI_Status Stat10;
	MPI_Status Stat11;
	MPI_Status Stat12;

	for(int l = 0 ; l < 3 ; l++ ) {

		/*****************************************invio dall'alto e ricevo dal basso*******************************************/

		// invio riga 2 alla penultima riga (Ly-2), 
		MPI_Send(&field[0][2][2][l], 1 , rowtype_M , up , 0+l , MPI_COMM_WORLD);
		MPI_Recv(&auxrow1 , Lz-4 , MPI_DOUBLE , down , 0+l , MPI_COMM_WORLD , &Stat1);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][Ly-2][k+2][l] = auxrow1[k];
		}

		// invio riga 3 all'ultima riga (Ly-1)
		MPI_Send(&field[0][3][2][l], 1 , rowtype_M , up , 1+l ,MPI_COMM_WORLD);
		MPI_Recv(&auxrow2 , Lz-4 , MPI_DOUBLE , down , 1+l , MPI_COMM_WORLD, &Stat2);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][Ly-1][k+2][l] = auxrow2[k];
		}		
	
		/*************************************************************************************************************************/

	
		/**************************************************invio dal basso ricevo dall'alto***************************************/
	
		//invio riga Ly-3 alla riga 1
		MPI_Send(&field[0][Ly-3][2][l] , 1 , rowtype_M , down , 2+l , MPI_COMM_WORLD);
		MPI_Recv(&auxrow3 , Lz-4 , MPI_DOUBLE , up , 2+l , MPI_COMM_WORLD , &Stat3);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][1][k+2][l] = auxrow3[k];
		} 

		//invio riga Ly-4 alla riga 0 
		MPI_Send(&field[0][Ly-4][2][l] , 1 , rowtype_M , down , 3+l ,MPI_COMM_WORLD);
		MPI_Recv(&auxrow4 , Lz-4 , MPI_DOUBLE , up , 3+l , MPI_COMM_WORLD , &Stat4);
		for (int k = 0 ; k < Lz-4 ; k++) {
			field[0][0][k+2][l] = auxrow4[k];
		} 
		/*************************************************************************************************************************/

	
		/***************************************************invio a sinista ricevo da destra**************************************/

		//invio la colonna 2 alla colonna Lz-2 
		MPI_Send(&field[0][2][2][l] , 1 , columntype_M , left , 4+l ,  MPI_COMM_WORLD);	
		MPI_Recv(&auxcolumn1 , Ly-4 , MPI_DOUBLE , right , 4+l , MPI_COMM_WORLD , &Stat5);
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][Lz-2][l] = auxcolumn1[j];
		} 
	
		//invio la colonna 3 alla colonna Lz-1 
		MPI_Send(&field[0][2][3][l] , 1 , columntype_M , left , 5+l ,MPI_COMM_WORLD);
		MPI_Recv(&auxcolumn2 , Ly-4 , MPI_DOUBLE , right , 5+l , MPI_COMM_WORLD , &Stat6);
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][Lz-1][l] = auxcolumn2[j];
		} 
		/*************************************************************************************************************************/


		/***************************************************invio a destra ricevo da sinistra**************************************/
	
		//invio la colonna Lz-3 nella colonna 1
		MPI_Send(&field[0][2][Lz-3][l] , 1 , columntype_M , right , 6+l , MPI_COMM_WORLD);
		MPI_Recv(&auxcolumn3 , Ly-4 , MPI_DOUBLE , left , 6+l , MPI_COMM_WORLD , &Stat7); 
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][1][l] = auxcolumn3[j];
		}

		//invio la colonna Lz-4 nella colonna 0 
		MPI_Send(&field[0][2][Lz-4][l] , 1 , columntype_M , right , 7+l ,  MPI_COMM_WORLD);	
		MPI_Recv(&auxcolumn4 , Ly-4 , MPI_DOUBLE , left , 7+l , MPI_COMM_WORLD , &Stat8); 
		for (int j = 0 ; j < Ly-4 ; j++) {
			field[0][j+2][0][l] = auxcolumn4[j];
		}
		/*************************************************************************************************************************/

		//invio angolo in alto a destra in angolo in basso a sinistra
		MPI_Send(&field[0][2][Lz-3][l], 1 , MPI_DOUBLE , upright , 8+l  , MPI_COMM_WORLD);	
		MPI_Recv(&field[0][Ly-2][1][l] , 1 , MPI_DOUBLE , downleft ,  8+l ,MPI_COMM_WORLD , &Stat9 );

		//invio angolo in alto a sx ricevo dall'angolo in basso a dx
		MPI_Send(&field[0][2][2][l], 1 , MPI_DOUBLE , upleft , 9+l , MPI_COMM_WORLD);		
		MPI_Recv(&field[0][Ly-2][Lz-2][l] , 1 , MPI_DOUBLE , downright , 9 +l , MPI_COMM_WORLD , &Stat10 );

		//invio angolo in basso a destra ricevo dall'angolo in alto a sinistra
		MPI_Send(&field[0][Ly-3][Lz-3][l], 1 , MPI_DOUBLE , downright , 10+l  ,	MPI_COMM_WORLD);
		MPI_Recv(&field[0][1][1][l] , 1 , MPI_DOUBLE , upleft , 10+l  , MPI_COMM_WORLD , &Stat11 );

		//invio angolo in basso a sx ricevo dall'angolo in alto a destra
		MPI_Send(&field[0][Ly-3][2][l], 1 , MPI_DOUBLE , downleft , 11+l  ,MPI_COMM_WORLD);	
		MPI_Recv(&field[0][1][Lz-2][l] , 1 , MPI_DOUBLE , upright ,  11+l , MPI_COMM_WORLD , &Stat12 );

		MPI_Barrier;
	}
}

	

void passing_fields(void)
{
	passing_field(phi); 		MPI_Barrier;
	passing_field(phiold);		MPI_Barrier;
	passing_3(h);			MPI_Barrier;
	passing_3(hold);		MPI_Barrier;
	passing_3(hmol);		MPI_Barrier;
	passing_field(h_phi);		MPI_Barrier;
	passing_field(h_phiold);	MPI_Barrier;
	passing_field(mu);  		MPI_Barrier;
	passing_field(Px);		MPI_Barrier;
        passing_field(Py);		MPI_Barrier;
        passing_field(Pz);		MPI_Barrier;
	passing_field(Pxold);		MPI_Barrier;
        passing_field(Pyold);		MPI_Barrier;
        passing_field(Pzold);		MPI_Barrier;
	passing_field(ux);		MPI_Barrier;
        passing_field(uy);		MPI_Barrier;
        passing_field(uz);		MPI_Barrier;
	passing_field(tauxy);		MPI_Barrier;
        passing_field(tauxz);		MPI_Barrier;
	passing_field(tauyz);		MPI_Barrier;
	passing_field(sxx);		MPI_Barrier;
        passing_field(sxy);		MPI_Barrier;
        passing_field(sxz);		MPI_Barrier;
	passing_field(syy);		MPI_Barrier;
        passing_field(syz);		MPI_Barrier;
	passing_field(szz);		MPI_Barrier;
	passing_15(f);			MPI_Barrier;
	passing_15(feq);		MPI_Barrier;
	passing_15(fold);		MPI_Barrier;
	passing_15(fpr);		MPI_Barrier;
	passing_15(C);			MPI_Barrier;
	passing_15(Cpr);		MPI_Barrier;
}


void passing_fields_P_function(void) 
{
	passing_3(h);			MPI_Barrier;
	passing_3(hmol);		MPI_Barrier;
	passing_field(tauxy);		MPI_Barrier;
        passing_field(tauxz);		MPI_Barrier;
	passing_field(tauyz);		MPI_Barrier;
	passing_field(sxx);		MPI_Barrier;
        passing_field(sxy);		MPI_Barrier;
        passing_field(sxz);		MPI_Barrier;
	passing_field(syy);		MPI_Barrier;
        passing_field(syz);		MPI_Barrier;
	passing_field(szz);		MPI_Barrier;
}		


void passing_to_MASTER (double ***field_M ,double field[Lx][Ly][Lz])
{
	MPI_Status Stat1;

	MPI_Type_vector(Ly-4, 1, 1, MPI_DOUBLE, &rowtype);
	MPI_Type_commit(&rowtype);
	
	int riga , colonna , lmuert;

	double aux[Lz-4];
	
	for (int j = 2 ; j < Ly-2 ; j++) {
		MPI_Send(&field[0][j][2], 1, rowtype , MASTER, j , MPI_COMM_WORLD);
	}	
	
	if (taskid == MASTER) {
		for ( int task = 0 ; task < numtasks ; task++) {
			for (int j = 2 ; j < Ly-2 ; j++) {
				riga = (int)(task/row)*(Ly-4) + j - 2;
				colonna = (task % row) * (Lz-4);
				MPI_Recv(&aux , Lz-4 , MPI_DOUBLE , task , j , MPI_COMM_WORLD , &Stat1);
				for (int k = 0 ; k < Lz-4 ; k++) {
					field_M[0][riga][colonna+k] = aux[k];
				}
			}
		}
	}

}


void passing_from_MASTER (double field_M[Lx_M][Ly_M][Lz_M] , double field[Lx][Ly][Lz] )
{
	MPI_Status Stat1;

	MPI_Type_vector(Ly-4, 1, 1, MPI_DOUBLE, &rowtype);
	MPI_Type_commit(&rowtype);
	
	int riga , colonna , lmuert;
	double aux[Lz-4];
	
		
	
	if (taskid == MASTER) {
		for ( int task = 0 ; task < numtasks ; task++) {
			for (int j = 0 ; j < Ly_M ; j++) {
				riga = j%(Ly-4) + 2 ;
				colonna = (task % row) * (Lz-4);
				MPI_Send(&field_M[0][j][colonna] , 1 , rowtype , task , riga , MPI_COMM_WORLD );
				
			}
		}
	}

	for (int j = 2 ; j < Ly-2 ; j++) {
		MPI_Recv(&aux, Lz-4, MPI_DOUBLE , MASTER, j , MPI_COMM_WORLD, &Stat1);
		for (int k = 2 ; k < Lz-4 ; k++ ) {
			field[0][j][k] = aux[k];
		}
	}
	
	passing_field(field);
	
}

/*
 * 	Lattice Boltzmann
 */

/* initialization */
void initialize(void)
{
	int i, j, k, l;

	e[0][0]= 0;
	e[0][1]= 0;
	e[0][2]= 0;

	e[1][0]= 1;
	e[1][1]= 0;
	e[1][2]= 0;

	e[2][0]= 0;
	e[2][1]= 1;
	e[2][2]= 0;

	e[3][0]= -1;
	e[3][1]= 0;
	e[3][2]= 0;

	e[4][0]= 0;
	e[4][1]= -1;
	e[4][2]= 0;

	e[5][0]= 0;
	e[5][1]= 0;
	e[5][2]= 1;

	e[6][0]= 0;
	e[6][1]= 0;
	e[6][2]= -1;

	e[7][0]= 1;
	e[7][1]= 1;
	e[7][2]= 1;

	e[8][0]= -1;
	e[8][1]= 1;
	e[8][2]= 1;

	e[9][0]= -1;
	e[9][1]= -1;
	e[9][2]= 1;

	e[10][0]= 1;
	e[10][1]= -1;
	e[10][2]= 1;

	e[11][0]= 1;
	e[11][1]= 1;
	e[11][2]= -1;

	e[12][0]= -1;
	e[12][1]= 1;
	e[12][2]= -1;

	e[13][0]= -1;
	e[13][1]= -1;
	e[13][2]= -1;

	e[14][0]= 1;
	e[14][1]= -1;
	e[14][2]= -1;

	for (i = 0; i < Lx; i++) {
		if (i == Lx-1) {iupa[i]  = 0;}		else {iupa[i]  = i+1;}
		if (i == 0)    {idwna[i] = Lx-1;}	else {idwna[i] = i-1;}
		for (j = 1; j < Ly-1; j++) {
			//if (j == Ly-1) {jupa[j]  = 0;}		else {jupa[j]  = j+1;}
			//if (j == 0)    {jdwna[j] = Ly-1;}	else {jdwna[j] = j-1;}
			jupa[j]  = j+1;
			jdwna[j] = j-1;			
			for (k = 1; k < Lz-1; k++) {
				//if (k == Lz-1) {kupa[k]  = 0;}		else {kupa[k]  = k+1;}
				//if (k == 0)    {kdwna[k] = Lz-1;}	else {kdwna[k] = k-1;}
				kupa[k]  = k+1;
				kdwna[k] = k-1;
				
				density[i][j][k] = densityinit;

				/* initialize f */
				for (l = 0; l < 15; l++) {
					f[i][j][k][l] = density[i][j][k]/15.0;
				}
			}
		}
	}

MPI_Barrier;
	
	passing_field(density);			

MPI_Barrier;					 				

	passing_15(f);				

MPI_Barrier;
}

/* predictor-corrector of the LB step */
void update0(double fpr[Lx][Ly][Lz][15], double fold[Lx][Ly][Lz][15])
{
	int i, j, k, l, imod, jmod, kmod;
	double rb;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				/* main LB equation */
				for (l = 0; l < 15; l++) {
					imod = (i - e[l][0] + Lx) % Lx;
					jmod = (j - e[l][1] + Ly) % Ly;
					kmod = (k - e[l][2] + Lz) % Lz;

					fpr[i][j][k][l] = fold[imod][jmod][kmod][l] + dt*C[imod][jmod][kmod][l];
				}
			}
		}
	}
}
void update(double fnew[Lx][Ly][Lz][15], double fpr[Lx][Ly][Lz][15], double fold[Lx][Ly][Lz][15])
{
	int i, j, k, l, imod, jmod, kmod;
	double rb;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				/* main LB equation */
				for (l = 0; l < 15; l++) {
					imod = (i - e[l][0] + Lx) % Lx;
					jmod = (j - e[l][1] + Ly) % Ly;
					kmod = (k - e[l][2] + Lz) % Lz;

					fnew[i][j][k][l] = fold[imod][jmod][kmod][l] + 0.5*dt*(C[imod][jmod][kmod][l] + Cpr[i][j][k][l]);
				}
			}
		}
	}
}

/* collision operator */
void collision(double feq[Lx][Ly][Lz][15], double f[Lx][Ly][Lz][15])
{
	int i, j, k, l;

	for (i = 0; i < Lx; i++)
		for (j = 0; j < Ly; j++)
			for (k = 0; k < Lz; k++)
				for (l = 0; l < 15; l++) {
					C[i][j][k][l] = (feq[i][j][k][l] - f[i][j][k][l])/tau1;
				}
}
void collisionpr(double feqpr[Lx][Ly][Lz][15], double fpr[Lx][Ly][Lz][15])
{
	int i, j, k, l;

	for (i = 0; i < Lx; i++)
		for (j = 0; j < Ly; j++)
			for (k = 0; k < Lz; k++)
				for (l = 0; l < 15; l++) {
					Cpr[i][j][k][l] = (feqpr[i][j][k][l] - fpr[i][j][k][l])/tau1;
				}
}

/* calculate equilibrium distribution for Q tensor*/
void equilibriumdistQ(void)
{
  double A0,A1,A2,B1,B2,C0,C1,C2,D1,D2;
  double G1xx,G2xx,G2xy,G2xz,G2yz,G1yy,G2yy,G1zz,G2zz;

  double rho,phil,Qxxl,Qxyl,Qyyl,Qxzl,Qyzl,Qzzl,usq,udote,omdote;
  double nnxxl,nnyyl;
  double Hxx,Hyy,Hxy,Hxz,Hyz,Qsqxx,Qsqxy,Qsqyy,Qsqzz,Qsqxz,Qsqyz,TrQ2;

  double sigxx,sigyy,sigxy,sigxz,sigyz,sigzz;

  double Force[3];

  double dbdtauxb,dbdtauyb,dbdtauzb;
  double hx,hy,hz;

  int i, j, k, l;
  int iup, jup, kup;
  int iup2, jup2, kup2;
  int idwn, jdwn, kdwn;
  int idwn2, jdwn2, kdwn2;

  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {
        int iup = iupa[i];
	int jup = jupa[j];
	int kup = kupa[k];
	      
	int iup2 = iupa[iupa[i]];
	int jup2 = jupa[jupa[j]];
        int kup2 = kupa[kupa[k]];
	      
	int idwn = idwna[i];
	int jdwn = jdwna[j];
	int kdwn = kdwna[k];
	      
	int idwn2 = idwna[idwna[i]];
	int jdwn2 = jdwna[jdwna[j]];
	int kdwn2 = kdwna[kdwna[k]];
	      
	rho = density[i][j][k];
	phil = phi[i][j][k];

	Qxxl = Qxx[i][j][k];
	Qxyl = Qxy[i][j][k];
	Qyyl = Qyy[i][j][k];
	Qxzl = Qxz[i][j][k];
	Qyzl = Qyz[i][j][k];
	Qzzl = -Qxxl - Qyyl;
	      
	nnxxl = Qxxl + 1.0/3.0;
	nnyyl = Qyyl + 1.0/3.0;

	Qsqxx = Qxxl + Qyyl;
	Qsqzz = Qsqxx*Qsqxx + Qxzl*Qxzl + Qyzl*Qyzl;
	Qsqxy = Qxyl*Qsqxx + Qxzl*Qyzl;
	Qsqxz = Qxxl*Qxzl - Qxzl*Qsqxx + Qxyl*Qyzl;
	Qsqyz = Qxyl*Qxzl - Qyzl*Qsqxx + Qyyl*Qyzl;
	Qsqxx = Qxxl*Qxxl + Qxyl*Qxyl + Qxzl*Qxzl;
	Qsqyy = Qyyl*Qyyl + Qxyl*Qxyl + Qyzl*Qyzl;
	TrQ2  = Qsqxx + Qsqyy + Qsqzz;
	      
	Hxx = DEH2xx[i][j][k];
	Hxy = DEH2xy[i][j][k];
	Hxz = DEH2xz[i][j][k];
	Hyy = DEH2yy[i][j][k];
	Hyz = DEH2yz[i][j][k];
	      
	/* there should be a minus sign in actual equation */
	sigxx=2.0/3.0*xi*((1.0+3.0*Qxxl)*(Hxx*(1.0-2.0*Qxxl-Qyyl)-Hyy*(Qxxl+2.0*Qyyl)-2.0*Hyz*Qyzl)+(Hxy*Qxyl+Hxz*Qxzl)*(1.0-6.0*Qxxl));	
	sigxy=xi*(Hxy*(2.0/3.0+Qxxl-4.0*Qxyl*Qxyl+Qyyl)-Hxx*Qxyl*(-1.0+4.0*Qxxl+2.0*Qyyl)-Hyy*Qxyl*(-1.0+4.0*Qyyl+2.0*Qxxl)+Hxz*(-4.0*Qxyl*Qxzl+Qyzl)+Hyz*(Qxzl-4.0*Qxyl*Qyzl));
	sigyy=2.0/3.0*xi*((1.0+3.0*Qyyl)*(Hyy*(1.0-Qxxl-2.0*Qyyl)-Hxx*(2.0*Qxxl+Qyyl)-2.0*Hxz*Qxzl)+(Hxy*Qxyl+Hyz*Qyzl)*(1.0-6.0*Qyyl));
	sigxz=xi*(Hxz*(2.0/3.0-4.0*Qxzl*Qxzl-Qyyl)-Hxx*Qxzl*(4.0*Qxxl+2.0*Qyyl)-Hyy*Qxzl*(1.0+4.0*Qyyl+2.0*Qxxl)+Hxy*(Qyzl-4.0*Qxyl*Qxzl)+Hyz*(Qxyl-4.0*Qxzl*Qyzl));
	sigyz=xi*(Hyz*(2.0/3.0-4.0*Qyzl*Qyzl-Qxxl)-Hyy*Qyzl*(4.0*Qyyl+2.0*Qxxl)-Hxx*Qyzl*(1.0+4.0*Qxxl+2.0*Qyyl)+Hxy*(Qxzl-4.0*Qxyl*Qyzl)+Hxz*(Qxyl-4.0*Qxzl*Qyzl));
	sigzz=-(sigxx+sigyy); // Davide to check this!
	
	sigxx += zeta*phi[i][j][k]*Qxxl;
	sigxy += zeta*phi[i][j][k]*Qxyl;
	sigxz += zeta*phi[i][j][k]*Qxzl;
	sigyy += zeta*phi[i][j][k]*Qyyl;
	sigyz += zeta*phi[i][j][k]*Qyzl;
	sigzz += zeta*phi[i][j][k]*Qzzl;
	      
	/* anti-symmetric part of stress tensor */
	dbdtauxb =  dy_(tauxy,i,j,k) + dz_(tauxz,i,j,k);
	dbdtauyb = -dx_(tauxy,i,j,k) + dz_(tauyz,i,j,k);
	dbdtauzb = -dy_(tauyz,i,j,k) - dx_(tauxz,i,j,k);

	/* thermodynamic stress go- */
	//Force[0] = dx_(DG2xx,i,j,k) + dy_(DG2xy,i,j,k) + dz_(DG2xz,i,j,k) - friction*ux[i][j][k];
	//Force[1] = dx_(DG2xy,i,j,k) + dy_(DG2yy,i,j,k) + dz_(DG2yz,i,j,k) - friction*uy[i][j][k];
	//Force[2] = dx_(DG2xz,i,j,k) + dy_(DG2yz,i,j,k) + dz_(DG2zz,i,j,k) - friction*uz[i][j][k];

	Force[0] = -phi[i][j][k]*dx_(mu,i,j,k) - globalP[0];
	Force[1] = -phi[i][j][k]*dy_(mu,i,j,k) - globalP[1];
	Force[2] = -phi[i][j][k]*dz_(mu,i,j,k) - globalP[2];

	A2 = (rho*temperature+0.0*phivr*pG[i][j][k])/10.0;
	A1 = A2;
	A0 = rho-14.0*A2;
	B2 = rho/24.0;
	B1 = 8.0*B2;
	C2 = -rho/24.0;
	C1 = 2.0*C2;
	C0 = -2.0*rho/3.0;
	D2 = rho/16.0;
	D1 = 8.0*D2;
	G2xx = phivr*sigxx/16.0;
	G2yy = phivr*sigyy/16.0;
	G2zz = phivr*sigzz/16.0;
	G2xy = phivr*sigxy/16.0;
	G2xz = phivr*sigxz/16.0;
	G2yz = phivr*sigyz/16.0;
	G1xx = 8.0*G2xx;
	G1yy = 8.0*G2yy;
	G1zz = 8.0*G2zz;
	      
	usq = ux[i][j][k]*ux[i][j][k] + uy[i][j][k]*uy[i][j][k] + uz[i][j][k]*uz[i][j][k];

	/* i = 0 */
	feq[i][j][k][0] = A0 + C0*usq;
	      
	/* i = 1, 2, 3, 4, 5, 6 */
	for (l = 1; l <= 6; l++) {
	  udote = ux[i][j][k]*e[l][0] + uy[i][j][k]*e[l][1] + uz[i][j][k]*e[l][2];

	  omdote  = Force[0]*e[l][0] + Force[1]*e[l][1] + Force[2]*e[l][2];
	  omdote += dbdtauxb*e[l][0] + dbdtauyb*e[l][1] + dbdtauzb*e[l][2];
	  omdote += -phivr*(Fh[i][j][k][0]*e[l][0] + Fh[i][j][k][1]*e[l][1] + Fh[i][j][k][2]*e[l][2]);
		
	  feq[i][j][k][l] = A1 + B1*udote + C1*usq + D1*udote*udote
		   	  + G1xx*e[l][0]*e[l][0] + G1yy*e[l][1]*e[l][1] + G1zz*e[l][2]*e[l][2]
			  + tau1*omdote/3.0;
	}

	/* i = 7, 8, 9, 10, 11, 12, 13, 14 */
	for (l = 7; l <= 14; l++) {
	  udote = ux[i][j][k]*e[l][0] + uy[i][j][k]*e[l][1] + uz[i][j][k]*e[l][2];
		
	  omdote  = Force[0]*e[l][0] + Force[1]*e[l][1] + Force[2]*e[l][2];
	  omdote += dbdtauxb*e[l][0] + dbdtauyb*e[l][1] + dbdtauzb*e[l][2];
	  omdote += -phivr*(Fh[i][j][k][0]*e[l][0] + Fh[i][j][k][1]*e[l][1] + Fh[i][j][k][2]*e[l][2]);
		
	  feq[i][j][k][l] = A2 + B2*udote + C2*usq + D2*udote*udote
		      	  + G2xx*e[l][0]*e[l][0] + G2yy*e[l][1]*e[l][1] + G2zz*e[l][2]*e[l][2] 
			  + 2.0*G2xy*e[l][0]*e[l][1] + 2.0*G2xz*e[l][0]*e[l][2] + 2.0*G2yz*e[l][1]*e[l][2]
		    	  + tau1*omdote/24.0;
	}
      }
    }
  }
}

/* calculate observables */
void observables(double f[Lx][Ly][Lz][15])
{
	int i, j, k, l;
	int iup,idwn;
	int jup,jdwn;
	int kup,kdwn;
	double phi_total;
	int R;

	/* calculate local density and velocity */
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				density[i][j][k] = 0.0;

				ux[i][j][k] = 0.0;
				uy[i][j][k] = 0.0;
				uz[i][j][k] = 0.0;

				for (l = 0; l < 15; l++) {
					density[i][j][k] += f[i][j][k][l];

					ux[i][j][k] += f[i][j][k][l]*e[l][0];
					uy[i][j][k] += f[i][j][k][l]*e[l][1];
					uz[i][j][k] += f[i][j][k][l]*e[l][2];
				}
				ux[i][j][k] = ux[i][j][k]/density[i][j][k];
				uy[i][j][k] = uy[i][j][k]/density[i][j][k];
				uz[i][j][k] = uz[i][j][k]/density[i][j][k];
			}
		}
	}
	phi_total = 0.0;

	Vx = 0.0;
	Vy = 0.0;
    	Vz = 0.0;

	R = 0;
	R1 = 0;
    	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				phi_total += phi[i][j][k]; 

	 	                Vx += phi[i][j][k]*ux[i][j][k];
	        	        Vy += phi[i][j][k]*uy[i][j][k];
	        	        Vz += phi[i][j][k]*uz[i][j][k];

				if (phi[i][j][k] > 1.0) {R += 1;}
				if (R > R1) {R1 = R;}
			}
			R = 0;
		}
    	}
	R1 = R1/2;

    	Vx = Vx/phi_total;
    	Vy = Vy/phi_total;
	Vz = Vz/phi_total;
}

/*
 *	Binary Fluid
 */

/* initialization */
void initializephi(void)
{
	int i, j, k;
	int R;
        double RR = 0.01;    
         double R1;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0 ; k < Lz; k++) {
				phiold[i][j][k] = 0.0;
#if DROPLET

                              R = sqrt((i-Lx/2)*(i-Lx/2) + (j-Ly/2)*(j-Ly/2) + (k-Lz/2)*(k-Lz/2));
                                if ((R < radius)/*||(R1 < radius)||(R2 < radius)||(R3 < radius)*/){
                                        phi[i][j][k] = 2.0;
                                }
                                else {
                                        phi[i][j][k] = 0.0;
                                }



				/*		                R1=drand48();
                                if (R1 < 0.1){		
			        phi[i][j][k] = 1.0+drand48();
                               }				
				else {
					phi[i][j][k] = drand48();
					}*/
#endif
			}
		}
	}
}

/* calculate chemical potential */
void phi_function(void)
{
  int i, j, k;
  int iup, jup, kup;
  int iup2, jup2, kup2;
  int idwn, jdwn, kdwn;
  int idwn2, jdwn2, kdwn2;

  double dphidx[Lx][Ly][Lz],dphidy[Lx][Ly][Lz],dphidz[Lx][Ly][Lz];
  double laplacianphi[Lx][Ly][Lz];
  double dphidxl,dphidyl,dphidzl;
  double d2phidxdx,d2phidxdy,d2phidxdz;
  double d2phidydy,d2phidydz,d2phidzdz;

  double dQxxdx,dQxxdy,dQxxdz,dQxydx,dQxydy,dQxydz,dQyydx,dQyydy,dQyydz;
  double dQxzdx,dQxzdy,dQxzdz,dQyzdx,dQyzdy,dQyzdz;
  double dQyxdx,dQyxdy,dQyxdz;
  double dQzxdx,dQzxdy,dQzxdz,dQzydx,dQzydy,dQzydz;
  double dQzzdx,dQzzdy,dQzzdz;
  double divQx,divQy,divQz;

  double gradphisq,phil,f;
  double vx, vy, vz;
  double vxdphidx, vydphidy, vzdphidz;

  double Qsqxx,Qsqxy,Qsqxz,Qsqyy,Qsqyz,Qsqzz;
  double Qxxl,Qxyl,Qxzl,Qyyl,Qyzl,Qzzl;
  double TrQ2;

  double Psquare,gradPsq;
  double divu;

  /* find chemical potential */
  phi_average = 0.0;
  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {

	dphidx[i][j][k]=dx_(phi,i,j,k);
	dphidy[i][j][k]=dy_(phi,i,j,k);
	dphidz[i][j][k]=dz_(phi,i,j,k);

      }
    }
  }


  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {

	laplacianphi[i][j][k]=laplacian_(phi,i,j,k);

      }
    }
  }




  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {
	
	Qxxl = Qxx[i][j][k];
	Qxyl = Qxy[i][j][k];
	Qyyl = Qyy[i][j][k];
	Qxzl = Qxz[i][j][k];
	Qyzl = Qyz[i][j][k];
	Qzzl = -Qxxl - Qyyl;
				
	Qsqxx = Qxxl + Qyyl;
	Qsqzz = Qsqxx*Qsqxx + Qxzl*Qxzl + Qyzl*Qyzl;
	Qsqxy = Qxyl*Qsqxx + Qxzl*Qyzl;
	Qsqxz = Qxxl*Qxzl - Qxzl*Qsqxx + Qxyl*Qyzl;
	Qsqyz = Qxyl*Qxzl - Qyzl*Qsqxx + Qyyl*Qyzl;
	Qsqxx = Qxxl*Qxxl + Qxyl*Qxyl + Qxzl*Qxzl;
	Qsqyy = Qyyl*Qyyl + Qxyl*Qxyl + Qyzl*Qyzl;
				
	TrQ2 = Qsqxx + Qsqyy + Qsqzz;

	d2phidxdx=dx_(dphidx,i,j,k);
	d2phidxdy=dx_(dphidy,i,j,k);
	d2phidxdz=dx_(dphidz,i,j,k);
	d2phidydy=dy_(dphidy,i,j,k);
	d2phidydz=dy_(dphidz,i,j,k);
	d2phidzdz=dz_(dphidz,i,j,k);

	dQxxdx = dx_(Qxx,i,j,k);
	dQxydx = dx_(Qxy,i,j,k);
	dQxzdx = dx_(Qxz,i,j,k);
	dQyxdx = dQxydx;
	dQyydx = dx_(Qyy,i,j,k);
	dQyzdx = dx_(Qyz,i,j,k);
	dQzxdx = dQxzdx;
	dQzydx = dQyzdx;
	dQzzdx = -(dQxxdx+dQyydx);

	dQxxdy = dy_(Qxx,i,j,k);
	dQxydy = dy_(Qxy,i,j,k);
	dQxzdy = dy_(Qxz,i,j,k);
	dQyxdy = dQxydy;
	dQyydy = dy_(Qyy,i,j,k);
	dQyzdy = dy_(Qyz,i,j,k);
	dQzxdy = dQxzdy;
	dQzydy = dQyzdy;
	dQzzdy = -(dQxxdy+dQyydy);

	dQxxdz = dz_(Qxx,i,j,k);
	dQxydz = dz_(Qxy,i,j,k);
	dQxzdz = dz_(Qxz,i,j,k);	     
	dQyxdz = dQxydz;
	dQyydz = dz_(Qyy,i,j,k);
	dQyzdz = dz_(Qyz,i,j,k);
	dQzxdz = dQxzdz;
	dQzydz = dQyzdz;
	dQzzdz = -(dQxxdz+dQyydz);

	divQx=dQxxdx+dQxydy+dQxzdz;
	divQy=dQxydx+dQyydy+dQyzdz;
	divQz=dQxzdx+dQyzdy+dQzzdz;

	mu[i][j][k] = alpha2phi*(phi[i][j][k]*phi[i][j][k]*phi[i][j][k] - 3.0*phi[i][j][k]*phi[i][j][k] + 2.0*phi[i][j][k])
	  - Kphi*laplacian_(phi,i,j,k)+ K2phi*laplacian_(laplacianphi,i,j,k)
		    + 0.25*(-Abulk/6.0*TrQ2
			    -Abulk/3.0*(Qsqxx*Qxxl+2.0*Qsqxy*Qxyl+2.0*Qsqxz*Qxzl+Qsqyy*Qyyl+2.0*Qsqyz*Qyzl+Qsqzz*Qzzl )
			    +Abulk/4.0*TrQ2*TrQ2);
	mu[i][j][k] -= W*(divQx*dphidx[i][j][k]+divQy*dphidy[i][j][k]+divQz*dphidz[i][j][k]);
	mu[i][j][k] -= W*(Qxxl*d2phidxdx+2.0*Qxyl*d2phidxdy+2.0*Qxzl*d2phidxdz+Qyyl*d2phidydy+2.0*Qyzl*d2phidydz+Qzzl*d2phidzdz);
  
	phi_average = phi_average + phi[i][j][k];
      }
    }
  }
  phi_average = phi_average/(Lx*Ly*Lz);

  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) { 				
	iup  = iupa[i];
	iup2 = iupa[iupa[i]];

	jup  = jupa[j];
	jup2 = jupa[jupa[j]];

	kup  = kupa[k];
	kup2 = kupa[kupa[k]];

	idwn  = idwna[i];
	idwn2 = idwna[idwna[i]];

	jdwn  = jdwna[j];
	jdwn2 = jdwna[jdwna[j]];

	kdwn  = kdwna[k];
	kdwn2 = kdwna[kdwna[k]];

	vx = ux[i][j][k];
	vy = uy[i][j][k];
	vz = uz[i][j][k];

	phil = phi[i][j][k];

	/* upwind scheme */
	vxdphidx = max(vx,0.0)*(2.0*phi[iup][j][k]  + 3.0*phi[i][j][k]   - 6.0*phi[idwn][j][k] +     phi[idwn2][j][k])/6.0 +
		   min(vx,0.0)*(   -phi[iup2][j][k] + 6.0*phi[iup][j][k] - 3.0*phi[i][j][k]    - 2.0*phi[idwn][j][k]) /6.0;
	vydphidy = max(vy,0.0)*(2.0*phi[i][jup][k]  + 3.0*phi[i][j][k]   - 6.0*phi[i][jdwn][k] +     phi[i][jdwn2][k])/6.0 +
		   min(vy,0.0)*(   -phi[i][jup2][k] + 6.0*phi[i][jup][k] - 3.0*phi[i][j][k]    - 2.0*phi[i][jdwn][k]) /6.0;
	vzdphidz = max(vz,0.0)*(2.0*phi[i][j][kup]  + 3.0*phi[i][j][k]   - 6.0*phi[i][j][kdwn] +     phi[i][j][kdwn2])/6.0 +
		   min(vz,0.0)*(   -phi[i][j][kup2] + 6.0*phi[i][j][kup] - 3.0*phi[i][j][k]    - 2.0*phi[i][j][kdwn]) /6.0;

	divu = dx_(ux,i,j,k) + dy_(uy,i,j,k) + dz_(uz,i,j,k);
	
	h_phi[i][j][k] = -vxdphidx - vydphidy - vzdphidz - phil*(divu) + M*laplacian_(mu,i,j,k);
      }
    }
  }
}

/* predictor-corrector for phi */
void updatephi0(double phipr[Lx][Ly][Lz], double phiold[Lx][Ly][Lz])
{
	int i, j, k;
	double dphidz;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				phipr[i][j][k] = phiold[i][j][k] + dt*h_phi[i][j][k];

				h_phiold[i][j][k] = h_phi[i][j][k];
			}
		}
	}
}
void updatephi(double phinew[Lx][Ly][Lz], double phiold[Lx][Ly][Lz])
{
	int i, j, k;
	double dphidz;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				phinew[i][j][k] = phiold[i][j][k] + 0.5*dt*(h_phi[i][j][k] + h_phiold[i][j][k]);
			}
		}
	}
}

/* Initialization */
void initializeQ(void)
{
  int i, j, k;
  int R;
  double Px[Lx][Ly][Lz],Py[Lx][Ly][Lz],Pz[Lx][Ly][Lz];

  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {
#if DROPLET
	R = sqrt((i-Lx/2)*(i-Lx/2) + (j-Ly/2)*(j-Ly/2) + (k-Lz/2)*(k-Lz/2));

	if (R < radius) {
/*	  Px[i][j][k] = 0.0;
	  Py[i][j][k] = 1.0;
	  Pz[i][j][k] = 0.0;

	  amplitude=0.0;
	  Qxx[i][j][k]=amplitude*(Px[i][j][k]*Px[i][j][k]-1.0/3.0);
	  Qxy[i][j][k]=amplitude*(Px[i][j][k]*Py[i][j][k]);
	  Qxz[i][j][k]=amplitude*(Px[i][j][k]*Pz[i][j][k]);
	  Qyy[i][j][k]=amplitude*(Py[i][j][k]*Py[i][j][k]-1.0/3.0);
	  Qyz[i][j][k]=amplitude*(Py[i][j][k]*Pz[i][j][k]);*/

         amplitude=(0.546-0.2723/2.0);

         Qxx[i][j][k]=0.2723/2.0+amplitude*cos(2.0*q0*j);
         Qxy[i][j][k]= 0.0;
         Qyy[i][j][k]= -0.2723;
         Qxz[i][j][k]= -amplitude*(sin(2.0*q0*j));
         Qyz[i][j][k]= 0.0;

	}

	else {

	  Qxx[i][j][k]=Qxy[i][j][k]=Qxz[i][j][k]=Qyy[i][j][k]=Qyz[i][j][k]=0.0;//0.0001;

	}
#endif
      }
    }
  }
}

/* calculate molecular field for Q tensor*/
void Q_function(void)
{
  int i, j, k;
  int iup, jup, kup;
  int iup2, jup2, kup2;
  int idwn, jdwn, kdwn;
  int idwn2, jdwn2, kdwn2;

  double dQxxdx,dQxxdy,dQxxdz,dQxydx,dQxydy,dQxydz,dQyydx,dQyydy,dQyydz;
  double dQxzdx,dQxzdy,dQxzdz,dQyzdx,dQyzdy,dQyzdz;
  double dQyxdx,dQyxdy,dQyxdz;
  double dQzxdx,dQzxdy,dQzxdz,dQzydx,dQzydy,dQzydz;
  double dQzzdx,dQzzdy,dQzzdz;
  double trt,trd2Q,TrE2;
  double txx,tyy,txy,tyx,tzz,txz,tzx,tyz,tzy;
  double d2Qxxdxdx,d2Qxxdydy,d2Qxxdxdy,d2Qxxdzdz,d2Qxxdxdz,d2Qxxdydz;
  double d2Qyydxdx,d2Qyydydy,d2Qyydxdy,d2Qyydzdz,d2Qyydxdz,d2Qyydydz;
  double d2Qxydxdx,d2Qxydydy,d2Qxydxdy,d2Qxydzdz,d2Qxydxdz,d2Qxydydz;
  double d2Qxzdxdx,d2Qxzdydy,d2Qxzdxdy,d2Qxzdzdz,d2Qxzdxdz,d2Qxzdydz;
  double d2Qyzdxdx,d2Qyzdydy,d2Qyzdxdy,d2Qyzdzdz,d2Qyzdxdz,d2Qyzdydz;
  double DGxx,DGyy,DGxy,DGyx,DGzz,DGxz,DGzx,DGyz,DGzy,TrG,divQx,divQy,divQz;
  double DGchol1xx,DGchol1yy,DGchol1zz;
	
  double Qsqxx,Qsqxy,Qsqxz,Qsqyy,Qsqyz,Qsqzz,Qxxl,Qxyl,Qxzl,Qyyl,Qyzl,Qzzl;
  double Hxx,Hyy,Hxy,Hxz,Hyz,TrQ2,TrDQI;
  double mDQ4xx,mDQ4xy,mDQ4yy,mDQ4xz,mDQ4yz,mDQ4zz,nnxxl,nnyyl;

  double dphidx, dphidy, dphidz;
  double phil, gradphisq, freeE;

  double duxdx, duxdy, duxdz;
  double duydx, duydy, duydz;
  double duzdx, duzdy, duzdz;
  double divu;

  double vx, vy, vz;
  double vxdQxxdx, vydQxxdy, vzdQxxdz;
  double vxdQxydx, vydQxydy, vzdQxydz;
  double vxdQxzdx, vydQxzdy, vzdQxzdz;
  double vxdQyydx, vydQyydy, vzdQyydz;
  double vxdQyzdx, vydQyzdy, vzdQyzdz;

  freeEQbulk=0.0;
  freeEphi=0.0;
  freeEel=0.0;
  freeEcol=0.0;
  freeEanc=0.0;
  freeEce=0.0;
  freeEtot=0.0; 
  sumx=0.0;
  sumy=0.0;
  sumz=0.0;
  sumphi=0.0;
  
  
  energy = 0.0;
  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {
	iup  = iupa[i];
	iup2 = iupa[iupa[i]];

	jup  = jupa[j];
	jup2 = jupa[jupa[j]];

	kup  = kupa[k];
	kup2 = kupa[kupa[k]];

	idwn  = idwna[i];
	idwn2 = idwna[idwna[i]];
	      
	jdwn  = jdwna[j];
	jdwn2 = jdwna[jdwna[j]];
	      
	kdwn  = kdwna[k];
	kdwn2 = kdwna[kdwna[k]];

	dQxxdx = dx_(Qxx,i,j,k);
	dQxydx = dx_(Qxy,i,j,k);
	dQxzdx = dx_(Qxz,i,j,k);
	dQyxdx = dQxydx;
	dQyydx = dx_(Qyy,i,j,k);
	dQyzdx = dx_(Qyz,i,j,k);
	dQzxdx = dQxzdx;
	dQzydx = dQyzdx;
	dQzzdx = -(dQxxdx+dQyydx);

	dQxxdy = dy_(Qxx,i,j,k);
	dQxydy = dy_(Qxy,i,j,k);
	dQxzdy = dy_(Qxz,i,j,k);
	dQyxdy = dQxydy;
	dQyydy = dy_(Qyy,i,j,k);
	dQyzdy = dy_(Qyz,i,j,k);
	dQzxdy = dQxzdy;
	dQzydy = dQyzdy;
	dQzzdy = -(dQxxdy+dQyydy);

	dQxxdz = dz_(Qxx,i,j,k);
	dQxydz = dz_(Qxy,i,j,k);
	dQxzdz = dz_(Qxz,i,j,k);	     
	dQyxdz = dQxydz;
	dQyydz = dz_(Qyy,i,j,k);
	dQyzdz = dz_(Qyz,i,j,k);
	dQzxdz = dQxzdz;
	dQzydz = dQyzdz;
	dQzzdz = -(dQxxdz+dQyydz);

	Qxxl = Qxx[i][j][k];
	Qxyl = Qxy[i][j][k];
	Qyyl = Qyy[i][j][k];
	Qxzl = Qxz[i][j][k];
	Qyzl = Qyz[i][j][k];
	Qzzl = -Qxxl-Qyyl;

	phil=phi[i][j][k];

	dphidx=dx_(phi,i,j,k);
	dphidy=dy_(phi,i,j,k);
	dphidz=dz_(phi,i,j,k);
	
	divQx=dQxxdx+dQxydy+dQxzdz;
	divQy=dQxydx+dQyydy+dQyzdz;
	divQz=dQxzdx+dQyzdy-dQxxdz-dQyydz;

	gradphisq=dphidx*dphidx+dphidy*dphidy+dphidz*dphidz;

	Qsqxx = Qxxl + Qyyl;
	Qsqzz = Qsqxx*Qsqxx + Qxzl*Qxzl + Qyzl*Qyzl;
	Qsqxy = Qxyl*Qsqxx + Qxzl*Qyzl;
	Qsqxz = Qxxl*Qxzl - Qxzl*Qsqxx + Qxyl*Qyzl;
	Qsqyz = Qxyl*Qxzl - Qyzl*Qsqxx + Qyyl*Qyzl;
	Qsqxx = Qxxl*Qxxl + Qxyl*Qxyl + Qxzl*Qxzl;
	Qsqyy = Qyyl*Qyyl + Qxyl*Qxyl + Qyzl*Qyzl;

	TrQ2 = Qsqxx + Qsqyy + Qsqzz;
	
	DGchol1xx=1.0/2.0*
        (2.0*dQxxdx*dQxxdx+2.0*(dQxydx*dQxydx+dQxxdy*dQxydx)+
         2.0*(dQxzdx*dQxzdx+dQxxdz*dQxzdx)+2.0*(dQxydy*dQyydx)+
         2.0*(dQxzdy*dQyzdx+dQxydz*dQyzdx)+2.0*(-dQxzdz*dQxxdx
         -dQxzdz*dQyydx));


        DGchol1yy=1.0/2.0*
        (2.0*(dQxydx*dQxxdy)+2.0*(dQyydx+dQxydy)*dQxydy
        +2.0*(dQyzdx+dQxydz)*dQxzdy+2.0*dQyydy*dQyydy
        +2.0*(dQyzdy+dQyydz)*dQyzdy+2.0*dQyzdz*(-dQxxdy-dQyydy));


        DGchol1zz=1.0/2.0*
        (2.0*dQxzdx*dQxxdz+2.0*(dQyzdx+dQxzdy)*dQxydz
        +2.0*(-dQxxdx-dQyydx+dQxzdz)*dQxzdz+2.0*(dQyzdy*dQyydz)
        +2.0*(-dQxxdy-dQyydy+dQyzdz)*dQyzdz+2.0*(-dQxxdz-dQyydz)*
        (-dQxxdz-dQyydz));


        DGxx=(dQxxdx*dQxxdx+2.0*dQxydx*dQxydx+dQyydx*dQyydx+
              2.0*dQxzdx*dQxzdx+2.0*dQyzdx*dQyzdx+
              (dQxxdx+dQyydx)*(dQxxdx+dQyydx));

        DGyy= (dQxxdy*dQxxdy+2.0*dQxydy*dQxydy+dQyydy*dQyydy+
               2.0*dQxzdy*dQxzdy+2.0*dQyzdy*dQyzdy+
               (dQxxdy+dQyydy)*(dQxxdy+dQyydy));

        DGzz= (dQxxdz*dQxxdz+2.0*dQxydz*dQxydz+dQyydz*dQyydz+
               2.0*dQxzdz*dQxzdz+2.0*dQyzdz*dQyzdz+
               (dQxxdz+dQyydz)*(dQxxdz+dQyydz));


	/*cholesteric contribution*/

#if CHOL

     txx = 4*q0*L1*(dQxydz-dQxzdy);
     txy = 4*q0*L1*(dQyydz-dQyzdy);
     txz = 4*q0*L1*(dQyzdz+dQxxdy+dQyydy);

     tyy = 4*q0*L1*(dQyzdx-dQxydz);
     tyx = 4*q0*L1*(dQxzdx-dQxxdz);
     tyz = 4*q0*L1*(dQzzdx-dQxzdz);

     tzz = 4*q0*L1*(dQxzdy-dQyzdx);
     tzy = 4*q0*L1*(dQxydy-dQyydx);
     tzx = 4*q0*L1*(dQxxdy-dQxydx);


     trt=txx+tyy+tzz;

     /*     trd2Q=d2Qxxdxdx+2.0*d2Qxydxdy+2.0*d2Qxzdxdz+2.0*d2Qyzdydz
	    +d2Qyydydy-d2Qxxdzdz-d2Qyydzdz;*/ // TO INCLUDE WITH DIFFERENT ELASTIC CONSTANTS!!

#endif


	/* second derivative contribution to molecular field Hab */
        DEHxx[i][j][k] = L1*laplacian_(Qxx,i,j,k) +  txx -(1.0/3.0)*trt;
	DEHxy[i][j][k] = L1*laplacian_(Qxy,i,j,k) + (txy + tyx)/2.0;
	DEHxz[i][j][k] = L1*laplacian_(Qxz,i,j,k) + (txz + tzx)/2.0;
	DEHyy[i][j][k] = L1*laplacian_(Qyy,i,j,k) + tyy -(1.0/3.0)*trt;
	DEHyz[i][j][k] = L1*laplacian_(Qyz,i,j,k) + (tyz + tzy)/2.0;
	      


	/* bulk contribution to molecular field Hab */
	gammac = gamma + 0.25*phi[i][j][k];

	DEHxx[i][j][k] += Abulk*(-(1.0-gammac/3.0)*Qxxl + gammac*(Qsqxx-TrQ2/3.0)-gammac*Qxxl*TrQ2);
	DEHxy[i][j][k] += Abulk*(-(1.0-gammac/3.0)*Qxyl + gammac*Qsqxy-gammac*Qxyl*TrQ2);
	DEHyy[i][j][k] += Abulk*(-(1.0-gammac/3.0)*Qyyl + gammac*(Qsqyy-TrQ2/3.0)-gammac*Qyyl*TrQ2);
	DEHxz[i][j][k] += Abulk*(-(1.0-gammac/3.0)*Qxzl + gammac*Qsqxz-gammac*Qxzl*TrQ2);
	DEHyz[i][j][k] += Abulk*(-(1.0-gammac/3.0)*Qyzl + gammac*Qsqyz-gammac*Qyzl*TrQ2);

	DEHxx[i][j][k] -= W*(dphidx*dphidx-gradphisq/3.0);
	DEHxy[i][j][k] -= W*(dphidx*dphidy);
	DEHxz[i][j][k] -= W*(dphidx*dphidz);
	DEHyy[i][j][k] -= W*(dphidy*dphidy-gradphisq/3.0);
	DEHyz[i][j][k] -= W*(dphidy*dphidz);


	DEHxx[i][j][k] += -aa*L1*q0*q0*Qxxl;
	DEHxy[i][j][k] += -aa*L1*q0*q0*Qxyl; 
	DEHxz[i][j][k] += -aa*L1*q0*q0*Qxzl;
	DEHyy[i][j][k] += -aa*L1*q0*q0*Qyyl;
	DEHyz[i][j][k] += -aa*L1*q0*q0*Qyzl;


          /*electric field contribution*/

        /* calculate E-field */
        Ex[i][j][k]= 0.0;//delVx/(-1.0+Lx);
        Ey[i][j][k]= delVy/(-1.0+Ly);//(delVy/(-1.0+Ly)*cos(theta2*Pi/180));
	Ez[i][j][k]= delVz/(-1.0+Lz);//(delVy/(-1.0+Lz)*sin(theta2*Pi/180));

        /* E-field contribution to H*/

        TrE2 = Inv12Pi*epsa*(Ex[i][j][k]*Ex[i][j][k]+Ey[i][j][k]*Ey[i][j][k]+
                             Ez[i][j][k]*Ez[i][j][k]);

        DEHxx[i][j][k] += Inv12Pi*epsa*Ex[i][j][k]*Ex[i][j][k]-TrE2/3.0
              -efb*(Ex[i][j][k]*(dQxydy+dQxzdz)-Ey[i][j][k]*dQxydx-Ez[i][j][k]*dQxzdx);
        DEHxy[i][j][k] += Inv12Pi*epsa*Ex[i][j][k]*Ey[i][j][k]
              -efb/2.0*(Ex[i][j][k]*(dQyydy+dQyzdz-dQxxdy)+Ey[i][j][k]*(dQxxdx+dQxzdz-dQyydx)-Ez[i][j][k]*(dQzydx+dQxzdy));
        DEHyy[i][j][k] += Inv12Pi*epsa*Ey[i][j][k]*Ey[i][j][k]-TrE2/3.0
              -efb*(Ey[i][j][k]*(dQxydx+dQyzdz)-Ex[i][j][k]*dQxydy-Ez[i][j][k]*dQzydy);
        DEHxz[i][j][k] += Inv12Pi*epsa*Ex[i][j][k]*Ez[i][j][k]
              -efb/2.0*(Ex[i][j][k]*(dQyzdy+dQzzdz-dQxxdz)-Ey[i][j][k]*(dQyzdx+dQxydz)+Ez[i][j][k]*(dQxxdx+dQxydy-dQzzdx));
        DEHyz[i][j][k] += Inv12Pi*epsa*Ey[i][j][k]*Ez[i][j][k]
              -efb/2.0*(Ey[i][j][k]*(dQxzdx+dQzzdz-dQyydz)-Ex[i][j][k]*(dQxzdy+dQxydz)+Ez[i][j][k]*(dQyydy+dQxydx-dQzzdy));


	Hxx = DEHxx[i][j][k];
	Hxy = DEHxy[i][j][k];
	Hxz = DEHxz[i][j][k];
	Hyy = DEHyy[i][j][k];
	Hyz = DEHyz[i][j][k];
	      
	DEH2xx[i][j][k] = DEHxx[i][j][k];
	DEH2xy[i][j][k] = DEHxy[i][j][k];
	DEH2xz[i][j][k] = DEHxz[i][j][k];
	DEH2yy[i][j][k] = DEHyy[i][j][k];
	DEH2yz[i][j][k] = DEHyz[i][j][k];

	/* elastic energy */
	energy += dQxxdx*dQxxdx + dQxxdy*dQxxdy + dQxxdz*dQxxdz;
	energy += dQxydx*dQxydx + dQxydy*dQxydy + dQxydz*dQxydz;
	energy += dQxzdx*dQxzdx + dQxzdy*dQxzdy + dQxzdz*dQxzdz;
	energy += dQyydx*dQyydx + dQyydy*dQyydy + dQyydz*dQyydz;
	energy += dQyzdx*dQyzdx + dQyzdy*dQyzdy + dQyzdz*dQyzdz;

	/* compute  stress tensor */ /*
	DGxx = (dQxxdx*dQxxdx + 2.0*dQxydx*dQxydx + dQyydx*dQyydx+
		2.0*dQxzdx*dQxzdx + 2.0*dQyzdx*dQyzdx + 
		(dQxxdx+dQyydx)*(dQxxdx+dQyydx));
	DGyy = (dQxxdy*dQxxdy + 2.0*dQxydy*dQxydy + dQyydy*dQyydy+
		2.0*dQxzdy*dQxzdy + 2.0*dQyzdy*dQyzdy +
		(dQxxdy+dQyydy)*(dQxxdy+dQyydy));
	DGzz = (dQxxdz*dQxxdz + 2.0*dQxydz*dQxydz + dQyydz*dQyydz+
		2.0*dQxzdz*dQxzdz + 2.0*dQyzdz*dQyzdz +
		(dQxxdz+dQyydz)*(dQxxdz+dQyydz));
	DGxy = (dQxxdx*dQxxdy + 2.0*dQxydx*dQxydy + dQyydx*dQyydy+
		2.0*dQxzdx*dQxzdy + 2.0*dQyzdx*dQyzdy +
		(dQxxdx+dQyydx)*(dQxxdy+dQyydy));
	DGxz = (dQxxdx*dQxxdz + 2.0*dQxydx*dQxydz + dQyydx*dQyydz+
		2.0*dQxzdx*dQxzdz + 2.0*dQyzdx*dQyzdz +
		(dQxxdx+dQyydx)*(dQxxdz+dQyydz));
	DGyz = (dQxxdy*dQxxdz + 2.0*dQxydy*dQxydz + dQyydy*dQyydz+
		2.0*dQxzdy*dQxzdz + 2.0*dQyzdy*dQyzdz +
		(dQxxdy+dQyydy)*(dQxxdz+dQyydz));

	gammac = gamma + 0.25*phil;

	freeE  = 0.25*alpha2phi*phil*phil*(phil-2.0)*(phil-2.0) + 0.5*Kphi*gradphisq;
	freeE += Abulk/2.0*(1.0-gammac/3.0)*TrQ2;
	freeE += Abulk*gammac/4.0*TrQ2*TrQ2;
	freeE += -Abulk*gammac/3.0*(Qsqxx*Qxxl + 2.0*Qsqxy*Qxyl + 2.0*Qsqxz*Qxzl + Qsqyy*Qyyl + 2.0*Qsqyz*Qyzl + Qsqzz*Qzzl);
	freeE += L1/2.0*(dQxxdx*dQxxdx + dQxxdy*dQxxdy + dQxxdz*dQxxdz);
	freeE += L1/2.0*(dQxydx*dQxydx + dQxydy*dQxydy + dQxydz*dQxydz);
	freeE += L1/2.0*(dQxzdx*dQxzdx + dQxzdy*dQxzdy + dQxzdz*dQxzdz);
	freeE += L1/2.0*(dQyydx*dQyydx + dQyydy*dQyydy + dQyydz*dQyydz);
	freeE += L1/2.0*(dQyzdx*dQyzdx + dQyzdy*dQyzdy + dQyzdz*dQyzdz);

	DG2xx[i][j][k] = -L1*DGxx;
	DG2yy[i][j][k] = -L1*DGyy;
	DG2zz[i][j][k] = -L1*DGzz;
	DG2xy[i][j][k] = -L1*DGxy;
	DG2xz[i][j][k] = -L1*DGxz;
	DG2yz[i][j][k] = -L1*DGyz;

	DGyx = -L1*DGxy;
	DGzx = -L1*DGxz;
	DGzy = -L1*DGyz;

	DG2xx[i][j][k] += -Kphi*dphidx*dphidx - phi[i][j][k]*mu[i][j][k] + freeE;
	DG2xy[i][j][k] += -Kphi*dphidx*dphidy;
	DG2xz[i][j][k] += -Kphi*dphidx*dphidz;
	DG2yy[i][j][k] += -Kphi*dphidy*dphidy - phi[i][j][k]*mu[i][j][k] + freeE;
	DG2yz[i][j][k] += -Kphi*dphidy*dphidz;
	DG2zz[i][j][k] += -Kphi*dphidz*dphidz - phi[i][j][k]*mu[i][j][k] + freeE; */
	
		//Energia libera del sistema:
              
    //1)   Bulk term of Q in the free-energy:
           
    freeEQbulk += Abulk/2.0*(1.0-gammac/3.0)*TrQ2;
	freeEQbulk += Abulk*gammac/4.0*TrQ2*TrQ2;
	freeEQbulk += -Abulk*gammac/3.0*(Qsqxx*Qxxl + 2.0*Qsqxy*Qxyl + 
    2.0*Qsqxz*Qxzl + Qsqyy*Qyyl + 2.0*Qsqyz*Qyzl + Qsqzz*Qzzl);       
         
    //2-3) Bulk term of phi and interface:
           
   freeEphi += 0.25*alpha2phi*phil*phil*(phil-2.0)*(phil-2.0) + 
   0.5*Kphi*gradphisq;      
         
    //4)   Elastic term: (d_betaQ_alpha_beta)^2
           
   freeEel += L1/2.0*(divQx*divQx+divQy*divQy+divQz*divQz);        
           
    //5)   Cholesteric term:
           
   //L1/2*4*q0*q0*(Q_alpha_beta)^2
           
   freeEcol += L1/2.0*(4*q0*q0*(Qxxl*Qxxl+Qyyl*Qyyl+Qzzl*Qzzl+
   2.0*Qxyl*Qxyl+2.0*Qxzl*Qxzl+2.0*Qyzl*Qyzl));
   
    //L1/2*4*q0*epsilon_alpha_zita_deltad_zitaQ_delta_beta*(Q_alpha_beta)
         
   freeEcol += L1/2.0*4.0*q0*(Qxxl*dQxzdy+Qxyl*dQyzdy+Qxzl*dQzzdy
   -Qxxl*dQxydz-Qxyl*dQyydz-Qxzl*dQyzdz-Qxyl*dQxzdx-Qyyl*dQyzdx-Qyzl*dQzzdx+
   Qxyl*dQxxdz+Qyyl*dQxydz+Qyzl*dQxzdz+Qxzl*dQxydx+Qyzl*dQyydx+Qzzl*dQyzdx-
   Qxzl*dQxxdy-Qyzl*dQxydy-Qzzl*dQxzdy);
   
   //L1/2*(epsilon_alpha_zita_deltad_zitaQ_delta_beta)^2
   
   freeEcol += L1/2.0*(-DGchol1xx+DGxx-DGchol1yy+DGyy-DGchol1zz+DGzz);                     
      
         
    //5) Anchoring term on the drop:
         
    freeEanc += W*(dphidx)*Qxxl*(dphidx)+W*(dphidy)*Qyyl*(dphidy)+
    W*(dphidz)*Qzzl*(dphidz)+2.0*W*(dphidx)*Qxyl*(dphidy)+
    2.0*W*(dphidx)*Qxzl*(dphidz)+2.0*W*(dphidy)*Qyzl*(dphidz);     
         
    //6) Electric field term: 
         
    freeEce += -Inv12Pi*epsa*(Ex[i][j][k]*Ex[i][j][k]*Qxxl+
    Ey[i][j][k]*Ey[i][j][k]*Qyyl+Ez[i][j][k]*Ez[i][j][k]*Qzzl+
    2.0*Ex[i][j][k]*Ey[i][j][k]*Qxyl+2.0*Ex[i][j][k]*Ez[i][j][k]*Qxzl+
    2.0*Ey[i][j][k]*Ez[i][j][k]*Qyzl);             
           
           
    /* Center of mass of the droplet */
    
    if(phi[i][j][k]>0){
    sumphi+=phi[i][j][k];
    sumx+=phi[i][j][k]*i;
    sumy+=phi[i][j][k]*j;
    sumz+=phi[i][j][k]*k;}      

	/* this will go to the bodyforce */		
	Fh[i][j][k][0]=Hxx*dQxxdx+2.0*Hxy*dQxydx
		+2.0*Hxz*dQxzdx+Hyy*dQyydx+2.0*Hyz*dQyzdx
		+(-Hyy-Hxx)*(-dQxxdx-dQyydx);
	Fh[i][j][k][1]=Hxx*dQxxdy+2.0*Hxy*dQxydy
		+2.0*Hxz*dQxzdy+Hyy*dQyydy+2.0*Hyz*dQyzdy
		+(-Hyy-Hxx)*(-dQxxdy-dQyydy);
	Fh[i][j][k][2]=Hxx*dQxxdz+2.0*Hxy*dQxydz
		+2.0*Hxz*dQxzdz+Hyy*dQyydz+2.0*Hyz*dQyzdz
		+(-Hyy-Hxx)*(-dQxxdz-dQyydz);

	/* work out antisymmetric part of the stress tensor, this will go to the bodyforce */
	tauxy[i][j][k]= -phivr*(Qxy[i][j][k]*(DEHxx[i][j][k]-DEHyy[i][j][k])-
		 		DEHxy[i][j][k]*(Qxx[i][j][k]-Qyy[i][j][k])+
		 		DEHxz[i][j][k]*Qyz[i][j][k]-DEHyz[i][j][k]*Qxz[i][j][k])
			+phivr*(DG2xy[i][j][k]-DGyx)/2.0;
	tauxz[i][j][k]= -phivr*(Qxz[i][j][k]*(2.0*DEHxx[i][j][k]+DEHyy[i][j][k])-
		 		DEHxz[i][j][k]*(2.0*Qxx[i][j][k]+Qyy[i][j][k])+
		 		DEHxy[i][j][k]*Qyz[i][j][k]-DEHyz[i][j][k]*Qxy[i][j][k])
			+phivr*(DG2xz[i][j][k]-DGzx)/2.0;
	tauyz[i][j][k]= -phivr*(Qyz[i][j][k]*(2.0*DEHyy[i][j][k]+DEHxx[i][j][k])-
		 		DEHyz[i][j][k]*(2.0*Qyy[i][j][k]+Qxx[i][j][k])+
		 		DEHxy[i][j][k]*Qxz[i][j][k]-DEHxz[i][j][k]*Qxy[i][j][k])
			+phivr*(DG2yz[i][j][k]-DGzy)/2.0;

	duxdx = dx_(ux,i,j,k);
	duxdy = dy_(ux,i,j,k);
	duxdz = dz_(ux,i,j,k);
				
	duydx = dx_(uy,i,j,k);
	duydy = dy_(uy,i,j,k);
	duydz = dz_(uy,i,j,k);

	duzdx = dx_(uz,i,j,k);
	duzdy = dy_(uz,i,j,k);
	duzdz = dz_(uz,i,j,k);

	/* -(Q+1/3) Tr(D.(Q+1/3)) term*/
	TrDQI = -(duxdx+duydy+duzdz)/3.0
	        -(Qxxl*duxdx+Qyyl*duydy-(Qxxl+Qyyl)*duzdz+Qxyl*(duxdy+duydx)+Qxzl*(duxdz+duzdx)+Qyzl*(duydz+duzdy));
	mDQ4xy = Qxyl*TrDQI;
	mDQ4yy = nnyyl*TrDQI;
	mDQ4xz = Qxzl*TrDQI;
	mDQ4yz = Qyzl*TrDQI;
	mDQ4zz = (-Qxxl-Qyyl+1.0/3.0)*TrDQI;
	mDQ4xx = nnxxl*TrDQI;

	Hxx=G*Hxx+duxdy*Qxyl*(1.0+xi)+duxdx*2.0*nnxxl*xi+
		duydx*Qxyl*(xi-1.0)+duxdz*Qxzl*(1.0+xi)+duzdx*Qxzl*(xi-1.0)+
		2.0*xi*mDQ4xx;
	Hxy=G*Hxy+0.5*(nnxxl*(xi-1.0)+nnyyl*(1.0+xi))*duxdy+
		Qxyl*xi*(duydy+duxdx)+0.5*(nnxxl*(1.0+xi)+nnyyl*(xi-1.0))*duydx+
		0.5*duxdz*Qyzl*(1.0+xi)+0.5*duzdx*Qyzl*(xi-1.0)+
		0.5*duydz*Qxzl*(1.0+xi)+0.5*duzdy*Qxzl*(xi-1.0)+2.0*xi*mDQ4xy;
	Hxz=G*Hxz+(-Qxxl-0.5*Qyyl*(1.0+xi)+xi/3.0)*duxdz+
		Qxzl*xi*(duxdx+duzdz)+duzdx*(Qxxl+0.5*Qyyl*(1.0-xi)+xi/3.0)+
		0.5*Qyzl*((1.0+xi)*duxdy+(xi-1.0)*duydx)+
		0.5*Qxyl*((xi-1.0)*duydz+(1.0+xi)*duzdy)+2.0*xi*mDQ4xz;
	Hyy=G*Hyy+duxdy*Qxyl*(xi-1.0)+duydy*2.0*nnyyl*xi+
		duydx*Qxyl*(xi+1.0)+duydz*Qyzl*(1.0+xi)+duzdy*Qyzl*(xi-1.0)+
		2.0*xi*mDQ4yy;
	Hyz=G*Hyz+(-Qyyl-0.5*Qxxl*(1.0+xi)+xi/3.0)*duydz+
		Qyzl*xi*(duydy+duzdz)+duzdy*(Qyyl+0.5*Qxxl*(1.0-xi)+xi/3.0)+
		0.5*Qxzl*((1.0+xi)*duydx+(xi-1.0)*duxdy)+
		0.5*Qxyl*((xi-1.0)*duxdz+(1.0+xi)*duzdx)+2.0*xi*mDQ4yz;
	      
	// ADD ADVECTION NOW      
	vx = ux[i][j][k];
	vy = uy[i][j][k];
	vz = uz[i][j][k];

	/* upwind scheme */
	vxdQxxdx = max(vx,0.0)*(2.0*Qxx[iup][j][k]  + 3.0*Qxx[i][j][k]   - 6.0*Qxx[idwn][j][k] +     Qxx[idwn2][j][k])/6.0 +
		   min(vx,0.0)*(   -Qxx[iup2][j][k] + 6.0*Qxx[iup][j][k] - 3.0*Qxx[i][j][k]    - 2.0*Qxx[idwn][j][k]) /6.0;
	vydQxxdy = max(vy,0.0)*(2.0*Qxx[i][jup][k]  + 3.0*Qxx[i][j][k]   - 6.0*Qxx[i][jdwn][k] +     Qxx[i][jdwn2][k])/6.0 +
		   min(vy,0.0)*(   -Qxx[i][jup2][k] + 6.0*Qxx[i][jup][k] - 3.0*Qxx[i][j][k]    - 2.0*Qxx[i][jdwn][k]) /6.0;
	vzdQxxdz = max(vz,0.0)*(2.0*Qxx[i][j][kup]  + 3.0*Qxx[i][j][k]   - 6.0*Qxx[i][j][kdwn] +     Qxx[i][j][kdwn2])/6.0 +
		   min(vz,0.0)*(   -Qxx[i][j][kup2] + 6.0*Qxx[i][j][kup] - 3.0*Qxx[i][j][k]    - 2.0*Qxx[i][j][kdwn]) /6.0;

	vxdQxydx = max(vx,0.0)*(2.0*Qxy[iup][j][k]  + 3.0*Qxy[i][j][k]   - 6.0*Qxy[idwn][j][k] +     Qxy[idwn2][j][k])/6.0 +
		   min(vx,0.0)*(   -Qxy[iup2][j][k] + 6.0*Qxy[iup][j][k] - 3.0*Qxy[i][j][k]    - 2.0*Qxy[idwn][j][k]) /6.0;
	vydQxydy = max(vy,0.0)*(2.0*Qxy[i][jup][k]  + 3.0*Qxy[i][j][k]   - 6.0*Qxy[i][jdwn][k] +     Qxy[i][jdwn2][k])/6.0 +
		   min(vy,0.0)*(   -Qxy[i][jup2][k] + 6.0*Qxy[i][jup][k] - 3.0*Qxy[i][j][k]    - 2.0*Qxy[i][jdwn][k]) /6.0;
	vzdQxydz = max(vz,0.0)*(2.0*Qxy[i][j][kup]  + 3.0*Qxy[i][j][k]   - 6.0*Qxy[i][j][kdwn] +     Qxy[i][j][kdwn2])/6.0 +
		   min(vz,0.0)*(   -Qxy[i][j][kup2] + 6.0*Qxy[i][j][kup] - 3.0*Qxy[i][j][k]    - 2.0*Qxy[i][j][kdwn]) /6.0;

	vxdQxzdx = max(vx,0.0)*(2.0*Qxz[iup][j][k]  + 3.0*Qxz[i][j][k]   - 6.0*Qxz[idwn][j][k] +     Qxz[idwn2][j][k])/6.0 +
		   min(vx,0.0)*(   -Qxz[iup2][j][k] + 6.0*Qxz[iup][j][k] - 3.0*Qxz[i][j][k]    - 2.0*Qxz[idwn][j][k]) /6.0;
	vydQxzdy = max(vy,0.0)*(2.0*Qxz[i][jup][k]  + 3.0*Qxz[i][j][k]   - 6.0*Qxz[i][jdwn][k] +     Qxz[i][jdwn2][k])/6.0 +
		   min(vy,0.0)*(   -Qxz[i][jup2][k] + 6.0*Qxz[i][jup][k] - 3.0*Qxz[i][j][k]    - 2.0*Qxz[i][jdwn][k]) /6.0;
	vzdQxzdz = max(vz,0.0)*(2.0*Qxz[i][j][kup]  + 3.0*Qxz[i][j][k]   - 6.0*Qxz[i][j][kdwn] +     Qxz[i][j][kdwn2])/6.0 +
		   min(vz,0.0)*(   -Qxz[i][j][kup2] + 6.0*Qxz[i][j][kup] - 3.0*Qxz[i][j][k]    - 2.0*Qxz[i][j][kdwn]) /6.0;
	      
	vxdQyydx = max(vx,0.0)*(2.0*Qyy[iup][j][k]  + 3.0*Qyy[i][j][k]   - 6.0*Qyy[idwn][j][k] +     Qyy[idwn2][j][k])/6.0 +
		   min(vx,0.0)*(   -Qyy[iup2][j][k] + 6.0*Qyy[iup][j][k] - 3.0*Qyy[i][j][k]    - 2.0*Qyy[idwn][j][k]) /6.0;
	vydQyydy = max(vy,0.0)*(2.0*Qyy[i][jup][k]  + 3.0*Qyy[i][j][k]   - 6.0*Qyy[i][jdwn][k] +     Qyy[i][jdwn2][k])/6.0 +
		   min(vy,0.0)*(   -Qyy[i][jup2][k] + 6.0*Qyy[i][jup][k] - 3.0*Qyy[i][j][k]    - 2.0*Qyy[i][jdwn][k]) /6.0;
	vzdQyydz = max(vz,0.0)*(2.0*Qyy[i][j][kup]  + 3.0*Qyy[i][j][k]   - 6.0*Qyy[i][j][kdwn] +     Qyy[i][j][kdwn2])/6.0 +
		   min(vz,0.0)*(   -Qyy[i][j][kup2] + 6.0*Qyy[i][j][kup] - 3.0*Qyy[i][j][k]    - 2.0*Qyy[i][j][kdwn]) /6.0;

	vxdQyzdx = max(vx,0.0)*(2.0*Qyz[iup][j][k]  + 3.0*Qyz[i][j][k]   - 6.0*Qyz[idwn][j][k] +     Qyz[idwn2][j][k])/6.0 +
		   min(vx,0.0)*(   -Qyz[iup2][j][k] + 6.0*Qyz[iup][j][k] - 3.0*Qyz[i][j][k]    - 2.0*Qyz[idwn][j][k]) /6.0;
	vydQyzdy = max(vy,0.0)*(2.0*Qyz[i][jup][k]  + 3.0*Qyz[i][j][k]   - 6.0*Qyz[i][jdwn][k] +     Qyz[i][jdwn2][k])/6.0 +
		   min(vy,0.0)*(   -Qyz[i][jup2][k] + 6.0*Qyz[i][jup][k] - 3.0*Qyz[i][j][k]    - 2.0*Qyz[i][jdwn][k]) /6.0;
	vzdQyzdz = max(vz,0.0)*(2.0*Qyz[i][j][kup]  + 3.0*Qyz[i][j][k]   - 6.0*Qyz[i][j][kdwn] +     Qyz[i][j][kdwn2])/6.0 +
		   min(vz,0.0)*(   -Qyz[i][j][kup2] + 6.0*Qyz[i][j][kup] - 3.0*Qyz[i][j][k]    - 2.0*Qyz[i][j][kdwn]) /6.0;

	Hxx-=vxdQxxdx+vydQxxdy+vzdQxxdz;
	Hxy-=vxdQxydx+vydQxydy+vzdQxydz;
	Hxz-=vxdQxzdx+vydQxzdy+vzdQxzdz;
	Hyy-=vxdQyydx+vydQyydy+vzdQyydz;
	Hyz-=vxdQyzdx+vydQyzdy+vzdQyzdz;

	DEHxx[i][j][k]=Hxx;
	DEHxy[i][j][k]=Hxy;
	DEHxz[i][j][k]=Hxz;
	DEHyy[i][j][k]=Hyy;
	DEHyz[i][j][k]=Hyz;
      }
    }
  }
  energy = energy/(Lx*Ly*Lz);
  freeEtot=freeEQbulk+freeEphi+freeEel+freeEcol+freeEanc+freeEce;
  xcm=sumx/sumphi;
  ycm=sumy/sumphi;
  zcm=sumz/sumphi;
}

/* predictor-corrector of time evolution of Q tensor */
void updateQ0(double Qxxpr[Lx][Ly][Lz],  double Qxypr[Lx][Ly][Lz],  double Qxzpr[Lx][Ly][Lz], double Qyypr[Lx][Ly][Lz],  double Qyzpr[Lx][Ly][Lz],
	      double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz])
{
	int i, j, k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				Qxxpr[i][j][k] = Qxxold[i][j][k] + dt*DEHxx[i][j][k];
				Qxypr[i][j][k] = Qxyold[i][j][k] + dt*DEHxy[i][j][k];
				Qxzpr[i][j][k] = Qxzold[i][j][k] + dt*DEHxz[i][j][k];
				Qyypr[i][j][k] = Qyyold[i][j][k] + dt*DEHyy[i][j][k];
				Qyzpr[i][j][k] = Qyzold[i][j][k] + dt*DEHyz[i][j][k];

				DEHxxold[i][j][k]=DEHxx[i][j][k];
				DEHxyold[i][j][k]=DEHxy[i][j][k];
				DEHxzold[i][j][k]=DEHxz[i][j][k];
				DEHyyold[i][j][k]=DEHyy[i][j][k];
				DEHyzold[i][j][k]=DEHyz[i][j][k];

			}
		}
	}
}
void updateQ(double Qxxnew[Lx][Ly][Lz], double Qxynew[Lx][Ly][Lz], double Qxznew[Lx][Ly][Lz], double Qyynew[Lx][Ly][Lz], double Qyznew[Lx][Ly][Lz],
	     double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz])
{
	int i, j, k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				Qxxnew[i][j][k] = Qxxold[i][j][k] + 0.5*dt*(DEHxxold[i][j][k] + DEHxx[i][j][k]);
				Qxynew[i][j][k] = Qxyold[i][j][k] + 0.5*dt*(DEHxyold[i][j][k] + DEHxy[i][j][k]);
				Qxznew[i][j][k] = Qxzold[i][j][k] + 0.5*dt*(DEHxzold[i][j][k] + DEHxz[i][j][k]);
				Qyynew[i][j][k] = Qyyold[i][j][k] + 0.5*dt*(DEHyyold[i][j][k] + DEHyy[i][j][k]);
				Qyznew[i][j][k] = Qyzold[i][j][k] + 0.5*dt*(DEHyzold[i][j][k] + DEHyz[i][j][k]);
			}
		}
	}
}

/*
 *	Plot to a file
 */

/* Q tensor output */
void plotQ(void)
{
  int i, j, k;

  int nrots,emax,enxt;
  double m[3][3],d[3],v[3][3];

  sprintf(filename1, "plot%f.txt", friction);

  output1 = fopen(filename1, "w");

  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {
        m[0][0]=Qxx[i][j][k];
	m[0][1]=Qxy[i][j][k];
	m[0][2]=Qxz[i][j][k];
	m[1][0]=Qxy[i][j][k];
	m[1][1]=Qyy[i][j][k];
	m[1][2]=Qyz[i][j][k];
	m[2][0]=Qxz[i][j][k];
	m[2][1]=Qyz[i][j][k];
	m[2][2]= -(m[0][0]+m[1][1]);
	jacobi(m,d,v,&nrots);

	if (d[0] > d[1]) {
	  emax=0;
	  enxt=1;
	}
	else {
	  emax=1;
	  enxt=0;
	}
	if (d[2] > d[emax]) {
	  emax=2;
	}
	else if (d[2] > d[enxt]) {
	  enxt=2;
	}
			  
	fprintf(output1, "%d %d %d %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E\n",
		i, j, k,Qxx[i][j][k],Qxy[i][j][k],Qxz[i][j][k],Qyy[i][j][k],Qyz[i][j][k],
		v[0][emax], v[1][emax], v[2][emax], d[emax],
		ux[i][j][k], uy[i][j][k], uz[i][j][k],
		mu[i][j][k],
		density[i][j][k],
		phi[i][j][k],
               f[i][j][k][0],f[i][j][k][1],f[i][j][k][2],f[i][j][k][3],f[i][j][k][4], 
               f[i][j][k][5],f[i][j][k][6],f[i][j][k][7],f[i][j][k][8],f[i][j][k][9],
               f[i][j][k][10],f[i][j][k][11],f[i][j][k][12],f[i][j][k][13],f[i][j][k][14]);
      }
      fprintf(output1, "\n");
    }
    fprintf(output1, "\n");
  }
  fclose(output1);
}


void multipleplot(int n)
{
	int i, j, k;
       int nrots,emax,enxt;
       double m[3][3],d[3],v[3][3];

	sprintf(filename1, "plot_%d.dat", n);
	sprintf(filename2, "freeE.dat");
	sprintf(filename3, "cm.dat");

	output1 = fopen(filename1, "w");
	output2 = fopen(filename2, "a");
	output3 = fopen(filename3, "a");
	
	 fprintf(output1,"%s %d %s %d %s %d %s\n","VARIABLES = \n \
\"X     \",\n \
\"Y    \",\n \
\"Z    \",\n \
\"Qxx   \",\n \
\"Qxy   \",\n \
\"Qxz   \",\n \
\"Qyy   \",\n \
\"Qyz   \",\n \
\"v0   \",\n \
\"v1   \",\n \
\"v2   \",\n \
\"dmax \",\n \
\"ux   \",\n \
\"uy   \",\n \
\"uz   \",\n \
\"mu \",\n \
\"density  \",\n \
\"phi   \",\n \
\"Ey   \",\n \
\"Ez   \",\n \
\"theta2   \",\n \
\"f0 \",\n \
\"f1 \",\n \
\"f2 \",\n \
\"f3 \",\n \
\"f4    \",\n \
\"f5    \",\n \
\"f6    \",\n \
\"f7   \",\n \
\"f8   \",\n \
\"f9   \",\n \
\"f10   \",\n \
\"f11   \",\n \
\"f12   \",\n \
\"f13   \",\n \
\"f14   \",\n \
ZONE T=\"   1\"   I=",    Lx,  "J=",Ly,  "K=",Lz,"F=POINT \n");



  for (i = 0; i < Lx; i++) {
    for (j = 0; j < Ly; j++) {
      for (k = 0; k < Lz; k++) {
        m[0][0]=Qxx[i][j][k];
	m[0][1]=Qxy[i][j][k];
	m[0][2]=Qxz[i][j][k];
	m[1][0]=Qxy[i][j][k];
	m[1][1]=Qyy[i][j][k];
	m[1][2]=Qyz[i][j][k];
	m[2][0]=Qxz[i][j][k];
	m[2][1]=Qyz[i][j][k];
	m[2][2]= -(m[0][0]+m[1][1]);
	jacobi(m,d,v,&nrots);

	if (d[0] > d[1]) {
	  emax=0;
	  enxt=1;
	}
	else {
	  emax=1;
	  enxt=0;
	}
	if (d[2] > d[emax]) {
	  emax=2;
	}
	else if (d[2] > d[enxt]) {
	  enxt=2;
	}
			  
	fprintf(output1, "%d %d %d %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E\n",
		i, j, k,Qxx[i][j][k],Qxy[i][j][k],Qxz[i][j][k],Qyy[i][j][k],Qyz[i][j][k],
		v[0][emax], v[1][emax], v[2][emax], d[emax],
		ux[i][j][k], uy[i][j][k], uz[i][j][k],
		mu[i][j][k],
		density[i][j][k],
		phi[i][j][k], Ey[i][j][k], Ez[i][j][k], theta2,
               f[i][j][k][0],f[i][j][k][1],f[i][j][k][2],f[i][j][k][3],f[i][j][k][4], 
               f[i][j][k][5],f[i][j][k][6],f[i][j][k][7],f[i][j][k][8],f[i][j][k][9],
               f[i][j][k][10],f[i][j][k][11],f[i][j][k][12],f[i][j][k][13],f[i][j][k][14]);
      }
      fprintf(output1, "\n");
    }
    fprintf(output1, "\n");
  }
  fclose(output1);
  
  fprintf(output2, "%d %f %f %f %f %f %f %f\n",n,freeEQbulk,freeEphi,freeEel,
  freeEcol,freeEanc,freeEce,freeEtot);
  fclose(output2); 
  
  fprintf(output3, "%d %f %f %f\n",n,xcm,ycm,zcm);
  fclose(output3);
  
}




void correction(int n)
{
	double momentum[3];
	int i,j,k;

	momentum[0] = 0.0;
	momentum[1] = 0.0;
	momentum[2] = 0.0;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				momentum[0] = momentum[0] + density[i][j][k]*ux[i][j][k];
				momentum[1] = momentum[1] + density[i][j][k]*uy[i][j][k];
				momentum[2] = momentum[2] + density[i][j][k]*uz[i][j][k];

				phi[i][j][k] += phi_average0 - phi_average;
				//if (phi[i][j][k] < 0.0) {phi[i][j][k] = 0.0;}
			}
		}
	}
	globalP[0] = momentum[0]/(Lx*Ly*Lz);
	globalP[1] = momentum[1]/(Lx*Ly*Lz);
	globalP[2] = momentum[2]/(Lx*Ly*Lz);

	if(n % stepskip == 0) {
		output1 = fopen("correction.txt", "a");
		fprintf(output1, "%d %E %E %E \n", n, momentum[0], momentum[1], momentum[2]);
		fclose(output1);
	}
}

/* add perturbation for Q tensor */

void add_Qperturbation(void)
{
	int i,j,k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				
				Qyy[i][j][k] += 0.01*(0.5 - drand48())*2.0;
				Qyz[i][j][k] += 0.01*(0.5 - drand48())*2.0;
			}
		}
	}
} 

/*
 *	mathematical functions
 */
double max(double x, double y)
{
	if (x > y) {
		return x;
	}
	else {
		return y;
	}
}
double min(double x, double y)
{
	if (x > y) {
		return y;
	}
	else {
		return x;
	}
}

/* 3D stencil */
double dx_(double phi[Lx][Ly][Lz], int i, int j, int k)
{
    	double dxphi;
    	int iup = iupa[i];
    	int jup = jupa[j];
    	int kup = kupa[k];

    	int idwn = idwna[i];
    	int jdwn = jdwna[j];
    	int kdwn = kdwna[k];

    	dxphi = (                         phi[iup][j][kup]
    	        + phi[iup][jdwn][k] + 2.0*phi[iup][j][k]    + phi[iup][jup][k]
    	                            +     phi[iup][j][kdwn]                   )/12.0
    	      - (                          phi[idwn][j][kup]
    	        + phi[idwn][jdwn][k] + 2.0*phi[idwn][j][k]    + phi[idwn][jup][k]
    	                             +     phi[idwn][j][kdwn]                    )/12.0;

    	return dxphi;
}
double dy_(double phi[Lx][Ly][Lz], int i, int j, int k)
{
    	double dyphi;
    	int iup = iupa[i];
    	int jup = jupa[j];
    	int kup = kupa[k];

    	int idwn = idwna[i];
    	int jdwn = jdwna[j];
    	int kdwn = kdwna[k];

    	dyphi = (                         phi[i][jup][kup]
    	        + phi[idwn][jup][k] + 2.0*phi[i][jup][k]    + phi[iup][jup][k]
    	                            +     phi[i][jup][kdwn]                   )/12.0
    	      - (                          phi[i][jdwn][kup]
    	        + phi[idwn][jdwn][k] + 2.0*phi[i][jdwn][k]    + phi[iup][jdwn][k]
    	                             +     phi[i][jdwn][kdwn]                    )/12.0;

    	return dyphi;
}
double dz_(double phi[Lx][Ly][Lz], int i, int j, int k)
{
    	double dzphi;
    	int iup = iupa[i];
    	int jup = jupa[j];
    	int kup = kupa[k];

    	int idwn = idwna[i];
    	int jdwn = jdwna[j];
    	int kdwn = kdwna[k];

    	dzphi = (                         phi[i][jup][kup]
    	        + phi[idwn][j][kup] + 2.0*phi[i][j][kup]    + phi[iup][j][kup]
    	                            +     phi[i][jdwn][kup]                   )/12.0
    	      - (                          phi[i][jup][kdwn]
    	        + phi[idwn][j][kdwn] + 2.0*phi[i][j][kdwn]    + phi[iup][j][kdwn]
    	                             +     phi[i][jdwn][kdwn]                    )/12.0;

	return dzphi;
}
double laplacian_(double phi[Lx][Ly][Lz], int i , int j , int k)
{
   	double laplacianphi;
   	int iup = iupa[i];
   	int jup = jupa[j];
   	int kup = kupa[k];
	
	int idwn = idwna[i];
	int jdwn = jdwna[j];
	int kdwn = kdwna[k];

    	laplacianphi = (                         phi[iup][j][kup]
        	       + phi[iup][jdwn][k] + 2.0*phi[iup][j][k]    + phi[iup][jup][k]
                       	                   +     phi[iup][j][kdwn]                   )/6.0

                     + (     phi[i][jdwn][kup]  +  2.0*phi[i][j][kup]  +     phi[i][jup][kup]
                       + 2.0*phi[i][jdwn][k]    - 24.0*phi[i][j][k]    + 2.0*phi[i][jup][k]
                       +     phi[i][jdwn][kdwn] +  2.0*phi[i][j][kdwn] +     phi[i][jup][kdwn])/6.0

                     + (                          phi[idwn][j][kup]
                       + phi[idwn][jdwn][k] + 2.0*phi[idwn][j][k]    + phi[idwn][jup][k]
                       	                    +     phi[idwn][j][kdwn]                    )/6.0;

	return laplacianphi;
}
void plotbifurcation(void)
{
	output1 = fopen("bifurcation.txt","a");

	fprintf(output1, "%E %E %E %E %d \n", friction, Vy, Vz, umax, R1);

	fclose(output1);
}
void plotV(int n)
{
	sprintf(filename1, "V%f.txt", friction);
	output1 = fopen(filename1, "a");

	fprintf(output1, "%d %E %E \n", n, Vy, Vz);

	fclose(output1);
}
void plotumax(int n)
{
	int i,j,k;

	double umag;
	umax = 0.0;
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				umag = sqrt(ux[i][j][k]*ux[i][j][k] + uy[i][j][k]*uy[i][j][k] + uz[i][j][k]*uz[i][j][k]);

				if (umag > umax) {umax = umag;}
			}
		}
	}
	sprintf(filename1, "umax%f.txt", friction);
	output1 = fopen(filename1, "a");

	fprintf(output1, "%d %E \n", n, umax);

	fclose(output1);
}
/*
 *	Stencils for 2DQ9 lattice vectors
 *
				dphidy = ( phi[i][jup][kup]  + 4.0*phi[i][jup][k]  + phi[i][jup][kdwn]
					   - phi[i][jdwn][kup] - 4.0*phi[i][jdwn][k] - phi[i][jdwn][kdwn])/12.0;

				dphidz = ( phi[i][jup][kup]  + 4.0*phi[i][j][kup]  + phi[i][jdwn][kup]
					   - phi[i][jup][kdwn] - 4.0*phi[i][j][kdwn] - phi[i][jdwn][kdwn])/12.0;

				laplacianphi = (     phi[i][jdwn][kup]  +  4.0*phi[i][j][kup]  +     phi[i][jup][kup]
					         + 4.0*phi[i][jdwn][k]    - 20.0*phi[i][j][k]    + 4.0*phi[i][jup][k]
					         +     phi[i][jdwn][kdwn] +  4.0*phi[i][j][kdwn] +     phi[i][jup][kdwn])/6.0;

 */


/* subroutine jacobi to find eigenvalues and eigenvectors from numerical recipes */


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
#define n 3
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c;
  double b[n],z[n];

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip< n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
	    && (fabs(d[iq])+g) == fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((fabs(h)+g) == fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=0;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=0;j<n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  printf("Too many iterations in routine jacobi");
  exit(0);
}
#undef n
#undef ROTATE

