#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"

typedef struct
{
	double Re,Im;
} dcomplex;

typedef unsigned int uint;
uint _output_file_index = 1;
int print_lo_shit = 0;

extern double p2nsma_(double*,int*);
extern double p2nspa_(double*,int*);
extern double p2nssa_(double*,int*);
extern double p2nsb_(double*,int*);
extern double p2nsmc_(double*,int*);
extern double p2nspc_(double*,int*);
extern double p2psa_(double*,int*);
extern double p2qga_(double*,int*);
extern double p2gqa_(double*,int*);
extern double p2gga_(double*,int*);
extern double p2ggb_(double*,int*);
extern double p2ggc_(double*,int*);
extern void hplog_(double*,int*,dcomplex*,dcomplex*,dcomplex*,dcomplex*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*);
extern void a02m_(double*,double*,double*,double*,int*,int*,int*,int*);
extern void struc_(double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,double*);

#define PI 3.1415926535897932384626433
const double Pi=PI;
const double Pi2=PI*PI;
const double Zeta2=PI*PI/6.;
const double Zeta3=1.2020569031595942854;
const double MZ=91.1876;
const double Cf=4./3.;
const double NC=3.;
const double TR=0.5;

double X[GRID_PTS];
double XG[NGP],WG[NGP];
int nf;
double Nf,beta0,beta1,beta2,log_mf2_mr2;
double alpha_pre[8]={0.},alpha_post[8]={0.};

static char const *output_prefix = "output/output-candiav1/";
typedef double (*SplitFunc)(int, double);
void OutputFileSplitFuncs(char const* filepath, SplitFunc funcs[], int num_funcs);
void OutputFileAlphaSBetas(char const* filepath);
void OutputFileInitialDistributions(char const* filepath);
void OutputFileGauleg(char const* filepath);
void OutputFileConvolution(char const* filepath, double vec[]);
void OutputFileInterpolation(char const* filepath, double vec[]);
void OutputLOCoefficients(double ***A);

void *Malloc(size_t size);
void *Calloc(size_t nmemb,size_t size);
void gauleg(double x1,double x2,double *x,double *w,int n);
double interp(double *A,double x);
double polint(double *xa,double *ya,int n,double x);
double convolution(int i,double kernel(int,double),double *A);
double RecRel_A(double *A,int k,double P0(int,double), int singlet);
double RecRel_B(double *A,double *B,int k,double P0(int,double),double P1(int,double));
double RecRel_C(double *A,double *B,double *C,int k,
                double P0(int,double),double P1(int,double),double P2(int,double));
double RecRel_Diag(double *C,int k,double P0(int,double), int nnlo);
double RecRel_Vert1(double *B,int k,double P0(int,double),double P1(int,double));
double RecRel_Vert2(double *C,int k,double P0(int,double),double P1(int,double),double P2(int,double));
double RecRel_Horiz(double *C,int k,double P0(int,double),double P1(int,double));
double Beta0(int f);
double Beta1(int f);
double Beta2(int f);
double Beta(int order,double alpha);
double alpha_rk(int order,double alpha0,double Qi,double Qf);
double post(double pre,int order);
double pre(double post,int order);
double Li2(double x);
double Li3(double x);
double S12(double x);
double fact(int n);
double xuv_1(int order,double x);
double xdv_1(int order,double x);
double xS_1(int order,double x);
double xg_1(int order,double x);
double xD_1(int order,double x);
double xuv_2(int order,double x);
double xus_2(int order,double x);
double xdv_2(int order,double x);
double xds_2(int order,double x);
double xss_2(int order,double x);
double xg_2(int order,double x);
double xuv_3(double x);
double xdv_3(double x);
double xg_3(double x);
double xdb_3(double x);
double xub_3(double x);
double xs_3(double x);
double P0NS(int i,double x);
double P0qq(int i,double x);
double P0qg(int i,double x);
double P0gq(int i,double x);
double P0gg(int i,double x);
double P1NSm(int i,double x);
double P1NSp(int i,double x);
double P1qq(int i,double x);
double P1qg(int i,double x);
double P1gq(int i,double x);
double P1gg(int i,double x);
double P2NSm(int i,double x);
double P2NSp(int i,double x);
double P2NSv(int i,double x);
double P2qq(int i,double x);
double P2qg(int i,double x);
double P2gq(int i,double x);
double P2gg(int i,double x);
double A2ns(int i,double z);
double A2gq(int i,double z);
double A2gg(int i,double z);
double A2hq(int i,double z);
double A2hg(int i,double z);




int main(int argc, char *argv[])
{	
	double lstep;
	const char *id[]={"g" ,"u" ,"d" ,"s" ,"c" ,"b" ,"t" ,
	                       "au","ad","as","ac","ab","at",
	                       "um","dm","sm","cm","bm","tm",
	                       "up","dp","sp","cp","bp","tp",
	                  "qm",     "dd","sd","cd","bd","td",
	                  "qp",     "ds","ss","cs","bs","ts"};
	double *Ntab;
	int *ntab;
	int i,j,k,l,s,t,n,nfi,nff,nf1,order,orderplus1,trunc,input,n_in=0,mode,one=1,npdf,npar,qct=0,vfns;
	double aux,aux2,beta,alpha0,alpha1,L1,L2,L3,q2,kr,mu0;
	double xuv,xdv,xub,xdb,xs,xc,xb,xg;
	double pdfs[10],dpdfs[10][17];
	double Q[9]={0.};
	double ****S;
	double **A[37];
	double ***B[37];
	double ****C[37];
	double *R[37];
	char extension[4],filename[20];
	FILE *fp;


	//User inputs and initialization
	if (argc!=7)
	{
		printf("\nUSAGE\n\n");
		printf("./candia.x <perturbative_order> <truncation_index> <input_model> <kr> <fns> <ext>\n\n");
		printf("<perturbative_order> can be 0, 1 or 2\n");
		printf("<truncation_index> cannot be less than <perturbative_order>\n");
		printf("<input_model>:\n");
		printf("   0= Les Houches toy model\n");
		printf("   1= MRST parametrization\n");
		printf("   2= MRST grid at 1.25 GeV^2 (minimal value in their grid)\n");
		printf("   3= Alekhin parametrization\n");
		printf("   4= Alekhin grid at 1 GeV^2\n");
		printf("<kr> is the ratio mu_r^2 / mu_f^2\n");
		printf("<fns>:\n");
		printf("   0= fixed flavor number scheme\n");
		printf("   1= variable flavor number scheme\n");
		printf("<ext> is the extension of the output files (max 3 characters allowed)\n\n");

		return 0;
	}

	//Perturbative order
	order=atoi(argv[1]);
	orderplus1=order+1;

	if (order<0 || order>2)
	{
		printf("\nInvalid argument for <order>\n");
		return 0;
	}

	kr = atof(argv[4]);

	if (kr==0.)
	{
		printf("\nInvalid argument for <kr>\n");
		return 0;
	}
	
	//Truncation index
	if (order>0)
	{
		trunc=atoi(argv[2]);

		if (trunc<order)
		{
			printf("\nInvalid argument for <truncation_index>\n");
			return 0;
		}
	}
	else
		trunc=0;

	Ntab=Malloc((sizeof(xtab)/sizeof(double))*sizeof(double));
	ntab=Calloc(sizeof(xtab)/sizeof(double),sizeof(int));


	
	S = malloc((trunc+1) * (sizeof *S));
	for (uint t=0; t<trunc+1; t++)
	{
		S[t] = malloc(2 * (sizeof *S[t]));
		for (uint i=0; i<2; i++)
		{
			S[t][i] = malloc(ITERATIONS * (sizeof *S[t][i]));
			for (uint n=0; n<ITERATIONS; n++)
				S[t][i][n] = Calloc(GRID_PTS, sizeof *S[t][i][n]);
		}
	}

	for (uint j=0; j<37; j++)
	{
		A[j] = malloc(ITERATIONS * (sizeof *A[j]));
		B[j] = malloc(ITERATIONS * (sizeof *B[j]));
		C[j] = malloc(ITERATIONS * (sizeof *C[j]));

		R[j] = calloc(GRID_PTS, sizeof *R[j]);

		for (uint n=0; n<ITERATIONS; n++)
		{
			A[j][n] = calloc(GRID_PTS, sizeof *A[j][n]);
			B[j][n] = malloc(ITERATIONS * (sizeof *B[j][n]));
			C[j][n] = malloc(ITERATIONS * (sizeof *C[j][n]));

			for (uint m=0; m<ITERATIONS; m++)
			{
				B[j][n][m] = calloc(GRID_PTS, sizeof *B[j][n][m]);
				C[j][n][m] = malloc(ITERATIONS * (sizeof *C[j][n][m]));

				for (uint l=0; l<ITERATIONS; l++)
					C[j][n][m][l] = calloc(GRID_PTS, sizeof *C[j][n][m][l]);
			}
		}
	}

	//Input model
	input=atoi(argv[3]);

	switch (input)
	{
		case 0: //Les Houches toy model

			nfi=3;
			Q[3]=sqrt(2.);
			Q[4]=sqrt(2.);
			Q[5]=4.5;
			Q[6]=175.;
			alpha1=0.35; //initial alpha_s at mu_f=sqrt(2.) GeV

			break;

		case 1: //MRST parametrization

			nfi=3;
			Q[3]=1.;
			Q[4]=1.43;
			Q[5]=4.3;
			Q[6]=175.;

			switch (order)
			{
				case 0:
					alpha1=0.13; //alpha_s(MZ)
					break;
				case 1:
					alpha1=0.119; //alpha_s(MZ)
					break;
				case 2:
					alpha1=0.1155; //alpha_s(MZ)
					break;
			}

			break;

		case 2: //MRST grid

			nfi=3;
			Q[3]=sqrt(1.25);
			Q[4]=1.43;
			Q[5]=4.3;
			Q[6]=175.;

			switch (order)
			{
				case 0:
					alpha1=0.13; //alpha_s(MZ)
					break;
				case 1:
					alpha1=0.119; //alpha_s(MZ)
					break;
				case 2:
					alpha1=0.1155; //alpha_s(MZ)
					break;
			}

			break;

		case 3: //Alekhin parametrization (BAD!!! charm missing!)

			nfi=3;
			Q[3]=3.;
			Q[4]=3.;
			Q[5]=4.5;
			Q[6]=180.;

			break;

		case 4: //Alekhin grid

			nfi=3;
			Q[3]=1.;
			Q[4]=1.5;
			Q[5]=4.5;
			Q[6]=180.;

			break;
	}
	//Flavor number scheme
	vfns=atoi(argv[5]);
	if (vfns<0 || vfns>1)
	{
		printf("\nInvalid argument for <fns>\n");
		return 0;
	}

	if (!vfns)
		for (i=nfi+1;i<=6;i++)
			if (Q[i]!=Q[nfi])
				Q[i]=2.*Qtab[sizeof(Qtab)/sizeof(double)-1];

	//Grid definition
	aux=-log10(xtab[0]);

	for (i=1;i<sizeof(xtab)/sizeof(double);i++)
		Ntab[i]=(double)(GRID_PTS-1)*log10(xtab[i]/xtab[i-1])/aux;

	ntab[0]=GRID_PTS-1;

	for (i=1;i<sizeof(xtab)/sizeof(double);i++)
	{
		ntab[i]=(int)Ntab[i];
		Ntab[i]-=(double)ntab[i];
		ntab[0]-=ntab[i];
	}

	for (i=1;i<sizeof(xtab)/sizeof(double);i++)
		if (ntab[i]==0)
		{
			ntab[i]=1;
			Ntab[i]-=1.;
			ntab[0]-=1;
		}

	for (;ntab[0]<0;ntab[0]++)
	{
		n=0;

		for (i=1;i<sizeof(xtab)/sizeof(double);i++)
			if (ntab[i]!=1)
				if ((n==0) || (Ntab[i]<=aux))
				{
					n=i;
					aux=Ntab[i];
				}

		ntab[n]--;
		Ntab[n]+=1.;
	}


	for (;ntab[0]>0;ntab[0]--)
	{
		n=0;

		for (i=1;i<sizeof(xtab)/sizeof(double);i++)
			if ((n==0) || (Ntab[i]>aux))
			{
				n=i;
				aux=Ntab[i];
			}

		ntab[n]++;
		Ntab[n]-=1.;
	}

	for (i=1;i<sizeof(xtab)/sizeof(double);i++)
		ntab[i]+=ntab[i-1];

	for (i=0;i<sizeof(xtab)/sizeof(double)-1;i++)
	{
		lstep=log10(xtab[i+1]/xtab[i])/(double)(ntab[i+1]-ntab[i]);

		for (j=ntab[i];j<ntab[i+1];j++)
			X[j]=xtab[i]*pow(10.,lstep*(double)(j-ntab[i]));
	}

	X[GRID_PTS-1]=1.;


	//Computation of Gauss-Legendre point and weights
	gauleg(0.,1.,XG,WG,NGP);

	if (Qtab[0]<Q[nfi])
	{
		printf("The final values of Q cannot be less than the initial Q (%.9g in this case)\n",Q[nfi]);
		printf("Please edit the Qtab array at the beginning of the main function\n");
		return 0;
	}

	aux=Qtab[sizeof(Qtab)/sizeof(double)-1];

	for (nff=6;aux<=Q[nff];nff--);

	if (aux>Q[6])
		i=7;
	else
		for (i=nfi+1;aux>Q[i];i++);

	Q[i]=aux;
	for (j=i+1;j<=8;j++)
		Q[j]=0.;

	sprintf(extension, "dat");

	//Writing the output summary file
	sprintf(filename,"Q0.%s",extension);
	fp=fopen(filename,"w");

	fprintf(fp,"file                     Q             Q^2\n");
	fprintf(fp,"------------------------------------------\n");

	for (i=1;i<=sizeof(Qtab)/sizeof(double);i++)
	{
		sprintf(filename,"Q%d.%s",i,extension);
		fprintf(fp,"%-10s %15.9g %15.9g\n",filename,Qtab[i-1],Qtab[i-1]*Qtab[i-1]);
	}

	fclose(fp);   


	//Setting PDF's initial values
	for (i=0;i<=GRID_PTS-2;i++)
	{
		aux=X[i];

		switch (input)
		{
			case 0: //Les Houches toy model

				//Singlet
				S[0][0][0][i]=xg_3(aux);
				S[0][1][0][i]=xuv_3(aux)+2.*xub_3(aux)+xdv_3(aux)+2.*xdb_3(aux)+2.*xs_3(aux);

				switch (order)
				{
					case 0:
						//LO nonsinglet
						A[7][0][i]=xub_3(aux);
						A[1][0][i]=xuv_3(aux)+A[7][0][i];
						A[8][0][i]=xdb_3(aux);
						A[2][0][i]=xdv_3(aux)+A[8][0][i];
						A[3][0][i]=A[9][0][i]=xs_3(aux);

						break;

					case 1:
						//NLO nonsinglet
						B[7][0][0][i]=xub_3(aux);
						B[1][0][0][i]=xuv_3(aux)+B[7][0][0][i];
						B[8][0][0][i]=xdb_3(aux);
						B[2][0][0][i]=xdv_3(aux)+B[8][0][0][i];
						B[3][0][0][i]=B[9][0][0][i]=xs_3(aux);

						break;

					case 2:
						//NNLO nonsinglet
						C[7][0][0][0][i]=xub_3(aux);
						C[1][0][0][0][i]=xuv_3(aux)+C[7][0][0][0][i];
						C[8][0][0][0][i]=xdb_3(aux);
						C[2][0][0][0][i]=xdv_3(aux)+C[8][0][0][0][i];
						C[3][0][0][0][i]=C[9][0][0][0][i]=xs_3(aux);

						break;
				}

				break;

			case 1: //MRST parametrization

				switch (order)
				{
					case 0:
						A[7][0][i]=0.2*xS_1(order,aux)-0.5*xD_1(order,aux);
						A[1][0][i]=xuv_1(order,aux)+A[7][0][i];
						A[8][0][i]=0.2*xS_1(order,aux)+0.5*xD_1(order,aux);
						A[2][0][i]=xdv_1(order,aux)+A[8][0][i];
						A[3][0][i]=A[9][0][i]=0.1*xS_1(order,aux);

						S[0][0][0][i]=xg_1(order,aux);
						S[0][1][0][i]=A[1][0][i]+A[7][0][i]+A[2][0][i]+A[8][0][i]+2.*A[3][0][i];

						break;

					case 1:
						B[7][0][0][i]=0.2*xS_1(order,aux)-0.5*xD_1(order,aux);
						B[1][0][0][i]=xuv_1(order,aux)+B[7][0][0][i];
						B[8][0][0][i]=0.2*xS_1(order,aux)+0.5*xD_1(order,aux);
						B[2][0][0][i]=xdv_1(order,aux)+B[8][0][0][i];
						B[3][0][0][i]=B[9][0][0][i]=0.1*xS_1(order,aux);

						S[0][0][0][i]=xg_1(order,aux);
						S[0][1][0][i]=B[1][0][0][i]+B[7][0][0][i]+B[2][0][0][i]+B[8][0][0][i]+2.*B[3][0][0][i];

						break;

					case 2:
						C[7][0][0][0][i]=0.2*xS_1(order,aux)-0.5*xD_1(order,aux);
						C[1][0][0][0][i]=xuv_1(order,aux)+C[7][0][0][0][i];
						C[8][0][0][0][i]=0.2*xS_1(order,aux)+0.5*xD_1(order,aux);
						C[2][0][0][0][i]=xdv_1(order,aux)+C[8][0][0][0][i];
						C[3][0][0][0][i]=C[9][0][0][0][i]=0.1*xS_1(order,aux);

						S[0][0][0][i]=xg_1(order,aux);
						S[0][1][0][i]=C[1][0][0][0][i]+C[7][0][0][0][i]+C[2][0][0][0][i]+C[8][0][0][0][i]+2.*C[3][0][0][0][i];

						break;
				}

				break;

			case 2: //MRST grid

				mode=-4*order+9;

				struc_(&aux,&Q[nfi],&mode,&xuv,&xdv,&xub,&xdb,&xs,&xc,&xb,&xg);

				//Singlet
				S[0][0][0][i]=xg;
				S[0][1][0][i]=xuv+xdv+2.*(xub+xdb+xs);

				switch (order)
				{
					case 0:
						//LO nonsinglet
						A[7][0][i]=xub;
						A[1][0][i]=xuv+xub;
						A[8][0][i]=xdb;
						A[2][0][i]=xdv+xdb;
						A[3][0][i]=A[9][0][i]=xs;

						break;

					case 1:
						//NLO nonsinglet
						B[7][0][0][i]=xub;
						B[1][0][0][i]=xuv+xub;
						B[8][0][0][i]=xdb;
						B[2][0][0][i]=xdv+xdb;
						B[3][0][0][i]=B[9][0][0][i]=xs;
						break;

					case 2:
						//NNLO nonsinglet
						C[7][0][0][0][i]=xub;
						C[1][0][0][0][i]=xuv+xub;
						C[8][0][0][0][i]=xdb;
						C[2][0][0][0][i]=xdv+xdb;
						C[3][0][0][0][i]=C[9][0][0][0][i]=xs;

						break;
				}

				break;

			case 3: //Alekhin parametrization

				q2=Q[nfi]*Q[nfi];

				a02m_(&aux,&q2,pdfs,*dpdfs,&npdf,&npar,&orderplus1,&one);

				//alpha_s
				alpha1=pdfs[0]; //initial alpha_s at mu_r=Q[nfi]

				switch (order)
				{
					case 0:
						A[7][0][i]=xus_2(order,aux);
						A[1][0][i]=xuv_2(order,aux)+A[7][0][i];
						A[8][0][i]=xds_2(order,aux);
						A[2][0][i]=xdv_2(order,aux)+A[8][0][i];
						A[3][0][i]=A[9][0][i]=xss_2(order,aux);

						S[0][0][0][i]=xg_2(order,aux);
						S[0][1][0][i]=A[1][0][i]+A[7][0][i]+A[2][0][i]+A[8][0][i]+2.*A[3][0][i];

						break;

					case 1:
						B[7][0][0][i]=xus_2(order,aux);
						B[1][0][0][i]=xuv_2(order,aux)+B[7][0][0][i];
						B[8][0][0][i]=xds_2(order,aux);
						B[2][0][0][i]=xdv_2(order,aux)+B[8][0][0][i];
						B[3][0][0][i]=B[9][0][0][i]=xss_2(order,aux);

						S[0][0][0][i]=xg_2(order,aux);
						S[0][1][0][i]=B[1][0][0][i]+B[7][0][0][i]+B[2][0][0][i]+B[8][0][0][i]+2.*B[3][0][0][i];

						break;

					case 2:
						C[7][0][0][0][i]=xus_2(order,aux);
						C[1][0][0][0][i]=xuv_2(order,aux)+C[7][0][0][0][i];
						C[8][0][0][0][i]=xds_2(order,aux);
						C[2][0][0][0][i]=xdv_2(order,aux)+C[8][0][0][0][i];
						C[3][0][0][0][i]=C[9][0][0][0][i]=xss_2(order,aux);

						S[0][0][0][i]=xg_2(order,aux);
						S[0][1][0][i]=C[1][0][0][0][i]+C[7][0][0][0][i]+C[2][0][0][0][i]+C[8][0][0][0][i]+2.*C[3][0][0][0][i];

						break;
				}

				break;

			case 4: //Alekhin grid

				q2=Q[nfi]*Q[nfi];

				a02m_(&aux,&q2,pdfs,*dpdfs,&npdf,&npar,&orderplus1,&one);

				//alpha_s
				alpha1=pdfs[0]; //initial alpha_s at mu_r=Q[nfi]

				//Singlet
				S[0][0][0][i]=pdfs[3];
				S[0][1][0][i]=pdfs[1]+pdfs[2]+2.*(pdfs[4]+pdfs[6]+pdfs[5]);

				switch (order)
				{
					case 0:
						//LO nonsinglet
						A[7][0][i]=pdfs[4];
						A[1][0][i]=pdfs[1]+pdfs[4];
						A[8][0][i]=pdfs[6];
						A[2][0][i]=pdfs[2]+pdfs[6];
						A[3][0][i]=A[9][0][i]=pdfs[5];

						break;

					case 1:
						//NLO nonsinglet
						B[7][0][0][i]=pdfs[4];
						B[1][0][0][i]=pdfs[1]+pdfs[4];
						B[8][0][0][i]=pdfs[6];
						B[2][0][0][i]=pdfs[2]+pdfs[6];
						B[3][0][0][i]=B[9][0][0][i]=pdfs[5];
						break;

					case 2:
						//NNLO nonsinglet
						C[7][0][0][0][i]=pdfs[4];
						C[1][0][0][0][i]=pdfs[1]+pdfs[4];
						C[8][0][0][0][i]=pdfs[6];
						C[2][0][0][0][i]=pdfs[2]+pdfs[6];
						C[3][0][0][0][i]=C[9][0][0][0][i]=pdfs[5];

						break;
				}

				break;
		}
	}


	// for (j=1; j<=9; j++)
	// {
	// 	for (k=0; k<GRID_PTS; k++)
	// 	{
	// 		fprintf(stderr, "%9.5lf\n", C[j][0][0][0][k]);
	// 	}
	// 	fprintf(stderr, "\n");
	// }
	// exit(0);


	//Computation of alpha_s at the quark thresholds
	mu0=( (input==1 || input==2)? MZ : Q[nfi] );

	for (nf1=nff;mu0<sqrt(kr)*Q[nf1];nf1--);
	if (nf1<nfi) nf1++;

	printf("nf1 = %d     alpha1 = %g\n",nf1,alpha1);

	beta0=Beta0(nf1);
	beta1=Beta1(nf1);
	beta2=Beta2(nf1);

	alpha_post[nf1]=alpha_rk(order,alpha1,mu0,sqrt(kr)*Q[nf1]);
	alpha_pre[nf1]=( HFT? pre(alpha_post[nf1],order) : alpha_post[nf1] );

	alpha_pre[nf1+1]=alpha_rk(order,alpha1,mu0,sqrt(kr)*Q[nf1+1]);
	alpha_post[nf1+1]=(HFT? post(alpha_pre[nf1+1],order) : alpha_pre[nf1+1] );

	for (nf=nf1-1;nf>=nfi;nf--)
	{
		beta0=Beta0(nf);
		beta1=Beta1(nf);
		beta2=Beta2(nf);

		alpha_post[nf]=alpha_rk(order,alpha_pre[nf+1],sqrt(kr)*Q[nf+1],sqrt(kr)*Q[nf]);
		alpha_pre[nf]=( HFT? pre(alpha_post[nf],order) : alpha_post[nf] );
	}

	for (nf=nf1+1;nf<=nff+1;nf++)
	{
		beta0=Beta0(nf-1);
		beta1=Beta1(nf-1);
		beta2=Beta2(nf-1);

		alpha_pre[nf]=alpha_rk(order,alpha_post[nf-1],sqrt(kr)*Q[nf-1],sqrt(kr)*Q[nf]);
		alpha_post[nf]=(HFT? post(alpha_pre[nf],order) : alpha_pre[nf] );
	}

	for (nf=nfi;nf<=nff+1;nf++)
		printf ("%d  %g  %15.9g  %15.9g\n",nf,Q[nf],alpha_pre[nf],alpha_post[nf]);


	/*
	int _nf = nf;
	nf = 4;
	Nf = 4.0;
	beta0 = Beta0(nf);

	SplitFunc LO[] = { P0NS, P0qq, P0qg, P0gq, P0gg };
	SplitFunc NLO[] = { P1NSp, P1NSm, P1qq, P1qg, P1gq, P1gg };
	SplitFunc NNLO[] = { P2NSp, P2NSm, P2NSv, P2qq, P2qg, P2gq, P2gg };
	OutputFileSplitFuncs("out-LO-splitfuncs.dat", LO, 5);
	OutputFileSplitFuncs("out-NLO-splitfuncs.dat", NLO, 6);
	OutputFileSplitFuncs("out-NNLO-splitfuncs.dat", NNLO, 7);
	
	nf = _nf;

	
	OutputFileAlphaSBetas("out-alphas_betas.dat");

	OutputFileInitialDistributions("out-init_dists.dat");

	OutputFileGauleg("out-gauleg.dat");

	double *vals = malloc(GRID_PTS * (sizeof *vals));
	for (uint i=0; i<GRID_PTS; i++)
		vals[i] = 1.0;

	double *out = malloc(GRID_PTS * (sizeof *out));
	for (uint i=0; i<GRID_PTS; i++)
		out[i] = convolution(i, P1NSp, vals);

	OutputFileConvolution("out-conv.dat", out);


	memset(vals, 0, sizeof vals);
	memset(out, 0, sizeof out);
	double x;
	for (uint i=0; i<GRID_PTS; i++)
	{
		x = X[i];
		vals[i] = sin(x);
	}
	for (uint i=0; i<10; i++)
	{
		x = ((double)i) / 67.0;
		out[i] = interp(vals, x);
	}
	OutputFileInterpolation("out-interp.dat", out);

	free(vals);
	free(out);

	*/
	/*
	for (j=1; j<=9; j++)
	{
		for (n=0; n<ITERATIONS; n++)
		{
			for (k=0; k<GRID_PTS; k++)
			{
				fprintf(stderr, "%8.5lf\n", A[j][n][k]);
			}
		}
		fprintf(stderr, "\n");
	}
	exit(0);
	*/


	//Evolution steps
	for (nf=nfi;;nf++)
	{
		printf("nf = %d\n",nf);

		Nf=(double)nf;

		//Setting PDF's values at the start of the current evolution step
		for (i=0;i<=GRID_PTS-1;i++)
		{
			switch (order)
			{
				case 0:
				{
					for (j=13;j<=18;j++) 
						A[j][0][i]=A[j-12][0][i]-A[j-6][0][i];

					A[25][0][i]=0.;
					for (j=13;j<=18;j++) 
						A[25][0][i]+=A[j][0][i];

					for (j=26;j<=30;j++) 
						A[j][0][i]=A[13][0][i]-A[j-12][0][i];

					for (j=19;j<=24;j++) 
						A[j][0][i]=A[j-18][0][i]+A[j-12][0][i];

					S[0][1][0][i]=0.;
					for (j=19;j<=24;j++) 
						S[0][1][0][i]+=A[j][0][i];

					for (j=32;j<=36;j++) 
						A[j][0][i]=A[19][0][i]-A[j-12][0][i];

				} break;

				case 1:
				{
					for (j=13;j<=18;j++) B[j][0][0][i]=B[j-12][0][0][i]-B[j-6][0][0][i];

					B[25][0][0][i]=0.;
					for (j=13;j<=18;j++) B[25][0][0][i]+=B[j][0][0][i];

					for (j=26;j<=30;j++) B[j][0][0][i]=B[13][0][0][i]-B[j-12][0][0][i];

					for (j=19;j<=24;j++) B[j][0][0][i]=B[j-18][0][0][i]+B[j-12][0][0][i];

					S[0][1][0][i]=0.;
					for (j=19;j<=24;j++) S[0][1][0][i]+=B[j][0][0][i];

					for (j=32;j<=36;j++) B[j][0][0][i]=B[19][0][0][i]-B[j-12][0][0][i];
				} break;
				case 2:
				{
					for (j=13;j<=18;j++) 
						C[j][0][0][0][i] = C[j-12][0][0][0][i]-C[j-6][0][0][0][i];

					C[25][0][0][0][i]=0.;
					for (j=13;j<=18;j++) 
						C[25][0][0][0][i] += C[j][0][0][0][i];

					for (j=26;j<=30;j++)
						C[j][0][0][0][i] = C[13][0][0][0][i] - C[j-12][0][0][0][i];

					for (j=19;j<=24;j++) 
						C[j][0][0][0][i] = C[j-18][0][0][0][i] + C[j-12][0][0][0][i];

					S[0][1][0][i]=0.;
					for (j=19;j<=24;j++)
						S[0][1][0][i] += C[j][0][0][0][i];

					for (j=32;j<=36;j++)
						C[j][0][0][0][i] = C[19][0][0][0][i] - C[j-12][0][0][0][i];

				} break;
			}
		}

		
		// for (k=0; k<GRID_PTS-1; k++)
		// {
		// 	fprintf(stderr, "%15.9lf\n", C[26][0][0][0][k]);
		// }
		// fprintf(stderr, "\n");
		// exit(0);
		

		if (Q[nf+1]==0.) break;

		beta0=Beta0(nf);
		beta1=Beta1(nf);
		beta2=Beta2(nf);
		fprintf(stderr, "betas: %.9lf %.9lf %.9lf\n", beta0, beta1, beta2);

		alpha0=alpha_post[nf];
		alpha1=alpha_pre[nf+1];
		printf("alpha: %.9g ---> %.9g\n",alpha0,alpha1);

		beta=Beta(order,alpha1);

		log_mf2_mr2=log(1./kr);

		//Solving recursively DGLAP equations for the current evolution step

		if (alpha0!=alpha1)
		{
			//Singlet
			for (i=0;i<=ITERATIONS-2;i++)
			{
				printf("Iteration %d (singlet)\n",i);

				//Evolution of singlet 0-coefficients
				for (k=0;k<=GRID_PTS-2;k++)
				{
					S[0][1][i+1][k]=RecRel_A(S[0][1][i],k,P0qq, 1)+RecRel_A(S[0][0][i],k,P0qg, 1);

					S[0][0][i+1][k]=RecRel_A(S[0][1][i],k,P0gq, 1)+RecRel_A(S[0][0][i],k,P0gg, 1);
				}

				if (order>0)
				{
					//Offset of singlet 1-coefficients
					for (k=0;k<=GRID_PTS-2;k++)
						for (j=0;j<=1;j++)
							S[1][j][i+1][k]=-S[0][j][i+1][k]*beta1/(4.*Pi*beta0)-S[1][j][i][k];

					//Evolution of singlet 1-coefficients
					for (k=0;k<=GRID_PTS-2;k++)
					{
						S[1][1][i+1][k]+=RecRel_B(S[0][1][i],S[1][1][i],k,P0qq,P1qq)+RecRel_B(S[0][0][i],S[1][0][i],k,P0qg,P1qg);

						S[1][0][i+1][k]+=RecRel_B(S[0][1][i],S[1][1][i],k,P0gq,P1gq)+RecRel_B(S[0][0][i],S[1][0][i],k,P0gg,P1gg);
					}
				}

				for (n=2;n<=trunc;n++)
				{
					//Offset of singlet n-coefficients
					for (k=0;k<=GRID_PTS-2;k++)
						for (j=0;j<=1;j++)
						{
							S[n][j][i+1][k]=-S[n-1][j][i+1][k]*beta1/(4.*Pi*beta0)
												-(double)n*S[n][j][i][k]-(double)(n-1)*S[n-1][j][i][k]*beta1/(4.*Pi*beta0);

							if (order==2)
								S[n][j][i+1][k]-=S[n-2][j][i+1][k]*beta2/(16.*Pi2*beta0)
													+(double)(n-2)*S[n-2][j][i][k]*beta2/(16.*Pi2*beta0);
						}

					if (order==1)
					{
						//Evolution of singlet NLO n-coefficients
						for (k=0;k<=GRID_PTS-2;k++)
						{
							S[n][1][i+1][k]+=RecRel_B(S[n-1][1][i],S[n][1][i],k,P0qq,P1qq)+RecRel_B(S[n-1][0][i],S[n][0][i],k,P0qg,P1qg);

							S[n][0][i+1][k]+=RecRel_B(S[n-1][1][i],S[n][1][i],k,P0gq,P1gq)+RecRel_B(S[n-1][0][i],S[n][0][i],k,P0gg,P1gg);
						}
					}

					if (order==2)
					{
						//Evolution of singlet NNLO n-coefficients
						for (k=0;k<=GRID_PTS-2;k++)
						{
							S[n][1][i+1][k]+=RecRel_C(S[n-2][1][i],S[n-1][1][i],S[n][1][i],k,P0qq,P1qq,P2qq)+
												RecRel_C(S[n-2][0][i],S[n-1][0][i],S[n][0][i],k,P0qg,P1qg,P2qg);

							S[n][0][i+1][k]+=RecRel_C(S[n-2][1][i],S[n-1][1][i],S[n][1][i],k,P0gq,P1gq,P2gq)+
												RecRel_C(S[n-2][0][i],S[n-1][0][i],S[n][0][i],k,P0gg,P1gg,P2gg);
						}
					}
				}
			}

			/*
			for (n=0; n<ITERATIONS; n++)
			{
				for (k=0; k<GRID_PTS-1; k++)
				{
					fprintf(stderr, "%8.4lf\n", A[13][n][k]);
				}
				fprintf(stderr, "\n");
			}
			*/

			//Nonsinglet
			switch (order)
			{
				case 0: //LO nonsinglet
					//print_lo_shit = 1;
					for (i=0;i<=ITERATIONS-2;i++)
					{
						printf("Iteration %d (LO nonsinglet)\n",i);

						for (k=0;k<=GRID_PTS-2;k++)
						{
							for (j=13;j<=12+nf;j++) {
								A[j][i+1][k]=RecRel_A(A[j][i],k,P0NS, 0);
							}

							for (j=32;j<=30+nf;j++)
								A[j][i+1][k]=RecRel_A(A[j][i],k,P0NS, 0);
						}
					}

					// OutputLOCoefficients(A);
				    
					break;

				case 1: //NLO nonsinglet
					for (i=1;i<=ITERATIONS-1;i++)
					{
						printf("Iteration %d (NLO nonsinglet)\n",i);

						for (k=0;k<=GRID_PTS-2;k++)
						{
							for (j=13;j<=12+nf;j++)
							{
								for (l=1;l<=i;l++)
									B[j][i][l][k]=RecRel_Diag(B[j][i-1][l-1],k,P0NS, 0);

								B[j][i][0][k]=-B[j][i][1][k]+RecRel_Vert1(B[j][i-1][0],k,P0NS,P1NSm);
							}

							for (j=32;j<=30+nf;j++)
							{
								for (l=1;l<=i;l++)
									B[j][i][l][k]=RecRel_Diag(B[j][i-1][l-1],k,P0NS, 0);

								B[j][i][0][k]=-B[j][i][1][k]+RecRel_Vert1(B[j][i-1][0],k,P0NS,P1NSp);
							}
						}
					}

					break;

				case 2: //NNLO nonsinglet
					for (s=1;s<=ITERATIONS-1;s++)
					{
						printf("Iteration %d (NNLO nonsinglet)\n",s);

						for (k=0;k<=GRID_PTS-2;k++)
						{
							for (j=26;j<=24+nf;j++)
							{
								for (t=1;t<=s;t++)
									for (n=1;n<=t;n++)
										C[j][s][t][n][k]=RecRel_Diag(C[j][s-1][t-1][n-1],k,P0NS, 1);

								C[j][s][s][0][k]=-0.5*C[j][s][s][1][k]+RecRel_Vert2(C[j][s-1][s-1][0],k,P0NS,P1NSm,P2NSm);

								for (t=s-1;t>=0;t--)
									C[j][s][t][0][k]=-2.*beta1*(C[j][s][t+1][0][k]+C[j][s][t+1][1][k])+RecRel_Horiz(C[j][s-1][t][0],k,P0NS,P1NSm);
							}

							for (j=32;j<=30+nf;j++)
							{
								for (t=1;t<=s;t++)
									for (n=1;n<=t;n++)
										C[j][s][t][n][k]=RecRel_Diag(C[j][s-1][t-1][n-1],k,P0NS, 1);

								C[j][s][s][0][k]=-0.5*C[j][s][s][1][k]+RecRel_Vert2(C[j][s-1][s-1][0],k,P0NS,P1NSp,P2NSp);

								for (t=s-1;t>=0;t--)
									C[j][s][t][0][k]=-2.*beta1*(C[j][s][t+1][0][k]+C[j][s][t+1][1][k])+RecRel_Horiz(C[j][s-1][t][0],k,P0NS,P1NSp);
							}

							for (t=1;t<=s;t++)
								for (n=1;n<=t;n++)
									C[25][s][t][n][k]=RecRel_Diag(C[25][s-1][t-1][n-1],k,P0NS, 1);

							C[25][s][s][0][k]=-0.5*C[25][s][s][1][k]+RecRel_Vert2(C[25][s-1][s-1][0],k,P0NS,P1NSm,P2NSv);

							for (t=s-1;t>=0;t--)
								C[25][s][t][0][k]=-2.*beta1*(C[25][s][t+1][0][k]+C[25][s][t+1][1][k])+RecRel_Horiz(C[25][s-1][t][0],k,P0NS,P1NSm);
						}
					}

					break;
			}

			//Summing the recursive series (for tabulated Q values)

			for (;Qtab[qct]<=Q[nf+1] && qct<sizeof(Qtab)/sizeof(double);qct++)
			{
				alpha1=alpha_rk(order,alpha0,Q[nf],Qtab[qct]);

				L1=log(alpha1/alpha0);

				if (order==1)
					L2=log((alpha1*beta1+4.*Pi*beta0)/(alpha0*beta1+4.*Pi*beta0));
				else
					L2=log((16.*Pi2*beta0+4.*Pi*alpha1*beta1+alpha1*alpha1*beta2)/(16.*Pi2*beta0+4.*Pi*alpha0*beta1+alpha0*alpha0*beta2));

				aux=4.*beta0*beta2-beta1*beta1;
				if (aux>=0)
					L3=atan(2.*Pi*(alpha1-alpha0)*sqrt(aux)/(2.*Pi*(8.*Pi*beta0+(alpha1+alpha0))+alpha1*alpha0*beta2))/sqrt(aux);
				else
					L3=atanh(2.*Pi*(alpha1-alpha0)*sqrt(-aux)/(2.*Pi*(8.*Pi*beta0+(alpha1+alpha0))+alpha1*alpha0*beta2))/sqrt(-aux);

				fprintf(stderr, "alpha0=%.9lf, alpha1=%.9lf\n", alpha0, alpha1);
				fprintf(stderr, "L1=%.9lf", L1);

				//Singlet
				for (j=0;j<=1;j++)
					for (k=0;k<=GRID_PTS-2;k++)
					{
						R[j*31][k]=S[0][j][0][k];

						for (i=1;i<=ITERATIONS-1;i++)
							for (n=0;n<=trunc;n++)
								R[j*31][k]+=S[n][j][i][k]*pow(alpha1,(double)n)*pow(L1,(double)i)/fact(i);
					}

				//Nonsinglet
				switch (order)
				{
					case 0:
						for (k=0;k<=GRID_PTS-2;k++)
						{
							for (j=13;j<=30+nf;j++)
							{
								R[j][k]=A[j][0][k];

								if (j==12+nf)
									j=31;
							}

							for (i=1;i<=ITERATIONS-1;i++)
							{
								for (j=13;j<=12+nf;j++)
									R[j][k]+=A[j][i][k]*pow(L1,(double)i)/fact(i);

								for (j=32;j<=30+nf;j++)
									R[j][k]+=A[j][i][k]*pow(L1,(double)i)/fact(i);
							}
						}

						break;

					case 1:
						for (k=0;k<=GRID_PTS-2;k++)
						{
							for (j=13;j<=30+nf;j++)
							{
								R[j][k]=B[j][0][0][k];

								if (j==12+nf)
									j=31;
							}

							for (i=1;i<=ITERATIONS-1;i++)
								for (l=0;l<=i;l++)
								{
									for (j=13;j<=12+nf;j++)
										R[j][k]+=B[j][i][l][k]/fact(l)/fact(i-l)*pow(L1,(double)l)*pow(L2,(double)(i-l));

									for (j=32;j<=30+nf;j++)
										R[j][k]+=B[j][i][l][k]/fact(l)/fact(i-l)*pow(L1,(double)l)*pow(L2,(double)(i-l));
								}
						}

						break;

					case 2:
						for (k=0;k<=GRID_PTS-2;k++)
						{
							for (j=25;j<=30+nf;j++)
							{
								R[j][k]=C[j][0][0][0][k];

								if (j==24+nf)
									j=31;
							}

							for (s=1;s<=ITERATIONS-1;s++)
								for (t=0;t<=s;t++)
									for (n=0;n<=t;n++)
									{
										for (j=25;j<=24+nf;j++)
											R[j][k]+=C[j][s][t][n][k]/fact(n)/fact(t-n)/fact(s-t)*pow(L1,(double)n)*pow(L2,(double)(t-n))*pow(L3,(double)(s-t));

										for (j=32;j<=30+nf;j++)
											R[j][k]+=C[j][s][t][n][k]/fact(n)/fact(t-n)/fact(s-t)*pow(L1,(double)n)*pow(L2,(double)(t-n))*pow(L3,(double)(s-t));
									}
						}

						break;
				}

				//Flavor reconstruction (for tabulated Q values)
				for (k=0;k<=GRID_PTS-2;k++)
				{
					if (order==2)
					{
						R[13][k]=R[25][k];
						for (j=26;j<=24+nf;j++) R[13][k]+=R[j][k];
						R[13][k]/=Nf;

						for (j=14;j<=12+nf;j++) R[j][k]=R[13][k]-R[j+12][k];
					}

					R[19][k]=R[31][k];
					for (j=32;j<=30+nf;j++) 
						R[19][k]+=R[j][k];
					R[19][k]/=Nf;

					for (j=20;j<=18+nf;j++) 
						R[j][k]=R[19][k]-R[j+12][k];

					for (j=1;j<=nf;j++)
					{
						R[j][k]  =0.5*(R[j+18][k]+R[j+12][k]);
						R[j+6][k]=0.5*(R[j+18][k]-R[j+12][k]);
					}

					if (order<2)
					{
						R[25][k]=0.;
						for (j=13;j<=12+nf;j++) 
							R[25][k]+=R[j][k];

						for (j=26;j<=24+nf;j++)
							R[j][k]=R[13][k]-R[j-12][k];
					}
				}

				//Writing the output
				sprintf(filename,"Q%d.%s",qct+1,extension);
				fp=fopen(filename,"w");

				for (k=0;k<=GRID_PTS-1;k++)
				{
					fprintf(fp,"%15.9g",X[k]);

					//Modify here if you want other PDFs to be displayed
					//
					//PDF indices
					//
					//0      gluons         g
					//1-6    quarks         u,d,s,c,b,t
					//7-12   antiquarks     au,ad,as,ac,ab,at
					//13-18  q_i^-          um,dm,sm,cm,bm,tm
					//19-24  q_i^+          up,dp,sp,cp,bp,tp
					//25     q^(-)
					//26-30  q_{NS,1i}^(-)  dd,sd,cd,bd,td
					//31     q^(+)
					//32-36  q_{NS,1i}^(+)  ds,ss,cs,bs,ts
					for (j=0;j<=12;j++)
						fprintf(fp,"   %15.9g",R[j][k]);

					fprintf(fp,"\n");
				}

				fclose(fp);
			}

			//Summing the recursive series (for threshold Q values)

			alpha1=alpha_pre[nf+1];

			L1=log(alpha1/alpha0);

			if (order==1)
				L2=log((alpha1*beta1+4.*Pi*beta0)/(alpha0*beta1+4.*Pi*beta0));
			else
				L2=log((16.*Pi2*beta0+4.*Pi*alpha1*beta1+alpha1*alpha1*beta2)/(16.*Pi2*beta0+4.*Pi*alpha0*beta1+alpha0*alpha0*beta2));

			aux=4.*beta0*beta2-beta1*beta1;
			if (aux>=0)
				L3=atan(2.*Pi*(alpha1-alpha0)*sqrt(aux)/(2.*Pi*(8.*Pi*beta0+(alpha1+alpha0))+alpha1*alpha0*beta2))/sqrt(aux);
			else
				L3=atanh(2.*Pi*(alpha1-alpha0)*sqrt(-aux)/(2.*Pi*(8.*Pi*beta0+(alpha1+alpha0))+alpha1*alpha0*beta2))/sqrt(-aux);

			//Singlet
			for (j=0;j<=1;j++)
				for (k=0;k<=GRID_PTS-2;k++)
					for (i=1;i<=ITERATIONS-1;i++)
						for (n=0;n<=trunc;n++)
							S[0][j][0][k]+=S[n][j][i][k]*pow(alpha1,(double)n)*pow(L1,(double)i)/fact(i);

			//Nonsinglet
			switch (order)
			{
				case 0:
					for (k=0;k<=GRID_PTS-2;k++)
						for (i=1;i<=ITERATIONS-1;i++)
						{
							for (j=13;j<=12+nf;j++)
								A[j][0][k]+=A[j][i][k]*pow(L1,(double)i)/fact(i);

							for (j=32;j<=30+nf;j++)
								A[j][0][k]+=A[j][i][k]*pow(L1,(double)i)/fact(i);
						}

					break;

				case 1:
					for (k=0;k<=GRID_PTS-2;k++)
						for (i=1;i<=ITERATIONS-1;i++)
							for (l=0;l<=i;l++)
							{
								for (j=13;j<=12+nf;j++)
									B[j][0][0][k]+=B[j][i][l][k]/fact(l)/fact(i-l)*pow(L1,(double)l)*pow(L2,(double)(i-l));

								for (j=32;j<=30+nf;j++)
									B[j][0][0][k]+=B[j][i][l][k]/fact(l)/fact(i-l)*pow(L1,(double)l)*pow(L2,(double)(i-l));
							}

					break;

				case 2:
					for (k=0;k<=GRID_PTS-2;k++)
						for (s=1;s<=ITERATIONS-1;s++)
							for (t=0;t<=s;t++)
								for (n=0;n<=t;n++)
								{
									for (j=25;j<=24+nf;j++)
										C[j][0][0][0][k]+=C[j][s][t][n][k]/fact(n)/fact(t-n)/fact(s-t)*pow(L1,(double)n)*pow(L2,(double)(t-n))*pow(L3,(double)(s-t));

									for (j=32;j<=30+nf;j++)
										C[j][0][0][0][k]+=C[j][s][t][n][k]/fact(n)/fact(t-n)/fact(s-t)*pow(L1,(double)n)*pow(L2,(double)(t-n))*pow(L3,(double)(s-t));
								}

					break;
			}

			//Flavor reconstruction (for threshold Q values)
			switch (order)
			{
				case 0:
					for (k=0;k<=GRID_PTS-2;k++)
					{
						A[19][0][k]=S[0][1][0][k];
						for (j=32;j<=30+nf;j++) A[19][0][k]+=A[j][0][k];
						A[19][0][k]/=Nf;

						for (j=20;j<=18+nf;j++) A[j][0][k]=A[19][0][k]-A[j+12][0][k];

						for (j=1;j<=nf;j++)
						{
							A[j][0][k]=0.5*(A[j+18][0][k]+A[j+12][0][k]);
							A[j+6][0][k]=0.5*(A[j+18][0][k]-A[j+12][0][k]);
						}
					}

					break;

				case 1:
					for (k=0;k<=GRID_PTS-2;k++)
					{
						B[19][0][0][k]=S[0][1][0][k];
						for (j=32;j<=30+nf;j++)
							B[19][0][0][k]+=B[j][0][0][k];
						B[19][0][0][k]/=Nf;

						for (j=20;j<=18+nf;j++)
							B[j][0][0][k]=B[19][0][0][k]-B[j+12][0][0][k];

						for (j=1;j<=nf;j++)
						{
							B[j][0][0][k]  =0.5*(B[j+18][0][0][k]+B[j+12][0][0][k]);
							B[j+6][0][0][k]=0.5*(B[j+18][0][0][k]-B[j+12][0][0][k]);
						}
					}

					break;

				case 2:
					for (k=0;k<=GRID_PTS-2;k++)
					{
						C[13][0][0][0][k]=C[25][0][0][0][k];
						for (j=26;j<=24+nf;j++) C[13][0][0][0][k]+=C[j][0][0][0][k];
						C[13][0][0][0][k]/=Nf;

						for (j=14;j<=12+nf;j++) C[j][0][0][0][k]=C[13][0][0][0][k]-C[j+12][0][0][0][k];

						C[19][0][0][0][k]=S[0][1][0][k];
						for (j=32;j<=30+nf;j++) C[19][0][0][0][k]+=C[j][0][0][0][k];
						C[19][0][0][0][k]/=Nf;

						for (j=20;j<=18+nf;j++) C[j][0][0][0][k]=C[19][0][0][0][k]-C[j+12][0][0][0][k];

						for (j=1;j<=nf;j++)
						{
							C[j][0][0][0][k]=0.5*(C[j+18][0][0][0][k]+C[j+12][0][0][0][k]);
							C[j+6][0][0][0][k]=0.5*(C[j+18][0][0][0][k]-C[j+12][0][0][0][k]);
						}
					}

					break;
			}
		}

		//Heavy flavors treatment
		if (HFT && vfns && order==2 && Q[nf+2]!=0.)
		{
			printf("%dth quark threshold (mass %g)\n",nf+1,Q[nf+1]);

			//Copy of pre-threshold distributions
			for (k=0;k<=GRID_PTS-2;k++)
			{
				for (j=0;j<=1;j++)
					S[0][j][1][k]=S[0][j][0][k];


				for (i=1;i<=nf;i++)
					for (j=i;j<=i+6;j+=6)
						C[j][1][0][0][k]=C[j][0][0][0][k];
			}

			aux=alpha_post[nf+1];
			fprintf(stderr, "value of alpha_s post threshold: %9.5lf\n", aux);

			//Computation of after-threshold distributions
			for (k=0;k<=GRID_PTS-2;k++)
			{
				for (i=1;i<=nf;i++)
					for (j=i;j<=i+6;j+=6)
						C[j][0][0][0][k] += pow(aux / (4.0*Pi), 2.0) * convolution(k, A2ns, C[j][1][0][0]);

				S[0][0][0][k] += pow(aux / (4.0*Pi), 2.0) * (convolution(k, A2gq, S[0][1][1])
							   + convolution(k, A2gg, S[0][0][1]));

				C[nf+1][0][0][0][k] = C[nf+7][0][0][0][k] = 0.5*pow(aux / (4.0*Pi), 2.0) * (convolution(k, A2hq, S[0][1][1])
                                                          + convolution(k, A2hg, S[0][0][1]));
			}
		}
	}


	

	for (uint t=0; t<trunc+1; t++)
	{
		for (uint i=0; i<2; i++)
		{
			for (uint n=0; n<ITERATIONS; n++)
				free(S[t][i][n]);

			free(S[t][i]);
		}
		free(S[t]);
	}
	free(S);


	for (uint j=0; j<37; j++)
	{
		for (uint n=0; n<ITERATIONS; n++)
		{
			for (uint m=0; m<ITERATIONS; m++)
			{
				for (uint l=0; l<ITERATIONS; l++)
					free(C[j][n][m][l]);

				free(B[j][n][m]);
				free(C[j][n][m]);
			}

			free(A[j][n]);
			free(B[j][n]);
			free(C[j][n]);
		}
		free(A[j]);
		free(B[j]);
		free(C[j]);
		free(R[j]);
	}
    
}






//Safe allocation functions

void *Malloc(size_t size)
{
	void *p=malloc(size);

	if (p==NULL)
	{
		printf("Allocation failed\n");
		exit(1);
	}

	return p;
}


void *Calloc(size_t nmemb,size_t size)
{
	void *p=calloc(nmemb,size);

	if (p==NULL)
	{
		printf("Allocation failed\n");
		exit(1);
	}

	return p;
}


//Functions needed to perform the numerical integrations

void gauleg(double x1,double x2,double *x,double *w,int n)
{
	//Computes the Gauss-Legendre abscissas and weights of order n in the interval [x1,x2]

	int m,j,i;
	double eps=1.e-14,z1,z,xm,xl,pp,p3,p2,p1;
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);

	for(i=1;i<=m;i++)
	{
		z=cos(Pi*((double)i-0.25)/((double)n+0.5));
		do
		{
			p1=1.;
			p2=0.;

			for(j=1;j<=n;j++)
			{
				p3=p2;
				p2=p1;
				p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j;
			}

			pp=n*(z*p1-p2)/(z*z-1.);
			z1=z;
			z=z1-p1/pp;
		}
		while(fabs(z-z1)>eps);

		x[i-1]=xm-xl*z;
		x[n-i]=xm+xl*z;
		w[i-1]=2.*xl/((1.-z*z)*pp*pp);
		w[n-i]=w[i-1];
	}
}


//Polynomial interpolation

double interp(double *A,double x)
{
	//Returns the polynomial interpolation of grade INTERP_PTS*2-1 of the array A at the point x

	double res,aux;
	int k;

	for (k=0;x>=X[k];k++);

	k-=INTERP_PTS;

	if (k<0) k=0;
	if (k>GRID_PTS-2*INTERP_PTS) k=GRID_PTS-2*INTERP_PTS;

	return polint(&X[k],&A[k],2*INTERP_PTS,x);
}


double polint(double *xa,double *ya,int n,double x)
{
	//Polynomial interpolation of grade n at point x using the pairs (xa[i],ya[i])

	int i,m,ns=0;
	double den,dif,dift,ho,hp,w,res;
	double *c = calloc(n, sizeof *c);
	double *d = calloc(n, sizeof *d);

	dif=fabs(x-xa[0]);
	for (i=0;i<=n-1;i++)
	{
		if ((dift=fabs(x-xa[i]))<dif)
		{
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}

	res=ya[ns--];

	for (m=1;m<n;m++)
	{
		for (i=0;i<=n-m-1;i++)
		{
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];

			if ((den=ho-hp)==0.0)
			{
				printf("Error in routine polint");
				exit(1);
			}

			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}

		res+=((2*ns<(n-1-m) ? c[ns+1] : d[ns--]));
	}

	free(c);
	free(d);

	return res;
}


double convolution(int i,double kernel(int,double),double *A)
{
	//Returns the convolution product between the kernel and A in X[i], multiplied by X[i]

	int j;
	double aux,x=X[i],y,w,a,b,res=(kernel(2,1.)*log(1.-x)+kernel(3,1.))*A[i];

	//Mapping
	for (j=0;j<=NGP-1;j++)
	{
		y=XG[j];
		w=WG[j];
		a=pow(x,1.-y);
		b=pow(x,y);
		res-=log(x)*a*kernel(1,a)*interp(A,b)*w;
		res-=log(x)*b*(kernel(2,b)*interp(A,a)-kernel(2,1.)*A[i])/(1.-b)*w;
	}

	return res;
}


//Recursion relations

double RecRel_A(double *A,int k,double P0(int,double), int singlet)
{
	if (singlet)
	{
		// for (uint _k=0; _k<GRID_PTS; _k++)
		// 	fprintf(stderr, "%15.9lf\n", A[_k]);
		
	}

	double conv = convolution(k,P0,A);
	return -2./beta0*conv;
}


double RecRel_B(double *A,double *B,int k,double P0(int,double),double P1(int,double))
{
	// for (uint _k=0; _k<GRID_PTS; _k++)
	// 	fprintf(stderr, "%15.9lf, %15.9lf\n", A[_k], B[_k]);

	double conv1 = convolution(k,P0,B);
	double conv2 = convolution(k,P1,A);
	double res=conv1*2./beta0;
	res+=conv2/(Pi*beta0);

	//Additional terms for mu_f!=mu_r
	if (log_mf2_mr2!=0.)
		res-=log_mf2_mr2*convolution(k,P0,A)/(2.*Pi);

	return -res;
}


double RecRel_C(double *A,double *B,double *C,int k,
                double P0(int,double),double P1(int,double),double P2(int,double))
{
	// for (uint _k=0; _k<GRID_PTS; _k++)
	// 	fprintf(stderr, "%15.9lf, %15.9lf, %15.9lf\n", A[_k], B[_k], C[_k]);

	double conv1 = convolution(k,P0,C);
	double conv2 = convolution(k,P1,B);
	double conv3 = convolution(k,P2,A);
	double res=conv1*2./beta0;
	res+=conv2/(Pi*beta0);
	res+=conv3/(2.*Pi2*beta0);

	//Additional terms for mu_f!=mu_r
	if (log_mf2_mr2!=0.)
	{
		res-=log_mf2_mr2*convolution(k,P0,B)/(2.*Pi);

		res+=log_mf2_mr2*convolution(k,P0,A)*(beta0*log_mf2_mr2-beta1/beta0)/(8.*Pi2);
		res-=log_mf2_mr2*convolution(k,P1,A)/(2.*Pi2);
	}

	return -res;
}


double RecRel_Diag(double *C,int k,double P0(int,double), int nnlo)
{
	if (nnlo)
	{
		// for (k=0; k<GRID_PTS; k++)
		// 	fprintf(stderr, "%15.9lf\n", C[k]);

	}

	double conv = convolution(k,P0,C);
	return -2./beta0*conv;
}


double RecRel_Vert1(double *B,int k,double P0(int,double),double P1(int,double))
{
	double res=-4./beta1*convolution(k,P1,B);

	//Additional terms for mu_f!=mu_r
	if (log_mf2_mr2!=0.)
		res+=log_mf2_mr2*2.*beta0/beta1*convolution(k,P0,B);

	return res;
}


double RecRel_Vert2(double *C,int k,double P0(int,double),double P1(int,double),double P2(int,double))
{
	// for (k=0; k<GRID_PTS; k++)
	// 	fprintf(stderr, "%15.9lf\n", C[k]);

	double conv = convolution(k,P2,C);
	double res=-4./beta2*conv;

	//Additional terms for mu_f!=mu_r
	if (log_mf2_mr2!=0.)
	{
		res-=log_mf2_mr2*(beta0*beta0*log_mf2_mr2-beta1)/beta2*convolution(k,P0,C);
		res+=log_mf2_mr2*4.*beta0/beta2*convolution(k,P1,C);
	}

	return res;
}


double RecRel_Horiz(double *C,int k,double P0(int,double),double P1(int,double))
{
	// for (k=0; k<GRID_PTS; k++)
	// 	fprintf(stderr, "%15.9lf\n", C[k]);

	double conv = convolution(k,P1,C);
	double res=-8.*conv;

	//Additional terms for mu_f!=mu_r
	if (log_mf2_mr2!=0.)
		res+=log_mf2_mr2*4.*beta0*convolution(k,P0,C);

	return res;
}


//Running coupling

double Beta0(int f)
{
	return 11./3.*NC-4./3.*TR*(double)f;
}


double Beta1(int f)
{
	return 34./3.*NC*NC-(4.*Cf+20./3.*NC)*TR*(double)f;
}


double Beta2(int f)
{
	double F=(double)f;

	return 2857./54.*NC*NC*NC+(2.*Cf*Cf-205./9.*Cf*NC-1415./27.*NC*NC)*TR*F+(44./9.*Cf+158./27.*NC)*TR*TR*F*F;
}


double Beta(int order,double alpha)
{
	double res=beta0;
	if (order>0) res+=beta1*alpha/(4.*Pi);
	if (order==2) res+=beta2*alpha*alpha/(16.*Pi2);
	res*=-alpha*alpha/(4.*Pi);

	return res;
}


double alpha_rk(int order,double alpha0,double Qi,double Qf)
{
	//Given alpha0=alpha_s(Qi), computes alpha_s(Qf) by a fourth order Runge-Kutta at NLO and NNLO and by an exact solution at LO

	const int STEPS=200;
	int i;
	double h=2.*log(Qf/Qi)/(double)STEPS;
	double k1,k2,k3,k4;
	double a=alpha0;

	if (Qi==Qf)
		return alpha0;

	if (order==0)
		return (2.*Pi*alpha0)/(2.*Pi+alpha0*beta0*log(Qf/Qi));

	for (i=1;i<=STEPS;i++)
	{
		k1=h*Beta(order,a);
		k2=h*Beta(order,a+k1/2.);
		k3=h*Beta(order,a+k2/2.);
		k4=h*Beta(order,a+k3);
		a+=k1/6.+k2/3.+k3/3.+k4/6.;
	}

	return a;
}


double post(double pre,int order)
{
	//Matching condition for alpha_s at quark mass thresholds: given alpha_s(nf) returns alpha_s(nf+1)

	double L=-log_mf2_mr2,res=pre;

	if (order==0)
		return pre;

	res+=pre*pre*L/(6.*Pi);

	if (order==2)
		res+=pow(pre,3.)*(14.+38.*L+4./3.*L*L)/(48.*Pi2);

	return res;
}


double pre(double post,int order)
{
	//Matching condition for alpha_s at quark mass thresholds: given alpha_s(nf+1) returns alpha_s(nf)

	double L=-log_mf2_mr2,res=post;

	if (order==0)
		return post;

	res-=post*post*L/(6.*Pi);

	if (order==2)
		res+=pow(post,3.)*(L*L/36.-19./24.*L-7./24.)/Pi2;

	return res;
}


//Special functions

double Li2(double x)
{
	int nw=2,n1=-1,n2=1;
	dcomplex HC1[3],HC2[3][3],HC3[3][3][3],HC4[3][3][3][3];
	double HR1[3],HR2[3][3],HR3[3][3][3],HR4[3][3][3][3];
	double HI1[3],HI2[3][3],HI3[3][3][3],HI4[3][3][3][3];

	hplog_(&x,&nw,HC1,*HC2,**HC3,***HC4,HR1,*HR2,**HR3,***HR4,HI1,*HI2,**HI3,***HI4,&n1,&n2);

	return HR2[2][1];
}


double Li3(double x)
{
	int nw=3,n1=-1,n2=1;
	dcomplex HC1[3],HC2[3][3],HC3[3][3][3],HC4[3][3][3][3];
	double HR1[3],HR2[3][3],HR3[3][3][3],HR4[3][3][3][3];
	double HI1[3],HI2[3][3],HI3[3][3][3],HI4[3][3][3][3];

	hplog_(&x,&nw,HC1,*HC2,**HC3,***HC4,HR1,*HR2,**HR3,***HR4,HI1,*HI2,**HI3,***HI4,&n1,&n2);

	return HR3[2][1][1];
}


double S12(double x)
{
	int nw=3,n1=-1,n2=1;
	dcomplex HC1[3],HC2[3][3],HC3[3][3][3],HC4[3][3][3][3];
	double HR1[3],HR2[3][3],HR3[3][3][3],HR4[3][3][3][3];
	double HI1[3],HI2[3][3],HI3[3][3][3],HI4[3][3][3][3];

	hplog_(&x,&nw,HC1,*HC2,**HC3,***HC4,HR1,*HR2,**HR3,***HR4,HI1,*HI2,**HI3,***HI4,&n1,&n2);

	return HR3[2][2][1];
}


//Factorial

double fact(int n)
{
	double res=1.;

	for (;n>1;n--) res*=(double)n;

	return res;
}


//MRST initial distributions (Phys.Lett.B 531 (2002) 216 for LO and NNLO, Eur.Phys.J.C 23 (2002) 73 for NLO; 1 GeV^2)

double xuv_1(int order,double x)
{
	switch (order)
	{
		case 0:
			return 0.474*pow(x,0.3)*pow(1.-x,3.12)*(1.-1.32*sqrt(x)+19.56*x);
		case 1:
			return 0.158*pow(x,0.25)*pow(1.-x,3.33)*(1.+5.61*sqrt(x)+55.49*x);
		case 2:
			return 0.262*pow(x,0.31)*pow(1.-x,3.5)*(1.+3.83*sqrt(x)+37.65*x);
	}
}


double xdv_1(int order,double x)
{
	switch (order)
	{
		case 0:
			return 0.668*pow(x,0.43)*pow(1.-x,4.03)*(1.-0.83*sqrt(x)+7.68*x);
		case 1:
			return 0.04*pow(x,0.27)*pow(1.-x,3.88)*(1.+52.73*sqrt(x)+30.65*x);
		case 2:
			return 0.061*pow(x,0.35)*pow(1.-x,4.03)*(1.+49.05*sqrt(x)+8.65*x);
	}
}


double xS_1(int order,double x)
{
	switch (order)
	{
		case 0:
			return 0.458*pow(x,-0.19)*pow(1.-x,7.51)*(1.+0.025*sqrt(x)+7.63*x);
		case 1:
			return 0.222*pow(x,-0.26)*pow(1.-x,7.1)*(1.+3.42*sqrt(x)+10.3*x);
		case 2:
			return 0.759*pow(x,-0.12)*pow(1.-x,7.66)*(1.-1.34*sqrt(x)+7.4*x);
	}
}


double xg_1(int order,double x)
{
	switch (order)
	{
		case 0:
			return 3.08*pow(x,0.1)*pow(1.-x,6.49)*(1.-2.96*sqrt(x)+9.26*x);
		case 1:
			return 1.9*pow(x,0.09)*pow(1.-x,3.7)*(1.+1.26*sqrt(x)-1.43*x)-0.21*pow(x,-0.33)*pow(1.-x,10.);
		case 2:
			return 0.669*pow(1.-x,3.96)*(1.+6.98*sqrt(x)-3.63*x)-0.23*pow(x,-0.27)*pow(1.-x,8.7);
	}
}


double xD_1(int order,double x)
{
	switch (order)
	{
		case 0:
			return 4.163*pow(x,1.76)*pow(1.-x,9.51)*(1.+7.2*x-24.8*x*x);
		case 1:
			return 1.195*pow(x,1.24)*pow(1.-x,9.1)*(1.+14.05*x-45.52*x*x);
		case 2:
			return 1.432*pow(x,1.24)*pow(1.-x,9.66)*(1.+9.86*x-29.04*x*x);
	}
}


//Alekhin initial distributions (Phys.Rev.D 68 (2003) 014002 ; 9 GeV^2)

double xuv_2(int order,double x)
{
	switch (order)
	{
		case 0:
			return 2.14745*pow(x,0.551)*pow(1.-x,3.672)*(1.+3.*x);
		case 1:
			return 4.03098*pow(x,0.7)*pow(1.-x,3.92)*(1.+1.14*x);
		case 2:
			return 4.42698*pow(x,0.725)*pow(1.-x,4.024)*(1.+1.05*x);
	}
}


double xus_2(int order,double x)
{
	switch (order)
	{
		case 0:
			return 0.16717*pow(x,-0.198)*pow(1.-x,9.2);
		case 1:
			return 0.166439*pow(x,-0.1968)*pow(1.-x,10.16);
		case 2:
			return 0.162806*pow(x,-0.2092)*pow(1.-x,10.49);
	}
}


double xdv_2(int order,double x)
{
	switch (order)
	{
		case 0:
			return 2.06778*pow(x,0.639)*pow(1.-x,4.48);
		case 1:
			return 2.81426*pow(x,0.722)*pow(1.-x,4.94);
		case 2:
			return 3.34753*pow(x,0.772)*pow(1.-x,5.14);
	}
}


double xds_2(int order,double x)
{
	switch (order)
	{
		case 0:
			return 0.117536*pow(x,-0.198)*pow(1.-x,3.8);
		case 1:
			return 0.143328*pow(x,-0.1968)*pow(1.-x,5.1);
		case 2:
			return 0.145752*pow(x,-0.2092)*pow(1.-x,5.6);
	}
}


double xss_2(int order,double x)
{
	switch (order)
	{
		case 0:
			return 0.0490347*pow(x,-0.198)*pow(1.-x,(9.2+3.8)/2.);
		case 1:
			return 0.060482*pow(x,-0.1968)*pow(1.-x,(10.16+5.1)/2.);
		case 2:
			return 0.0616588*pow(x,-0.2092)*pow(1.-x,(10.49+5.6)/2.);
	}
}


double xg_2(int order,double x)
{
	switch (order)
	{
		case 0:
			return 1.60302*pow(x,-0.302)*pow(1.-x,5.3)*(1.-1.94*sqrt(x)+2.8*x);
		case 1:
			return 3.91645*pow(x,-0.146)*pow(1.-x,8.2)*(1.-3.76*sqrt(x)+7.7*x);
		case 2:
			return 4.15733*pow(x,-0.128)*pow(1.-x,9.4)*(1.-3.84*sqrt(x)+8.6*x);
	}
}


//Les Houches benchmark initial distributions (hep-ph/0204316 ; 2 GeV^2)

double xuv_3(double x)
{
	return 5.1072*pow(x,0.8)*pow(1.-x,3.);
}

double xdv_3(double x)
{
	return 3.06432*pow(x,0.8)*pow(1.-x,4.);
}

double xg_3(double x)
{
	return 1.7*pow(x,-0.1)*pow(1.-x,5.);
}

double xdb_3(double x)
{
	return 0.1939875*pow(x,-0.1)*pow(1.-x,6.);
}

double xub_3(double x)
{
	return xdb_3(x)*(1.-x);
}

double xs_3(double x)
{
	return 0.2*(xub_3(x)+xdb_3(x));
}


//Kernels
//P(x) = P1(x) + P2(x)/(1-x)_+ + P3 delta(1-x)

//Structure of the generic kernel:
// double P(int i,double x)
// {
// 	switch (i)
// 	{
// 		case 1:
// 			return P1(x); (regular part)
// 		case 2:
// 			return P2(x); (plus distribution part)
// 		case 3:
// 			return P3; (delta function part)
// 	}
// }

//LO kernels

double P0NS(int i,double x)
{
	switch (i)
	{
		case 1:
			return Cf*(-1.-x);
		case 2:
			return 2.*Cf;
		case 3:
			return 3./2.*Cf;
	}
}


double P0qq(int i,double x)
{
	return P0NS(i,x);
}


double P0qg(int i,double x)
{
	switch (i)
	{
		case 1:
			return 2.*TR*Nf*(2.*x*x-2.*x+1.);
		case 2:
		case 3:
			return 0.;
	}
}


double P0gq(int i,double x)
{
	switch (i)
	{
		case 1:
			return Cf*(x-2.+2./x);
		case 2:
		case 3:
			return 0.;
	}
}


double P0gg(int i,double x)
{
	switch (i)
	{
		case 1:
			return 2.*NC*(1./x-2.+x-x*x);
		case 2:
			return 2.*NC;
		case 3:
			return beta0/2.;
	}
}


//NLO kernels

double P1NSm(int i,double x)
{
	switch (i)
	{
		case 1:
			return
				(Cf*(4.*Nf*TR*(1.+x)*(-1.+11.*x)+NC*(89.+(-134.+6.*Pi2-223.*x)*x)+
				 6.*Cf*(-27.+Pi2+(27.+ Pi2)*x*x)))/(18.*(1.+x))
				+Li2(-x)*(2.*Cf*(2.*Cf-NC)*(1.+x*x))/(1.+x)
				+log(x)*(Cf*(30.*Cf-23.*NC+4.*Nf*TR+12.*Cf*x+(-24.*Cf+NC+4.*Nf*TR)*x*x))/(6.*(-1.+x))
				-pow(log(x),2.)*(Cf*(2.*NC*(1.+x*x)+Cf*(-1.+x)*(3.+x*(2.+3.*x))))/(2.*(-1.+x*x))
				+log(x)*log(1.-x)*(2.*Cf*Cf*(1.+x*x))/(-1.+x)
				+log(x)*log(1.+x)*(2.*Cf*(2.*Cf-NC)*(1.+x*x))/(1.+x);
		case 2:
			return -(Cf*(NC*(-67.+3.*Pi2)+20.*Nf*TR))/9.;
		case 3:
			return (Cf*(-4.*Nf*(3.+4.*Pi2)*TR+NC*(51.+44.*Pi2-216.*Zeta3)+9.*Cf*(3.-4.*Pi2+48.*Zeta3)))/72.;
	}
}


double P1NSp(int i,double x)
{
	switch(i)
	{
		case 1:
			return
				(Cf*(4.*Nf*TR*(1.+x)*(-1.+11.*x)-6.*Cf*(3.+Pi2+(-3.+Pi2)*x*x)+
				 NC*(-((1.+x)*(-17.+151.*x))+6.*Pi2*(1+x+x*x))))/(18.*(1+x))
				+Li2(-x)*(-2.*Cf*(2.*Cf-NC)*(1.+x*x))/(1.+x)
				+log(x)*(Cf*(6.*Cf*(1.+2.*x)-(11.*NC-4.*Nf*TR)*(1.+x*x)))/(6.*(-1.+x))
				+pow(log(x),2.)*(Cf*(Cf*pow(-1.+x,3.)-2.*NC*x*(1.+x*x)))/(2.*(-1.+x*x))
				+log(x)*log(1.-x)*(2.*Cf*Cf*(1.+x*x))/(-1.+x)
				+log(x)*log(1.+x)*(-2.*Cf*(2.*Cf-NC)*(1.+x*x))/(1.+x);
		case 2:
			return -(Cf*(NC*(-67.+3.*Pi2)+20.*Nf*TR))/9.;
		case 3:
			return (Cf*(-4.*Nf*(3.+4.*Pi2)*TR+NC*(51.+44.*Pi2-216.*Zeta3)+9.*Cf*(3.-4.*Pi2+48.*Zeta3)))/72.;
	}
}


double P1qq(int i,double x)
{
	switch(i)
	{
		case 1:
			return
				(Cf*(4.*Nf*TR*(20.+x+46.*x*x+9.*pow(x,3.)-56.*pow(x,4.))+
				     x*(-6.*Cf*(3.+Pi2+(-3.+Pi2)*x*x)+NC*(-((1.+x)*(-17.+151.*x))+6.*Pi2*(1.+x+x*x)))))/(18.*x*(1.+x))
				+Li2(-x)*(-2.*Cf*(2.*Cf-NC)*(1.+x*x))/(1.+x)
				+log(x)*(Cf*(6.*Cf*(1.+2.*x)-11.*NC*(1.+x*x)+8.*Nf*TR*(-1.+2.*x*(-3.+2.*x*(1.+x)))))/(6.*(-1.+x))
				+pow(log(x),2.)*(Cf*(Cf*pow(-1.+x,3.)-2.*(2.*Nf*TR*(-1.+x)*pow(1.+x,2.)+NC*x*(1.+x*x))))/(2.*(-1.+x*x))
				+log(x)*log(1.-x)*(2.*Cf*Cf*(1.+x*x))/(-1.+x)
				+log(x)*log(1.+x)*(-2.*Cf*(2.*Cf-NC)*(1.+x*x))/(1.+x);
		case 2:
			return -(Cf*(NC*(-67.+3.*Pi2)+20.*Nf*TR))/9.;
		case 3:
			return (Cf*(-4.*Nf*(3.+4.*Pi2)*TR+NC*(51.+44.*Pi2-216.*Zeta3)+9.*Cf*(3.-4.*Pi2+48.*Zeta3)))/72.;
	}
}


double P1qg(int i,double x)
{
	switch (i)
	{
		case 1:
			return
				(Nf*TR*(3.*Cf*x*(42.-87.*x+60.*x*x+Pi2*(-2.-4.*(-1.+x)*x))-
				        2.*NC*(-20.+x*(18.+x*(-225.+6.*Pi2+218.*x)))))/(9.*x)
				-Li2(-x)*4.*NC*Nf*TR*(1.+2.*x*(1.+x))
				+log(x)*(Nf*TR*(6.*NC+8.*NC*x*(6.+11.*x)+3.*Cf*(3.-4.*x+8.*x*x)))/3.
				-log(1.-x)*8.*(Cf-NC)*Nf*TR*(-1.+x)*x
				+pow(log(x),2.)*Nf*TR*(Cf-2.*NC-2.*(Cf+2.*NC)*x+4.*Cf*x*x)
				+pow(log(1.-x),2.)*2.*(Cf-NC)*Nf*TR*(1.+2.*(-1.+x)*x)
				-log(x)*log(1.-x)*4.*Cf*Nf*TR*(1.+2.*(-1.+x)*x)
				-log(x)*log(1.+x)*4.*NC*Nf*TR*(1.+2.*x*(1.+x));
		case 2:
		case 3:
			return 0.;
	}
}


double P1gq(int i,double x)
{
	switch (i)
	{
		case 1:
			return
				(Cf*(-9.*Cf*x*(5.+7.*x)-16.*Nf*TR*(5.+x*(-5.+4.*x))+2.*NC*(9.+x*(19.+6.*Pi2+x*(37.+44.*x)))))/(18.*x)
				+Li2(-x)*(2.*Cf*NC*(2.+x*(2.+x)))/x
				+log(x)*(Cf*(3.*Cf*(4.+7.*x)-2.*NC*(36.+x*(15.+8.*x))))/6.
				+log(1.-x)*(Cf*(-4.*Nf*TR*(2.+(-2.+x)*x)-3.*Cf*(6.+x*(-6.+5.*x))+NC*(22.+x*(-22.+17.*x))))/(3.*x)
				+pow(log(x),2.)*(Cf*(Cf*(-2.+x)+2.*NC*(2.+x)))/2.
				-pow(log(1.-x),2.)*((Cf*(Cf-NC)*(2.+(-2.+x)*x))/x)
				+log(x)*log(1.-x)*(-2.*Cf*NC*(2.+(-2.+x)*x))/x
				+log(x)*log(1.+x)*(2.*Cf*NC*(2.+x*(2.+x)))/x;
		case 2:
		case 3:
			return 0.;
	}
}


double P1gg(int i,double x)
{
	switch(i)
	{
		case 1:
			return
				(24.*Cf*Nf*TR*(-1.+x)*(1.+x)*(-1.+x*(11.+5.*x))+
				 NC*(NC*x*(-((1.+x)*(25.+109.*x))+6.*Pi2*(3.+2.*x*(2.+x+x*x)))+
				     4.*Nf*TR*(-23.+x*(6.+x*(10.+x*(4.+23.*x))))))/(18.*x*(1.+x))
				+Li2(-x)*(4.*NC*NC*pow(1.+x+x*x,2.))/(x*(1.+x))
				+log(x)*(-4.*NC*Nf*TR*(1.+x)-6.*Cf*Nf*TR*(3.+5.*x)+NC*NC*(-25.+11.*(1.-4.*x)*x))/3.
				+pow(log(x),2.)*(-2.*(Cf*Nf*TR*(-1.+x)*pow(1+x,2.)+NC*NC*pow(-1.+(-1.+x)*x,2.)))/(-1.+x*x)
				+log(x)*log(1.-x)*(4.*NC*NC*pow(1.+(-1.+x)*x,2.))/((-1.+x)*x)
				+log(x)*log(1.+x)*(4.*NC*NC*pow(1.+x+x*x,2.))/(x*(1.+x));
		case 2:
			return -(NC*(NC*(-67.+3.*Pi2)+20.*Nf*TR))/9.;
		case 3:
			return -(Cf*Nf*TR)+(NC*(-4.*Nf*TR+NC*(8.+9.*Zeta3)))/3.;
	}
}


//NNLO kernels

double P2NSm(int i,double x)
{
	switch(i)
	{
		case 1:
			return p2nsma_(&x,&nf)/8.;
		case 2:
			return p2nsb_(&x,&nf)/8.;
		case 3:
			return p2nsmc_(&x,&nf)/8.;
	}
}


double P2NSp(int i,double x)
{
	switch(i)
	{
		case 1:
			return p2nspa_(&x,&nf)/8.;
		case 2:
			return p2nsb_(&x,&nf)/8.;
		case 3:
			return p2nspc_(&x,&nf)/8.;
	}
}


double P2NSv(int i,double x)
{
	switch(i)
	{
		case 1:
			return (p2nsma_(&x,&nf)+p2nssa_(&x,&nf))/8.;
		case 2:
			return p2nsb_(&x,&nf)/8.;
		case 3:
			return p2nsmc_(&x,&nf)/8.;
	}
}


double P2qq(int i,double x)
{
	switch(i)
	{
		case 1:
			return (p2nspa_(&x,&nf)+p2psa_(&x,&nf))/8.;
		case 2:
			return p2nsb_(&x,&nf)/8.;
		case 3:
			return p2nspc_(&x,&nf)/8.;
	}
}


double P2qg(int i,double x)
{
	switch(i)
	{
		case 1:
			return p2qga_(&x,&nf)/8.;
		case 2:
		case 3:
			return 0.;
	}
}


double P2gq(int i,double x)
{
	switch(i)
	{
		case 1:
			return p2gqa_(&x,&nf)/8.;
		case 2:
		case 3:
			return 0.;
	}
}


double P2gg(int i,double x)
{
	switch(i)
	{
		case 1:
			return p2gga_(&x,&nf)/8.;
		case 2:
			return p2ggb_(&x,&nf)/8.;
		case 3:
			return p2ggc_(&x,&nf)/8.;
	}
}


//Threshold kernels

double A2ns(int i,double z)
{
	double L=log(z);

	switch(i)
	{
		case 1:
			return Cf*TR*((1.+z*z)/(1.-z)*(2./3.*L*L+20./9.*L)+8./3.*(1.-z)*L+44./27.-268./27.*z);
		case 2:
			return Cf*TR*224./27.;
		case 3:
			return Cf*TR*(-8./3.*Zeta3+40./9.*Zeta2+73./18.);
	}
}


double A2gq(int i,double z)
{
	double M=log(1.-z);

	switch(i)
	{
		case 1:
			return Cf*TR*(4./3.*(2./z-2.+z)*M*M+8./9.*(10./z-10.+8.*z)*M+(448./z-448.+344.*z)/27.);
		case 2:
			return 0.;
		case 3:
			return 0.;
	}
}


double A2gg(int i,double z)
{
	double L=log(z),M=log(1.-z);

	switch(i)
	{
		case 1:
			return Cf*TR*(4./3.*(1.+z)*L*L*L+(6.+10.*z)*L*L+(32.+48.*z)*L-8./z+80.-48.*z-24.*z*z)+NC*TR*(4./3.*(1.+z)*L*L+(52.+88.*z)*L/9.-4./3.*z*M+(556./z-628.+548.*z-700.*z*z)/27.);
		case 2:
			return NC*TR*224./27.;
		case 3:
			return -Cf*TR*15.+NC*TR*10./(9.*27.);
	}
}


double A2hq(int i,double z)
{
	double L=log(z);

	switch(i)
	{
		case 1:
			return Cf*TR*((1.+z)*(32.*S12(1.-z)+16.*L*Li2(1.-z)-16.*Zeta2*L-4./3.*L*L*L)+(32./(3.*z)+8-8.*z-32./3.*z*z)*Li2(1.-z)+(-32./(3.*z)-8.+8.*z+32./3.*z*z)*Zeta2+(2.+10.*z+16./3.*z*z)*L*L-(56./3.+88./3.*z+448./9.*z*z)*L-448./(27.*z)-4./3.-124./3.*z+1600./27.*z*z);
		case 2:
			return 0.;
		case 3:
			return 0.;
	}
}


double A2hg(int i,double z)
{
	double L=log(z),M=log(1.-z),P=log(1.+z),S=S12(1.-z);

	switch(i)
	{
		case 1:
			return Cf*TR*((1.-2.*z+2.*z*z)*(8.*Zeta3+4./3.*M*M*M-8.*M*Li2(1.-z)+8.*Zeta2*L-4.*L*M*M+2./3.*L*L*L-8.*L*Li2(1.-z)+8.*Li3(1.-z)-24.*S)+z*z*(-16.*Zeta2*L+4./3.*L*L*L+16.*L*Li2(1.-z)+32.*S)-(4.+96.*z-64.*z*z)*Li2(1.-z)-(4.-48.*z+40.*z*z)*Zeta2-(8.+48.*z-24.*z*z)*L*M+(4.+8.*z-12.*z*z)*M*M-(1.+12.*z-20.*z*z)*L*L-(52.*z-48.*z*z)*M-(16.+18.*z+48.*z*z)*L+26.-82.*z+80.*z*z) + NC*TR*((1.-2.*z+2.*z*z)*(-4./3.*M*M*M+8.*M*Li2(1.-z)-8.*Li3(1.-z))+(1.+2.*z+2.*z*z)*(-8.*Zeta2*P-16.*P*Li2(-z)-8.*L*P*P+4.*L*L*P+8.*L*Li2(-z)-8.*Li3(-z)-16.*S12(-z))+(16.+64.*z)*(2.*S+L*Li2(1.-z))-(4./3.+8./3.*z)*L*L*L+(8.-32.*z+16.*z*z)*Zeta3-(16.+64.*z)*Zeta2*L+(16.*z+16.*z*z)*(Li2(-z)+L*P)+(32./(3.*z)+12.+64.*z-272./3.*z*z)*Li2(1.-z)-(12.+48.*z-260./3.*z*z+32./(3.*z))*Zeta2-4.*z*z*L*M-(2.+8.*z-10.*z*z)*M*M+(2.+8.*z+46./3.*z*z)*L*L+(4.+16.*z-16.*z*z)*M-(56./3.+172./3.*z+1600./9.*z*z)*L-448./(27.*z)-4./3.-628./3.*z+6352./27.*z*z);
		case 2:
			return 0.;
		case 3:
			return 0.;
	}
}



void OutputFileSplitFuncs(char const* filepath, SplitFunc funcs[], int num_funcs)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix, filepath);
	FILE * f = fopen(filename_buf, "w");

	SplitFunc func;
	double x;
	for (int k=0; k<GRID_PTS; k++)
	{
		x = X[k];
		fprintf(f, "%15.9g\t", x);
		for (int i=0; i<num_funcs; i++)
		{
			func = funcs[i];
			fprintf(f, "%8.4lf\t", (*func)(1, x));
			fprintf(f, "%8.4lf\t", (*func)(2, x));
			fprintf(f, "%8.4lf\t", (*func)(3, x));
		}
		fprintf(f, "\n");
	}
	
	fclose(f);
}


void OutputFileAlphaSBetas(char const* filepath)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix, filepath);
	FILE * f = fopen(filename_buf, "w");

	for (uint _nf=3; _nf<=5; _nf++)
	{
		double _beta0 = Beta0(_nf);
		double _beta1 = Beta1(_nf);
		double _beta2 = Beta2(_nf);
		fprintf(f, "%d\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\n",
				_nf, _beta0, _beta1, _beta2, alpha_pre[_nf], alpha_post[_nf]);
	}
	
	fclose(f);
}


void OutputFileInitialDistributions(char const* filepath)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix, filepath);
	FILE * f = fopen(filename_buf, "w");

	double x;
	for (uint k=0; k<GRID_PTS-1; k++)
	{
		x = X[k];
		fprintf(f, "%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\n",
				x, xuv_3(x), xdv_3(x), xg_3(x), xdb_3(x), xub_3(x), xs_3(x));
	}
	
	fclose(f);
}

void OutputFileGauleg(char const* filepath)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix, filepath);
	FILE * f = fopen(filename_buf, "w");

	for (uint i=0; i<NGP; i++)
	{
		fprintf(f, "%8.4lf\t%8.4lf\n", XG[i], WG[i]);
	}
	
	fclose(f);
}

void OutputFileConvolution(char const* filepath, double vec[])
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix, filepath);
	FILE * f = fopen(filename_buf, "w");

	for (uint i=0; i<GRID_PTS; i++)
	{
		fprintf(f, "%8.4lf\n", vec[i]);
	}
	
	fclose(f);
}

void OutputFileInterpolation(char const* filepath, double vec[])
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix, filepath);
	FILE * f = fopen(filename_buf, "w");

	for (uint i=0; i<67; i++)
	{
		fprintf(f, "%8.4lf\n", vec[i]);
	}
	
	fclose(f);
}



void OutputLOCoefficients(double ***A)
{
	fprintf(stderr, "printing out coefficients...\n");
	char filepath_buf[128] = {0};
	sprintf(filepath_buf, "./output/output-candiav1/%d.dat", _output_file_index++);
	FILE * fp = fopen(filepath_buf, "w");
	if (fp == NULL)
	{
		fprintf(stderr, "could not open file '%s'\n", filepath_buf);
		exit(1);
	}

	for (uint j=13; j<=12+nf; j++)
	{
		for (uint n=0; n<ITERATIONS; n++)
		{
			for (uint k=0; k<GRID_PTS; k++)
				fprintf(fp, "%15.9lf ", A[j][n][k]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}

	for (uint j=32; j<=30+nf; j++)
	{
		for (uint n=0; n<ITERATIONS; n++)
		{
			for (uint k=0; k<GRID_PTS; k++)
				fprintf(fp, "%15.9lf ", A[j][n][k]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}
	
	fclose(fp);

	fprintf(stderr, "done printing out coefficients\n");
}
