/* smooth 21cm density field
    
    function:
        delta(x)=integ dx^3 S(x-x') delta(x')
        S(x-x')=exp(-r*r/(2R*R))

    smoothing scale: default=1.25 Mpc/h
	foreground substract scale: 
			1. high k cut: kc=0.5 h/Mpc
			2. foreground substract: Rpar=15,60 Mpc/h;
    unit:
        x=Mpc/h
        
    OUTPUT:
        delta(k): smooth_deltak.dat      */
//=============================================================
//DOngzi Li 2015.Dec, follow Hongming Zhu's program


#include<stdio.h>
#include<iostream>
using namespace std;
#include<math.h>
#include<stdlib.h>
#include<fftw3.h>

#include"tide.h"
//Global variable for whole project:
/*
double PI=3.1415926;
double nc=1024.;    //number of cells per dimension
double nck=(nc/2.+1); //number of cells for fftw output
double box=1200.;    //N-body simulation box size, 1.2Gpc/h
double dk=2.*PI/box;//fundamental frequency
*/

int smooth_fs()
{
//for smooth_fs
double kc=0.5;  //high frequency cut off h/Mpc
double Rpar=15;//foreground substract, higher signal better
double R=1.25; // smooth scale

char inPath[80]={"/home/zhm/tides00/1.000den00.bin"};
char outPath[80]={"/home/zhm/dongzi/tidez1/3dtidenew/00smooth1.25_z1_rpar15.bin"};
char outPath2[80]={"/home/zhm/dongzi/tidez1/3dtidenew/00slice_smooth1.25_z1_rpar15.bin"};


//variables and fftw setting
//=================================
    double *deltar;       //input field
    double *deltars;      //output field
    fftw_complex *deltak; //Fourier field
    fftw_plan plan,planback;       //fftw setting variable

    int nn=int(nc);
    int nnc=int(nck);
    int ntot=nn*nn*nn;
    int nktot=nn*nn*nnc;
    deltak = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    deltar = (double *) fftw_malloc(ntot * sizeof(double));
    deltars = (double *) fftw_malloc(ntot * sizeof(double));

    plan = fftw_plan_dft_r2c_3d(nn,nn,nn,deltar,deltak,FFTW_ESTIMATE);
    planback = fftw_plan_dft_c2r_3d(nn,nn,nn,deltak,deltars,FFTW_ESTIMATE);

//reading data to deltar    
//=====================================
    cout<<"reading  "<<inPath<<endl;
    
    float intemp,check;
    FILE *in=fopen(inPath,"rb");
    for(int i=0;i<ntot;i++)
    {
        check=fread(&intemp,4,1,in);
        if (check!=1)
        {cout<<"reading problem. check="<<check<<endl;
         return -1;
        }

        deltar[i]=intemp-1.;

            }
    fclose(in);
//========================================
cout<<"fftw forward"<<endl;
    fftw_execute(plan);
    fftw_destroy_plan(plan);

//========================================
cout<<"smooth in fourier space"<<endl;

    int no=0;//for matrix index
    double ksquare,kx,ky,kz;
	double sincx,sincy,sincz;
    double smooth;
	double winfs;

	double kc2=kc/dk;//to compare with ksquare
	kc2*=kc2;
	double Rpar2=Rpar*Rpar;

    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nnc;k++)
            {
              
                kx=double(k);
                ky=double(j<nnc?j:(j-nn));
                kz=double(i<nnc?i:(i-nn));
                ksquare=kx*kx+ky*ky+kz*kz;
				//high energy cut off
				if (ksquare > kc2) 
				{
					deltak[no][0]=0.;
					deltak[no][1]=0.;
					no++;
					continue;
				}

				//foreground substraction information loss
				winfs=1.-exp(-kz*kz*dk*dk*Rpar2/2.);

				//deconvolve the filter in N-body simulation
				if (kx==0) sincx=1.; else sincx=sin(PI*kx/nc)/(PI*kx/nc);
				if (ky==0) sincy=1.; else sincy=sin(PI*ky/nc)/(PI*ky/nc);
				if (kz==0) sincz=1.; else sincz=sin(PI*kz/nc)/(PI*kz/nc);

                smooth=exp(-R*R/2.*ksquare*dk*dk)/sincx/sincy/sincz;

                deltak[no][0]*=(smooth*winfs);
                deltak[no][1]*=(smooth*winfs);
              
                no++;
            }
        }
    }
//======================================
cout<<"fftw back"<<endl;    
    fftw_execute(planback);
//======================================
//output and gaussianize

cout<<"save to:  "<<outPath<<endl;
cout<<"save to:  "<<outPath2<<endl;

    FILE *out=fopen(outPath,"wb");
	FILE *out2=fopen(outPath2,"wb");

    double tot=double(ntot);
    float outemp;
	no=0;
    for(int i=0;i<nn;i++)
    {   
		for (int j=0;j<nn;j++)
		{
			for (int k=0;k<nn;k++)
			{
			  if (isnan(deltars[no])!=0)
			  {   
				 cout<<"nan appears, no: "<<i<<endl;
				 return -1;
				 }
				 deltars[no]/=tot;//renormalization. fftw doesn't compose /L^3

        //semi gaussianization, set deltars[no]=log(deltars[no]+1)
			  if (deltars[no]>0) 
			   {
				 deltars[no]=log(1.+deltars[no]);
				 }

				 outemp=float(deltars[no]+1.);
				 fwrite(&outemp,4,1,out);
				
				 //for slice iamge
				 if (i==512) fwrite(&outemp,4,1,out2);
                 no++;
				 }
		}
	}
				 fclose(out);
				 fclose(out2);
//=====================================
    fftw_destroy_plan(planback);
    fftw_free(deltak);
    fftw_free(deltar);
    fftw_free(deltars);  
    cout<<"end smooth"<<endl;
return 0;
}
  
