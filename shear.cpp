/* calculate the gamma1,2 shear field, and then reconstruct the 3d noisy density field
    
    function:
        delta(x)=integ dx^3 S(x-x') delta(x')
        S(x-x')=exp(-r*r/(2R*R))

    smoothing scale: default=1.25 Mpc/h
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
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_errno.h>

#include"tide.h"
//Global variable for whole project:
/*
double PI=3.1415926;
double nc=40.;     //number of cells per dimension
double nck=(nc/2.+1); //number of cells for fftw output
double box=1200.;    //N-body simulation box size, 1.2Gpc/h
double dk=2.*PI/box;//fundamental frequency
double R=1.25;      //smoothing scales
*/
int shear()
{
cout<<"start shear field calculation"<<endl;
/*    
char outPath1[]={"/home/zhm/dongzi/tidez1/pipdata/shear1.bin"};
char outPath2[]={"/home/zhm/dongzi/tidez1/pipdata/shear2.bin"};
char inPath[]={"/home/zhm/dongzi/tidez1/pipdata/smooth1.25lpk.bin"};
char inPath_pl[]={"/home/zhm/dongzi/tidez1/pipdata/coeffz1.dat"};
*/

char outPath[]={"/project/zhm/ksz/z1/tidekperp/00shear"};
char outname[];

char outPath2[]={"/project/zhm/ksz/z1/tidekperp/00k3d_noisy_z1_rpar15_gaz0.5.bin"};
char inPath[]={"/project/zhm/ksz/z1/tidekperp/00smooth1.25_z1_rpar15.bin"};
char inPath_pl[]={"/home/zhm/dongzi/tidez1/coeffz1.dat"};


//double alpha=0.9709,beta=2.1645;//for a=1,z=0
//variables and fftw setting
//=================================
    int nn=int(nc);
    int nnc=int(nck);
    int ntot=nn*nn*nn;
    int nktot=nn*nn*nnc;

    double *deltar;       //input field
    double *deltars1;
    double *deltars2;
    double *deltars3;      //output field

    fftw_complex *deltak1; //Fourier field g,w1
    fftw_complex *deltak2; //Fourier field g,w2
    fftw_complex *deltak3; //Fourier field g,w2
    fftw_plan plan,planback1,planback2,planback3;       //fftw setting variable
    deltak1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    deltak2 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    deltak3 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));

    deltar = (double *) fftw_malloc(ntot * sizeof(double));
    deltars1 = (double *) fftw_malloc(ntot * sizeof(double));
    deltars2 = (double *) fftw_malloc(ntot * sizeof(double));
    deltars3 = (double *) fftw_malloc(ntot * sizeof(double));
    

    plan = fftw_plan_dft_r2c_3d(nn,nn,nn,deltar,deltak1,FFTW_ESTIMATE);
    planback1 = fftw_plan_dft_c2r_3d(nn,nn,nn,deltak1,deltars1,FFTW_ESTIMATE);
    planback2 = fftw_plan_dft_c2r_3d(nn,nn,nn,deltak2,deltars2,FFTW_ESTIMATE);
    planback3 = fftw_plan_dft_c2r_3d(nn,nn,nn,deltak3,deltars3,FFTW_ESTIMATE);

//reading data to deltar    
//=====================================
    cout<<"reading  "<<inPath<<endl;
    
    float intemp,check;
    FILE *in=fopen(inPath,"rb");
    for(int i=0;i<ntot;i++)
    {
        check=fread(&intemp,4,1,in);
        if (check!=1) return -1;
        deltar[i]=intemp-1.;
    }
    fclose(in);
//========================================
cout<<"fftw forward"<<endl;
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(deltar);
//========================================

//======================
//for window calculation, related with linear power spectrum
//==============================================================
//reading data: k, P1s    
cout<<"reading   "<<inPath_pl<<endl;
    int nw=500;//number of linear power spectrum
    double kl[500],pl[500];
    double ktemp,ptemp;
    in=fopen(inPath_pl,"r");
    for(int i=0;i<nw;i++)
    {
        check=fscanf(in,"%lf %lf \n",&ktemp,&ptemp);
        if (check!=2)
        {
            cout<<"error reading "<<inPath_pl<<endl;
            return -1;
        }

       //prepare for dlogp/dlogk 
        kl[i]=ktemp;
        pl[i]=ptemp;
    }
    fclose(in);
//=============================
cout<<"gsl set up"<<endl;
gsl_interp_accel *acc = gsl_interp_accel_alloc();
gsl_spline *spline=gsl_spline_alloc (gsl_interp_cspline ,nw);
gsl_spline_init(spline,kl,pl,nw);

    
 //=====================================================

cout<<"convolve in fourier space"<<endl;

    int no=0;//for matrix index
    double kabs,kx,ky,kz,ipart,rpart;
    double win;
    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nnc;k++)
            {
              if(no<nktot)
              {
                kx=double(k);
                ky=double(j<nnc?j:(j-nn));
                kz=double(i<nnc?i:(i-nn));
                kabs=kx*kx+ky*ky+kz*kz;
                kabs=sqrt(kabs)*dk;
//==============================================================
//calculate window
//========================================
           if (kabs>50 ||kabs<3e-3)
           {
				win=1.;                
		   }
		   else
		   {
			   win=gsl_spline_eval (spline,kabs,acc);
		   }
//=================================================================
				//multiply i
				rpart=deltak1[no][0];
				ipart=deltak1[no][1];

                deltak1[no][0]=-ipart*win*kx*dk;
                deltak1[no][1]=rpart*win*kx*dk;
                deltak2[no][0]=-ipart*win*ky*dk;
                deltak2[no][1]=rpart*win*ky*dk;
                deltak3[no][0]=-ipart*win*kz*dk;
                deltak3[no][1]=rpart*win*kz*dk;
						no++;
              }
            }
        }
    }
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

//=============================
/*
FILE *wtf=fopen("/home/zhm/dongzi/tidez1/pipdata/sheark.dat","w");
  for (int i=0;i<nktot;i++)
  {fprintf(wtf,"%e  %e  %e  %e \n",deltak1[i][0],deltak1[i][1],deltak1[i][0],deltak2[i][1]);
  }
fclose(wtf);
*/
//======================================
cout<<"fftw back"<<endl;    
    fftw_execute(planback1);
    fftw_execute(planback2);
    fftw_execute(planback3);
    fftw_destroy_plan(planback1);
    fftw_destroy_plan(planback2);
    fftw_destroy_plan(planback3);
    fftw_free(deltak1);
    fftw_free(deltak2);
    fftw_free(deltak3);

//======================================
cout<<"calculate shear1,2,x,y,z"<<endl;
    double *shear1;
    double *shear2;   
    double *shearx;
    double *sheary;
    double *shearz;

    shear1 = (double *) fftw_malloc(ntot * sizeof(double));
    shear2 = (double *) fftw_malloc(ntot * sizeof(double));
    shearx = (double *) fftw_malloc(ntot * sizeof(double));
    sheary = (double *) fftw_malloc(ntot * sizeof(double));
    shearz = (double *) fftw_malloc(ntot * sizeof(double));
    
    double tot=double(ntot);
    for(int i=0;i<ntot;i++)
    {
        //renormaliza
      //  cout<<deltars1[i];   
        deltars1[i]/=tot;
        deltars2[i]/=tot;
        deltars3[i]/=tot;
//cout<<deltars2[i];
        shear1[i]=pow(deltars1[i],2)-pow(deltars2[i],2);
        shear2[i]=2.*deltars1[i]*deltars2[i];
        shearx[i]=2.*deltars1[i]*deltars3[i];
        sheary[i]=2.*deltars2[i]*deltars3[i];
        shearz[i]=2.*pow(deltars3[i],2)-pow(deltars1[i],2)-pow(deltars2[i],2); //correct the paper
    }
//===================================================================
    fftw_free(deltars1);
    fftw_free(deltars2);  
    fftw_free(deltars3);  
//================================
    cout<<"save 5 shear field"<<endl;

FILE *outshear[5];
int shearno=5;
    for(int i=1;i<=shearno;i++)
    {
        sprintf(outname,"%s%d%s",outPath,i,".bin");
        outshear[i]=fopen(outname,"wb");
    }

cout<<outname<<endl;

    float outempshear;
    for(int i=0;i<ntot;i++)
    {
            outempshear=float(shear1[i]);
            fwrite(&outempshear,4,1,outshear[0]);
            if (isnan(outempshear)!=0) return -1;
    
            outempshear=float(shear2[i]);
            fwrite(&outempshear,4,1,outshear[1]);
            if (isnan(outempshear)!=0) return -1;

            outempshear=float(shearx[i]);
            fwrite(&outempshear,4,1,outshear[2]);
            if (isnan(outempshear)!=0) return -1;

            outempshear=float(sheary[i]);
            fwrite(&outempshear,4,1,outshear[3]);
            if (isnan(outempshear)!=0) return -1;

            outempshear=float(shearz[i]);
            fwrite(&outempshear,4,1,outshear[4]);
            if (isnan(outempshear)!=0) return -1;
    }






//======================================================================
cout<<"start to calculate noisy field"<<endl;
	double *k3dr;
	fftw_complex *k3d;
	fftw_complex *sheark1,*sheark2;
    fftw_complex *shearkx,*shearky,*shearkz;

    sheark1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    sheark2 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    shearkx = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    shearky = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    shearkz = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
   

    k3d = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    k3dr = (double *) fftw_malloc(ntot * sizeof(double));

	//transform again to calculate the noisy field
    fftw_plan plan1,plan2,planx,plany,planz;
    fftw_plan planbackk;
    plan1 = fftw_plan_dft_r2c_3d(nn,nn,nn,shear1,sheark1,FFTW_ESTIMATE);
    plan2 = fftw_plan_dft_r2c_3d(nn,nn,nn,shear2,sheark2,FFTW_ESTIMATE);
    planx = fftw_plan_dft_r2c_3d(nn,nn,nn,shearx,shearkx,FFTW_ESTIMATE);
    plany = fftw_plan_dft_r2c_3d(nn,nn,nn,sheary,shearky,FFTW_ESTIMATE);
    planz = fftw_plan_dft_r2c_3d(nn,nn,nn,shearz,shearkz,FFTW_ESTIMATE);

    planbackk = fftw_plan_dft_c2r_3d(nn,nn,nn,k3d,k3dr,FFTW_ESTIMATE);
//========================================
cout<<"fftw forward"<<endl;
    fftw_execute(plan1);
    fftw_execute(plan2);
    fftw_execute(planx);
    fftw_execute(plany);
    fftw_execute(planz);
    fftw_destroy_plan(plan1);
    fftw_destroy_plan(plan2);
    fftw_destroy_plan(planx);
    fftw_destroy_plan(plany);
    fftw_destroy_plan(planz);
    fftw_free(shear1);
    fftw_free(shear2);
    fftw_free(shearx);
    fftw_free(sheary);
    fftw_free(shearz);
//========================================
cout<<"calculate noisy k3d"<<endl;

    no=0;//for matrix index
    double kx2,ky2,kz2,k2;

    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nnc;k++)
            {
              if(no<nktot)
              {
                kx=double(k);
                ky=double(j<nnc?j:(j-nn));
                kz=double(i<nnc?i:(i-nn));
                if (kx==0&&ky==0&&kz==0)
                {
                    k3d[no][0]=0.;
                    k3d[no][1]=0.;
                    no++;
                    continue;
                }

                kx2=kx*kx;
                ky2=ky*ky;
                kz2=kz*kz;
                k2=kx2+ky2+kz2;
                k3d[no][0]=1./k2*(
                           (kx2-ky2)*sheark1[no][0]+2.*kx*ky*sheark2[no][0]
                           +(2.*kx*kz*shearkx[no][0]+2.*ky*kz*shearky[no][0]
                           +(2.*kz2-kx2-ky2)*shearkz[no][0])/2.);
                k3d[no][1]=1./k2*(
                           (kx2-ky2)*sheark1[no][1]+2.*kx*ky*sheark2[no][1]
                           +(2.*kx*kz*shearkx[no][1]+2.*ky*kz*shearky[no][1]
                           +(2.*kz2-kx2-ky2)*shearkz[no][1])/2.);
//cout<<sheark1[no][0];
                no++;
			  }
            }
        }
    }
    
     fftw_free(sheark1);
     fftw_free(sheark2);
     fftw_free(shearkx);
     fftw_free(shearky);
     fftw_free(shearkz);
//=====================================    
    cout<<"fftw back"<<endl;    
    fftw_execute(planbackk);
	fftw_destroy_plan(planbackk);
     fftw_free(k3d);
//======================================
//output

cout<<"save to:  "<<outPath2<<endl;
    FILE *out=fopen(outPath2,"wb");

    float outemp;
    for(int i=0;i<ntot;i++)
    {
        outemp=float(k3dr[i]/tot);
        fwrite(&outemp,4,1,out);
        if (isnan(outemp)!=0) return -1;

       // cout<<outemp;
    }
       fclose(out);

//=====================================
     fftw_free(k3dr);
	
    cout<<"end noisy field calculation"<<endl;
return 0;
}
  
