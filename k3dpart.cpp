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
cout<<"start seperate field calculation"<<endl;
char inPath[]={"/home/zhm/dongzi/tidez1/3dtidenew/00shear"};
char inname[];

char backupPath[]={"/home/zhm/dongzi/tidez1/3dtidenew/00r_part_rpar15_kc0.6_bin10.dat"};
char inPath1[]={"/home/zhm/tides00/1.000den00.bin"};


//double alpha=0.9709,beta=2.1645;//for a=1,z=0
//variables and fftw setting
//=================================
    int nn=int(nc);
    int nnc=int(nck);
    int ntot=nn*nn*nn;
    int nktot=nn*nn*nnc;
	double renorm;//after transformation, should renormalize powerspectrum
	renorm=box/nc/nc;
	renorm=pow(renorm,3);

    int kbin=10; //number of bins for power spectrum
    //change it you need to change pk[10] as well
    double kstart,kend,kgap;
    kstart=1.;
    kend=double(nnc);
    kgap=log10(kend/kstart)/kbin;
    cout<<"kbin, kend, kgap"<<kbin<<" "<<kend<<" "<<kgap<<endl;
    
    double pk[10][10],pd[10][10],pkd[10][10],bb[10][10],pn[10][10],window[10][10];
    double kcount[10][10];
    
    for (int i=0;i<kbin;i++)
		for (int j=0;j<kbin;j++)
		{
    {
        pk[i][j]=0.;
        pd[i][j]=0.;
        pkd[i][j]=0.;
        kcount[i][j]=0.;
		bb[i][j]=0.;
		pn[i][j]=0.;
		window[i][j]=0.;
    }
		}
//======================================
cout<<"reading shear field"<<endl;
    double *shear1;
    double *shear2;   
    double *shearx;
    double *shearz;

    shear1 = (double *) fftw_malloc(ntot * sizeof(double));
    shear2 = (double *) fftw_malloc(ntot * sizeof(double));
    shearx = (double *) fftw_malloc(ntot * sizeof(double));
    shearz = (double *) fftw_malloc(ntot * sizeof(double));

    FILE *inshear[5];
    int shearno=5;
    for(int i=1;i<=shearno;i++)
    {
        sprintf(inname,"%s%d%s",inPath,i,".bin");
        inshear[i]=fopen(inname,"rb");
    }

cout<<inname<<endl;

    float intempshear,check;
    for(int i=0;i<ntot;i++)
    {
        check=fread(&intempshear,4,1,inshear[0]);
        if (check!=1) return -1;
        shear1[i]=intempshear;

        check=fread(&intempshear,4,1,inshear[1]);
        if (check!=1) return -1;
        shear2[i]=intempshear;

        check=fread(&intempshear,4,1,inshear[2]);
        if (check!=1) return -1;
        shearx[i]=intempshear;

        check=fread(&intempshear,4,1,inshear[4]);
        if (check!=1) return -1;
        shearz[i]=intempshear;
    }


//======================================================================
cout<<"start to calculate noisy field"<<endl;
	double *delta;
	fftw_complex *deltak1;
    //for original delta field

	fftw_complex *sheark1,*sheark2;
    fftw_complex *shearkx,*shearkz;

    sheark1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    sheark2 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    shearkx = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    shearkz = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
   

    deltak1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    delta = (double *) fftw_malloc(ntot * sizeof(double));

	//transform again to calculate the noisy field
    fftw_plan plan1,plan2,planx,plany,planz;
    fftw_plan plandelta;
    plan1 = fftw_plan_dft_r2c_3d(nn,nn,nn,shear1,sheark1,FFTW_ESTIMATE);
    plan2 = fftw_plan_dft_r2c_3d(nn,nn,nn,shear2,sheark2,FFTW_ESTIMATE);
    planx = fftw_plan_dft_r2c_3d(nn,nn,nn,shearx,shearkx,FFTW_ESTIMATE);
    planz = fftw_plan_dft_r2c_3d(nn,nn,nn,shearz,shearkz,FFTW_ESTIMATE);

    plandelta = fftw_plan_dft_r2c_3d(nn,nn,nn,delta,deltak1,FFTW_ESTIMATE);
//========================================
cout<<"fftw forward"<<endl;

    fftw_execute(plan1);
    fftw_execute(plan2);
    fftw_execute(planx);
    fftw_execute(planz);
    fftw_destroy_plan(plan1);
    fftw_destroy_plan(plan2);
    fftw_destroy_plan(planx);
    fftw_destroy_plan(planz);
    fftw_free(shear1);
    fftw_free(shear2);
    fftw_free(shearx);
    fftw_free(shearz);
//========================================
    //reading original delta
     cout<<"reading  "<<inPath2<<endl;   
    float intemp,check;
    FILE *in1=fopen(inPath1,"rb");

    for(int i=0;i<ntot;i++)
    {
        check=fread(&intemp,4,1,in1);
        if (check!=1) return -1;
        deltar[i]=intemp-1.;
    }
    fclose(in1);

    fftw_execute(plandelta);
    fftw_destroy_plan(plandelta);
    fftw_free(delta);
 
//========================================

cout<<"calculate noisy k3d in different parts, and correlation"<<endl;
	int negco=0; //to count how many negative correlation

    int no=0;//for matrix index
    int index1,index2; //for grouping

    double pkk[4],pdd[4],pdk[4];//,bk;
    double kx,ky,kz,korth,kpar;

	double sincx,sincy,sincz,rev_tophat;
    no=0;//for matrix index
    double kx2,ky2,kz2,k2;
    double k3d[4][2];

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
                    no++;
                    continue;
                }

                kx2=kx*kx;
                ky2=ky*ky;
                kz2=kz*kz;
                k2=kx2+ky2+kz2;
              
                k3d[0][0]=(kx2-ky2)*sheark1[no][0]/k2;
                k3d[1][0]=2.*kx*ky*sheark2[no][0]/k2;
                k3d[2][0]=2.*kx*kz*shearkx[no][0]/k2;
                k3d[3][0]=(2.*kz2-kx2-ky2)*shearkz[no][0]/k2;

                k3d[0][1]=(kx2-ky2)*sheark1[no][1]/k2;
                k3d[1][1]=2.*kx*ky*sheark2[no][1]/k2;
                k3d[2][1]=2.*kx*kz*shearkx[no][1]/k2;
                k3d[3][1]=(2.*kz2-kx2-ky2)*shearkz[no][1]/k2;
                //===============================================

				if (kx==0) sincx=1.; else sincx=sin(PI*kx/nc)/(PI*kx/nc);
				if (ky==0) sincy=1.; else sincy=sin(PI*ky/nc)/(PI*ky/nc);
				if (kz==0) sincz=1.; else sincz=sin(PI*kz/nc)/(PI*kz/nc);

                rev_tophat=1./sincx/sincy/sincz;

                deltak1[no][0]*=rev_tophat;
                deltak1[no][1]*=rev_tophat;
//==================================================================
                korth=sqrt(kx2+ky2);
                kz=fabs(kz);
                
                if (korth==0)korth++;
                if(kz==0)kz++;

                   korth=log10(korth);
                   kz=log10(kz);

                index1=int(floor(korth/kgap));
				index2=int(floor(kz/kgap));


                if (index1 < kbin && index2 < kbin)// && pdk > 0)
                {

                for(int shearno=0;shearno<4;shearno++)
                {
                pkk[shearno]=pow(k3d[shearno][0],2)+pow(k3d[shearno][1],2);
                pdd[shearno]=pow(deltak1[no][0],2)+pow(deltak1[no][1],2);
                pdk[shearno]=deltak1[no][0]*k3d[shearno][0]+deltak1[no][1]*k3d[shearno][1];               
                    pk[shearno][index1][index2]+=pkk[shearno];
                    pd[shearno][index1][index2]+=pdd[shearno];
                    pkd[shearno][index1][index2]+=pdk[shearno];
                    kcount[shearno][index1][index2]++;
                }
                }
                no++;
			  }
            }
        }
    }
    
     fftw_free(sheark1);
     fftw_free(sheark2);
     fftw_free(shearkx);
     fftw_free(shearkz);
     fftw_free(deltak1);
//======================================
//output

//calculate b(k//,k orth), Pn(k//,k orth)

	 for (int i=0;i<kbin;i++)
     {
		for (int j=0;j<kbin;j++)
		{
				if (pd[i][j]>0)
				{
				bb[i][j]=pkd[i][j]/pd[i][j];
				pn[i][j]=pk[i][j]-pow(bb[i][j],2)*pd[i][j];

				window[i][j]=pd[i][j]/(pd[i][j]+pn[i][j]/pow(bb[i][j],2));
               // if (bb[i][j]<0) bb[i][j]=0.;// drop these mode
				}	 
			 	}
		}

cout<<"backup save to: "<<backupPath<<endl;

	FILE *out2=fopen(backupPath,"w");
		fprintf(out2," negative b:  %i \n",negco);
		fprintf(out2," b, window, pn, pd,pk,pkd, count \n");

    for (int i=0;i<kbin;i++)
    {
		for (int j=0;j<kbin;j++)
		{
            
        if (kcount[i][j]!=0)
        {
            pn[i][j]/=kcount[i][j];
            pd[i][j]/=kcount[i][j];
            pkd[i][j]/=kcount[i][j];
            pk[i][j]/=kcount[i][j];
        }
        
	fprintf(out2,"%e  %e  %e  %e  %e  %e  %e \n",bb[i][j],window[i][j],pn[i][j]*renorm,pd[i][j]*renorm,pk[i][j]*renorm,pkd[i][j]*renorm,kcount[i][j]);
    }
	}
	fclose(out2);
return 0;
}
  
