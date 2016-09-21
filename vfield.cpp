/* calculate velocity field from density field and compress it into 2D
    
    function:
		v(k)=-i*f*H*delta(k)*k_vec/ksquare          
		f=dlnD/dlna
		
    unit:
        x=Mpc/h
		v=km/s
    
	note: set h=0.678

    OUTPUT:
        3D velocity field, 2D velocity field      */
//=============================================================
//DOngzi Li 2015.Dec, follow Hongming Zhu's program


#include<stdio.h>
#include<iostream>
using namespace std;
#include<math.h>
#include<stdlib.h>
#include<fftw3.h>

#include"tide.h"

int vfield()
{
double h=3;
double fH=0.96*h*100;
cout<<"start vfield"<<endl;

char inPath[]={"/project/zhm/ksz/z1/tidekperp/00k3d_clean_z2_rpar10_kc0.5_l300.bin"};
char outPath1[]={"/project/zhm/ksz/z1/tidekperp/00v_3d_z2_rpar10_kc0.5_l300.bin"};

//variables and fftw setting
//=====================================
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
//========================================
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
cout<<"calculate kz velocity field"<<endl;

    int no=0;//for matrix index
    double ksquare,kx,ky,kz;
    double ipart,rpart,paraz;
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
                ksquare=kx*kx+ky*ky+kz*kz;
				if (ksquare==0) paraz=0; //drop this mode
				else paraz=fH*kz/ksquare/dk;

                rpart=deltak[no][0];
                ipart=deltak[no][1];

				deltak[no][0]=-ipart*paraz;
				deltak[no][1]=rpart*paraz;
              }
                no++;
            }
        }
    }
//======================================
cout<<"fftw back"<<endl;    
    fftw_execute(planback);
//======================================
cout<<"save to:  "<<outPath1<<endl;
    FILE *out1=fopen(outPath1,"wb");

    double tot=double(ntot);
    float outemp;
    for(int i=0;i<ntot;i++)
    {   
        if (isnan(deltars[i])!=0)
        {   
            cout<<"nan appears, no: "<<i<<endl;
            return -1;
        }
        deltars[i]/=tot;//renormalization. fftw doesn't compose /L^3

        outemp=float(deltars[i]);
        fwrite(&outemp,4,1,out1);
    }
       fclose(out1);
//========================================
    fftw_destroy_plan(planback);
    fftw_free(deltak);
    fftw_free(deltar);
    fftw_free(deltars);  
    cout<<"end velocity field calculation"<<endl;
return 0;
}
  
