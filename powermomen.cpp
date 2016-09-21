/* calculate the power spectrums, of reconstructed velocity and ksz
        2D

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


int powermomen()
{
cout<<"start powermomen"<<endl;
//char inPath1[]={"/project/zhm/ksz/z1/tidepart/v_2d_z1_origin.bin"};
char inPath1[]={"/project/zhm/ksz/z1/result/00momen_2d_z1_origin.bin"};
char inPath2[]={"/project/zhm/ksz/z1/tidepart/00ns_tide_momen_2d_z1_rpar15_kc1.3_l300.bin"};
char plotPath2[]={"/project/zhm/ksz/z1/tidepart/00ns_tide_powermomen_z1_rpar15_kc1.3_l300.dat"};
//variables and fftw setting
//=================================
    int nn=int(nc);
    int nnc=int(nck);
    int ntot=nn*nn;
    int nktot=nn*nnc;

	double renorm;//after transformation, should renormalize powerspectrum
	renorm=box/nc/nc;
	renorm=pow(renorm,2);
  //  double pn[nktot];   
    //=========================to group power spectra

    int kbin=10; //number of bins for power spectrum
    //change it you need to change pk[10] as well
    double kstart,kend,kgap;
    kstart=1.;
    kend=double(nnc);
    kgap=log10(kend/kstart)/kbin;

//=================================================

    double pk[10],pd[10],pkd[10];
    double kcount[10];
    
    for (int i=0;i<kbin;i++)
    {
        pk[i]=0.;
        pd[i]=0.;
        pkd[i]=0.;
        kcount[i]=0.;
    }

    double *deltar1,*deltar2;       //input field
    fftw_complex *deltak1, *deltak2; //Fourier field 
    fftw_plan plan1,plan2;       //fftw setting variable


    deltak1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    deltak2 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));

    deltar1 = (double *) fftw_malloc(ntot * sizeof(double));
    deltar2 = (double *) fftw_malloc(ntot * sizeof(double));

    plan1 = fftw_plan_dft_r2c_2d(nn,nn,deltar1,deltak1,FFTW_ESTIMATE);
    plan2 = fftw_plan_dft_r2c_2d(nn,nn,deltar2,deltak2,FFTW_ESTIMATE);
//======================================================================
//reading data: k, P1s    

//reading data to deltar    
//=====================================
    cout<<"reading  "<<inPath1<<endl;
     cout<<"reading  "<<inPath2<<endl;   
    FILE *in1=fopen(inPath1,"rb");
    FILE *in2=fopen(inPath2,"rb");
float intemp,check;
    for(int i=0;i<ntot;i++)
    {
        check=fread(&intemp,4,1,in1);
        if (check!=1) return -1;
        deltar1[i]=intemp;
        check=fread(&intemp,4,1,in2);
        if (check!=1) return -1;
        deltar2[i]=intemp;
    }
    fclose(in1);
    fclose(in2);
//========================================
cout<<"fftw forward"<<endl;
    fftw_execute(plan1);
    fftw_destroy_plan(plan1);
    fftw_execute(plan2);
    fftw_destroy_plan(plan2);

//========================================
cout<<"calculate power spectra"<<endl;

    int no=0;//for matrix index
    int index; //for grouping

    double pkk,pdd,pdk;//,bk;

    double kx,ky,kmo;
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nnc;k++)
            {
    //           bk=pdk/pdd;
    //           pn[no]=pkk-bk*bk*pdd;//noise

//==============================================================
//group power spectra according to abs(k)
                kx=double(k);
                ky=double(j<nnc?j:(j-nn));
                kmo=sqrt(kx*kx+ky*ky);
                if (kmo==0) 
                {
                    no++;
                    continue;
                }
                pkk=pow(deltak2[no][0],2)+pow(deltak2[no][1],2);
                pdd=pow(deltak1[no][0],2)+pow(deltak1[no][1],2);
                pdk=deltak1[no][0]*deltak2[no][0]+deltak1[no][1]*deltak2[no][1];               
                kmo=log10(kmo);
                index=int(floor(kmo/kgap));
                if (index < kbin)
                {
                    pk[index]+=pkk;
                    pd[index]+=pdd;
                    pkd[index]+=pdk;
                    kcount[index]++;
                }
//=======================================================
              
                no++;
            }
        }
    
//=============================================================
cout<<"save to:  "<<plotPath2<<endl;

    FILE *out1=fopen(plotPath2,"w");

    for (int i=0;i<kbin;i++)
    {
        if (kcount[i]!=0)
        {
            pk[i]/=kcount[i];
            pd[i]/=kcount[i];
            pkd[i]/=kcount[i];
if (isnan(pk[i]*pd[i]*pkd[i])!=0) return -1;

        }
        double kout=pow(10,(i+0.5)*kgap)*dk;
    fprintf(out1,"%e  %e   %e  %e  %e \n",kout,pd[i]*renorm,pk[i]*renorm,pkd[i]*renorm,kcount[i]);
    }
    fclose(out1);


//=======================================
//output
/*
cout<<"save to:  "<<outPath1<<endl;
cout<<"save to:  "<<outPath2<<endl;

    FILE *out1=fopen(outPath1,"wb");
    FILE *out2=fopen(outPath2,"wb");   

float outemp;
    for(int i=0;i<ntot;i++)
    {
        outemp=float(shear1[i]);
        fwrite(&outemp,4,1,out1);
        outemp=float(shear2[i]);   
        fwrite(&outemp,4,1,out2);
    }
    fclose(out1);
    fclose(out2);
*/
//=====================================
cout<<"end powerspectrum"<<endl;
    fftw_free(deltak1);
    fftw_free(deltak2);
    fftw_free(deltar1);
    fftw_free(deltar2);
return 0;
}
  
