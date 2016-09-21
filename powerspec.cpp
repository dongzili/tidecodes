/* calculate the power spectrums, P_k3d, P_delta, P_k3d,delta
    
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


int powerspec()
{


char outPath[]={"/project/zhm/ksz/z1/tidekperp/00k3d_clean_z1_rpar15_kc0.6_l300.bin"};
//char plotPath1[]={"/project/zhm/ksz/z1/tidekperp/noise.dat"};
char plotPath2[]={"/project/zhm/ksz/z1/tidekperp/00powerspec_rpar15_kc0.6_l300.dat"};
char inPath1[]={"/home/zhm/tidesData/tides00/1.000den00.bin"};
char inPath2[]={"/project/zhm/ksz/z1/tidekperp/00k3d_noisy_z1_rpar15_kc0.6_l300.bin"};
char inPath_fil[]={"/project/zhm/ksz/z1/tidekperp/00noisebias_rpar15_kc0.6_l300.dat"};
//char inPath_fil[]={"/project/zhm/ksz/z1/tidekperp/10den_noisebias_z1_rpar15.dat"};
//variables and fftw setting
//=================================
//    double kcut=kc/dk; //for reconstructed momentum higher than this, cut!

    int nn=int(nc);
    int nnc=int(nck);
    int ntot=nn*nn*nn;
    int nktot=nn*nn*nnc;

	double renorm;//after transformation, should renormalize powerspectrum
	renorm=box/nc/nc;
	renorm=pow(renorm,3);
    double renorm_weigh=0.;//for weight average of the field
  //  double pn[nktot];   
    //=========================to group power spectra

    int kbin=10; //number of bins for power spectrum
    //change it you need to change pk[10/[10] as well
    double kstart,kend,kgap;
    kstart=1.;
    kend=double(nnc);
    kgap=log10(kend/kstart)/kbin;

//================================================

    double pk[10][10],pd[10][10],pkd[10][10];
    double count[10][10],kcount[10][10],pcount[10][10];//one for kappa, one for power_kappa renorm
    
    for (int i=0;i<kbin;i++)
    {
        for (int j=0;j<kbin;j++)
        {
        pk[i][j]=0.;
        pd[i][j]=0.;
        pkd[i][j]=0.;
        count[i][j]=0.;
        kcount[i][j]=0.;
        pcount[i][j]=0.;
        }
    }

    double *deltar1,*deltar2,*deltars;       //input field
    fftw_complex *deltak1, *deltak2; //Fourier field 
    fftw_plan plan1,plan2,planback;       //fftw setting variable


    deltak1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    deltak2 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));

    deltar1 = (double *) fftw_malloc(ntot * sizeof(double));
    deltar2 = (double *) fftw_malloc(ntot * sizeof(double));
    deltars = (double *) fftw_malloc(ntot * sizeof(double));

    plan1 = fftw_plan_dft_r2c_3d(nn,nn,nn,deltar1,deltak1,FFTW_ESTIMATE);
    plan2 = fftw_plan_dft_r2c_3d(nn,nn,nn,deltar2,deltak2,FFTW_ESTIMATE);
	planback = fftw_plan_dft_c2r_3d(nn,nn,nn,deltak2,deltars,FFTW_ESTIMATE);
//======================================================================
//reading data: k, P1s    
cout<<"reading   "<<inPath_fil<<endl;
    float intemp,check;
    int nbin=10;//number of linear power spectrum
    double bb[10][10],window[10][10];
    double ktemp,ptemp;
FILE *in=fopen(inPath_fil,"r");
    for(int i=0;i<nbin;i++)
		for (int j=0;j<nbin;j++)
	{
		 {
        check=fscanf(in,"%lf %lf \n",&ktemp,&ptemp);
        if (check!=2)
        {
            cout<<"error reading "<<inPath_fil<<endl;
            return -1;
        }

       //prepare for dlogp/dlogk 
        bb[i][j]=ktemp;
        window[i][j]=ptemp;
		}    
	}
    fclose(in);
//=============================

//reading data to deltar    
//=====================================
    cout<<"reading  "<<inPath1<<endl;
     cout<<"reading  "<<inPath2<<endl;   
    FILE *in1=fopen(inPath1,"rb");
    FILE *in2=fopen(inPath2,"rb");

    for(int i=0;i<ntot;i++)
    {
        check=fread(&intemp,4,1,in1);
        if (check!=1) return -1;
        deltar1[i]=intemp-1.;

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
    fftw_free(deltar1);
    fftw_free(deltar2);

//========================================
cout<<"calculate power spectra"<<endl;

    int no=0;//for matrix index
    int index1,index2; //for grouping
    double rev_tophat,sincx,sincy,sincz;
    double pkk,pdd,pdk;//,bk;

    double kx,ky,kz,kmo,korth,kpar;
    for(int i=0;i<nn;i++)
    {
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
                kz=double(i<nnc?i:(i-nn));
                kmo=sqrt(kx*kx+ky*ky+kz*kz);

				korth=sqrt(kx*kx+ky*ky);
				kpar=fabs(kz);

                if (korth==0 && kpar==0) //korth==0,infinite noise, kpar=0 don't contribute to vz
                {
					deltak2[no][0]=0.;
					deltak2[no][1]=0.;
                    no++;
                    continue;
                }

                if (korth==0) korth++;
                if (kpar==0) kpar++;
                korth=log10(korth);
				kpar=log10(kpar);
                
				index1=int(floor(korth/kgap));
				index2=int(floor(kpar/kgap));
/*
                pdk=deltak1[no][0]*deltak2[no][0]+deltak1[no][1]*deltak2[no][1];               
				if(kmo>kcut)//drop large modes
				{
					deltak2[no][0]=0.;
					deltak2[no][1]=0.;
					no++;
					continue;
				}
*/
//====================================================================
//for real effective modes
				if (index1<kbin && index2<kbin )
				{
                    if (bb[index1][index2]>0)
                    {
				deltak2[no][0]=deltak2[no][0]*window[index1][index2]/bb[index1][index2];
				deltak2[no][1]=deltak2[no][1]*window[index1][index2]/bb[index1][index2];
					}
                else
				{
					deltak2[no][0]=0.;
					deltak2[no][1]=0.;
				}

				//deconvolve the filter in N-body simulation
				if (kx==0) sincx=1.; else sincx=sin(PI*kx/nc)/(PI*kx/nc);
				if (ky==0) sincy=1.; else sincy=sin(PI*ky/nc)/(PI*ky/nc);
				if (kz==0) sincz=1.; else sincz=sin(PI*kz/nc)/(PI*kz/nc);

                rev_tophat=1./sincx/sincy/sincz;

                deltak1[no][0]*=rev_tophat;
                deltak1[no][1]*=rev_tophat;


				pkk=pow(deltak2[no][0],2)+pow(deltak2[no][1],2);
                pdd=pow(deltak1[no][0],2)+pow(deltak1[no][1],2);
                pdk=deltak1[no][0]*deltak2[no][0]+deltak1[no][1]*deltak2[no][1];               
                    pk[index1][index2]+=pkk;
                    pd[index1][index2]+=pdd;
                    pkd[index1][index2]+=pdk;
                    count[index1][index2]++;
                    kcount[index1][index2]+=window[index1][index2];//weight average
                    pcount[index1][index2]+=pow(window[index1][index2],2);//weight sum for p_kappa
                }
				else
				{
					deltak2[no][0]=0.;
					deltak2[no][1]=0.;
				}
//=======================================================
              
                no++;
            }
        }
    }
//=============================================================
cout<<"save to:  "<<plotPath2<<endl;

    FILE *out1=fopen(plotPath2,"w");
	
    for (int i=0;i<kbin;i++)
    {
        for (int j=0;j<kbin;j++)
        {
        if (kcount[i][j]!=0)
        {
			renorm_weigh+=kcount[i][j];
            pk[i][j]/=pcount[i][j];
            pd[i][j]/=count[i][j];
            pkd[i][j]/=kcount[i][j];
            if (isnan(pk[i][j]*pd[i][j]*pkd[i][j])!=0) return -1;

        }
        double kout1=pow(10,i*kgap)*dk;
        double kout2=pow(10,j*kgap)*dk;
    fprintf(out1,"%e  %e %e  %e  %e  %e %e %e \n",kout1,kout2,pd[i][j]*renorm,pk[i][j]*renorm,pkd[i][j]*renorm,count[i][j],kcount[i][j],pcount[i][j]);
    }
}
    fclose(out1);

cout<<"renorm-weigh: "<<renorm_weigh<<endl;
//=======================================
cout<<"fftw back"<<endl;
fftw_execute(planback);
//==============================================
//output

cout<<"save to:  "<<outPath<<endl;

    out1=fopen(outPath,"wb");

	float outemp;
    double tot=double(ntot);
    for(int i=0;i<ntot;i++)
    {
		//deltars[i]/=renorm_weigh;
		deltars[i]/=tot;
        outemp=float(deltars[i]+1.);//to keep in the same form as original field
        fwrite(&outemp,4,1,out1);
    }
    fclose(out1);
//=====================================
cout<<"end powerspectrum"<<endl;

	fftw_free(deltars);
    fftw_free(deltak1);
    fftw_free(deltak2);
return 0;
}
  
