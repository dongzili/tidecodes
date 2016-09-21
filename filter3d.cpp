/* calculate the power spectrums, P_k3d, P_delta, P_k3d,delta
    
       unit:
        x=Mpc/h
        
    OUTPUT:
        delta(k): smooth_deltak.dat      */
//=============================================================
//DOngzi Li 1015.Dec, follow Hongming Zhu's program


#include<stdio.h>
#include<iostream>
using namespace std;
#include<math.h>
#include<stdlib.h>
#include<fftw3.h>
#include"tide.h"


int filter3d()
{


//char outPath1[]={"/home/zhm/dongzi/tide1.17/tide0/k3d_clean.bin"};
//char plotPath1[]={"/home/zhm/dongzi/tide1.17/tide0/noise.dat"};
char backupPath[]={"/project/zhm/ksz/z1/tidekperp/3d00noisefilter_rpar15_kc0.6_bin10.dat"};
char inPath1[]={"/home/zhm/tidesData/tides00/1.000den00.bin"};
char inPath2[]={"/project/zhm/ksz/z1/tidekperp/00k3d_noisy_z1_rpar15.bin"};
//variables and fftw setting
//=================================
    int nn=int(nc);
    int nnc=int(nck);
    int ntot=nn*nn*nn;
    int nktot=nn*nn*nnc;
  //  double pn[nktot];   
    //=========================to group power spectra
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
    
    double pk[10][10][10],pd[10][10][10],pkd[10][10][10],bb[10][10][10],pn[10][10][10],window[10][10][10];
    double kcount[10][10][10];
    
    for (int i=0;i<kbin;i++)
		for (int j=0;j<kbin;j++)
            for (int k=0;k<kbin;k++)
		{
    {
        {
        pk[i][j][k]=0.;
        pd[i][j][k]=0.;
        pkd[i][j][k]=0.;
        kcount[i][j][k]=0.;
		bb[i][j][k]=0.;
		pn[i][j][k]=0.;
		window[i][j][k]=0.;
    }
		}
        }

    double *deltar1,*deltar2;       //input field
    fftw_complex *deltak1, *deltak2; //Fourier field 
    fftw_plan plan1,plan2;       //fftw setting variable


    deltak1 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));
    deltak2 = (fftw_complex*) fftw_malloc(nktot * sizeof(fftw_complex));

    deltar1 = (double *) fftw_malloc(ntot * sizeof(double));
    deltar2 = (double *) fftw_malloc(ntot * sizeof(double));

    plan1 = fftw_plan_dft_r2c_3d(nn,nn,nn,deltar1,deltak1,FFTW_ESTIMATE);
    plan2 = fftw_plan_dft_r2c_3d(nn,nn,nn,deltar2,deltak2,FFTW_ESTIMATE);


//reading data to deltar    
//=====================================
    cout<<"reading  "<<inPath1<<endl;
     cout<<"reading  "<<inPath2<<endl;   
    float intemp,check;
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
		if (intemp>1e2) cout<<intemp<<"  "<<i<<endl;
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
	
cout<<"negative correlation"<<endl;
	int negco=0; //to count how many negative correlation

    int no=0;//for matrix index
    int index1,index2,index3; //for grouping

    double pkk,pdd,pdk;//,bk;
    double kx,ky,kz;

	double sincx,sincy,sincz,rev_tophat;

    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nnc;k++)
            {
				kx=double(k);
                ky=double(j<nnc?j:(j-nn));
                kz=double(i<nnc?i:(i-nn));
                ky=fabs(ky);
				kz=fabs(kz);

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

    //           bk=pdk/pdd;
    //           pn[no]=pkk-bk*bk*pdd;//noise

//==============================================================
                if (kx==0 && ky==0 && kz==0) 
                {
                    no++;
                    continue;
                }
                if (kx==0) kx++;
                if (ky==0) ky++;
                if (kz==0) kz++;

                   kx=log10(kx);
                   ky=log10(ky);
                   kz=log10(kz);

                index1=int(floor(kx/kgap));
                index2=int(floor(ky/kgap));
				index3=int(floor(kz/kgap));


                if (index1 < kbin && index2 < kbin && index3<kbin)// && pdk > 0)
                {
					if (pdk<=0)
					{
						negco++;
					}
                    pk[index1][index2][index3]+=pkk;
                    pd[index1][index2][index3]+=pdd;
                    pkd[index1][index2][index3]+=pdk;
                    kcount[index1][index2][index3]++;
                }
//=======================================================
              
                no++;
            }
        }
    }
//=============================================================
//calculate b(k//,k orth), Pn(k//,k orth)

	 for (int i=0;i<kbin;i++)
     {
		for (int j=0;j<kbin;j++)
		{
            for (int k=0;k<kbin;k++)
            {
				if (pd[i][j][k]>0)
				{
				bb[i][j][k]=pkd[i][j][k]/pd[i][j][k];
				pn[i][j][k]=pk[i][j][k]-pow(bb[i][j][k],2)*pd[i][j][k];

				window[i][j][k]=pd[i][j][k]/(pd[i][j][k]+pn[i][j][k]/pow(bb[i][j][k],2));
               // if (bb[i][j][k]<0) bb[i][j][k]=0.;// drop these mode
				}	 
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
            for(int k=0;k<kbin;k++)
            {
            
        if (kcount[i][j][k]!=0)
        {
            pn[i][j][k]/=kcount[i][j][k];
            pd[i][j][k]/=kcount[i][j][k];
            pkd[i][j][k]/=kcount[i][j][k];
            pk[i][j][k]/=kcount[i][j][k];
        }
        
	fprintf(out2,"%e  %e  %e  %e  %e  %e  %e \n",bb[i][j][k],window[i][j][k],pn[i][j][k]*renorm,pd[i][j][k]*renorm,pk[i][j][k]*renorm,pkd[i][j][k]*renorm,kcount[i][j][k]);
    }
	}
    }
	fclose(out2);


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
    fftw_free(deltak1);
    fftw_free(deltak2);
    fftw_free(deltar1);
    fftw_free(deltar2);
return 0;
}
  
