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


int filterv()
{


//char outPath1[]={"/home/zhm/dongzi/tide1.17/tide0/k3d_clean.bin"};
//char plotPath1[]={"/home/zhm/dongzi/tide1.17/tide0/noise.dat"};
char backupPath[]={"/project/zhm/ksz/z1/tidekperp/00v_noisefilter_z2_rpar10_kc0.5_l300.dat"};
char filterPath[]={"/project/zhm/ksz/z1/tidekperp/00v_noisebias_z2_rpar10_kc0.5_l300.dat"};
char inPath1[]={"/home/zhm/tidesData/tides00/2.000velz00.bin"};
char inPath2[]={"/project/zhm/ksz/z1/tidekperp/00v_3d_z2_rpar10_kc0.5_l300.bin"};
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
    fftw_free(deltar1);
    fftw_free(deltar2);
//========================================
	
cout<<"calculate vfilter"<<endl;
	int negco=0; //to count how many negative correlation

    int no=0;//for matrix index
    int index1,index2; //for grouping

    double pkk,pdd,pdk;//,bk;

    double kx,ky,kz,korth;
    for(int i=0;i<nn;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nnc;k++)
            {
                pkk=pow(deltak2[no][0],2)+pow(deltak2[no][1],2);
                pdd=pow(deltak1[no][0],2)+pow(deltak1[no][1],2);
                pdk=deltak1[no][0]*deltak2[no][0]+deltak1[no][1]*deltak2[no][1];               

    //           bk=pdk/pdd;
    //           pn[no]=pkk-bk*bk*pdd;//noise

//==============================================================
//group power spectra according to abs(k)
                kx=double(k);
                ky=double(j<nnc?j:(j-nn));
                kz=double(i<nnc?i:(i-nn));
                korth=sqrt(kx*kx+ky*ky);
				kz=fabs(kz);
                if (korth==0 && kz==0) 
                {
                    no++;
                    continue;
                }
                if (korth==0)korth++;
                if (kz==0)kz++;

                korth=log10(korth);
				kz=log10(kz);

                index1=int(floor(korth/kgap));
				index2=int(floor(kz/kgap));
                if (index1 < kbin && index2 < kbin)// && pdk > 0)
                {
					if (pdk<=0)
					{
						negco++;
					}
                    pk[index1][index2]+=pkk;
                    pd[index1][index2]+=pdd;
                    pkd[index1][index2]+=pdk;
                    kcount[index1][index2]++;
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
    FILE *out= fopen(filterPath,"w");
		fprintf(out2," negative b:  %i \n",negco);
		fprintf(out2," b, window, pn \n");

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
    fprintf(out,"%e %e \n",bb[i][j],window[i][j]);
    }
	}
	fclose(out2);
    fclose(out);


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
return 0;
}
  
