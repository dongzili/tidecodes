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

double PI=3.1415926;
double nc=1024.;    //number of cells per dimension
double nck=(nc/2.+1); //number of cells for fftw output
double box=1200.;    //N-body simulation box size, 1.2Gpc/h
double dk=2.*PI/box;//fundamental frequency

int main()
{
cout<<"start powerv"<<endl;

/*
char plotPath2[]={"/home/zhm/dongzi/tidez1/3dtide_data/00den_2dpowerspec_rpar15.dat"};
char inPath1[]={"/home/zhm/tidesData/tides00/2.000den00.bin"};
char inPath2[]={"/home/zhm/dongzi/tidez1/3dtide_data/00k3d_noisy_z1_rpar15.bin"};
char inPath_fil[]={"/home/zhm/dongzi/tidez1/3dtide_data/10den_noisebias_z1_rpar15.dat"};
*/
char plotPath2[]={"/home/zhm/dongzi/tidez1/3dtide_data/00v_2dpowerspec_rpar15.dat"};
char inPath1[]={"/home/zhm/tidesData/tides00/2.000velz00.bin"};
char inPath2[]={"/home/zhm/dongzi/tidez1/3dtide_data/00v_3d_z1_rpar15.bin"};
char inPath_fil[]={"/home/zhm/dongzi/tidez1/3dtide_data/10v_noisebias_z2_rpar10_kc0.5.dat"};
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
//=================================================
//filter grouping information
	int kbin_fil=10;
	double kstart_fil,kend_fil,kgap_fil;
	kstart_fil=1.;
	kend_fil=double(nnc);
	kgap_fil=log10(kend_fil/kstart_fil)/kbin_fil;
//================================================

    double pk[10][10],pd[10][10],pkd[10][10];
    double count[10][10],kcount[10][10],pcount[10][10];//one for kappa, one for power_kappa renorm
    
    for (int i=0;i<kbin_fil;i++)
	{
		for (int j=0;j<kbin_fil;j++)
    {
        pk[i][j]=0.;
        pd[i][j]=0.;
        pkd[i][j]=0.;
        count[i][j]=0.;
        kcount[i][j]=0.;
        pcount[i][j]=0.;
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
//======================================================================
//reading data: k, P1s    
cout<<"reading   "<<inPath_fil<<endl;
    float intemp,check;
    int nbin=10;//number of linear power spectrum
    double bb[10][10],window[10][10];
    double ktemp,ptemp;
FILE *in=fopen(inPath_fil,"r");
    for(int i=0;i<nbin;i++)
	{
		for (int j=0;j<nbin;j++)
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
    int index1,index2; //for grouping

    double pkk,pdd,pdk;//,bk;

    double kx,ky,kz,korth,kpar;
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

				korth=sqrt(kx*kx+ky*ky);
				kpar=fabs(kz);

                if (korth==0 || kpar==0) //korth==0,infinite noise, kpar=0 don't contribute to vz
                {
					deltak2[no][0]=0.;
					deltak2[no][1]=0.;
                    no++;
                    continue;
                }

                korth=log10(korth);
				kpar=log10(kpar);
                
				index1=int(floor(korth/kgap_fil));
				index2=int(floor(kpar/kgap_fil));
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
				if (index1<kbin_fil && index2<kbin_fil )
				{
                    if (bb[index1][index2]!=0)
                    {
                        if(bb[index1][index2]<0) bb[index1][index2]=1.;
				deltak2[no][0]=deltak2[no][0]*window[index1][index2]/bb[index1][index2];
				deltak2[no][1]=deltak2[no][1]*window[index1][index2]/bb[index1][index2];
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
					}
                
//=======================================================
              
                no++;
            }
        }
    }
//=============================================================
cout<<"save to:  "<<plotPath2<<endl;

    FILE *out1=fopen(plotPath2,"w");
	double kout1,kout2;

    for (int i=0;i<kbin_fil;i++)
	{
		for (int j=0;j<kbin_fil;j++)
    {
        if (kcount[i][j]!=0)
        {
			renorm_weigh+=kcount[i][j];
            pk[i][j]/=pcount[i][j];
            pd[i][j]/=count[i][j];
            pkd[i][j]/=kcount[i][j];
if (isnan(pk[i][j]*pd[i][j]*pkd[i][j])!=0) return -1;

        }
        kout1=pow(10,(i+0.5)*kgap_fil)*dk;
		kout2=pow(10,(j+0.5)*kgap_fil)*dk;
    fprintf(out1,"%e  %e  %e  %e  %e %e %e %e \n",kout1,kout2,pd[i][j]*renorm,pk[i][j]*renorm,pkd[i][j]*renorm,count[i][j],kcount[i][j],pcount[i][j]);
    }
	}
    fclose(out1);

cout<<"renorm-weigh: "<<renorm_weigh<<endl;
//=======================================

//=====================================
cout<<"end powerv"<<endl;
    fftw_free(deltak1);
    fftw_free(deltak2);
    fftw_free(deltar1);
    fftw_free(deltar2);
return 0;
}
  
