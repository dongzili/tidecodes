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

int momenfield()
{
    cout<<"start momenfield"<<endl;
char inPath1[]={"/home/zhm/tidesData/tides00/2.000den00.bin"};
//char inPath1[]={"/home/zhm/tidesData/tides00/2.000den00.bin"};
char inPath2[]={"/project/zhm/ksz/z1/tidepart/00v_3d_clean_z2_rpar10_kc0.5.bin"};
char outPath[]={"/project/zhm/ksz/z1/photoz/00tide_photoz_momen_2d_z2_rpar10_kc0.5_dchi200.bin"};

    /*
char inPath2[]={"/home/zhm/tidesData/tides00/2.000velz00.bin"};

char outPath[]={"/project/zhm/ksz/z1/photoz/00momen_2d_z1_origin.bin"};

*/
//variables and fftw setting
//=====================================
    double *deltar1;       //original field
    double *deltar2;      //reconstructed velocity field

    int nn=int(nc);
    int ntot=nn*nn*nn;
    deltar1 = (double *) fftw_malloc(ntot * sizeof(double));
    deltar2 = (double *) fftw_malloc(ntot * sizeof(double));

//========================================
//reading data to deltar    
//=====================================
    
    cout<<"reading  "<<inPath1<<endl;
    cout<<"reading  "<<inPath2<<endl;

    float intemp,check;
    FILE *in1=fopen(inPath1,"rb");
    FILE *in2=fopen(inPath2,"rb");
    for(int i=0;i<ntot;i++)
    {
		//reading original density field
        check=fread(&intemp,4,1,in1);
        if (check!=1)
        {cout<<"reading problem. check="<<check<<endl;
         return -1;
        }
        deltar1[i]=intemp-1.;

		//reading velocity field
		check=fread(&intemp,4,1,in2);
        if (check!=1)
        {cout<<"reading problem. check="<<check<<endl;
         return -1;
        }
        deltar2[i]=intemp;


            }
    fclose(in1);
	fclose(in2);

//=================================
//bin it into photoz

    double binz=6;
    double zslice=200;
    double *den2d, *v2d;       //original field

    int slicetot=nn*nn*binz;
    den2d = (double *) fftw_malloc(slicetot * sizeof(double));
    v2d = (double *) fftw_malloc(slicetot * sizeof(double));
    double count=ntot/binz;
    int no=0;
	for(int i=0;i<nn;i++)
	{
		for(int j=0;j<nn;j++)
		{
            for (int k=0;k<binz;k++)
            {
			    den2d[no]=0.;
                v2d[no]=0.;
                no++;
            }
		}
	}	
    //den2d[slice][j][k]
//================================
no=0;
int noslice=0;
   for(int i=0;i<nn;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nn;k++)
            {
              if(no<ntot)
              {
                  int slice=int(floor(i/zslice));
                  noslice=(slice*nn+j)*nn+k;
				  den2d[noslice]+=(1.+deltar1[no]);
                  v2d[noslice]+=deltar2[no];
                //  if(momen==0)cout<<"i,j,k,momen"<<i<<" "<<j<<" "<<k<<endl;
             }
                no++;
            }
        }
    }
    fftw_free(deltar1);
    fftw_free(deltar2);  
//=====================================
//matrix for 2d field
	double p2d[1024][1024];
	for(int i=0;i<nn;i++)
	{
		for(int j=0;j<nn;j++)
		{
			p2d[i][j]=0.;
		}
	}	
//=================================
cout<<"calculate 2D velocity field"<<endl;
double winsz=1;
double momen;
noslice=0;
   for(int i=0;i<binz;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nn;k++)
            {
				  momen=den2d[noslice]/count*v2d[noslice];
				  winsz=1.;
				  p2d[j][k]+=(momen*winsz);
                  noslice++;
                //  if(momen==0)cout<<"i,j,k,momen"<<i<<" "<<j<<" "<<k<<endl;
            }
        }
    }
   fftw_free(den2d);
   fftw_free(v2d);
//======================================
cout<<"save to:  "<<outPath<<endl;
    FILE *out2=fopen(outPath,"wb");
    float outemp2;
    for(int j=0;j<nn;j++)
	{
		for (int k=0;k<nn;k++)
		{
        if (isnan(p2d[j][k])!=0)
        {   
            cout<<"nan appears, no: "<<j<<", "<<k<<endl;
            return -1;
        }
        outemp2=float(p2d[j][k]);
        fwrite(&outemp2,4,1,out2);
		}
		}
       fclose(out2);
//========================================
    cout<<"end momentum field calculation"<<endl;
return 0;
}
  
