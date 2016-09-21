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

int momenorigin()
{
    cout<<"start momen origin"<<endl;
    /*
char inPath1[80]={"/home/zhm/dongzi/tidez1/3dtidenew/00smooth1.25_z1_rpar15.bin"};
char inPath2[80]={"/home/zhm/dongzi/tidez1/3dtidenew/00v_3d_z1_rpar15_filtered.bin"};
char outPath[80]={"/home/zhm/dongzi/tidez1/3dtidenew/00momen_2d_z1_rpar15.bin"};
*/
char inPath1[80]={"/home/zhm/tides00/1.000den00.bin"};
char inPath2[80]={"/home/zhm/tides00/1.000velz00.bin"};

char outPath[80]={"/home/zhm/dongzi/tidez1/3dtidenew/00momen_2d_z1_origin.bin"};

//variables and fftw setting
//=====================================
    double *deltar1;       //original field
    double *deltar2;      //reconstructed velocity field

    int nn=int(nc);
    int ntot=nn*nn*nn;
    deltar1 = (double *) fftw_malloc(ntot * sizeof(double));
    deltar2 = (double *) fftw_malloc(ntot * sizeof(double));

//=================================
//matrix for 2d field
	double v2d[1024][1024];
	for(int i=0;i<nn;i++)
	{
		for(int j=0;j<nn;j++)
		{
			v2d[i][j]=0.;
		}
	}	
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

//=====================================
cout<<"calculate 2D velocity field"<<endl;
double winsz=1;
double momen;
int no=0;
   for(int i=0;i<nn;i++)
    {
        for(int j=0;j<nn;j++)
        {
            for(int k=0;k<nn;k++)
            {
              if(no<ntot)
              {
				  momen=(1.+deltar1[no])*deltar2[no];
				  winsz=1.;
				  v2d[j][k]+=(momen*winsz);
                //  if(momen==0)cout<<"i,j,k,momen"<<i<<" "<<j<<" "<<k<<endl;
             }
                no++;
            }
        }
    }
//======================================
cout<<"save to:  "<<outPath<<endl;
    FILE *out2=fopen(outPath,"wb");
    float outemp2;
    for(int j=0;j<nn;j++)
	{
		for (int k=0;k<nn;k++)
		{
        if (isnan(v2d[j][k])!=0)
        {   
            cout<<"nan appears, no: "<<j<<", "<<k<<endl;
            return -1;
        }
        outemp2=float(v2d[j][k]);
        fwrite(&outemp2,4,1,out2);
		}
		}
       fclose(out2);

//========================================
    fftw_free(deltar1);
    fftw_free(deltar2);  
    cout<<"end momentum field calculation"<<endl;
return 0;
}
  
