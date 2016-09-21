#include<vector>
#include<algorithm>
#include<utility>
#include<map>

#include<stdio.h>
#include<iostream>
using namespace std;
#include<math.h>
#include<stdlib.h>
#include<fftw3.h>

#include<time.h>
//for sort
namespace Solution3
{
  template<class T>
  struct CompareDeref
  {
    bool operator()( const T& a, const T& b ) const
      { return *a < *b; }
  };


  template<class T, class U>
  struct Pair2nd
  {
    const U& operator()( const std::pair<T,U>& a ) const
      { return a.second; }
  };


  template<class IterIn, class IterOut>
  void sort_idxtbl( IterIn first, IterIn last, IterOut out )
  {
    std::multimap<IterIn, int, CompareDeref<IterIn> > v;
    for( int i=0; first != last; ++i, ++first )
      v.insert( std::make_pair( first, i ) );
    std::transform( v.begin(), v.end(), out,
                    Pair2nd<IterIn const,int>() );
  }
}

//=========================================================
int main()
{
   clock_t t=clock();
char inPath[]={"/home/zhm/dongzi/tidez1/rpar15/00smooth1.25_z1_rpar15.bin"};

char outPath[]={"/home/zhm/dongzi/tidez1/gauss_3d_r15kc0.5/00gausmooth1.25_z1_rpar15.bin"};
char outPath2[]={"/home/zhm/dongzi/tidez1/gauss_3d_r15kc0.5/00slice_smooth1.25_z1_rpar15.bin"};
char outPath_sort[]={"/home/zhm/dongzi/tidez1/gauss_3d_r15kc0.5/00sequence_smooth1.25_z1_rpar15.bin"};
//variables and fftw setting
//=================================
    double *deltar;       //input field
    int nn=1024;
    int ntot=nn*nn*nn;
    deltar = (double *) fftw_malloc(ntot * sizeof(double));
//=================================
/*test
for (int i=0;i<ntot;i++)
{
    deltar[i]=30-i;
}
*/
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

        deltar[i]=intemp;

            }
    fclose(in);
//========================================
cout<<"sort"<<endl;
std::vector<int>id(ntot);
Solution3::sort_idxtbl(deltar,deltar+ntot,id.begin());
//===============================================
//gaussianize
double sigma=2.,pmin=-4.,pmax=4.,pbin=8000.,prenorm=0.;
double ppointer=pmin,pgap=(pmax-pmin)/pbin, sigma2=2*sigma*sigma;
double bin[8000];
for(int i=0;i<pbin;i++)
{
    bin[i]=exp(-ppointer*ppointer/sigma2);
    ppointer+=pgap;
    prenorm+=bin[i];
}
ppointer=pmin;
double bincount;
int no=0;
for(int i=0;i<pbin;i++)
{
    bincount=round(bin[i]/prenorm*ntot);
    while(bincount>0 && no<ntot)
    {
        deltar[id[no]]=ppointer;
        bincount--;
        no++;
    }
    ppointer+=pgap;
}
//=========================================

cout<<"save to:  "<<outPath<<endl;
cout<<"save to:  "<<outPath2<<endl;

    FILE *out=fopen(outPath,"wb");
	FILE *out2=fopen(outPath2,"wb");
    FILE *out3=fopen(outPath_sort,"wb");
    float outemp;
    int outemp2;
	no=0;
    for(int i=0;i<nn;i++)
    {   
		for (int j=0;j<nn;j++)
		{
			for (int k=0;k<nn;k++)
			{
			  if (isnan(deltar[no])!=0)
			  {   
				 cout<<"nan appears, no: "<<i<<endl;
				 return -1;
				 }
				 outemp=float(deltar[no]+1.);
				 fwrite(&outemp,4,1,out);
				
				 //for slice iamge
				 if (i==512) fwrite(&outemp,4,1,out2);

                 //for sequence
                 outemp2=int(id[no]);
                 fwrite(&outemp2,sizeof(int),1,out3);

                 no++;
				 }
		}
	}
				 fclose(out);
				 fclose(out2);
//========================================
/*test
  for( int i=0; i<ntot; ++i )
  std::cout << "i=" << i
            << ", id[i]=" << id[i]
            << ", deltar[id[i]]=" << deltar[id[i]]
            <<", deltar[i]="<<deltar[i]
            << std::endl;
  std::cout << "#################" << std::endl;
*/
//=====================================
    fftw_free(deltar);
    t=clock()-t;
    cout<<"time = "<<float(t)/CLOCKS_PER_SEC/60.<<"min"<<endl;
return 0;
}
    

