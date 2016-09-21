#include<stdio.h>
#include<iostream>
using namespace std;
#include"smooth_fs.h"
#include"recon.h"
#include"powerspec.h"
#include"filter_neg.h"
#include"vfield.h"
#include"filterv.h"
#include"powerv.h"
#include"momenfield.h"
#include"powermomen.h"
#include"momenorigin.h"
#include"filter3d.h"
#include"shear.h"

//Global variable for whole project:
double PI=3.1415926;
double nc=1024.;    //number of cells per dimension
double nck=(nc/2.+1); //number of cells for fftw output
double box=1200.;    //N-body simulation box size, 1.2Gpc/h
double dk=2.*PI/box;//fundamental frequency
//========================================================

int main()
{
    int check=0;
   // check=smooth_fs();
//  check=shear(); 
  //  check=filter3d();
    /*
    if (check!=0) cout<<"error smooth"<<endl;
    check=recon();
   if (check!=0) cout<<"error k3d_noisy"<<endl;
        check=filter_neg();
    if (check!=0) cout<<"error filter_neg"<<endl;
  */ 
    /*
    check=powerspec();
    if (check!=0) cout<<"error powerspec"<<endl;
     check=vfield();
    if (check!=0) cout<<"error vfield"<<endl; 
 check=filterv();
    if (check!=0) cout<<"error filterv"<<endl;

  check=powerv();
    if (check!=0) cout<<"error powerv"<<endl;
   */ 
check=momenfield();
    if (check!=0) cout<<"error momenfield"<<endl;
 check=powermomen();
    if (check!=0) cout<<"error powermomen"<<endl;

// check=momenorigin();
  //  if (check!=0) cout<<"error momenorigin"<<endl;

cout<<"end whole program"<<endl;
    return 0;
}

