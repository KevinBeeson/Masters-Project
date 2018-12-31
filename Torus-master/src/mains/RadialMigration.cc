/***************************************************************************//**
\file testfalPot.cc
\brief Agian, should replace with something general.


*                                                                              *
*  testfalPot.cc                                                               *
*                                                                              *
*  C++ code written by Walter Dehnen, 1995-96,                                 *
*                      Paul McMillan, 2007-08,                                 *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/

#include <fstream>
#include <iostream>
#include "falPot.h"
#include <cmath>
#include <vector>

using std::cout;
using std::cerr;
using std::ifstream;
using namespace std;

int main(int argc,char *argv[])
{
    ifstream file;
    int    iso=1;
    Frequencies KNO;
    double PhiSperical;
    double currentPoint=1;
    double saved, step, rot,total,initialPhi,angPhi ;
    long double zsmall1,diffPer=1e+10,zsmall2,zlarge1,zlarge2,rsmall1,rsmall2,Rsmall1,Rsmall2,rlarge1,rlarge2,Rlarge1,Rlarge2;
    vector<double> q(3),p(3),Dv(3),Dp(3);

    cout<<"How many rotations do you want?";
    cin>>rot;
    //Only changes the resolution of the intergrator not the differenciator
    cout<<"whats the resolution that you want?";
    cin>>step;
    cout<<"How many saved points do you want?";
    cin>>saved;

    //sees if you have enough input potentials
    if(argc !=3) {
      cerr << "This test needs an input potential and innitial arguements, e.g. pot/DB97Mod1.Tpot and  starting.txt\n";
      exit(1);
    }

    cout<< "Reading input potential and starting points described in file " << argv[1] << '\n';
    file.open(argv[1]);
    GalaxyPotential Phi(file);
    //Phi(file);
    file.close();

    //Opens the starting coordinates and momenta with the structure shown bellow
    // r phi q[2]
    // pr pphi pz
    file.open(argv[2]);
    file >>q[0] >> q[1]>> q[2] >>p[0]>>p[1]>>p[2];
    file.close();

    //changes the rotattions to radians, the program will save depending on the change of phi saved will determine after how much change in phi it will save a point
    rot=2*3.1452*rot;
    saved=(rot-q[1])/saved;

    //This commented out section of code is to look at the Phi speed to be in circular orbit in that point
    // rsmall1= q[0]-q[0]/diffPer-1/diffPer;
    // rlarge1=q[0]+q[0]/diffPer+1/diffPer;
    // Rlarge1=pow(rlarge1*rlarge1+q[2]*q[2],0.5);
    // Rsmall1=pow(rsmall1*rsmall1+q[2]*q[2],0.5);
    // rsmall2= q[0]-2*(q[0]/diffPer+1/diffPer);
    // rlarge2=q[0]+2*(q[0]/diffPer+1/diffPer);
    // Rlarge2=pow(rlarge2*rlarge2+q[2]*q[2],0.5);
    // Rsmall2=pow(rsmall2*rsmall2+q[2]*q[2],0.5);
    // Dp[0]=(-1*Phi(Rlarge2,q[2])+8*Phi(Rlarge1,q[2])-8*Phi(Rsmall1,q[2])+Phi(Rsmall2,q[2]))*diffPer/(12*(q[0]+1));
    //
    // zsmall1= q[2] -q[2]/diffPer-1/diffPer;
    // zlarge1=q[2]+q[2]/diffPer+1/diffPer;
    // Rlarge1=pow(q[0]*q[0]+zlarge1*zlarge1,0.5);
    // Rsmall1=pow(q[0]*q[0]+zsmall1*zsmall1,0.5);
    // zsmall2= q[2] -2*(q[2]/diffPer+1/diffPer);
    // zlarge2=q[2]+2*(q[2]/diffPer+1/diffPer);
    // Rlarge2=pow(q[0]*q[0]+zlarge2*zlarge2,0.5);
    // Rsmall2=pow(q[0]*q[0]+zsmall2*zsmall2,0.5);
    // Dp[2]=(-Phi(Rlarge2,zlarge2)+8*Phi(Rlarge1,zlarge1)-8*Phi(Rsmall1,zsmall1)+Phi(Rsmall2,zlarge2))*diffPer/(12*(q[2]+1));
    // PhiSperical= pow(Dp[0]/q[0],0.5);
    // cout<<setprecision(16)<<"PhiSperical"<<PhiSperical<<"/n";


    //Angular momentum and the initial phi.
    initialPhi=q[1];
    angPhi=q[0]*q[0]*p[1];


    //File used to save the elements
    ofstream myfile;
    myfile.open("R.txt");
    myfile << "r \t pr \t phi \t pphi \t z \t pz \t t \n";

    //t is to take account of time
    int t=0;
    //programs stops after so many rotations
    while (q[1]-initialPhi<rot){

      t+=1;
      //  Gets the derivative at that point in space using Runge-Kutta 4
      rsmall1= q[0]-q[0]/diffPer-1/diffPer;
      rlarge1=q[0]+q[0]/diffPer+1/diffPer;
      Rlarge1=pow(rlarge1*rlarge1+q[2]*q[2],0.5);
      Rsmall1=pow(rsmall1*rsmall1+q[2]*q[2],0.5);
      rsmall2= q[0]-2*(q[0]/diffPer+1/diffPer);
      rlarge2=q[0]+2*(q[0]/diffPer+1/diffPer);
      Rlarge2=pow(rlarge2*rlarge2+q[2]*q[2],0.5);
      Rsmall2=pow(rsmall2*rsmall2+q[2]*q[2],0.5);
      Dp[0]=(-1*Phi(Rlarge2,q[2])+8*Phi(Rlarge1,q[2])-8*Phi(Rsmall1,q[2])+Phi(Rsmall2,q[2]))*diffPer/(12*(q[0]+1));

      zsmall1= q[2] -q[2]/diffPer-1/diffPer;
      zlarge1=q[2]+q[2]/diffPer+1/diffPer;
      Rlarge1=pow(q[0]*q[0]+zlarge1*zlarge1,0.5);
      Rsmall1=pow(q[0]*q[0]+zsmall1*zsmall1,0.5);
      zsmall2= q[2] -2*(q[2]/diffPer+1/diffPer);
      zlarge2=q[2]+2*(q[2]/diffPer+1/diffPer);
      Rlarge2=pow(q[0]*q[0]+zlarge2*zlarge2,0.5);
      Rsmall2=pow(q[0]*q[0]+zsmall2*zsmall2,0.5);
      Dp[2]=(-1*Phi(Rlarge2,zlarge2)+8*Phi(Rlarge1,zlarge1)-8*Phi(Rsmall1,zsmall1)+Phi(Rsmall2,zlarge2))*diffPer/(12*(q[2]+1));

      // The change in p using the E-L equation
      Dv[0]=-Dp[0]+p[1]*p[1]*q[0];
      Dv[1]=0;
      Dv[2]=-Dp[2];

      //registers the position history and momenta history and prints it in a text file  after moving a certain Phi.
      if ((q[1]-initialPhi)/saved>currentPoint){
        currentPoint+=1;

        // Moves a certain amount for the momentum
        p[0] = p[0] + 1 / step* Dv[0];
        q[0] = q[0] + 1 / step*p[0];
        myfile << q[0] << ","<< p[0] << ",";

        //as the angular momentum is constant I will base p[1] i.e. the dphi/dt on it
        p[1] = angPhi/(q[0]*q[0]);
        q[1] = q[1] + 1 / step*p[1];
        myfile << q[1] << ","<< p[1] << ",";
        p[2] = p[2] + 1 / step* Dv[2];
        q[2] = q[2] + 1 / step*p[2];
        myfile << q[2] << ","<< p[2] << ",";
        // cout << r[i]<< "\t"<< p[i]<<"\n";
        myfile <<t/step << "\n";
      }
      else{
        //same as above but doesnt save it
        p[0] = p[0] + 1 / step* Dv[0];
        q[0] = q[0] + 1 / step*p[0];
        p[1] = angPhi/(q[0]*q[0]);
        q[1] = q[1] + 1 / step*p[1];
        p[2] = p[2] + 1 / step* Dv[2];
        q[2] = q[2] + 1 / step*p[2];
        }

    }
    myfile.close();
    return 0;
}
