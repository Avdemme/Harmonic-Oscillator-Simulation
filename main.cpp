#include "stdio.h"
#include "stdlib.h"
#include <fstream>
#include <iostream>
#include "math.h"
#include "string.h"

#define ROW 35
#define COL 91

double rand_uniform ( );
void ClearTheScreen ( );
using namespace std;

int main ( int argc, char** argv ) {
  
  double r1[2], v1[2], R1[2], V1[2], vNorm;
  double dt, max; double a0[2], a1[2], A0[2], A1[2];
  double r2[2], v2[2], a2[2], R2[2], V2[2], A2[2];
  bool grid[ROW][COL]; int x, y, X, Y;
  
  cout << "time_step [s] = ";
  cin >> dt;
  cout << "stop-time [s] = ";
  cin >> max;
  

  double equil = 300 ;
  double G_Newton = 1e-6;
 
  double powerLaw = 1.;
  double M = 0.2 ;
  double m = 1. ;
  double k = G_Newton * M * m;
  double r0[2] = { 300 , 0. } ;
  double v0[2] = { 0.6 , 0. };
  double R0[2] = { 0. , 0. };

  double V0[2] = {-(m/M)*v0[0],-(m/M)*v0[1]};
  

  
  double radius2 = sqrt(pow(R0[0]-r0[0],2.)+pow(R0[1]-r0[1],2.)) - equil;
  double sign = fabs(radius2)/radius2; if ( radius2 == 0. ) sign = 1.;
  radius2 = pow(radius2,2.); double min = radius2;
  double alpha = pow(radius2,(-powerLaw/2.)+0.5)*sign; //directional b/c of spring case

  a0[0] = (k / m ) * (R0[0]-r0[0])/alpha;
 


  A0[0] = (k / M ) * (r0[0]-R0[0])/alpha;
  

   
 
  
  a0[1] = (k / m ) * (R0[1]-r0[1])/alpha;
 
  A0[1] = (k / M ) * (r0[1]-R0[1])/alpha;
  
  double time = 0.;
  while ( time < max ) {
    
    for ( int i = 0; i < ROW; i++ ) {
      for ( int j = 0; j < COL; j++ ) {
	grid[i][j] = false;
      }
    } //end of the initial grid setup
    
    if ( !time ) {
      r1[0] = r0[0] + 0.5 * v0[0] * dt + 0.25 * a0[0] * pow(dt,2.);
      r1[1] = r0[1] + 0.5 * v0[1] * dt + 0.25 * a0[1] * pow(dt,2.);
      R1[0] = R0[0] + 0.5 * V0[0] * dt + 0.25 * A0[0] * pow(dt,2.);
      R1[1] = R0[1] + 0.5 * V0[1] * dt + 0.25 * A0[1] * pow(dt,2.);
    }
    else {
      r1[0] = r0[0] + 0.5 * v0[0] * dt;
      r1[1] = r0[1] + 0.5 * v0[1] * dt;
      R1[0] = R0[0] + 0.5 * V0[0] * dt;
      R1[1] = R0[1] + 0.5 * V0[1] * dt;
    }
    
    radius2 = sqrt(pow(R0[0]-r0[0],2.)+pow(R0[1]-r0[1],2.)) - equil;
    sign = fabs(radius2)/radius2; if ( radius2 == 0. ) sign = 1.;
    radius2 = pow(radius2,2.);
    alpha = pow(radius2,(-powerLaw/2.)+0.5)*sign;
    if ( !alpha || pow(R1[0]-r1[0],2.)+pow(R1[1]-r1[1],2.) < 1e-2*min ) break; // CRASH!

    
    a1[0] = ( k / m ) * (R1[0]-r1[0])/alpha;
    a1[1] = ( k / m ) * (R1[1]-r1[1])/alpha;
    A1[0] = ( k / M ) * (r1[0]-R1[0])/alpha;
    A1[1] = ( k / M ) * (r1[1]-R1[1])/alpha;
    
    v2[0] = v0[0] + a1[0] * dt;
    v2[1] = v0[1] + a1[1] * dt;
    V2[0] = V0[0] + A1[0] * dt;
    V2[1] = V0[1] + A1[1] * dt;
    
    r2[0] = r1[0] + 0.5 * v2[0] * dt;
    r2[1] = r1[1] + 0.5 * v2[1] * dt;
    R2[0] = R1[0] + 0.5 * V2[0] * dt;
    R2[1] = R1[1] + 0.5 * V2[1] * dt;
        
    X = int(floor((R0[0]/100)+0.5)) + int(floor(double(COL)/2.));
    Y =-int(floor((R0[1]/200)+0.5)) + int(floor(double(ROW)/2.));
    x = int(floor((r0[0]/100)+0.5)) + int(floor(double(COL)/2.));
    y =-int(floor((r0[1]/200)+0.5)) + int(floor(double(ROW)/2.));
    if ( x < COL && y < ROW && x >= 0 && y >= 0 )
      grid[y][x] = true;
 
    ClearTheScreen();
    for ( int i = 0; i < ROW; i++ ) {
      for ( int j = 0; j < COL; j++ ) {
	if ( grid[i][j] ) cout << "m";
	else if ( j == X && i == Y ) cout << "M";
	else cout << " ";
      }
      cout << endl;
    } //draw the new grid to the screen by row
    
    time += dt;
    r0[0] = r2[0]; r0[1] = r2[1];
    v0[0] = v2[0]; v0[1] = v2[1];
    R0[0] = R2[0]; R0[1] = R2[1];
    V0[0] = V2[0]; V0[1] = V2[1];
    
   
      
    
    if ( (r0[0] - R0[0]) < 5 ){
      V0[0] = - V0[0];
      v0[0] = - v0[0];
    }
      

  }
  
  return -1;
  
}

double rand_uniform ( ) {
  
  return (double)rand() / (double)RAND_MAX;
  
}

void ClearTheScreen ( ) {
  
  printf("\033[%d;%dH",1,1);
  return;
  
} //hopefully this func works for all OSes!!
