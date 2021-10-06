#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h> 
#include <cstring>
#include <algorithm> // for function min
#include <random>
#include <cmath>
#include <chrono>
using namespace std;

// ------- Problem data:
int n=10; // number of aircraft
int mode=6; // to indicate the scenario type, default=6 (pseudo random in 2D)
float vmin=400; float vmax=400; // random speeds will be generated in [vmin,vmax], default=[400,400]
float* x_0;  float* y_0;    float* z_0;// vectors of initial positions (x0,y0,z0)
float* x_t;  float* y_t;    float* z_t;// vectors of final position (only for drawing scenario in a .tex file)
float* hat_v; float* hat_theta; float* hat_phi; // vectors of nominal speed and heading angles (hat_teta wrt x, hat_phi wrt z)
float* vx; float* vy;   float* vz; // vectors of nominal components of velocity (vx,vy,vz)=hat_v(sin hat_phi*cos hat_theta,sin hat_phi*sin hat_theta,cos hat_phi)
char* outputFile=0;
float D=5; // safety distance
int HP; int VP;  // polyhedral problem: number of horizontal and vertical planes
const int hp=3; //  default number of horizontal planes
const int vp=2; //  default number of vertical planes

// ------- Other parameters:
//int random_seed=std::chrono::system_clock::now().time_since_epoch().count();;
int random_seed=14;
bool verbose=0;
bool is3D=false;

// --------------------------------------------------------------------------
// -----                  AUXILIARY FUNCTIONS                        -------
// --------------------------------------------------------------------------

// 1. randomFloat(float a, float b): returns a random float in the interval [a,b]
float randomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}
// 2. randomInt(int a, int b): returns a random integer in the interval [a,b]
int randomInt(int a, int b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    int r = round(random * diff);
    return a + r;
}
// 3. areSepatered(int i, int j): returns true if the distance between the initial 
//       positions of i and j is greater than or equal to D
bool areSepatered(int i, int j){
	float xr0=x_0[i]-x_0[j];
	float yr0=y_0[i]-y_0[j];
    float zr0=z_0[i]-z_0[j];
	bool areSepatered= (xr0*xr0+yr0*yr0+zr0*zr0 >= D*D);
	return areSepatered;
}
// 4. conflict(int i, int j): returns true if there is a conflict between i and j, 
//       i.e., if the separation constraints with safety distance D are not satisfied by their trajectories
bool conflict(int i, int j){
	double xr0=x_0[i]-x_0[j];
	double yr0=y_0[i]-y_0[j];
    double zr0=z_0[i]-z_0[j];
	double vrx=vx[i]-vx[j];
	double vry=vy[i]-vy[j];
	double vrz=vz[i]-vz[j];
	float negative_zero=-1/std::numeric_limits<float>::infinity();
	double a = (vrx*xr0+vry*yr0+vrz*zr0);
    double dmin = (xr0*xr0+yr0*yr0+zr0*zr0) - (a*a)/(vrx*vrx+vry*vry+vrz*vrz); 	//minimum distance among i and j
	if (dmin-D*D >=0 || dmin-D*D==negative_zero) {return false;}
	else {return true;}
}
// 5. numConflicts(int i,bool* explored): returns the number of conflicts between aircraft i and
//     the rest of aircraft j such that explored[j]=true
int numConflicts(int i,bool* explored){
	int num_conflicts=0;
	for (int j = 0; j < n; ++j){
		if(explored[j] && conflict(i,j))
			num_conflicts++;
	}
	return num_conflicts;
}
// 6. crossing(float pos1, float vel1, float pos2): calculates horizontal or vertical crossing time
float crossing(float pos1, float vel1, float pos2){
	return (pos2-pos1)/vel1;
}
// 7. getInstantBoundary(int i,float height,float width,float altitude): used in random and pseudo-random scenarios.
//  Get instants in wich 'i' will cross boundary of the space window height x width. 
//  This serves to calculate final position (x_t,y_t) in which i leaves the observed airspace, which is returned by this method.
float getInstantBoundary(int i,float height,float width,float altitude){
	float t[6];
	for (int k = 0; k < 6; ++k) t[k]=-1;
	float instant_boundary=std::numeric_limits<float>::infinity();
	if(vy[i]!=0){
		t[0]=crossing(y_0[i],vy[i],height); // north
		t[1]=crossing(y_0[i],vy[i],0); // south
	}
	if(vx[i]!=0){
		t[2]=crossing(x_0[i],vx[i],0); // west
		t[3]=crossing(x_0[i],vx[i],width); // east
	}
	if(vz[i]!=0){
		t[4]=crossing(z_0[i],vz[i],0); // down
		t[5]=crossing(z_0[i],vz[i],altitude); // up
	}
	for (int k = 0; k < 6; ++k){
		if(t[k]>0 && t[k]<instant_boundary) //we should exclude t[k]=0 (we are looking for 'exit' time, not 'entry' time)
			instant_boundary=t[k];
	}
	return instant_boundary;
}
// 8. correctAngle(float angle): returns an equivalent angle in the interval [-pi,pi]
float correctAngle(float angle){
	while( !( angle>= -M_PI && angle<= M_PI) ){
		if(angle > M_PI)
			angle=angle-2*M_PI;
		if(angle <-M_PI)
			angle=angle+2*M_PI;
	}
	return angle;
}
// 9. Norm of a vector
double norm(float* vector, int dimension){
	double accum = 0.;
	double norm = 0;
	for (int j = 0; j < dimension; ++j) 
		{
        accum += (double)(vector[j]) * (double)(vector[j]);
    	}
    norm = sqrt((double) accum);
	return norm;
}
//--------------------------------------------------------------------------
//                             I/O functions
//--------------------------------------------------------------------------
void printBenchmarkInfo(){
	printf("\n--------------------------------------------------------------------------------------------------\n");
	printf("Scenario generated: Instance with %i aircraft, characterized by initial positions and velocities:\n",n);
	printf("--------------------------------------------------------------------------------------------------\n");
    //---------------------benchmark output
    if(outputFile){
    	printf("Instance saved in file %s\n",outputFile);
   		freopen(outputFile, "a", stdout);
    }
    printf("p0={\n"); //initial position
    for (int i = 0; i < n; ++i){
    	printf("%f \t %f",x_0[i],y_0[i]);
    	if (is3D) printf("\t %f",z_0[i]);
    	printf("\n");
    }
	printf("}\n");
	// printf("p1={\n"); //final position
    // for (int i = 0; i < n; ++i){
    // 	printf("%f \t %f",x_t[i],y_t[i]);
    // 	if (is3D) printf("\t %f",z_t[i]);
    // 	printf("\n");
    // }
    // printf("}\n");
    printf("V_polar=(v,theta"); //polar velocity i.e. (hat_v,hat_theta,hat_phi)
    if (is3D) printf(",phi");
    printf(")={\n");
    for (int i = 0; i < n; ++i){
    	printf("%f \t %f",hat_v[i],correctAngle(hat_theta[i]));
    	if (is3D) printf("\t %f",correctAngle(hat_phi[i]));
    	printf("\n");
    }
    printf("}\n");
    //---------------------complementary output
    printf("(Vx,Vy,Vz)={\n"); //components of the velocity
    for (int i = 0; i < n; ++i){
    	printf("%f \t %f",vx[i],vy[i]);
    	if (is3D) printf("\t %f",vz[i]);
    	printf("\n");
    }
	printf("}\n");
    if(outputFile)
		freopen("/dev/tty", "w", stdout); /*restore output*/
    printf("--------------------------------------------------------------------------------------------------\n");
	printf("Some additional information about the instance you just generated:\n");
	printf("--------------------------------------------------------------------------------------------------\n");
    printf("Pairs in conflict:\n");
    int num_conflicts=0;
    int conf_aircraft[n];
    for (int i = 0; i < n; ++i) conf_aircraft[i]=0;
    for (int i = 0; i < n; ++i) {
    	for (int j = i+1; j < n; ++j){
    		if(conflict(i,j)){
    			printf("(%i, %i)\n",i+1,j+1);
    			conf_aircraft[i]++;
    			conf_aircraft[j]++;
    			num_conflicts++;
    		}
    	}
    }
    printf("Nº conflicts per aircraft:\n");
    int aircraft_conflict=0;
    for (int i = 0; i < n; ++i){
    	printf("%i, ", conf_aircraft[i]);
    	if(conf_aircraft[i]!=0) aircraft_conflict++;
    }
    printf("\nProportion of aircraft with one conflict or more: %f\nTotal pairs in conflict: %i\n",aircraft_conflict/(n*1.0),num_conflicts);
	printf("Proportion of pairs of conflicts: %f\n",(2.0*num_conflicts)/(n*(n-1)));
	printf("Graphical representation saved in Figure.tex\n");

    // Generate Latex file with figure 
    string colors[18]={ "brown", "cyan", "blue","green", 
    				   "lightgray", "lime", "magenta", "olive", "gray",  "orange", 
    				   "pink", "purple", "red", "teal", "violet", "teal", "yellow","darkgray"};
	int color_index=0;
	ofstream f("Figure.tex");
	if (!f.is_open()){ printf("Error when trying to write latex file\n");return;}
	f<<"\\documentclass[border=10pt,varwidth]{standalone}\n";
	f<<"\\usepackage{tikz,tikz-3dplot}\n";
    f<<"\\usetikzlibrary{shapes,snakes}\n";
	f<<"\\begin{document}\n";
	float scale;
	if (mode==6 ||mode==7 || mode==14 || mode==1 || mode==0 || mode==8 || mode==9 || mode==15) {scale=0.02;}
	else {scale=0.05;}
	if (is3D){
		f<<"\\tdplotsetmaincoords{70}{145}\n";
		f<<"\\begin{tikzpicture} [scale="<<scale<<", tdplot_main_coords, axis/.style={->,black,thick}, vector/.style={-stealth,black,very thick},vector guide/.style={dotted,black,thick},]\n";
		f<<"\\coordinate (O) at (0,0,0);\n";
		f<<"\\pgfmathsetmacro{\\ax}{1}\n";
		f<<"\\pgfmathsetmacro{\\ay}{-1}\n";
		f<<"\\pgfmathsetmacro{\\az}{0.5}\n";
		if (mode==10||mode==11|| mode==12 || mode==13){
			f<<"\\draw[axis] (0,0,0) -- (150,0,0) node[anchor=north east]{$x$};\n";  
			f<<"\\draw[axis] (0,0,0) -- (0,150,0) node[anchor=south]{$y$};\n"; 
			f<<"\\draw[axis] (0,0,0) -- (0,0,150) node[anchor=south]{$z$};\n";
		}
		else{
			f<<"\\draw[axis] (0,0,0) -- (400,0,0) node[anchor=north east]{$x$};\n";  
			f<<"\\draw[axis] (0,0,0) -- (0,400,0) node[anchor=south]{$y$};\n"; 
			f<<"\\draw[axis] (0,0,0) -- (0,0,400) node[anchor=south]{$z$};\n";
		}
        //f<<"\\node[star, fill,  minimum size=3pt, inner sep=0pt,label=right:alt_step] at (0,0,"<<altitude_step<<") {};\n";
	}
	else{f<<"\\begin{tikzpicture}[scale="<<scale<<"]\n";}

    for(int i=0;i<n;i++){
      f<<"\\coordinate ("<<i+1<<") at ("<<x_0[i]<<", "<<y_0[i];
      if(is3D) f<<", "<<z_0[i];
      f<<");\n";
      f<<"\\draw[color="<<colors[color_index]<<",thick, ->] ("<<i+1<<") -- ("<<x_t[i]<<", "<<y_t[i];
      if(is3D) f<<", "<<z_t[i];
      f<<");\n";
      if(x_0[i]>=0.5) f<<"\\node[circle, fill,  minimum size=3pt, inner sep=0pt,label=right:"<<i+1<<"] at ("<<i+1<<") {};\n";
      else if(x_0[i]<=-0.5) f<<"\\node[circle, fill,  minimum size=3pt, inner sep=0pt,label=left:"<<i+1<<"] at ("<<i+1<<") {};\n";
      else {
      	if(y_0[i]>=0) f<<"\\node[circle, fill,  minimum size=4pt, inner sep=0pt,label=above:"<<i+1<<"] at ("<<i+1<<") {};\n";
		else f<<"\\node[circle, fill,  minimum size=3pt, inner sep=0pt,label=below:"<<i+1<<"] at ("<<i+1<<") {};\n";
	  }
	  if(color_index<18) color_index++;
	  if(color_index==18) color_index=0;
    } 

    f<<"\\end{tikzpicture}\n";
    f<<"\\end{document}\n";
    f.close();
}

// --------------------------------------------------------------------------
// -----                  SCENARIO GENERATORS                         -------
// --------------------------------------------------------------------------
/*
There are 10 scenarios in total, each one has an associated code in the variable "mode"
- 2D Scenarios:
	*Predefined:
	    + circle 		--> mode= 0
	    + randomCircle  --> mode= 1
	    + rombo  		--> mode= 2
		+ randomRombo   --> mode= 3
		+ grid  		--> mode= 4
		+ randomGrid    --> mode= 5
	*Random:
		+ pseudorandom  --> mode= 6
		+ random 		--> mode= 7
- 3D Scenarios:
	*Predefined:
		+ sphere			--> mode= 8
	    + randomsphere 		--> mode= 9
	    + Polyhedral    	--> mode= 10
	    + randomPolyhedralP --> mode= 11
	    + Grid3d        	--> mode= 12
	    + randomGrid3d  	--> mode= 13
	*Random:
		+ pseudorandom  	--> mode= 14
		+ random 			--> mode= 15  
*/

// --------------------------------------------------------------------------
// -----                 PREDEFINED SCENARIOS                        -------
// --------------------------------------------------------------------------
//**************   Predefined 2D scenarios   *******************
// ####         1. CIRCLE SCENARIOS     ###
// randomCircleP: set initial positions (x_0,y_0), vectors of velocity (vx,vy)
//      and final positions (x_t,y_t) of the random circle problem.
//		Heading angles are deviated from those of circle problem a random amount in the interval [hamin,hamax]
//      sector_ini=0 and sector_angle=2*M_PI gives the classical circle configuration
// 		sector_ini=M_PI and sector_angle=M_PI/2 gives the "quarter" circle configuration in S.  Cafieri &  N.  Durand. (2014):   
// 		"Aircraft  deconfliction  with  speed  regulation:  new  modelsfrom mixed-integer optimization", JOGO.
//		Otherwise, these parameters allow the user to select where the first aircraft is placed (sector_ini) and the
//		total angle covered by the sector between the first and the last aircraft (sector_angle)
void randomCircleP(float radius,float sector_ini_x,double sector_angle_x,float hamin,float hamax){
	float angle_step;
	if (sector_angle_x==2*M_PI)
		angle_step=sector_angle_x/n;
	else
		angle_step=sector_angle_x/(n-1);
	float v,dev;
	for (int i = 0; i < n; ++i){
		x_0[i]=radius*cos(i*angle_step+sector_ini_x);
		y_0[i]=radius*sin(i*angle_step+sector_ini_x);
		v=randomFloat(vmin,vmax);
		dev=randomFloat(hamin,hamax);
		vx[i]=-v*cos(i*angle_step+sector_ini_x+dev);
		vy[i]=-v*sin(i*angle_step+sector_ini_x+dev);
		hat_v[i]=v;
		hat_theta[i]=i*angle_step+dev+sector_ini_x;
		x_t[i]=x_0[i]-2*radius*cos(i*angle_step+sector_ini_x+dev);
		y_t[i]=y_0[i]-2*radius*sin(i*angle_step+sector_ini_x+dev);
		float final_position[2]= { x_t[i],y_t[i] };
		float final_position_norm=norm(final_position,2);
		if(abs(x_t[i])>radius || abs(y_t[i])>radius)
			{x_t[i]=radius*x_t[i]/final_position_norm; y_t[i]=radius*y_t[i]/final_position_norm;}
	}
    if (hamin==0 and hamax==0){
        printf("\n--------------------------------------------------------------------------------------------------\n");
	    printf("Circle Problem Generator: circle with radius: %.2f\n",radius);
        printf("Aircraft generated with heading angles in the range: [%.2f,%.2f]",sector_ini_x,sector_angle_x);
    }
    else{
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Random-circle Problem Generator: circle with radius: %.2f\n",radius);
        printf("Aircraft generated with heading angles in the range: [%.2f,%.2f]",sector_ini_x,sector_angle_x);
    }
}
// CircleP(float radius): scenario generated with the function randomCircleP and 0 angle deviation
void CircleP(float radius, float sector_ini_x,double sector_angle_x){
	randomCircleP(radius,sector_ini_x,sector_angle_x,0,0);
}


// ####      2.   RHOMBOIDAL and GRID SCENARIOS     ###
// randomRomboP: set initial positions (x_0,y_0),  vectors of velocity (vx,vy) and 
//     final positions (x_t,y_t) of the random rombo problem. Heading angles are deviated from 
//     those of rombo problem a random amount in the interval [hamin,hamax]
// 	   horiz_trails (resp. verti_trails)= number of horizontal (resp. vertical) trails
//     nx (resp. ny)= number of aircraft per horizontal (resp. vertical) trail
//     dx (resp. dy)= separation between horizontal (resp. vertical) trails
void randomRomboP(double alpha_r,int horiz_trails,int verti_trails,int nx,int ny,float dx,float dy,float d_aircraft,float hamin,float hamax){
	float v,dev;
	float first_height=sin(alpha_r)*(ny-1)*d_aircraft+dx; // vertical position of the first horizontal trail (from lower to upper)
	float first_width=(nx-1)*d_aircraft+dy; // horizontal position of the first vertical trail (from left to right)
	float last_width=first_width+(verti_trails-1)*dy; // end of horizontal trails (for drawing)
	// Loop to create trajectories at horizontal trails:
	for (int i = 0; i < horiz_trails; ++i){  
		for (int k = 0; k < nx; ++k){			// * *------------------------------- i= horiz_trails-1
			x_0[i*nx+k]=k*d_aircraft;			// * *------------------------------- i=2
			y_0[i*nx+k]=first_height+i*dx;		// * *------------------------------- i=1
			v=randomFloat(vmin,vmax);			// * *------------------------------- i=0
			dev=randomFloat(hamin,hamax);		// *=aircraft, in this example nx=2 aircraft per trail
			vx[i*nx+k]=v*cos(dev);
			vy[i*nx+k]=v*sin(dev);
			hat_v[i*nx+k]=v;
			hat_theta[i*nx+k]=dev;
			x_t[i*nx+k]=x_0[i*nx+k]+last_width*cos(dev);
			y_t[i*nx+k]=y_0[i*nx+k]+last_width*sin(dev);
		}
	}
	int cont=horiz_trails*nx; // number of aircraft for which we already generated trajectories
	float ini_width;
	int sign;
	if(alpha_r<=M_PI/2){ // corresponds with this slope: "/"
		ini_width=first_width;
		sign=1;
	}
	else{ // corresponds with this slope: "\"
		ini_width=last_width;
		sign=-1;
	}
	// Loop to create trajectories at vertical trails
	for (int i = 0; i < verti_trails; ++i){     
		for (int k = 0; k < ny; ++k){       				        			//----/---/---/---/---/---/---- 
			x_0[cont+i*ny+k]=ini_width+sign*i*dy+k*d_aircraft*cos(alpha_r); 	//----/---/---/---/---/---/-----
			y_0[cont+i*ny+k]=k*d_aircraft*sin(alpha_r);		        			//----/---/---/---/---/---/-----								
			v=randomFloat(vmin,vmax);						    				//---/---/---/---/---/---/-----
			dev=randomFloat(hamin,hamax);				       					// i=0  i=1  i=2     i=verti_trails-1
			vx[cont+i*ny+k]=v*cos(alpha_r+dev);
			vy[cont+i*ny+k]=v*sin(alpha_r+dev);
			hat_v[cont+i*ny+k]=v;
			hat_theta[cont+i*ny+k]=alpha_r+dev;
			x_t[cont+i*ny+k]=x_0[cont+i*ny+k]+(first_height+dx*horiz_trails)*cos(alpha_r+dev);
			y_t[cont+i*ny+k]=y_0[cont+i*ny+k]+(first_height+dx*horiz_trails)*sin(alpha_r+dev);
		}
	}

    if (hamin==0 and hamax==0 and alpha_r==M_PI/2){
        printf("\n--------------------------------------------------------------------------------------------------\n");
	    printf("Grid Problem Generator");    }
    else if (hamin==0 and hamax==0) {
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Rombo Problem Generator");
    }
    else if (alpha_r==M_PI/2){
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Random-grid Problem Generator");
    }
    else{
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Random-rombo Problem Generator");
    }
}
// RomboP: scenario generated with the function randomRomboP and 0 angle deviation
void RomboP(float alpha_r,int horiz_trails,int verti_trails,int nx,int ny,float dx,float dy,float d_aircraft){
	randomRomboP(alpha_r,horiz_trails,verti_trails,nx,ny,dx,dy,d_aircraft,0,0);
}
// GridP: scenario generated with the function randomRomboP, alpha_r=pi/2 and 0 angle deviation
void GridP(int horiz_trails,int verti_trails,int nx,int ny,float dx,float dy,float d_aircraft){
	randomRomboP(M_PI/2,horiz_trails,verti_trails,nx,ny,dx,dy,d_aircraft,0,0);
}
// GridP: scenario generated with the function randomRomboP and alpha_r=pi/2
void randomGridP(int horiz_trails,int verti_trails,int nx,int ny,float dx,float dy,float d_aircraft,float hamin,float hamax){
	randomRomboP(M_PI/2,horiz_trails,verti_trails,nx,ny,dx,dy,d_aircraft,hamin,hamax);
}



//************** Predefined 3D scenarios   *******************
// ####      3.   SPHERE SCENARIOS     ###
// randomSphereP: set initial positions (x_0,y_0,z0), vectors of velocity (vx,vy,vz) 
// and final positions (x_t,y_t,z_t) of the sphere problem, with the possibility to be in a certain portion of the sphere
// Heading angles are deviated from those of sphere problem a random amount in the interval [hamin,hamax]
void randomSphereP(float radius,float sector_ini_x,double sector_angle_x,float sector_ini_z,double sector_angle_z,float hamin,float hamax){
	double theta[n]; double phi[n]; // starting angles of the 3d configuration
	float v,dev1,dev2;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> uniform01(0.0, 1.0); // uniform distribution
	for (int i = 0; i < n; ++i){
	   	theta[i] = sector_ini_x+sector_angle_x*uniform01(gen); // we generate uniformly distributed points on the sphere 
		phi[i] = abs(acos(cos((double)(sector_ini_z)) - (cos((double)(sector_ini_z))-cos((double)(sector_ini_z)+sector_angle_z)) * uniform01(gen))); //phi must be in [sector_ini_z,sector_ini_z+sector_angle_z]
		x_0[i] = radius*sin(phi[i]) * cos(theta[i]);
		y_0[i] = radius*sin(phi[i]) * sin(theta[i]);
		z_0[i] = radius*cos(phi[i]);
		v=randomFloat(vmin,vmax);
		dev1=randomFloat(hamin,hamax);
        dev2=randomFloat(hamin,hamax);
		if (phi[i]+dev2<sector_ini_z) {dev2=sector_ini_z-phi[i];}
		else if (phi[i]+dev2>sector_ini_z+sector_angle_z) {dev2=sector_ini_z+sector_angle_z-phi[i];}
		vx[i]=-v*cos(theta[i]+dev1)*sin(phi[i]+dev2);
		vy[i]=-v*sin(theta[i]+dev1)*sin(phi[i]+dev2);
        vz[i]=-v*cos(phi[i]+dev2);
		hat_v[i]=v;
		hat_theta[i]=theta[i]+dev1;
        hat_phi[i]=phi[i]+dev2;
		x_t[i]=x_0[i]-2*radius*cos(theta[i]+dev1)*sin(phi[i]+dev2);
		y_t[i]=y_0[i]-2*radius*sin(theta[i]+dev1)*sin(phi[i]+dev2);
        z_t[i]=z_0[i]-2*radius*cos(phi[i]+dev2);
		// we want the final positions of the aircraft being ON the surface of the sphere
		float final_position[3]= { x_t[i],y_t[i],z_t[i] };
		float final_position_norm=norm(final_position,3);
		if(final_position_norm>radius) //points outside the sphere
			{x_t[i]=radius*x_t[i]/final_position_norm;  y_t[i]=radius*y_t[i]/final_position_norm; z_t[i]=radius*z_t[i]/final_position_norm;}
	}
    if (hamin==0 and hamax==0){
        printf("\n--------------------------------------------------------------------------------------------------\n");
	    printf("Sphere Problem Generator: sphere with radius: %.2f\n",radius);
        printf("Aircraft generated with theta in the range [%.2f,%.2f], and phi in the range [%.2f,%.2f]",sector_ini_x,sector_ini_x+sector_angle_x,sector_ini_z,sector_ini_z+sector_angle_z);
    }
    else{
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Random-sphere Problem Generator: sphere with radius: %.2f\n",radius);
        printf("Aircraft generated with theta in the range [%.2f,%.2f], and phi in the range [%.2f,%.2f]",sector_ini_x,sector_ini_x+sector_angle_x,sector_ini_z,sector_ini_z+sector_angle_z);
    }
}
// SphereP: scenario generated with the function randomSphereP and 0 angle deviation
void SphereP(float radius,float sector_ini_x,double sector_angle_x,float sector_ini_z,double sector_angle_z){
	randomSphereP(radius,sector_ini_x,sector_angle_x,sector_ini_z,sector_angle_z,0,0);
}

// ####       4. POLYHEDRAL SCENARIOS     ###
// ----------------------------------------
// HP = horizontal planes = planes parallel to xy plane
// VP = vertical planes = planes parallel to xz plane
// mx = horizontal trails on horizontal planes = trails parallel to x_axis and perpendicular to y_axis
// my = "vertical" trails on horizontal planes = trails forming an alpha angle with the HP_HT (and with x_axis)
// mz = "vertical" trails on vertical planes = trails forming a beta angle with the VP_HT (and with x axis)
// ----------------------------------------
// randomPolyhedralP: set initial positions (x_0,y_0,z_0), vectors of velocity (vx,vy,vz) and 
//     final positions (x_t,y_t,z_t) of the random polyhedral problem. Heading angles are deviated from 
//     those of polyhedral problem of a random amount in the interval [hamin,hamax]
void randomPolyhedralP(double* alpha,double* beta,int* mx,int* my,int* mz,int nx,int ny,int nz,float d_HP,float d_VP,float dx,float dy,float dz,float d_aircraft,float hamin,float hamax){
	float v,dev;
	float ini_my;   float ini_mz;  int sign; //necessary since alpha and beta could be > or < 90°
	
    float first_VP=d_VP;   // horizontal position of the first vertical plane  (from left to right)
    // since there is a different beta[j] for each vertical plane j, we use the following approach to compute the vertical position of the first horizontal plane 
	// (similar to the first position of the first horizontal trail on each horizontal plane, i.e. first_mx[j])
	float a=0;
	for (int j=0; j < VP; ++j){
		if(sin(beta[j])*(nz-1)*d_aircraft>a)
			a=sin(beta[j])*(nz-1)*d_aircraft;
	}
    float first_HP=a+d_HP; // vertical position of the first horizontal plane (from lower to upper) 
	

	// We focus on horizontal planes:
    float first_mx[HP]; // vector of vertical positions of the first horizontal trail (from lower to upper) for each horizontal plane
	float first_my=(nx-1)*d_aircraft+dy; // horizontal position of the first sloping trail (from left to right), the same for each horizontal plane
	float last_my[HP]; //  horizontal position of the last sloping trail for each horizontal plane
	int cont=0;	 // number of aircraft for which we already generated trajectories
	for (int j=0; j < HP; ++j){
		// For each j, we consider the j-th HP:
		first_mx[j]=sin(alpha[j])*(ny-1)*d_aircraft+dx; // vertical positions of the first horizontal trail on the j-th horizontal plane
    	last_my[j]=first_my+dy*(my[j]-1); //  horizontal position of the last sloping trail on the j-th horizontal plane
		// Loop to create trajectories at horizontal trails:
        for (int i = 0; i < mx[j]; ++i){
            for (int k = 0; k < nx; ++k){			    // * *------------------------------- i= mx[j]-1
                x_0[cont+i*nx+k]=k*d_aircraft;			// * *------------------------------- i=2
                y_0[cont+i*nx+k]=first_mx[j]+i*dx;		// * *------------------------------- i=1
                z_0[cont+i*nx+k]=first_HP+j*d_HP;		// * *------------------------------- i=0
                v=randomFloat(vmin,vmax);				// *=aircraft, in this example nx=2 aircraft per trail     
                dev=randomFloat(hamin,hamax);		            
                vx[cont+i*nx+k]=v*cos(dev);
                vy[cont+i*nx+k]=v*sin(dev);
                vz[cont+i*nx+k]=0;
                hat_v[cont+i*nx+k]=v;
                hat_theta[cont+i*nx+k]=dev;
                x_t[cont+i*nx+k]=x_0[cont+i*nx+k]+last_my[j]*cos(dev);
                y_t[cont+i*nx+k]=y_0[cont+i*nx+k]+last_my[j]*sin(dev);
                z_t[cont+i*nx+k]=z_0[cont+i*nx+k];
            }
        }
		cont=cont+mx[j]*nx;

		if(alpha[j]<=M_PI/2){ // corresponds with this slope: "/"
			ini_my=first_my;
			sign=1;}
		else{// corresponds with this slope: "\"
			ini_my=last_my[j];
			sign=-1;}
		// Loop to create trajectories at sloping trails
        for (int i = 0; i < my[j]; ++i){  
            for (int k = 0; k < ny; ++k){       				        		//----/---/---/---/---/---/---- 
                x_0[cont+i*ny+k]=ini_my+sign*i*dy+k*d_aircraft*cos(alpha[j]); 	//----/---/---/---/---/---/-----
                y_0[cont+i*ny+k]=k*d_aircraft*sin(alpha[j]);		        	//----/---/---/---/---/---/------									
                z_0[cont+i*ny+k]=first_HP+j*d_HP;								//----/---/---/---/---/---/------
                v=randomFloat(vmin,vmax);						    			// i=0  i=1  i=2     	  i=my-1
                dev=randomFloat(hamin,hamax);				       				
                vx[cont+i*ny+k]=v*cos(alpha[j]+dev);
                vy[cont+i*ny+k]=v*sin(alpha[j]+dev);
                vz[cont+i*ny+k]=0;
                hat_v[cont+i*ny+k]=v;
                hat_theta[cont+i*ny+k]=alpha[j]+dev;
                x_t[cont+i*ny+k]=x_0[cont+i*ny+k]+(first_mx[j]+dx*mx[j])*cos(alpha[j]+dev);
                y_t[cont+i*ny+k]=y_0[cont+i*ny+k]+(first_mx[j]+dx*mx[j])*sin(alpha[j]+dev);
                z_t[cont+i*ny+k]=z_0[cont+i*ny+k];
            }
        }
		cont=cont+my[j]*ny;
    }

	// We focus on vertical planes, where only sloping (parallel to (x,z)-plane and perpendicular to (x,y)-plane) trajectories are generated:
    float first_mz=dz; // horizontal position of the first sloping trail (from left to right), the same for each vertical plane
	float last_mz[VP]; // vector of horizontal positions of the last sloping trail (from left to right) for each vertical plane
	for (int j=0; j < VP; ++j){
		// For each j, we consider the j-th VP:
		last_mz[j]=first_mz+dz*(mz[j]-1);
		if(beta[j]<=M_PI/2){ // corresponds with this slope: "/"
			ini_mz=first_mz;
			sign=1;}
		else{// corresponds with this slope: "\"
			ini_mz=last_mz[j];
			sign=-1;}
    	// Loop to create trajectories at sloping trails (assuming aircraft start flying from the bottom):
        for (int i = 0; i < mz[j]; ++i){  
            for (int k = 0; k < nz; ++k){			            
                x_0[cont+i*nz+k]=ini_mz+sign*i*dz+k*d_aircraft*cos(beta[j]); 	//----/---/---/---/---/---/-----
                y_0[cont+i*nz+k]=first_VP+j*d_VP;	        					//----/---/---/---/---/---/------									
                z_0[cont+i*nz+k]=k*d_aircraft*sin(beta[j]);						//----/---/---/---/---/---/------
                v=randomFloat(vmin,vmax);						    			// i=0  i=1  i=2     	  i=mz[j]-1
                dev=randomFloat(hamin,hamax);				       				
                vx[cont+i*nz+k]=v*cos(beta[j]+dev);
                vy[cont+i*nz+k]=0;
                vz[cont+i*nz+k]=v*sin(beta[j]+dev);
                hat_v[cont+i*nz+k]=v;
                hat_theta[cont+i*nz+k]=beta[j]+dev;
                x_t[cont+i*nz+k]=x_0[cont+i*nz+k]+(first_HP+d_HP*HP)*cos(beta[j]+dev);
                y_t[cont+i*nz+k]=y_0[cont+i*nz+k];
                z_t[cont+i*nz+k]=z_0[cont+i*nz+k]+(first_HP+d_HP*HP)*sin(beta[j]+dev);
            }
        }
		cont=cont+mz[j]*nz;
    }
    // output message:
    bool grid = true;
    for(int i=0;i<HP;++i){
        if (alpha[i]!=M_PI/2){ grid = false;}
        if (beta[i]!=M_PI/2) { grid =false;}
    }
    if (hamin==0 and hamax==0 and grid){
        printf("\n--------------------------------------------------------------------------------------------------\n");
	    printf("Cubic Problem Generator");    
    }
    else if (hamin==0 and hamax==0) {
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Polyhedral Problem Generator");
    }
    else if (grid){
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Random-cubic Problem Generator");
    }
    else{
        printf("\n--------------------------------------------------------------------------------------------------\n");
        printf("Random-polyhedral Problem Generator");
    }
}
// PolyhedralP: scenario generated with the function randomPolyhedralP and 0 angle deviation
void PolyhedralP(double* alpha,double* beta,int* mx,int* my,int* mz,int nx,int ny,int nz,float d_HP,float d_VP,float dx,float dy,float dz,float d_aircraft){
    randomPolyhedralP(alpha,beta,mx,my,mz,nx,ny,nz,d_HP,d_VP,dx,dy,dz,d_aircraft,0,0);
}
// Grid3dP: scenario generated with the function randomPolyhedralP and alpha=beta=pi/2 and 0 angle deviation
void Grid3dP(int* mx,int* my,int* mz,int nx,int ny,int nz,float d_HP,float d_VP,float dx,float dy,float dz,float d_aircraft){
	double grid3dalpha[HP];
	double grid3dbeta[VP];
	for(int i=0;i<HP;++i){
		grid3dalpha[i]=M_PI/2;
	}
	for(int i=0;i<VP;++i){
		grid3dbeta[i]=M_PI/2;
	}
	randomPolyhedralP(grid3dalpha,grid3dbeta,mx,my,mz,nx,ny,nz,d_HP,d_VP,dx,dy,dz,d_aircraft,0,0);
}
// randomGrid3dP: scenario generated with the function randomPolyhedral and alpha=beta=pi/2
void randomGrid3dP(int* mx,int* my,int* mz,int nx,int ny,int nz,float d_HP,float d_VP,float dx,float dy,float dz,float d_aircraft,float hamin,float hamax){
	double grid3dalpha[HP];
	double grid3dbeta[VP];
	for(int i=0;i<HP;++i){
		grid3dalpha[i]=M_PI/2;
	}
	for(int i=0;i<VP;++i){
		grid3dbeta[i]=M_PI/2;
	}
	randomPolyhedralP(grid3dalpha,grid3dbeta,mx,my,mz,nx,ny,nz,d_HP,d_VP,dx,dy,dz,d_aircraft,hamin,hamax);
}


//************** Random scenarios   *******************
// This auxiliary method generates the initial position of aircraft for a pseudo-random scenario
void pseudoRandomP_iniPos(char* air_config,float height,float width, float altitude, float* hamin,float* hamax, float* phimin, float* phimax){
	// 1. -------- Partition airspace into a 2d or  3d matrix of size (m x l) or (m x l x a) respectively s.t. m/l=height/width m/a=height/altitude l/a=width/altitude
	// This is done to create cells in the airspace and place one aircraft in each cell, instead of randomly generating the initial positions in the full airspace. 
    //  This guarantees that the aircraft are spread.
	float m,l,a; // dimensions of the matrix. Each aircraft is placed on one cell of the resulting matrix
	
    int new_n; // new_n is the total number of cells   
	// The value of new_n depends on configuration:
	// ***** 2D cases *****
	// "W-N": aircraft placed at left and top of R, new_n=m+l-1
	// "N-S": aircraft placed at bottom and top of R, new_n=2l
	// "W-E": aircraft placed at left and right of R, new_n=2m
	// "all": aircraft placed at all sides of R, new_n=2*m+2*l-4
	// ***** 3D cases ***** //D=bottom, U=top
	// "W-U": aircraft placed at (x=0,y,z) and top faces of P, new_n=m*a+m*l-m
	// "N-U": aircraft placed at (x,y=0,z) and top faces of P, new_n=l*a+m*l-l
	// "U-D": aircraft placed at bottom and top faces of P, new_n=2*m*l
    // "W-N": aircraft placed at (x=0,y,z) and (x,y=0,z) faces of P, new_n=m*a+l*a-a
	// "N-S": aircraft placed at (x,y=0,z) and (x,y=height,z) faces of P, new_n=2*l*a
    // "W-E": aircraft placed at (x=0,y,z) and (x=width,y,z) faces of P, new_n=2*m*a
	// "all": aircraft placed at all faces of P, new_n=(2*m+2*l-4)*a+(l-2)*(m-2)*2

	// new_n can be different from n, we will leave empty cells if needed

	// In the following, we generate m,l,a and new_n depending on the selected airspace configuration:
	int all=0; // will be set to 1 if the selected configuration is 'all'													
	if(is3D==false){
        a=1; //it should be 0, but we set it to 1 in order to avoid an undefined altitude_step = altitude/a in the following
		if(strcmp(air_config,"W-N")==0){					// Scenario W-N: (in the example, m=5 (rows), l=4 (columns))	
			l=(n+1)/(height/width+1);						//  x  x  x  x  	
			m=ceil(height*l/width);							//  x
			l=ceil(l);										//  x
			new_n=m+l-1; // number of "x" in the example:8	//  x								
		}													//  x
		else if(strcmp(air_config,"N-S")==0){				// Scenario N-S:
			l=n/2.0;										//  x  x  x  x  x
			m=ceil(height*l/width); 
			l=ceil(l);
			new_n=2*l;										//  x  x  x  x  x	
		}													// Scenario W-E:
		else if(strcmp(air_config,"W-E")==0){				//  x           x
			m=n/2.0;										//  x           x
			l=ceil(width*m/height);							//  x           x
			m=ceil(m);										//  x           x
			new_n=2*m;										//  x           x
		}										 			// Scenario all:
		else{												//  x  x  x  x  x
			l=(n+4)/(2*(height/width+1));					//  x           x
			m=ceil(height*l/width);							//  x           x
			l=ceil(l);										//  x           x
			new_n=2*m+2*l-4;								//  x  x  x  x  x
			all=1;
			if(strcmp(air_config,"all")!=0 and verbose)
				printf("Pseudo-random Problem Generator: unvalid or unexisting airspace configuration for random scenario, using default\n");
		}
		printf("\n----------------------------------------------------------------------------------------------\n");
		printf("Pseudo-random Problem Generator: grid with %i rows, %i columns, %i cells, %i of them not used\n",(int)m,(int)l,new_n,new_n-n);
	}
	else{
        if(strcmp(air_config,"W-U")==0){
			float m_true;  					
            m_true=(1+sqrt(1+4*n*(altitude+width)/height))/(2*(altitude+width)/height);		
            // we approximate al the parameters to the previous integer
            l=floor(m_true*width/height);					
            a=floor(m_true*altitude/height);
            m=floor(m_true);
            new_n=m*a+m*l-m;
			if(new_n<n){
                // we approximate the maximum among m, l, and a to the next integer
				if (max(max(m,l),a)==m){m=ceil(m_true);l=floor(m_true*width/height);	a=floor(m_true*altitude/height);}
				else if (max(max(m,l),a)==l){m=floor(m_true);  l=ceil(m_true*width/height);	a=floor(m_true*altitude/height);}
				else{m=floor(m_true);l=floor(m_true*width/height);	a=ceil(m_true*altitude/height);}
				new_n=m*a+m*l-m;
			}
			if(new_n<n){
                // we approximate the maximum among the remaining parameters to the next integer
				if (m==ceil(m_true)&&max(l,a)==l){l=ceil(m_true*width/height);a=floor(m_true*altitude/height);}
				else if (m==ceil(m_true)&&max(l,a)==a){l=floor(m_true*width/height);a=ceil(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==m){m=ceil(m_true);a=floor(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==a){m=floor(m_true);a=ceil(m_true*altitude/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==m){m=ceil(m_true);l=floor(m_true*width/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==l){m=floor(m_true);l=ceil(m_true*width/height);}
				new_n=m*a+m*l-m;
			}
			if(new_n<n){
                // we approximate all the parameters to the next integer
				l=ceil(m_true*width/height);					
            	a=ceil(m_true*altitude/height);
            	m=ceil(m_true);
            	new_n=m*a+m*l-m;
			}
        }													
        else if(strcmp(air_config,"S-U")==0){
			float l_true;				
            l_true=(1+sqrt(1+4*n*(altitude+height)/width))/(2*(altitude+height)/width);						
            m=floor(l_true*height/width);					
            a=floor(l_true*altitude/width);
            l=floor(l_true);
            new_n=l*a+m*l-l;
			if(new_n<n){
				if (max(max(m,l),a)==m){m=ceil(l_true*height/width);l=floor(l_true);a=floor(l_true*altitude/width);}
				else if (max(max(m,l),a)==l){m=floor(l_true*height/width);l=ceil(l_true);a=floor(l_true*altitude/width);}
				else{m=floor(l_true*height/width);l=floor(l_true);a=ceil(l_true*altitude/width);}
            	new_n=l*a+m*l-l;
			}
			if(new_n<n){
				if (m==ceil(l_true*height/width)&&max(l,a)==l){l=ceil(l_true);a=floor(l_true*altitude/width);}
				else if (m==ceil(l_true*height/width)&&max(l,a)==l){l=floor(l_true);a=ceil(l_true*altitude/width);}
				else if (l==ceil(l_true)&&max(m,a)==m){m=ceil(l_true*height/width);a=floor(l_true*altitude/width);}
				else if (l==ceil(l_true)&&max(m,a)==a){m=floor(l_true*height/width);a=ceil(l_true*altitude/width);}
				else if (a==ceil(l_true*altitude/width)&&max(m,l)==m){m=ceil(l_true*height/width);l=floor(l_true);}
				else if (a==ceil(l_true*altitude/width)&&max(m,l)==l){m=floor(l_true*height/width);l=ceil(l_true);}
            	new_n=l*a+m*l-l;
			}
			if(new_n<n){
				m=ceil(l_true*height/width);					
            	a=ceil(l_true*altitude/width);
            	l=ceil(l_true);
            	new_n=l*a+m*l-l;
			}
        }
        else if(strcmp(air_config,"U-D")==0){
			float m_true;				
            m_true=sqrt(n*height/(width*2));						
            l=floor(m_true*width/height);					
            a=floor(m_true*altitude/height);
            m=floor(m_true);
            new_n=2*m*l;
			if(new_n<n){
				if (max(max(m,l),a)==m){m=ceil(m_true);l=floor(m_true*width/height);	a=floor(m_true*altitude/height);}
				else if (max(max(m,l),a)==l){m=floor(m_true);  l=ceil(m_true*width/height);	a=floor(m_true*altitude/height);}
				else{m=floor(m_true);l=floor(m_true*width/height);	a=ceil(m_true*altitude/height);}
				new_n=2*m*l;
			}
			if(new_n<n){
				if (m==ceil(m_true)&&max(l,a)==l){l=ceil(m_true*width/height);a=floor(m_true*altitude/height);}
				else if (m==ceil(m_true)&&max(l,a)==a){l=floor(m_true*width/height);a=ceil(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==m){m=ceil(m_true);a=floor(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==a){m=floor(m_true);a=ceil(m_true*altitude/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==m){m=ceil(m_true);l=floor(m_true*width/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==l){m=floor(m_true);l=ceil(m_true*width/height);}
				new_n=2*m*l;
			}
			if(new_n<n){
				l=ceil(m_true*width/height);					
            	a=ceil(m_true*altitude/height);
            	m=ceil(m_true);
            	new_n=2*m*l;
			}							
        }
        else if(strcmp(air_config,"W-N")==0){
			float a_true;			
            a_true=(1+sqrt(1+4*n*(width+height)/altitude))/(2*(width+height)/altitude);						
            m=floor(a_true*height/altitude);					
            l=floor(a_true*width/altitude);
            a=floor(a_true);
            new_n=l*a+m*a-a;
			if(new_n<n){
				if (max(max(m,l),a)==m){m=ceil(a_true*height/altitude);l=floor(a_true*width/altitude);a=floor(a_true);}
				else if (max(max(m,l),a)==l){m=floor(a_true*height/altitude);l=ceil(a_true*width/altitude);a=floor(a_true);}
				else{m=floor(a_true*height/altitude);l=floor(a_true*width/altitude);a=ceil(a_true);}
				new_n=l*a+m*a-a;

			}
			if(new_n<n){
				if (m==ceil(a_true*height/altitude)&&max(l,a)==l){l=ceil(a_true*width/altitude);a=floor(a_true);}
				else if (m==ceil(a_true*height/altitude)&&max(l,a)==a){l=floor(a_true*width/altitude);a=ceil(a_true);}
				else if (l==ceil(a_true*width/altitude)&&max(m,a)==m){m=ceil(a_true*height/altitude);a=floor(a_true);}
				else if (l==ceil(a_true*width/altitude)&&max(m,a)==a){m=floor(a_true*height/altitude);a=ceil(a_true);}
				else if (a==ceil(a_true)&&max(m,l)==m){m=ceil(a_true*height/altitude);l=floor(a_true*width/altitude);}
				else if (a==ceil(a_true)&&max(m,l)==l){m=floor(a_true*height/altitude);l=ceil(a_true*width/altitude);}
            	new_n=l*a+m*a-a;
			}
			if(new_n<n){
				m=ceil(a_true*height/altitude);					
            	l=ceil(a_true*width/altitude);
            	a=ceil(a_true);
            	new_n=l*a+m*a-a;
			}								
        }													
        else if(strcmp(air_config,"N-S")==0){
			float l_true;				
            l_true=sqrt(n*width/(altitude*2));						
            m=floor(l_true*height/width);					
            a=floor(l_true*altitude/width);
            l=floor(l_true);
            new_n=2*l*a;
			if(new_n<n){
				if (max(max(m,l),a)==m){m=ceil(l_true*height/width);l=floor(l_true);a=floor(l_true*altitude/width);}
				else if (max(max(m,l),a)==l){m=floor(l_true*height/width);l=ceil(l_true);a=floor(l_true*altitude/width);}
				else{m=floor(l_true*height/width);l=floor(l_true);a=ceil(l_true*altitude/width);}
            	new_n=2*l*a;
			}
			if(new_n<n){
				if (m==ceil(l_true*height/width)&&max(l,a)==l){l=ceil(l_true);a=floor(l_true*altitude/width);}
				else if (m==ceil(l_true*height/width)&&max(l,a)==l){l=floor(l_true);a=ceil(l_true*altitude/width);}
				else if (l==ceil(l_true)&&max(m,a)==m){m=ceil(l_true*height/width);a=floor(l_true*altitude/width);}
				else if (l==ceil(l_true)&&max(m,a)==a){m=floor(l_true*height/width);a=ceil(l_true*altitude/width);}
				else if (a==ceil(l_true*altitude/width)&&max(m,l)==m){m=ceil(l_true*height/width);l=floor(l_true);}
				else if (a==ceil(l_true*altitude/width)&&max(m,l)==l){m=floor(l_true*height/width);l=ceil(l_true);}
            	new_n=2*l*a;
			}
			if(new_n<n){
				m=ceil(l_true*height/width);					
            	a=ceil(l_true*altitude/width);
            	l=ceil(l_true);
            	new_n=2*l*a;
			}									
        }	
		else if(strcmp(air_config,"W-E")==0){				
            float m_true;
			m_true=sqrt(n*height/(altitude*2));						
            l=floor(m_true*width/height);					
            a=floor(m_true*altitude/height);
            m=floor(m_true);
            new_n=2*m*a;
			if(new_n<n){
				if (max(max(m,l),a)==m){m=ceil(m_true);l=floor(m_true*width/height);	a=floor(m_true*altitude/height);}
				else if (max(max(m,l),a)==l){m=floor(m_true);  l=ceil(m_true*width/height);	a=floor(m_true*altitude/height);}
				else{m=floor(m_true);l=floor(m_true*width/height);	a=ceil(m_true*altitude/height);}
				new_n=2*m*a;
			}
			if(new_n<n){
				if (m==ceil(m_true)&&max(l,a)==l){l=ceil(m_true*width/height);a=floor(m_true*altitude/height);}
				else if (m==ceil(m_true)&&max(l,a)==a){l=floor(m_true*width/height);a=ceil(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==m){m=ceil(m_true);a=floor(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==a){m=floor(m_true);a=ceil(m_true*altitude/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==m){m=ceil(m_true);l=floor(m_true*width/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==l){m=floor(m_true);l=ceil(m_true*width/height);}
				new_n=2*m*a;
			}
			if(new_n<n){
				l=ceil(m_true*width/height);					
            	a=ceil(m_true*altitude/height);
            	m=ceil(m_true);
            	new_n=2*m*a;
			}									
        }										 			
        else{
			all=1;
			float m_true;												
            m_true=(2*((altitude+width+height)/height)+sqrt(4*pow((altitude+width+height)/height,2.0)-2*(8-n)*(altitude*height+width*height+width*altitude)/(pow(height, 2.0))))/(2*(altitude*height+width*height+width*altitude)/(pow(height, 2.0)));
            l=floor(m_true*width/height);					
            a=floor(m_true*altitude/height);
            m=floor(m_true);
			if (l==1){new_n=a*m;}
			else if (m==1){new_n=a*l;}
			else if (a==1){new_n=m*l;}
            else {new_n=(2*m+2*l-4)*a+(l-2)*(m-2)*2;}
			if(new_n<n){
				if (max(max(m,l),a)==m){m=ceil(m_true);l=floor(m_true*width/height);	a=floor(m_true*altitude/height);}
				else if (max(max(m,l),a)==l){m=floor(m_true);  l=ceil(m_true*width/height);	a=floor(m_true*altitude/height);}
				else{m=floor(m_true);l=floor(m_true*width/height);	a=ceil(m_true*altitude/height);}
				if (l==1){new_n=a*m;}
				else if (m==1){new_n=a*l;}
				else if (a==1){new_n=m*l;}
            	else {new_n=(2*m+2*l-4)*a+(l-2)*(m-2)*2;}
			}
			if(new_n<n){
				if (m==ceil(m_true)&&max(l,a)==l){l=ceil(m_true*width/height);a=floor(m_true*altitude/height);}
				else if (m==ceil(m_true)&&max(l,a)==a){l=floor(m_true*width/height);a=ceil(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==m){m=ceil(m_true);a=floor(m_true*altitude/height);}
				else if (l==ceil(m_true*width/height)&&max(m,a)==a){m=floor(m_true);a=ceil(m_true*altitude/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==m){m=ceil(m_true);l=floor(m_true*width/height);}
				else if (a==ceil(m_true*altitude/height)&&max(m,l)==l){m=floor(m_true);l=ceil(m_true*width/height);}
				if (l==1){new_n=a*m;}
				else if (m==1){new_n=a*l;}
				else if (a==1){new_n=m*l;}
            	else {new_n=(2*m+2*l-4)*a+(l-2)*(m-2)*2;}
			}
			if(new_n<n){
				l=ceil(m_true*width/height);					
            	a=ceil(m_true*altitude/height);
            	m=ceil(m_true);
            	if (l==1){new_n=a*m;}
				else if (m==1){new_n=a*l;}
				else if (a==1){new_n=m*l;}
            	else {new_n=(2*m+2*l-4)*a+(l-2)*(m-2)*2;}
			}									
        }
        if(strcmp(air_config,"all")!=0 and verbose)     
            printf("random Problem Generator: unvalid or unspecified airspace configuration for random scenario, using default (all configuration).\n");
        
        printf("\n--------------------------------------------------------------------------------------------------------------------\n");
	    printf("Pseudo-random Problem Generator: three-dimensional matrix with dimensions (%i, %i, %i) and %i cells, %i of them not used\n",(int)m,(int)l,(int)a,new_n,new_n-n);
	}
	// 2.---------- Generate initial positions (x_0,y_0,z_0), one in each cell
	float horizontal_step=width/l;
	float vertical_step=height/m;
	float altitude_step=altitude/a;  //if 2D, altitude_step = 0/1 = 0
	if( (horizontal_step< D && strcmp(air_config,"W-E")!=0)||
	    (vertical_step< D) && strcmp(air_config,"N-S")!=0 ||
        ((altitude_step< D) && strcmp(air_config,"U-D")!=0 && is3D==true) ){
            if(is3D==false)
		        {D=min(horizontal_step,vertical_step);}
            else
                {D=min(horizontal_step,min(vertical_step,altitude_step));} 
		printf("random Problem Generator: safety distance D too small, using %f\n",D);
	}
	int cont=0; // counts the number of aircraft for which we have alrady generated (x_0,y_0,z_0)
	int excess_n=new_n-n; // if more cells than n some of them will be skipped
	int skip_xN[(int)a];
    int skip_xS[(int)a];
    int skip_yW[(int)a];
    int skip_yE[(int)a];
	int skip_xU[(int)l];
    int skip_xD[(int)l];
    // randomly generated the cell that will be skipped if needed for each altitude step
    for (int j=0; j<a; ++j){
        skip_xN[j]=randomInt(0,l-1);
        skip_xS[j]=randomInt(0,l-1);
        skip_yW[j]=randomInt(0,m-1);
        skip_yE[j]=randomInt(0,m-1);
    }
	for (int j=0; j<l; ++j){
    	skip_xU[j]=randomInt(0,m-1);
    	skip_xD[j]=randomInt(0,m-1);
	}
	//in order to record if we generate N,W,S,E faces
    int generate_northface=0;
    int generate_southface=0;
    int generate_westface=0;
    int generate_eastface=0;

	bool sep_ok;
	int u;
	int north_final=a; //row of the north face from which no more cells will be skipped
	int west_final=a; //row of the west face from which no more cells will be skipped
	int south_final=a; //row of the south face from which no more cells will be skipped

	
    if(strcmp(air_config,"W-E")!=0&&strcmp(air_config,"W-U")!=0&&strcmp(air_config,"U-D")!=0){ // north face
		u=0;
		generate_northface=1;
		if (excess_n==0){north_final=0;}    // we won't skip any other cell
        for (int j = 0; j < a; ++j){ 		// we place aircraft starting from the bottom of the face
            for (int i = 0; i < l; ++i){	// we place aircraft starting from the left of the face
                if(excess_n!=0 && skip_xN[j]==i){ // in row j, we skip column i
                	excess_n--;
                    continue;
                }
				if(cont>=n){excess_n=0; continue;} // we already generated n aircraft
				if(excess_n==0){u+=1; if(u==1 && north_final>0){north_final = j+1;}} // u=1 at the first row j s.t. no cell must be skipped (will be needed for the separation below)
				
                sep_ok=false;
                while(!sep_ok){// make sure that this aircraft and the previous are separated
                    x_0[cont]=randomFloat( horizontal_step*i, horizontal_step*(i+1) );
                    y_0[cont]=randomFloat( height-vertical_step, height );
                    z_0[cont]=randomFloat( j*altitude_step, altitude_step*(j+1) );
					if (is3D==false)
						{if(i==0) sep_ok=true;
						else sep_ok=areSepatered(cont,cont-1);}
					else{ // 3D configuration
						if(i==0) { //first column of the face
							if(j==0){sep_ok=true;}  //first row of the face (from bottom to up)
							else { //check if the aircraft is separated with the one in the cell BELOW (no aircraft on the PREVIOUS cell since i==0)
								if(u>=2){sep_ok=areSepatered(cont,cont-l);} //no cell has been skipped at the previous j
								else {sep_ok=areSepatered(cont,cont-l+1);} //a cell has been skipped at the previous j
							} 
						}
						else {
							if(j==0){sep_ok=areSepatered(cont,cont-1);} //check if the aircraft is separated with the one in the PREVIOUS cell
							else { //check if the aircraft is separated with the one in the PREVIOUS cell AND the one in the cell BELOW 
								//we must distinguish among different cases depending on the position of the cell
								if(u>=2 || u==1 && i>= skip_xN[j-1]|| i<skip_xN[j] && i>= skip_xN[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l));}
								else if (u==1 && i < skip_xN[j-1]|| i>skip_xN[j] && i>= skip_xN[j-1] || i<skip_xN[j] && i<skip_xN[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l+1));}
								else if (i>skip_xN[j] && i<= skip_xN[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l+2));}
							} 
						}
					}
                }
                // heading angle bounds depends on the cell
                if(i==0){ // left corner of the north face
                    hamin[cont]=-M_PI/2;
                    hamax[cont]=0;
                    if(j==0){ //left bottom
                        phimin[cont]=0;
                        phimax[cont]=M_PI/2;}
                    else if (j==a-1){ //left up
                        phimin[cont]=M_PI/2;
                        phimax[cont]=M_PI;}
                    else{ //left middle
                        phimin[cont]=0;
                        phimax[cont]=M_PI;}            
                }
                else if(i==l-1){// right corner of the north face
                    hamin[cont]=-M_PI;
                    hamax[cont]=-M_PI/2;
                    if(j==0){ //right bottom
                        phimin[cont]=0;
                        phimax[cont]=M_PI/2;}
                    else if (j==a-1){ //right up
                        phimin[cont]=M_PI/2;
                        phimax[cont]=M_PI;}
                    else{  //right middle
                        phimin[cont]=0;
                        phimax[cont]=M_PI;}
                }
                else{// middle cell of the north face
                    hamin[cont]=-M_PI;
                    hamax[cont]=0;
                    if(j==0){ // bottom
                        phimin[cont]=0;
                        phimax[cont]=M_PI/2;}
                    else if (j==a-1){ // up
                        phimin[cont]=M_PI/2;
                        phimax[cont]=M_PI;}
                    else{  // middle
                        phimin[cont]=0;
                        phimax[cont]=M_PI;}
                }
            cont++;
            }
        }
		if (verbose){cout<<"Number of aircraft on the north face: "<<cont<<"\n";}
    }
    if(strcmp(air_config,"N-S")!=0&&strcmp(air_config,"N-U")!=0&&strcmp(air_config,"U-D")!=0){ //west face
        int k=cont;
        u=0;
        generate_westface=1;
        if (excess_n==0){west_final=0;}
        for (int j = 0; j < a; ++j){ 
            for (int i = 0; i < m; ++i){
                if(l==1&&generate_northface==1){continue;}
                if(excess_n!=0 && skip_yW[j]==i){ 
                    excess_n--;
                    continue;
                }
                if(cont>=n){excess_n=0; continue;}
                if(excess_n==0){u+=1; if(u==1 && west_final>0){west_final = j+1;}} 
                if (i==0 && generate_northface==1) continue; // first column already generated

                sep_ok=false;
                while(!sep_ok){
                    x_0[cont]=randomFloat( 0, horizontal_step );
                    y_0[cont]=randomFloat( vertical_step*(m-i-1), vertical_step*(m-i) );
                    z_0[cont]=randomFloat( j*altitude_step, altitude_step*(j+1) );
                    
                    if (is3D==false){
                        if(i==0) sep_ok=true;
                        else if(i==1 && cont>1)	sep_ok=areSepatered(cont,0);
                        else sep_ok=areSepatered(cont,cont-1);
                    }
                    else{
                        if (generate_northface==0){ //first face to be generated
                            if(i==0){ 
                                if(j==0){sep_ok=true;}
                                else {
                                    if(u>=2){sep_ok=areSepatered(cont,cont-m);} //no cell has been skipped at the previous j
                                    else {sep_ok=areSepatered(cont,cont-m+1);} // a cell has been skipped at the previous j
                                }
                            }
                            else{ //i>0
                                if (j==0){sep_ok=areSepatered(cont,cont-1);} //check the PREVIOUS cell
                                else{   //check the PREVIOUS cell and the cell BELOW
                                    if(u>=2 || (u==1 && i>= skip_yW[j-1]) || (i<skip_yW[j] && i>= skip_yW[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m));}
                                    else if ((u==1 && i< skip_yW[j-1]) || (i>skip_yW[j] && i>= skip_yW[j-1]) || (i<skip_yW[j] && i<=skip_yW[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+1));}
                                    else if (i>skip_yW[j] && i<=skip_yW[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+2));}
                                }
                            }
                        }
                        else{ // generate_northface==1
                            if (m<=2){  
                                if (j==0){sep_ok=(areSepatered(cont,0));}
                                else{sep_ok=(areSepatered(cont,cont-m+1)&&areSepatered(cont,(j+1)*l-1-min(north_final,j+1)));}
                            }
                            else{
                                if (i==1){ //check the cell at the same altitude on the left corner of north face and the cell BELOW
                                    if (j==0){sep_ok=(areSepatered(cont,0));}
                                    else{
                                        if(u>=2 || ( u==1 && i>= skip_yW[j-1] ) || ( i<skip_yW[j] && i>= skip_yW[j-1])){sep_ok=(areSepatered(cont,(j+1)*l-1-min(north_final,j+1))&&areSepatered(cont,cont-m+1));}
                                        else if ((u==1 && i< skip_yW[j-1] ) || ( i>skip_yW[j] && i>= skip_yW[j-1] ) || ( i<skip_yW[j] && i<skip_yW[j-1])){ sep_ok=(areSepatered(cont,(j+1)*l-1-min(north_final,j+1))&&areSepatered(cont,cont-m+2));}
                                        else if (i>skip_yW[j] && i <= skip_yW[j-1]){sep_ok=(areSepatered(cont,(j+1)*l-1-min(north_final,j+1))&&areSepatered(cont,cont-m+3));}
                                    
                                    } 
                                }
                                else{ //i>1 
                                    if (j==0){sep_ok=areSepatered(cont,cont-1);}
                                    else{ //check the PREVIOUS cell and the cell BELOW
                                        if(u>=2 || (u==1 && i>= skip_yW[j-1] ) || ( i<skip_yW[j] && i>= skip_yW[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+1));}
                                        else if ((u==1 && i< skip_yW[j-1] ) || ( i>skip_yW[j] && i>= skip_yW[j-1] )|| (i<skip_yW[j] && i<skip_yW[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+2));}
                                        else if (i>skip_yW[j] && i<= skip_yW[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+3));}
                                    }
                                }
                            }
                        }
                    }
                }
                // heading angle bounds depends on the cell
                if(i==0){ // north corner of west face
                    hamin[cont]=-M_PI/2;
                    hamax[cont]=0;
                    if(j==0){ //north bottom
                        phimin[cont]=0;
                        phimax[cont]=M_PI/2;}
                    else if (j==a-1){ //north up
                        phimin[cont]=M_PI/2;
                        phimax[cont]=M_PI;}
                    else{ //north middle
                        phimin[cont]=0;
                        phimax[cont]=M_PI;}
                }
                else if(i==m-1){// south corner of the west face
                    hamin[cont]=0;
                    hamax[cont]=M_PI/2;
                    if(j==0){ //south bottom
                        phimin[cont]=0;
                        phimax[cont]=M_PI/2;}
                    else if (j==a-1){ //south up
                        phimin[cont]=M_PI/2;
                        phimax[cont]=M_PI;}
                    else{  //south middle
                        phimin[cont]=0;
                        phimax[cont]=M_PI;}
                }
                else{// middle cell of the west face
                    hamin[cont]=-M_PI/2;
                    hamax[cont]=M_PI/2;
                    if(j==0){ // bottom
                        phimin[cont]=0;
                        phimax[cont]=M_PI/2;}
                    else if (j==a-1){ // up
                        phimin[cont]=M_PI/2;
                        phimax[cont]=M_PI;}
                    else{  // middle
                        phimin[cont]=0;
                        phimax[cont]=M_PI;}
                }
                cont++;
            }
        }
		if (verbose){cout<<"Number of aircraft on the west face: "<<cont-k<<"\n";}
    }
    if(strcmp(air_config,"W-E")!=0 && strcmp(air_config,"W-N")!=0 && strcmp(air_config,"N-U")!=0 && strcmp(air_config,"W-U")!=0&&strcmp(air_config,"U-D")!=0){// south face
		int k=cont;
		u=0;
		generate_southface=1;
		if (excess_n==0){south_final=0;}
		for (int j = 0; j < a; ++j){ 
			for (int i = 0; i < l; ++i){
				if(m==1&&generate_westface==1){continue;}
				if(excess_n!=0 && skip_xS[j]==i){
					excess_n--;
					continue;
				}
				if (cont>=n){excess_n=0; continue;}				
				if(excess_n==0){u+=1; if(u==1&&south_final>0) {south_final=j+1;}}
				if (i==0 && generate_westface==1 ) continue; 
				
                sep_ok=false;
				while(!sep_ok){
					x_0[cont]=randomFloat( horizontal_step*i, horizontal_step*(i+1) );
					y_0[cont]=randomFloat( 0, vertical_step );
					z_0[cont]=randomFloat( j*altitude_step, altitude_step*(j+1) );
					if(is3D==false){
                        if(i==0) sep_ok=true;
						else sep_ok=areSepatered(cont,cont-1);
					}
					else{
						if(generate_westface==0 && generate_northface==1){ // just north face generated, we don't ignore any column - > 'N-S' configuration
							if(i==0){
								if(j==0){sep_ok=true;}
								else {//check if the aircraft is separated with the one in the cell BELOW
									if(u>=2){sep_ok=areSepatered(cont,cont-l);} //no cell has been skipped at the previous j
									else {sep_ok=areSepatered(cont,cont-l+1);} // a cell has been skipped at the previous j
								}
							}
							else{
								if (j==0){sep_ok=areSepatered(cont,cont-1);}
								else {//check if the aircraft is separated with the one in the cell BELOW *and* the one in the PREVIOUS cell
									if(u>=2 || (u==1 && i>= skip_xS[j-1]) || ( i<skip_xS[j] && i>= skip_xS[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l));}
									else if ((u==1 && i < skip_xS[j-1]) || ( i>skip_xS[j] && i>= skip_xS[j-1] ) || ( i<skip_xS[j] && i<skip_xS[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l+1));}
									else if (i>skip_xS[j] && i<= skip_xS[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l+2));}
								}
							}
						}
						else if(generate_westface==1 && generate_northface==1){  //'all' configuration
							//index of the last aircraft placed on the north face = l*a-north_final-1 (until j<north_final one cell is skipped for each j)
							if (l<=2){if (j==0){sep_ok=areSepatered(cont,(l*a-north_final-1)+(m-1)-min(west_final,j+1));}
										else {sep_ok=(areSepatered(cont,(l*a-north_final-1)+(j+1)*(m-1)-min(west_final,j+1))&&areSepatered(cont,cont-l+1));}
									}
							else{							
								if (i==1){ //check the cell at the same altitude on the south corner of west face and the cell BELOW
									if (j==0){ sep_ok=(areSepatered(cont,(l*a-north_final-1)+(m-1)-min(west_final,j+1)));}
									else{ 
										if(u>=2 || (u==1 && i>= skip_xS[j-1] ) || ( i<skip_xS[j] && i>= skip_xS[j-1])){  sep_ok=(areSepatered(cont,(l*a-1-north_final)+(j+1)*(m-1)-min(west_final,j+1))&&areSepatered(cont,cont-l));}
										else if ((u==1 && i< skip_xS[j-1] ) || ( i>skip_xS[j] && i>= skip_xS[j-1] ) || ( i<skip_xS[j] && i<skip_xS[j-1])){ sep_ok=(areSepatered(cont,(l*a-1-north_final)+(j+1)*(m-1)-min(west_final,j+1))&&areSepatered(cont,cont-l+1));}
										else if (i>skip_xS[j] && i<= skip_xS[j-1]){ sep_ok=(areSepatered(cont,(l*a-1-north_final)+(j+1)*(m-1)-min(west_final,j+1))&&areSepatered(cont,cont-l+2));}
									} //cont+1-l*a-m
								}
								else
								{
									if (j==0){sep_ok=areSepatered(cont,cont-1);}
									else {//check if the aircraft is separated with the one in the PREVIOUS cell AND the one in the cell BELOW 
										if(u>=2 || (u==1 && i>= skip_xS[j-1]) || ( i<skip_xS[j] && i>= skip_xS[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l+1));}
										else if ((u==1 && i < skip_xS[j-1]) || ( i>skip_xS[j] && i>= skip_xS[j-1] ) || ( i<skip_xS[j] && i<skip_xS[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l+2));}
										else if (i>skip_xS[j] && i<=skip_xS[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l+3));}
									}
								}
							}
						}
					}
				}

				// heading angle bounds depends on the cell
				if(i==0){// left corner of south face
					hamin[cont]=0;
					hamax[cont]=M_PI/2;
					if(j==0){ //left bottom
						phimin[cont]=0;
						phimax[cont]=M_PI/2;}
					else if (j==a-1){ //left up
						phimin[cont]=M_PI/2;
						phimax[cont]=M_PI;}
					else{ //left middle
						phimin[cont]=0;
						phimax[cont]=M_PI;}
				}
				else if(i==l-1){// right corner of south face
					hamin[cont]=M_PI/2;
					hamax[cont]=M_PI;
					if(j==0){ //right bottom
						phimin[cont]=0;
						phimax[cont]=M_PI/2;}
					else if (j==a-1){ //right up
						phimin[cont]=M_PI/2;
						phimax[cont]=M_PI;}
					else{ //right middle
						phimin[cont]=0;
						phimax[cont]=M_PI;}
				}
				else{// middle cell of south face
					hamin[cont]=0;
					hamax[cont]=M_PI;
					if(j==0){ // bottom
						phimin[cont]=0;
						phimax[cont]=M_PI/2;}
					else if (j==a-1){ // up
						phimin[cont]=M_PI/2;
						phimax[cont]=M_PI;}
					else{ //middle
						phimin[cont]=0;
						phimax[cont]=M_PI;}
				}
				cont++;
			}
		}
		if (verbose){cout<<"Number of aircraft on the south face: "<<cont-k<<"\n";}
    } 
	if(strcmp(air_config,"W-N")!=0 && strcmp(air_config,"N-S")!=0&&strcmp(air_config,"W-U")!=0 && strcmp(air_config,"N-U")!=0&&strcmp(air_config,"U-D")!=0){// east face
		int k=cont;
		generate_eastface=1;
		u=0;
		for (int j = 0; j < a; ++j){
			for (int i = 0; i < m; ++i){//
				if (m<=2&&generate_northface==1){ continue; }
				else{                                        
					if(excess_n!=0 && skip_yE[j]==i){
						excess_n--;//                                                       
						continue;
					}
					if (cont>=n){excess_n=0; continue;}
					if(excess_n==0){u+=1;}
					if (i==0 && generate_southface==1 ) continue; 
					if (i==m-1 && generate_northface==1 ) continue; 

					sep_ok=false;
					while(!sep_ok){
						x_0[cont]=randomFloat( width-horizontal_step, width );
						y_0[cont]=randomFloat( vertical_step*i, vertical_step*(i+1) );
						z_0[cont]=randomFloat( altitude_step*j, altitude_step*(j+1) );
						if (is3D==false){
							if(i==0) sep_ok=true;
							else if (generate_northface==1&&i==m-2) {sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,l-1));}
							else {sep_ok=areSepatered(cont,cont-1);}
						}
						else{
                            if(generate_southface==0 && generate_northface==0){  //only west has been already generated, we don't ignore any row -> 'W-E' configuration
                                if(i==0){
                                    if(j==0){sep_ok=true;}
                                    else {//check if the aircraft is separated with the one in the cell BELOW
                                        if(u>=2){sep_ok=areSepatered(cont,cont-m);} //no cell has been skipped at the previous j
                                        else {sep_ok=areSepatered(cont,cont-m+1);} // a cell has been skipped at the previous j
                                    }
                                }
                                else{
                                    if (j==0){sep_ok=areSepatered(cont,cont-1);}
                                    else {//check if the aircraft is separated with the one in the cell BELOW *and* the one in the PREVIOUS cell
                                        if(u>=2 || (u==1 && i>= skip_yE[j-1]) || ( i<skip_yE[j] && i>= skip_yE[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m));}
                                        else if ((u==1 && i < skip_yE[j-1]) || ( i>skip_yE[j] && i>= skip_yE[j-1] ) || ( i<skip_yE[j] && i<skip_yE[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+1));}
                                        else if (i>skip_yE[j] && i<= skip_yE[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+2));}
                                    }
                                }
                            }
                            else if(generate_westface==1 && generate_northface==1 && generate_southface==1){  //'all' configuration
                                if (m==3) { //we have just one column for i=1, and we must check separation among the aircraft and BOTH previous AND next cell (and the cell below)
                                    if (j==0){sep_ok=(areSepatered(cont,l*a-north_final-1+(m-1)*a-west_final+(j+1)*(l-1)-min(south_final,j+1))&&areSepatered(cont,(j+1)*l-1-min(north_final,j+1)));}
                                    else{sep_ok=(areSepatered(cont,l*a-north_final-1+(m-1)*a-west_final+(j+1)*(l-1)-min(south_final,j+1))&&areSepatered(cont,(j+1)*l-1-min(north_final,j+1))&&areSepatered(cont,cont-1));}
                                }
                                if (i==1 ){ //check the cell at the same altitude on the right corner of south face and the cell BELOW
                                //index of the last aircraft placed on the west face= (l*a-1-north_final)+(m-1)*a-west_final 
                                    if (j==0){sep_ok=(areSepatered(cont,(l*a-1-north_final+(m-1)*a-west_final)+(j+1)*(l-1)-min(south_final,j+1)));}
                                    else{
                                        if(u>=2 || (u==1 && i>= skip_yE[j-1] ) || ( i<skip_yE[j] && i>= skip_yE[j-1])){sep_ok=(areSepatered(cont,(l*a+(m-1)*a-1-north_final-west_final)+(j+1)*(l-1)-min(south_final,j+1))&&areSepatered(cont,cont-m+2));}
                                        else if ((u==1 && i< skip_yE[j-1] ) || ( i>skip_yE[j] && i>= skip_yE[j-1] ) || ( i<skip_yE[j] && i<skip_yE[j-1])){sep_ok=(areSepatered(cont,(l*a+(m-1)*a-1-north_final-west_final)+(j+1)*(l-1)-min(south_final,j+1))&&areSepatered(cont,cont-m+3));}
                                        else if (i>skip_yE[j] && i<= skip_yE[j-1]){sep_ok=(areSepatered(cont,(l*a+(m-1)*a-1-north_final-west_final )+(j+1)*(l-1)-min(south_final,j+1))&&areSepatered(cont,cont-m+4));}
                                    } //cont+1-l*a-m
                                }
                                else if (1<i<m-2)
                                {
                                    if (j==0){sep_ok=areSepatered(cont,cont-1);}
                                    else {//check if the aircraft is separated with the one in the PREVIOUS cell AND the one in the cell BELOW 
                                        if(u>=2 || (u==1 && i>= skip_yE[j-1]) || ( i<skip_yE[j] && i>= skip_yE[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+2));}
                                        else if ((u==1 && i < skip_yE[j-1]) || ( i>skip_yE[j] && i>= skip_yE[j-1] ) || ( i<skip_yE[j] && i<skip_yE[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+3));}
                                        else if (i>skip_yE[j] && i<= skip_yE[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-m+4));}
                                    }
                                }
                                else { //i=m-1
                                    if (j==0){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,(j+1)*l-1-min(north_final,j+1)));}
                                    else {//check if the aircraft is separated with the one in the PREVIOUS cell, the one in the NEXT cell, AND the one in the cell BELOW 
                                        if(u>=2 || (u==1 && i>= skip_yE[j-1]) || ( i<skip_yE[j] && i>= skip_yE[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,(j+1)*l-1-min(north_final,j+1))&&areSepatered(cont,cont-m+2));}
                                        else if ((u==1 && i < skip_yE[j-1]) || ( i>skip_yE[j] && i>= skip_yE[j-1] ) || ( i<skip_yE[j] && i<skip_yE[j-1])){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,(j+1)*l-1-min(north_final,j+1))&&areSepatered(cont,cont-m+3));}
                                        else if (i>skip_yE[j] && i<= skip_yE[j-1]){sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,(j+1)*l-1-min(north_final,j+1))&&areSepatered(cont,cont-m+4));}
                                    }
                                }
                                
                            }
						}
					}
					
					// heading angle bounds depends on the cell
					if(i==0){ //  north corner of east face
						hamin[cont]=M_PI/2;
						hamax[cont]=M_PI;
						if(j==0){ //north bottom
							phimin[cont]=0;
							phimax[cont]=M_PI/2;}
						else if (j==a-1){ //north up
							phimin[cont]=M_PI/2;
							phimax[cont]=M_PI;}
						else{ //north middle
							phimin[cont]=0;
							phimax[cont]=M_PI;}
					}
					else if(i==m-1){// south corner of east face
						hamin[cont]=-M_PI;
						hamax[cont]=-M_PI/2;
						if(j==0){ //south bottom
							phimin[cont]=0;
							phimax[cont]=M_PI/2;}
						else if (j==a-1){ //south up
							phimin[cont]=M_PI/2;
							phimax[cont]=M_PI;}
						else{  //south middle
							phimin[cont]=0;
							phimax[cont]=M_PI;}
					}
					else{// middle cell of east face
						hamin[cont]=M_PI/2;
						hamax[cont]=3*M_PI/2;
						if(j==0){ // bottom
							phimin[cont]=0;
							phimax[cont]=M_PI/2;}
						else if (j==a-1){ // up
							phimin[cont]=M_PI/2;
							phimax[cont]=M_PI;}
						else{ // middle
							phimin[cont]=0;
							phimax[cont]=M_PI;}
					}
					cont++;
				}
			}
		}
		if (verbose){cout<<"Number of aircraft on the east face: "<<cont-k<<"\n";}
    }
    if(is3D==true&&strcmp(air_config,"W-N")!=0 && strcmp(air_config,"N-S")!=0&&strcmp(air_config,"W-E")!=0){ //upper face
		int k=cont;
		for (int j = 0; j < m; ++j){
			for (int i = 0; i < l; ++i){
				if ((m<=2 && generate_northface==1 && generate_southface==1) || (l<=2 && generate_westface==1 && generate_eastface==1)){continue;} //no space for other aircraft
				else{                                          
					if(excess_n!=0 && skip_xU[j]==i){//                                         
						excess_n--;//                                                       
						continue;
					}
					if (cont>=n){excess_n=0; continue;}
					if(i==0 && generate_westface==1 ) continue; //&& skip_yW[(int)(a-1)!=j]
					if(i==l-1 && generate_eastface==1) continue; // && skip_yW[(int)(a-1)!=i]
					if(j==0 && generate_southface==1 ) continue;
					if(j==m-1 && generate_northface==1 ) continue;
					
					sep_ok=false;
					while(!sep_ok){
						x_0[cont]=randomFloat(horizontal_step*i, horizontal_step*(i+1) );
						y_0[cont]=randomFloat( vertical_step*j, vertical_step*(j+1) );
						z_0[cont]=randomFloat(altitude-altitude_step,altitude );
						//******todo add the separations for the corners
						if(i==0 && j==0) sep_ok=true;
						else sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l));
					} 
					// heading angle bounds depends on the cell
					phimin[cont]=M_PI/2;
					phimax[cont]=M_PI;
					if(i==0){ //  south corner of upper face
						if(j==0){ // south left
							hamin[cont]=0;
							hamax[cont]=M_PI/2;}
						else if (j==l-1){ //south right 
							hamin[cont]=M_PI/2;
							hamax[cont]=M_PI;}
						else{ //south middle
							hamin[cont]=0;
							hamax[cont]=M_PI;}
					}
					else if(i==m-1){// north corner of upper face
						if(j==0){ // north left
							hamin[cont]=-M_PI/2;
							hamax[cont]=0;}
						else if (j==l-1){ //north right 
							hamin[cont]=-M_PI;
							hamax[cont]=-M_PI/2;}
						else{ //north middle
							hamin[cont]=-M_PI;
							hamax[cont]=0;}
					}
					else{// middle cell of upper face
						if(j==0){ // left
							hamin[cont]=-M_PI/2;
							hamax[cont]=M_PI/2;}
						else if (j==l-1){ // right
							hamin[cont]=M_PI/2;
							hamax[cont]=3*M_PI/2;}
						else{ // middle
							hamin[cont]=0;
							hamax[cont]=2*M_PI;}
					}
					cont++;
				}
			}
		}
		if (verbose){cout<<"Number of aircraft on the upper face: "<<cont-k<<"\n";}
    }
    if(is3D==true&&strcmp(air_config,"W-N")!=0 && strcmp(air_config,"N-S")!=0&&strcmp(air_config,"W-E")!=0&& strcmp(air_config,"N-U")!=0&& strcmp(air_config,"W-U")!=0){ //bottom face    
		int k=cont;
		for (int j = 0; j < m; ++j){
            for (int i = 0; i < l; ++i){
				if ((m<=2 && generate_northface==1 && generate_southface==1) || (l<=2 && generate_westface==1 && generate_eastface==1)){continue;} //no space for other aircraft				                                
                else{     
					if(excess_n!=0 && skip_xD[j]==i){
                        excess_n--;                                       
                        continue;
					}
					if (cont>=n){excess_n=0; continue; }
					if(i==0 && generate_westface==1 ) continue; 
					if(i==l-1 && generate_eastface==1 ) continue; 
					if(j==0 && generate_southface==1 ) continue; 
					if(j==m-1 && generate_northface==1 ) continue; 

                    sep_ok=false;
                    while(!sep_ok){
                        x_0[cont]=randomFloat(horizontal_step*i, horizontal_step*(i+1) );
                        y_0[cont]=randomFloat(vertical_step*j, vertical_step*(j+1) );
                        z_0[cont]=randomFloat(0,altitude_step);
                        //******todo add the separations for the corners
						if(i==0 &&j==0) sep_ok=true;
                        else sep_ok=(areSepatered(cont,cont-1)&&areSepatered(cont,cont-l));
                    }

					
                    // heading angle bounds depends on the cell
                    phimin[cont]=0;
                    phimax[cont]=M_PI/2;
                    if(i==0){ //  south corner of bottom face
                        if(j==0){ // south left
                            hamin[cont]=0;
                            hamax[cont]=M_PI/2;}
                        else if (j==l-1){ //south right 
                            hamin[cont]=M_PI/2;
                            hamax[cont]=M_PI;}
                        else{ //south middle
                            hamin[cont]=0;
                            hamax[cont]=M_PI;}
                    }
                    else if(i==m-1){// north corner of bottom face
                        if(j==0){ // north left
                            hamin[cont]=-M_PI/2;
                            hamax[cont]=0;}
                        else if (j==l-1){ //north right 
                            hamin[cont]=-M_PI;
                            hamax[cont]=-M_PI/2;}
                        else{ //north middle
                            hamin[cont]=-M_PI;
                            hamax[cont]=0;}
                    }
                    else{// middle cell of bottom face
                        if(j==0){ // left
                            hamin[cont]=-M_PI/2;
                            hamax[cont]=M_PI/2;}
                        else if (j==l-1){ // right
                            hamin[cont]=M_PI/2;
                            hamax[cont]=3*M_PI/2;}
                        else{ // middle
                            hamin[cont]=0;
                            hamax[cont]=2*M_PI;}
                    }
                    cont++;
                }
            }
		}
		if (verbose){cout<<"Number of aircraft on the bottom face: "<<cont-k<<"\n";}
    }
	if (verbose){
        if(is3D==false){
		    printf("coordinates generated in a rectangle %f x %f\n",width,height);
		    for (int i = 0; i < n; ++i)
			    printf("p_%i=(%f,%f)\n",i,x_0[i],y_0[i]);
        }
        else{
            printf("coordinates generated in a parallelepiped %f x %f x %f\n",width,height,altitude);
		    for (int i = 0; i < n; ++i)
			    printf("p_%i=(%f,%f,%f)\n",i,x_0[i],y_0[i],z_0[i]);
        }

	}	
}
// This auxiliary method generates the vectors of velocity of aircraft for a random scenario
void pseudoRandomP_vectors(int nc, float pc,int maxc, float* hamin,float* hamax, float* phimin, float* phimax){
	float random;
	int total_conflicts=0;
	bool* explored=new bool[n];// vector to indicate if velocity was generated for an aircraft or not
	for (int i = 0; i < n; ++i)
		explored[i]=false;
	bool* tried=new bool[n];
	for (int i = 0; i < n; ++i)
		tried[i]=false;
	int total_explored=0;
	int i;
	int conflicts_i;
	float vi, thetai, phii;
	int num_trials, num_loops, max_conflicts; 
	int max_trials;
	if (is3D==true){max_trials=1000;}
    else {max_trials=1000;}
	bool meet_requested; // true if we are close to the end and still do not have nc
	bool zeroIniConf; // true if the first value generated for conflicts_i is 0
	
	while(total_explored < n){
		// A. randomly select an unexplored aircraft i:
		i=randomInt(0,n-1);
		while(explored[i]) i=randomInt(0,n-1);
		// B. randomly chose velocity and heading angles for i:
		vi=randomFloat(vmin,vmax);
		thetai=randomFloat(hamin[i],hamax[i]);
		if (is3D==true)
			{phii=randomFloat(phimin[i],phimax[i]);}
		else {phii=M_PI/2;}
		vx[i]=vi*cos(thetai)*sin(phii);
		vy[i]=vi*sin(thetai)*sin(phii);
		vz[i]=vi*cos(phii);
		// C. i will have at least 1 conflict with probability pc (if i is not the first for which we generate (vx,vy,vz) and
		// we have not reached the maximum number of conflicts yet):
		max_conflicts=min(min(maxc,total_explored),max(0,nc-total_conflicts)); // max num of conflicts we can generate
		random = ((float) rand()) / (float) RAND_MAX;
		meet_requested=ceil( (1.0*(nc-total_conflicts))/ ((maxc+1)/2) )>=n-total_explored; 

		if(total_explored>0 && !meet_requested && total_conflicts<nc && random < pc) {
			conflicts_i=randomInt(1,max_conflicts);
			zeroIniConf=false;
			if (verbose)
				printf("conflicts_i generated randomly\n");
		}
		else if (meet_requested){  // we are close to the end and still do not have nc
			conflicts_i=max_conflicts;
			zeroIniConf=false;
			if (verbose)
				printf("conflicts_i generated to meet requested\n");
		}
		else {
			conflicts_i=0;
			zeroIniConf=true;
			if (verbose){
				printf("conflicts_i set to 0\n");
				printf("%f  < %f, nc=%i \n",random, pc,nc );
			}
		}
		// Loop: generate random trajectories until the desired number of conflicts are yielded
		if (verbose){
			if(is3D==true){
				printf("conflicts already generated: %i, current aircraft i= %i conflicts_i=%i, p0=(%f,%f,%f), v=(%f,%f,%f)\n",
				total_conflicts, i,conflicts_i,x_0[i],y_0[i],z_0[i],vx[i],vy[i],vz[i]);}
			else{
				printf("conflicts already generated: %i, current aircraft i= %i conflicts_i=%i, p0=(%f,%f), v=(%f,%f)\n",
				total_conflicts, i,conflicts_i,x_0[i],y_0[i],vx[i],vy[i]);}
		}
		num_trials=0; // number of times we try to generate a vector with specific numConflicts
		num_loops=0; // number of different values we try for numConflicts
		while(num_trials==0){
			if (verbose){
				if(is3D==true){
					printf("\t trying with num_conflicts= %i,  v=(%f,%f,%f)...\n",conflicts_i,vx[i],vy[i],vz[i] );}
				else{
					printf("\t trying with num_conflicts= %i,  v=(%f,%f)...\n",conflicts_i,vx[i],vy[i]);}
			}
			while((num_trials< max_trials) && (numConflicts(i,explored)!=conflicts_i)){
				vi= randomFloat(vmin,vmax);
				thetai=randomFloat(hamin[i],hamax[i]);
				if (is3D==true)
					{phii=randomFloat(phimin[i],phimax[i]);}
				else {phii=M_PI/2;}
				vx[i]=vi*cos(thetai)*sin(phii);
				vy[i]=vi*sin(thetai)*sin(phii);
				vz[i]=vi*cos(phii);
				num_trials++;
			}
			if(numConflicts(i,explored)!=conflicts_i){ // try with another number of conflicts
				num_trials=0;
				tried[conflicts_i]=true;
				if(num_loops<max_conflicts){
					if(meet_requested && conflicts_i>0)
						conflicts_i--;
					else if(!zeroIniConf && num_loops==max_conflicts-1)
						conflicts_i=0;
					else{
						conflicts_i=randomInt(1,max_conflicts);
						while(tried[conflicts_i])
							conflicts_i=randomInt(1,max_conflicts);
					}
					num_loops++;
				}
				else{
					num_loops++;
					conflicts_i=num_loops;
				}
			}
			else{ // trajectory with conflicts_i generated successfully
				num_trials++;
				for (int k = 0; k < n; ++k)
					tried[k]=false;
			}
		}
		// D. Set speed and angle of i, and update loop variables
		hat_v[i]=vi;
		hat_theta[i]=thetai;
		hat_phi[i]=phii;
		total_conflicts +=conflicts_i;
		total_explored++;
		explored[i]=true;
		if(total_explored!=n)
			pc=(4.0*(nc-total_conflicts))/((n-total_explored)*(1+maxc));
	}
	if(nc!=total_conflicts){
		printf("----------------------------------------------------------------------------------------------\n");
		printf("Random Problem Generator: WARNING deviation from requested number of conflicts of %i\n",total_conflicts-nc );
		if (nc < total_conflicts){
			printf("We could not generate a scenario with the requested congestion, you can try to:\n");
			printf("\ti)increase nc or maxc\n");
			printf("\tii) increase window size\n");
			printf("\tiii) decrease n\n");
			printf("\tiv) change random seed by using -seed option\n");
		}
		else{
			printf("We could not generate a scenario with the requested congestion, you can try to:\n");
			printf("\ti)decrease nc or maxc\n");
			printf("\tii) decrease window size\n");
			printf("\tiii) increase n\n");
			printf("\tiv) change random seed by using -seed option\n");
		}
	}
	
	delete[] explored;
	delete[] tried;
}

// ####       5. PSEUDO-RANDOM SCENARIOS     ###
// pseudoRandomP: generates scenario with initial positions near the borders of a parallelepiped P of size height x width x altitude
// heading angles are generated randomly, depending on the initial positions of aircraft (in order to make sure that
// aircraft are flying crossing the parallelepiped P). There are different configurations for the initial positions,
// since aircraft can be on the 6 faces of P or only on 2 of them.
void pseudoRandomP(char* air_config,float height,float width,float altitude,int nc, float pc,int maxc){
	float* hamin=new float[n]; // vectors to keep min and max heading angle theta for each aircraft depending on its initial position
	float* hamax=new float[n];
	float* phimin=new float[n]; // vectors to keep min and max heading angle phi for each aircraft depending on its initial position
	float* phimax=new float[n];

	pseudoRandomP_iniPos(air_config,height,width,altitude,hamin,hamax,phimin,phimax); // generates initial positions for aircraft
	pseudoRandomP_vectors(nc,pc,maxc,hamin,hamax,phimin,phimax);// generate velocity vectors V_i i=1...n 
	// calculate ending points of the trajectories:
	float t;
	for (int i = 0; i < n; ++i){
		t=getInstantBoundary(i,height,width,altitude); // get instant in wich i will cross boundary (space window)
        x_t[i]=min(x_0[i]+t*vx[i],width);
		y_t[i]=min(y_0[i]+t*vy[i],height);
		z_t[i]=min(z_0[i]+t*vz[i],altitude);
	}
	
	delete[] hamin;
	delete[] hamax;
}

// ####       6. RANDOM SCENARIOS     ###
// randomP: generates random initial positions in 3D parallelepiped (if z_max !=0); vectors also generated randomly
void randomP(float x_max,float y_max,float z_max){	
	bool sep_ok;
	for (int i = 0; i < n; ++i){
		sep_ok=false;
		while(!sep_ok){
			x_0[i]=randomFloat(0,x_max);
			y_0[i]=randomFloat(0,y_max);
			z_0[i]=randomFloat(0,z_max);
			sep_ok=true;
			for (int j = 0; j < i; ++j){
				sep_ok=areSepatered(i,j);
				if(!sep_ok) break;
			}
		}
	}
	for (int i = 0; i < n; ++i){
		hat_v[i]= randomFloat(vmin,vmax);
		hat_theta[i]=randomFloat(-M_PI,M_PI);
		vx[i]=hat_v[i]*cos(hat_theta[i]);
		vy[i]=hat_v[i]*sin(hat_theta[i]);
		vz[i]=0;
		if(is3D){
			hat_phi[i]=randomFloat(0,M_PI);
			vx[i]=vx[i]*sin(hat_phi[i]);
			vy[i]=vy[i]*sin(hat_phi[i]);
			vz[i]=hat_v[i]*cos(hat_phi[i]);
		}	
	}
	float t;
	for (int i = 0; i < n; ++i){
		t=getInstantBoundary(i,x_max,y_max,z_max);// get instant in wich i will cross boundary (space window)
		x_t[i]=min(x_0[i]+t*vx[i],x_max);
		y_t[i]=min(y_0[i]+t*vy[i],y_max);
		z_t[i]=min(z_0[i]+t*vz[i],z_max);
	}
    printf("\n--------------------------------------------------------------------------------------------------\n");
    printf("Random Problem Generator");    
}


// --------------------------------------------------------------------------
// -----                       MAIN PROGRAM                           -------
// --------------------------------------------------------------------------
int main(int argc, char **argv) {
	char* air_config=new char[5]; // random problem: airspace configuration (all,N-N,N-S,W-E,W-N)
	char* char_mode=new char[5]; // when the mode of the generator is described by a code instead of a number
	srand(random_seed);
	// circle/sphere problem:
	float radius=200; // circle/sphere problem: radius 
	float sector_ini_x=0; double sector_angle_x=2*M_PI; //circle/sphere problem: sector where aircraft will be placed
	float sector_ini_z=0; double sector_angle_z=M_PI;
	// rhomboidal problem:
	float alpha_r=M_PI/6; // crossing angle 
    int horiz=2; int verti=2; //  number of horizontal/sloping trails
    // polyhedral problem:
	double* alpha; double* beta;
	int* mx; int* my; int* mz;
    float d_HP=3*D; float d_VP=3*D;
    // rhomboidal AND polyhedral
	int nx=1; int ny=1; int nz=1;// number of aircraft at each trail
	float dx=3*D; float dy=3*D; float dz=3*D; //  distance between each trail
	float d_aircraft=2*D; //  distance between aircraft on the same trail
	// random circle and rhomboidal problems:
	float hamin=-M_PI/6; float hamax=M_PI/6; // angle incr limit 
	// random and pseudo-random scenarios:
	int nc=-1; //  total number of conflicts between different pairs of aircraft
	float pc=-1; // probability that one aircraft have a conflict with at least one other aircraft
	int maxc=-1; //  max number of conflicts that a fixed aircraft can have with others
	float height=400; float width=400; float altitude=400; //space window size


	// 1..............Read input options
	int i=1;
    int j;
	int insert_HP=0;
	int insert_VP=0;
	int insert_alpha=0;
	int insert_beta=0;
	int insert_mx=0;
	int insert_my=0;
	int insert_mz=0;
	int insert_width=0;
	int insert_height=0;
	int insert_altitude=0;
	char* key;
    while(i<argc){
        j=2;
		key=argv[i];
      	if(strcmp(key,"-V")!=0 && i+1>=argc) {printf("Invalid syntax\n");return -1;}
      	// general input parameters:
	    if(strcmp(key,"-n")==0){
	    	n=atoi(argv[i+1]);
	    }
	    else if(strcmp(key,"-mode")==0) {
	    	char_mode=argv[i+1];
	    	if(strcmp(char_mode,"CP")==0)// circle
	    		mode=0;
	    	else if(strcmp(char_mode,"RCP")==0)// random circle
	    		mode=1;
	    	else if(strcmp(char_mode,"RP")==0)// rhomboidal
	    		mode=2;
	    	else if(strcmp(char_mode,"RRP")==0)// random rhomboidal
	    		mode=3;
	    	else if(strcmp(char_mode,"GP")==0) // grid
	    		mode=4;
	    	else if(strcmp(char_mode,"RGP")==0)// random grid
	    		mode=5;
	    	else if(strcmp(char_mode,"PR2")==0)// pseudo-random 2D (num conflicts according to user)
	    		mode=6;
			else if(strcmp(char_mode,"PR3")==0){// pseudo-random 3D (num conflicts according to user)
	    		mode=14;
	    		is3D=true;
	    	}
	    	else if(strcmp(char_mode,"R2")==0)// random2D
	    		mode=7;
			else if(strcmp(char_mode,"R3")==0){// random3D
	    		mode=15;
	    		is3D=true;
	    	}
	    	else if(strcmp(char_mode,"SP")==0){// sphere
	    		mode=8;
	    		is3D=true;
	    	}
	    	else if(strcmp(char_mode,"RSP")==0){// random sphere
	    		mode=9;
	    		is3D=true;
	    	}
            else if(strcmp(char_mode,"PL")==0){// polyhedral
	    		mode=10;
	    		is3D=true;
	    	}
			else if(strcmp(char_mode,"RPL")==0){// random polyhedral
	    		mode=11;
	    		is3D=true;
	    	}
			else if(strcmp(char_mode,"QP")==0){// cubic
	    		mode=12;
	    		is3D=true;
	    	}
			else if(strcmp(char_mode,"RQP")==0){// random cubic
	    		mode=13;
	    		is3D=true;
	    	}
	    	else{// other scenarios
	    		mode=atoi(char_mode);
	    		if (mode>=8) 
	    			is3D=true;
	    	}
	    }
	    else if(strcmp(key,"-seed")==0) {
	    	srand(atoi(argv[i+1])); 
	    }
	    else if(strcmp(key,"-vmin")==0){
	    	vmin=atof(argv[i+1]); 
	    }
	    else if(strcmp(key,"-vmax")==0){
	    	vmax=atof(argv[i+1]); 
	    }
	    else if(strcmp(key,"-vdefault")==0){
	    	vmin=atof(argv[i+1]); 
	    	vmax=vmin;
	    }
	    else if(strcmp(key,"-D")==0){
	    	D=atof(argv[i+1]);
	    	d_HP=3*D;  
			d_VP=3*D;
   			dx=3*D;
   			dy=3*D; 
   			dz=3*D; 
			d_aircraft=2*D; 
	    }
	    else if(strcmp(key,"-V")==0){
	    	verbose=true; 
	    	i=i-1;
	    }
	    else if(strcmp(key,"-f")==0){
	    	outputFile=argv[i+1];
	    }
	    // random problem :
	    else if(strcmp(key,"-nc")==0)
	    	nc=atoi(argv[i+1]); 
	    else if(strcmp(key,"-pc")==0) 
	    	pc=atof(argv[i+1]); 
		else if(strcmp(key,"-maxc")==0) 
	    	maxc=atoi(argv[i+1]); 
	    else if(strcmp(key,"-airconfig")==0) 
			air_config=argv[i+1]; 
	    else if(strcmp(key,"-h")==0) 
	     	{height=atof(argv[i+1]);
			insert_height=1;}
	    else if(strcmp(key,"-w")==0) 
	     	{width=atof(argv[i+1]);
			insert_width=1;}
		else if(strcmp(key,"-a")==0) 
	     	{altitude=atof(argv[i+1]);
			insert_altitude=1;}
	    //polyhedral problems
        else if(strcmp(key,"-HP")==0){
	    	HP=atoi(argv[i+1]);
			alpha = new double[HP]; // crossing angle diagonal trails on horizontal planes
    		mx = new int[HP];  my = new int[HP]; //number of trails on each horizontal plane (mx,my)
			for (int h=0;h<HP;++h)
				{mx[h]=2; my[h]=2; alpha[h]=M_PI/6;}
			insert_HP=1; //to record that the user is inserting HP 
		}
	    else if(strcmp(key,"-VP")==0){
	    	VP=atoi(argv[i+1]);
			beta = new double[VP]; // crossing angle diagonal trails on vertical planes
			mz= new int[VP];      //number of trails on each vertical plane 
			for (int v=0;v<VP;++v)
				{mz[v]=2; beta[v]=2*M_PI/3;}
			insert_VP=1; //to record that the user is inserting VP
		}
		else if(strcmp(key,"-dHP")==0) 
	    	d_HP=atoi(argv[i+1]);  
		else if(strcmp(key,"-dVP")==0) 
	    	d_VP=atoi(argv[i+1]);
		else if(strcmp(key,"-beta")==0){
			insert_beta=1; //to record that the user is inserting beta
			if (insert_VP==0){VP=vp; beta = new double[VP];} //if the user does not insert VP
			int k=i+1;
			int b=0;
			while(k<argc && argv[k][0]!='-'){ //until the user does not insert a new parameter OR stop giving inputs
				if (b<VP)
					{beta[b]=atof(argv[k]);
					k=k+1;
					b=b+1;}
				else  //if the user insert more than VP values for beta, we must store the number of values we must ignore
					{k=k+1;
					b=b+1;}
			}

			if (b<VP){ //if the user insert less than VP beta values we set the others to the default value
				for(int v=b;v<VP;++v)
					beta[v]=2*M_PI/3;					
				printf("\nInserted less than VP=%i values for beta. The last %i are set to default=%f.\n",VP,VP-b,beta[b]);
			}
			else if(b>VP){  //if the user insert more than VP beta values, we ignore the last ones
				printf("\nInserted more than VP=%i values for beta. The last %i values are ignored.\n",VP,b-VP);}
			j=k-i; //when updating i, we will set i=i+j =i+k-i =k
		}
		else if(strcmp(key,"-betad")==0) { 
			insert_beta=1; //to record that the user is inserting beta
			if (insert_VP==0){VP=vp; beta = new double[VP];} //if the user does not insert VP
			int k=i+1;
			int b=0;
			while(k<argc && argv[k][0]!='-'){ //until the user does not insert a new parameter OR stop giving inputs
				if (b<VP)
					{beta[b]=atof(argv[k]);
					beta[b]=M_PI*beta[b]/180;
					k=k+1;
					b=b+1;}
				else
					{k=k+1;
					b=b+1;}
			}	
			
			if (b<VP){ //if the user insert less than VP beta values we set the others to the default value
				for(int v=b;v<VP;++v)
					beta[v]=2*M_PI/3;
				printf("\nInserted less than VP=%i values for beta. The last %i are set to default=%f.\n",VP,VP-b,beta[b]);
			}
			else if(b>VP){ //if the user insert more than VP beta values, we ignore the last ones
				printf("\nInserted more than VP=%i values for beta. The last %i values are ignored.\n",VP,b-VP);}
			j=k-i;
	    }
		else if(strcmp(key,"-mz")==0){
			insert_mz=1; //to record that the user is inserting mz
			if (insert_VP==0){VP=vp;  mz = new int[VP];} //if the user does not insert VP
			int k=i+1;
			int b=0;
			while(k<argc && argv[k][0]!='-'){
				if (b<VP)
					{mz[b]=atoi(argv[k]);
					k=k+1;
					b=b+1;}
				else
					{k=k+1;
					b=b+1;}
			}
			if (b<VP){ //if the user insert less than VP mz values we set the others to the default value
				for(int v=b;v<VP;++v)
					mz[v]=2;
				printf("\nInserted the number of vertical trails only for %i vertical planes. The last %i are set to default=%i.\n",b,VP-b,mz[b]);
			}
			else if(b>VP){ //if the user insert more than VP mz values, we ignore the last ones
					printf("\nInserted %i vertical trails, but there are only VP=%i vertical planes. The last %i values are ignored.\n",b,VP,b-VP);
			}

			j=k-i;
		}
        // rhomboidal AND polyhedral problems:
	    else if(strcmp(key,"-nx")==0) 
	    	nx=atoi(argv[i+1]);  
	    else if(strcmp(key,"-ny")==0) 
	    	ny=atoi(argv[i+1]);  
        else if(strcmp(key,"-nz")==0) 
	    	nz=atoi(argv[i+1]);  
	    else if(strcmp(key,"-dx")==0) 
	    	dx=atof(argv[i+1]);  
	    else if(strcmp(key,"-dy")==0) 
	    	dy=atof(argv[i+1]);
        else if(strcmp(key,"-dz")==0) 
	    	dz=atof(argv[i+1]);
		else if(strcmp(key,"-d_aircraft")==0) 
	    	d_aircraft=atof(argv[i+1]);  
	    else if(strcmp(key,"-alpha")==0){
			insert_alpha=1; //to record that the user is inserting alpha
			if (insert_HP==0){HP=hp; alpha = new double[HP];} //the user does not insert HP
			if (mode==2 || mode==3 || mode ==4 || mode==5) //rhomboidal instance
				alpha_r=atoi(argv[i+1]); 
			else{ //pholyhedral instance
				int k=i+1;
				int a=0;
				while(k<argc && argv[k][0]!='-'){ //until the user does not insert any other parameter OR stop giving inputs
	    			if (a<HP)
						{alpha[a]=atof(argv[k]);
						k=k+1;
						a=a+1;}
					else //the user could insert more than HP values for alpha: we must remember how many values to ignore 
						{k=k+1;
						a=a+1;}
				}
				if (a<HP){ //if the user insert less than HP alpha values we set the others to the default value
					for(int h=a;h<HP;++h)
						alpha[h]=M_PI/6;
					printf("\nInserted less than HP=%i values for alpha. The last %i are set to default=%f.\n",HP,HP-a,alpha[a]);
				}
				else if(a>HP){ //if the user insert more than HP alpha values, we ignore the last ones
					printf("\nInserted more than HP=%i values for alpha. The last %i values are ignored.\n",HP,a-HP);
				}
				
				j=k-i;
			}
		}
	    else if(strcmp(key,"-alphad")==0){
			insert_alpha=1; //to record that the user is inserting alpha
			if (insert_HP==0){HP=hp; alpha = new double[HP];}//the user does not insert HP
			if (mode==2 || mode==3 || mode ==4 || mode==5) //rhomboidal instance
				{alpha_r=atoi(argv[i+1]); 
				alpha_r=M_PI*alpha_r/180;}
			else{
				int k=i+1;
				int a=0;
				while(k<argc && argv[k][0]!='-' ){ //until the user does not insert any other parameter OR stop giving inputs
	    			if (a<HP)
						{alpha[a]=atof(argv[k]);
						alpha[a]=M_PI*alpha[a]/180;
						k=k+1;
						a=a+1;}
					else
						{k=k+1;
						a=a+1;}
				}

				if (a<HP){ //if the user insert less than HP alpha values we set the others to the default value
					for(int h=a;h<HP;++h)
						alpha[h]=M_PI/6;
					printf("\nInserted less than HP=%i values for alpha. The last %i are set to default=%f.\n",HP,HP-a,alpha[a]);
				}
				else if(a>HP){ //if the user insert more than HP alpha values, we ignore the last ones
					printf("\nInserted more than HP=%i values for alpha. The last %i values are ignored.\n",HP,a-HP);
				}
				
				j=k-i;
			}
	    }
        else if(strcmp(key,"-mx")==0){
			insert_mx=1; //to record that the user is inserting mx
			if (insert_HP==0){HP=hp; mx = new int[HP];} //if the user does not insert HP
			if (mode==2 || mode==3 || mode ==4 || mode==5) //rhomboidal instance
				horiz=atoi(argv[i+1]);
			else{		//polyhedral instance
				int k=i+1;
				int a=0;
				while(k<argc && argv[k][0]!='-'){
				if (a<HP)
					{mx[a]=atoi(argv[k]);
					k=k+1;
					a=a+1;}
				else
					{k=k+1;
					a=a+1;}
				}

				if (a<HP){ //if the user insert less than HP mx values we set the others to the default value
					for(int h=a;h<HP;++h)
						mx[h]=2;
					printf("\nInserted the number of horizontal trails only for %i horizontal planes. The last %i are set to default=%i.\n",a,HP-a,mx[a]);
				}
				else if(a>HP){ //if the user insert more than HP mx values, we ignore the last ones
					printf("\nInserted %i horizontal trails, but there are only HP=%i horizontal planes. The last %i values are ignored.\n",a,HP,a-HP);
				}

				j=k-i;
			}
		}
	    else if(strcmp(key,"-my")==0){
			insert_my=1; //to record that the user is inserting my
			if (insert_HP==0){HP=hp;  my = new int[HP];} //if the user does not insert HP
			if (mode==2 || mode==3 || mode ==4 || mode==5) //rhomboidal instance
				verti=atoi(argv[i+1]);
			else{   // polyhedral instance
				int k=i+1;
				int a=0;
				while(k<argc && argv[k][0]!='-'){
				if (a<HP)
					{my[a]=atoi(argv[k]);
					k=k+1;
					a=a+1;}
				else
					{k=k+1;
					a=a+1;}
				}

				if (a<HP){ //if the user insert less than HP my values we set the others to the default value
					for(int h=a;h<HP;++h)
						my[h]=2;
					printf("\nInserted the number of vertical trails only for %i horizontal planes. The last %i are set to default=%i.\n",a,HP-a,my[a]);
				}
				else if(a>HP){ //if the user insert more than HP my values, we ignore the last ones
					printf("\nInserted %i vertical trails, but there are only HP=%i horizontal planes. The last %i values are ignored.\n",a,HP,a-HP);
				}

				j=k-i;
			}
		}
	    // circle/sphere problem:
	    else if(strcmp(key,"-r")==0) 
	     	radius=atof(argv[i+1]);
	    else if(strcmp(key,"-secInix")==0) 
	     	sector_ini_x=atof(argv[i+1]);
	    else if(strcmp(key,"-secAngx")==0) 
	     	sector_angle_x=atof(argv[i+1]);
	    else if(strcmp(key,"-secInixd")==0){
	    	sector_ini_x=atof(argv[i+1]);
	    	sector_ini_x=M_PI*sector_ini_x/180;
	    } 
	    else if(strcmp(key,"-secAngxd")==0){
	    	sector_angle_x=atof(argv[i+1]);
	    	sector_angle_x=M_PI*sector_angle_x/180;
	    }
	    else if(strcmp(key,"-secIniz")==0) 
	     	sector_ini_z=atof(argv[i+1]);
	    else if(strcmp(key,"-secAngz")==0) 
	     	sector_angle_z=atof(argv[i+1]);
	    else if(strcmp(key,"-secInizd")==0){
	    	sector_ini_z=atof(argv[i+1]);
	    	sector_ini_z=M_PI*sector_ini_z/180;
	    } 
	    else if(strcmp(key,"-secAngzd")==0){
	    	sector_angle_z=atof(argv[i+1]);
	    	sector_angle_z=M_PI*sector_angle_z/180;
	    } 
	    // random circle, sphere, rhomboidal, grid, polyhedral, and cubic:
	    else if(strcmp(key,"-hamin")==0)
	    	hamin=atof(argv[i+1]); 
	    else if(strcmp(key,"-hamax")==0)
	    	hamax=atof(argv[i+1]); 
	    else
	        printf("Skipping unrecognized option: %s\n", key);
      	
		i=i+j; //either i+2 or k
    }

	// if the user does not insert HP and/or alpha and/or mx and/or my
	if (insert_HP==0){
		HP=hp;
		if(insert_alpha==0){
			alpha = new double[HP];
			for (int i=0;i<HP;++i)
				{alpha[i]=M_PI/6;}
		}
		if (insert_mx==0){
			mx = new int[HP];
			for (int i=0;i<HP;++i)
				mx[i]=2;
		}
		if (insert_my==0){
			my = new int[HP];
			for (int i=0;i<HP;++i)
				my[i]=2;
		}
	}
	// if the user does not insert VP and/or beta and/or mz
	if (insert_VP==0){
		VP=vp;
		if(insert_beta==0){
			beta = new double[VP];
			for (int i=0;i<VP;++i)
				{beta[i]=2*M_PI/3;}
		}
		if (insert_mz==0){
			mz = new int[VP];
			for (int i=0;i<VP;++i)
				mz[i]=2;
		}
	}
	if(mode==14 && insert_width==0){width=100;}
	if(mode==14 && insert_height==0){height=100;}
	if(mode==14 && insert_altitude==0){altitude=100;}
	
    // 2............Print alerts
    if(vmin<=0 || vmax <=0 || vmax < vmin){
   		printf("WARNING: At least one unvalid bound for velocity: %f, %f",vmin, vmax);
   		vmin=400;
   		vmax=400;
   		printf(", set to default: vmin=vmax=%f\n",vmin);
   		
    }
	if  ((mode==8 || mode==9)&& sector_ini_z+sector_angle_z>M_PI){
		printf("WARNING: the angle Phi must be in [0,PI]. ");
		sector_angle_z=M_PI;
		sector_ini_z=0;
		printf("The sector boundaries are set to default: Sector_init_phi=0, Sector_angle_phi=PI.");
	}
    if (mode==6 || mode==14){
    	// detect unvalid parameters:
    	if (nc>n*(n-1)/2.0) nc=-1;
    	if (pc>1) pc=-1;
    	if (maxc>n) maxc=-1;
    	// nc and pc missing:
    	if(nc==-1 && pc==-1){ 
    		pc=0.5;
    		printf("random Problem Generator: missing probability of conflicts set to 0.5 \n");  
    	}
    	// maxc missing:
    	if(maxc==-1 && (nc==-1 || pc==-1)){ 

    		if (nc==-1) maxc=n-1;
    		else if (pc==-1) maxc=round(8.0*nc/n-1); // this will yield pc=0.5
    		printf("random Problem Generator: missing max conflicts per aircraft set to %i\n",maxc);  
    	}
    	// at this point at least 2 parameters among nc, pc and maxc are not -1, we calculate the remaining one:
    	if(nc==-1)
    		nc=round(pc*(n/2.0)*((1+maxc)/2.0));
    	else if(pc==-1){
    		pc=nc*(2.0/n)*(2.0/(1+maxc));
    		if (pc>1){
    			pc=1;
				printf("WARNING random Problem Generator: saturated scenario, all aircraft having more than (maxc+1)/2 conflicts\n");
    		}
    	}
    	else if(maxc==-1)
    		maxc=round((4.0*nc)/(n*pc)-1);
    	else{
    		float new_pc=nc*(2.0/n)*(2.0/(1+maxc));
    		if (new_pc>1){
				pc=1;
				printf("WARNING random Problem Generator: saturated scenario, all aircraft having more than (maxc+1)/2 conflicts\n");
	    	}
    		else if(new_pc!=pc){
    			pc=new_pc;
				printf("WARNING random Problem Generator: unvalid probability of conflicts set to %f\n",pc);	
    		}
    	}
    	if (verbose)	
    		printf("random Problem Generator: congestion measures nc=%i, maxc=%i, pc=%f\n",nc,maxc,pc);

    }
    else if(mode==2 ||mode==3 ||mode==4 ||mode==5 ){// rhomboidal or grid problems
        n=nx*horiz+ny*verti;
        if( mode==2 || mode==3 ){// rhomboidal problem
            alpha_r=correctAngle(alpha_r);
            if(alpha_r < 0){
                alpha_r=M_PI/6;
                printf("Rhomboidal Problem generator: unvalid crossing angle, set to default=%f\n",alpha_r);
            }
        } 
    }
    else if(mode==10 ||mode==11 ||mode==12 ||mode==13){// polyhedral or cubic problems
		n=0;
		for(int i=0; i<HP; ++i){
        	n=n+nx*mx[i]+ny*my[i];}
		for(int i=0; i<VP; ++i){
			n=n+nz*mz[i];}
		
		if( mode==10 || mode==11 ){ //polyhedral
			for (int i=0; i<HP; ++i){
				alpha[i]=correctAngle(alpha[i]);
				if(alpha[i] < 0){
					alpha[i]=M_PI/6;
					printf("Unvalid crossing angle alpha[%i], set to default=%f\n",i,alpha[i]);
				}
			}
			for (int i=0; i<VP;++i){
				beta[i]=correctAngle(beta[i]);
				if(beta[i] < 0){
					beta[i]=2*M_PI/3;
					printf("Polyhderal Problem generator: unvalid crossing angle beta[%i], set to default=%f\n",i,beta[i]);
				}
			}
        }
    }

	
	// 3...............Create scenarios
	x_0=new float[n];
	y_0=new float[n];
	z_0=new float[n];
	x_t=new float[n];
	y_t=new float[n];
	z_t=new float[n];
	vx=new float[n];
	vy=new float[n];
	vz=new float[n];
	hat_v =new float[n];
	hat_theta=new float[n];
	hat_phi=new float[n];

	if(!is3D){
		for (int i = 0; i < n; ++i){
			z_0[i]=0;
			z_t[i]=0;
			vz[i]=0;
			hat_phi[i]=M_PI/2;
		}
	}

    switch(mode){// CP,RCP,RP,RRP,GP,RGP,PR2,R2,SP,RSP,PL,RPL,QP,RQP,PR3,R3
    	case 0:{CircleP(radius,sector_ini_x,sector_angle_x);break;} //CP
    	case 1:{randomCircleP(radius,sector_ini_x,sector_angle_x,hamin,hamax);break;} //RCP
    	case 2:{RomboP(alpha_r,horiz,verti,nx,ny,dx,dy,d_aircraft);break;} //RP
    	case 3:{randomRomboP(alpha_r,horiz,verti,nx,ny,dx,dy,d_aircraft,hamin,hamax);break;} //RRP
    	case 4:{GridP(horiz,verti,nx,ny,dx,dy,d_aircraft);break;} //GP
    	case 5:{randomGridP(horiz,verti,nx,ny,dx,dy,d_aircraft,hamin,hamax);break;} //RGP
    	case 6:{pseudoRandomP(air_config,height,width,0,nc,pc,maxc);break;} //PR2
		case 7:{randomP(height,width,0);break;} //R2
    	case 8:{SphereP(radius,sector_ini_x,sector_angle_x,sector_ini_z,sector_angle_z);break;} //SP
		case 9:{randomSphereP(radius,sector_ini_x,sector_angle_x,sector_ini_z,sector_angle_z,hamin,hamax);break;} //RSP
        case 10:{PolyhedralP(alpha,beta,mx,my,mz,nx,ny,nz,d_HP,d_VP,dx,dy,dz,d_aircraft);break;} //PL
		case 11:{randomPolyhedralP(alpha,beta,mx,my,mz,nx,ny,nz,d_HP,d_VP,dx,dy,dz,d_aircraft,hamin,hamax);break;} //RPL
		case 12:{Grid3dP(mx,my,mz,nx,ny,nz,d_HP,d_VP,dx,dy,dz,d_aircraft);break;}  //QP
		case 13:{randomGrid3dP(mx,my,mz,nx,ny,nz,d_HP,d_VP,dx,dy,dz,d_aircraft,hamin,hamax);break;}  //RQP
		case 14:{pseudoRandomP(air_config,height,width,altitude,nc,pc,maxc);break;} //PR3
		case 15:{randomP(height,width,altitude);break;} //R3
    }


    // 4............Print output
	printBenchmarkInfo();
    
    return 0;
}  
