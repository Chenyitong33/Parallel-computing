#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <omp.h>
#include <ctime>
#include <iostream>
using namespace std;

#define THREADS 16

// paricle
struct body {
	double m; //mass
	double x, y; //x,y position
	double vx, vy; //x,y velocity
	double ax, ay; //x,y accel
};

// barnes-hut tree
struct node {
	double totalmass;	
	double centerx, centery;
	double xmin, xmax;
	double ymin, ymax;
	double diag;
	struct body * bodyp;
	struct node * q1;
	struct node * q2;
	struct node * q3;
	struct node * q4;
};

enum quadrant { q1, q2, q3, q4 };

enum quadrant getquadrant(double x, double y, double xmin, double xmax, double ymin, double ymax)
//makes a rectangle with bounds of xmin,xmax,ymin,ymax, and returns the quadrant that (x,y) is in
{
	double midx, midy;
	
	midx = xmin + 0.5*(xmax - xmin);
	midy = ymin + 0.5*(ymax - ymin);
	
	if(y > midy)
	{
		if(x > midx)
		{
			return q1;
		} else
		{
			return q2;
		}
	} else {
		if(x > midx)
		{
			return q4;
		} else
		{
			return q3;
		}		
	}
	
}
//***************************** BH tree **********************************
//************************************************************************
struct node * createnode(struct body * bodyp, double xmin, double xmax, double ymin, double ymax)
//creates a leaf node to insert into the tree
{
	struct node * rootnode;
  
	if(!(rootnode=(node*)malloc(sizeof(struct node))))
	{
		//printf("Allocation failed, exit");
		cerr << "Allocation failed, exit!" << endl;
		return 0;
	}
	 
  //rootnode=malloc(sizeof(struct node));
	rootnode->totalmass = bodyp->m;
	rootnode->centerx = bodyp->x;
	rootnode->centery = bodyp->y;
	rootnode->xmin = xmin;
	rootnode->xmax = xmax;
	rootnode->ymin = ymin;
	rootnode->ymax = ymax;
	// calculate the diagonal
  rootnode->diag = sqrt(( pow(xmax - xmin, 2) + pow(ymax - ymin, 2) ));
		
	rootnode->bodyp = bodyp;
	rootnode->q1 = NULL;
	rootnode->q2 = NULL;
	rootnode->q3 = NULL;
	rootnode->q4 = NULL;
	
	return rootnode;
}

void updatecenterofmass(struct node * nodep, struct body * bodyp)
//updates the center of mass after inserting a point into a branch
{
	nodep->centerx = (nodep->totalmass*nodep->centerx + bodyp->m*bodyp->x)/(nodep->totalmass + bodyp->m);
	nodep->centery = (nodep->totalmass*nodep->centery + bodyp->m*bodyp->y)/(nodep->totalmass + bodyp->m);
	nodep->totalmass += bodyp->m;
	return;
}

void insertbody(struct body * insbody, struct node * nodep)
//inserts a body into the tree, converting leaf nodes into branches if necessary
{
	enum quadrant existingquad, newquad;
	double xmid, ymid;
	
	xmid = nodep->xmin + 0.5*(nodep->xmax - nodep->xmin);
	ymid = nodep->ymin + 0.5*(nodep->ymax - nodep->ymin);
		
	if(nodep->bodyp != NULL) 
	//if the node is a leaf convert to a branch by inserting the leaf point into one of its subquadrants
	{
		existingquad = getquadrant(nodep->bodyp->x, nodep->bodyp->y, nodep->xmin, nodep->xmax, nodep->ymin, nodep->ymax);
			
		switch (existingquad)
		{
			case q1:
				nodep->q1 = createnode(nodep->bodyp, xmid, nodep->xmax, ymid, nodep->ymax);
				break;
			case q2:
				nodep->q2 = createnode(nodep->bodyp, nodep->xmin, xmid, ymid, nodep->ymax);
				break;
			case q3:
				nodep->q3 = createnode(nodep->bodyp, nodep->xmin, xmid, nodep->ymin, ymid);
				break;
			case q4:
				nodep->q4 = createnode(nodep->bodyp, xmid, nodep->xmax, nodep->ymin, ymid);
				break;
		}
		nodep->bodyp = NULL;
	}
	
	newquad = getquadrant(insbody->x, insbody->y, nodep->xmin, nodep->xmax, nodep->ymin, nodep->ymax);
	
	updatecenterofmass(nodep,insbody);
	
	switch (newquad) 
  //insert the new point into one of the quadrants if empty, otherwise recurse deeper into tree
	{
	case q1:
		if(nodep->q1 == NULL)
		{		
			nodep->q1 = createnode(insbody, xmid, nodep->xmax, ymid, nodep->ymax);		
		} else {
			insertbody(insbody,nodep->q1);
		}
		break;
	case q2:
		if(nodep->q2 == NULL)
		{			
			nodep->q2 = createnode(insbody, nodep->xmin, xmid, ymid, nodep->ymax);
		} else {			
			insertbody(insbody,nodep->q2);
		}
		break;
	case q3:
		if(nodep->q3 == NULL)
		{			
			nodep->q3 = createnode(insbody, nodep->xmin, xmid, nodep->ymin, ymid);
		} else {			
			insertbody(insbody,nodep->q3);
		}
		break;
	case q4:
		if(nodep->q4 == NULL)
		{			
			nodep->q4 = createnode(insbody, xmid, nodep->xmax, nodep->ymin, ymid);
		} else {			
			insertbody(insbody,nodep->q4);
		}
		break;
	}
		
}

//sum the forces on body bodyp from points in tree with root node nodep
void sumforce(struct node * nodep, struct body * bodyp, double G, double theta)

{
	double dx, dy, r, rsqr; //x distance, y distance, distance, distance^2
	double accel;
	double a_over_r;
		
	dx = nodep->centerx - bodyp->x;
	dy = nodep->centery - bodyp->y;
	
	rsqr = pow(dx,2) + pow(dy,2);
	r = sqrt(rsqr);
	
	// r > diag/theta, where 0 < theta < 1 is an opening angle,
	if( (((r/nodep->diag) > 1/theta) || (nodep->bodyp))&&(nodep->bodyp!=bodyp) )
	{
		accel = (G * nodep->totalmass) / rsqr;  //acceleration
		a_over_r = accel/r;  
		
		bodyp->ax += a_over_r*dx;
		bodyp->ay += a_over_r*dy;		
		
	} else {
		//If not then the algorithm is recursively applied to the children, summing up the forces obtained to obtain a net force.
	  if(nodep->q1) { sumforce(nodep->q1, bodyp, G, theta); }
		if(nodep->q2) { sumforce(nodep->q2, bodyp, G, theta); }
		if(nodep->q3) { sumforce(nodep->q3, bodyp, G, theta); }
		if(nodep->q4) { sumforce(nodep->q4, bodyp, G, theta); }
	}
	return;
}

//recursively free subnodes, then free self
void destroytree(struct node * nodep)

{
	if(nodep != NULL)
	{
		if(nodep->q1 != NULL) 
		{ 
			destroytree(nodep->q1); 
		}
		if(nodep->q2 != NULL) 
		{ 
			destroytree(nodep->q2); 
		}
		if(nodep->q3 != NULL) 
		{ 
			destroytree(nodep->q3); 
		}
		if(nodep->q4 != NULL) 
		{ 
			destroytree(nodep->q4); 
		}
		free(nodep);
	}
}

//***************************** N-body problem **********************************
//*******************************************************************************

#include <mpi.h>
int  numtasks, rank; 
//********************* broadcast functions **************************************
int broadcastbody(struct body * bodies, int srcrank)
{
	MPI_Bcast(&(bodies->m),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->x),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->y),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->vx),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->vy),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
}

void broadcastbodies(int * nbodies, struct body ** bodies)
{
	MPI_Bcast(nbodies,1,MPI_INT,0,MPI_COMM_WORLD);
	if(rank!=0)
	{
		if(!(*bodies = (body*)malloc( (*nbodies)*sizeof(struct body))))
		{
			//printf("Error: rank %d failed to allocate memory for %d bodies, Abort\n", rank, *nbodies);
      cout << "Error: rank " << rank << "failed to allocate memory for " << *nbodies << endl;
			MPI_Abort(MPI_COMM_WORLD,0);
		}
	}
	for(int i = 0; i < *nbodies; i++)
	{
		broadcastbody((*bodies)+i,0);
	}
	return;
}

int broadcastaccel(struct body * bodies, int srcrank)
{
	MPI_Bcast(&(bodies->ax),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->ay),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
}

//********************* Random Initialization of the system***************************************
struct body * initbodies(const int nbodies, const double mass, const double mzero)
//randomly initialize bodies in r and theta, and give them a spin around mzero 
{
	struct body * bodies;
	
	srand(time(NULL));   // should only be called once
  
  if(!(bodies = (body*)malloc(nbodies*sizeof(struct body))))  //allocate mem
	{
		//printf("Error: failed to allocate memory for %d bodies, exit\n",nbodies);
		cerr << "Error: failed to allocate memory for " << nbodies << endl;
		return NULL;
	}

	//bodies = malloc(nbodies*sizeof(struct body));
  
  for(int i = 0; i < nbodies; i++)  //initialize with random position
	{
		bodies[i].m = mass;
		bodies[i].vx = 0.0;
		bodies[i].vy = 0.0;
		bodies[i].ax = 0.0;
		bodies[i].ay = 0.0;
		
		
		bodies[i].x = ((float)rand())/RAND_MAX * 1000.0 - 500.0;
		bodies[i].y = ((float)rand())/RAND_MAX * 1000.0 - 500.0;
		
		
		bodies[i].vx = ((float)rand())/RAND_MAX * 10.0 - 5.0 ;
		bodies[i].vy = ((float)rand())/RAND_MAX * 10.0 - 5.0 ;
		
		
		//printf("body %d: x=%f y=%f\n",i,bodies[i].x,bodies[i].y);
	}
	
	bodies[0].m = mzero;  //put bodies[0] at the center and give it a different mass
	bodies[0].x = 0.0;
	bodies[0].y = 0.0;
	bodies[0].vx = 0.0;
	bodies[0].vy = 0.0;
	
	return bodies;
	
}

//*********************************** I/ O ***************************************

int outputbodies(const char * outfilename, const struct body * bodies, const int nbodies, const double timestep)
{
	FILE * outfile; 
	remove(outfilename);
	if(!(outfile = fopen(outfilename, "w")))
	{
		return 0;
	}
	fprintf(outfile,"time step = %f\nnumber of bodies = %d\n",timestep, nbodies);
	for(int i = 0; i < nbodies; i++)
	{
		fprintf(outfile,"mass = %f pos_x = %f pos_ y = %f vel_x = %f vel_y = %f\n", bodies[i].m, bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy);
	}
	fflush(outfile);
	fclose(outfile);
	return 1;
}


//********************************************************************************
// Recalculate forces on bodies, then integrate force and velocity using 'leapfrog'.
int evolve(struct body * bodies, const int nbodies, const double timestep, const double G, const double treeratio)
{
	struct node * rootnode;
	double xmin, xmax, ymin, ymax;
	
	xmin = 0.0;
	xmax = 0.0;
	ymin = 0.0;
	ymax = 0.0;
	
	for(int i = 0; i < nbodies; i++)  //reset accel to 0
	{
		bodies[i].ax = 0.0;
		bodies[i].ay = 0.0;	
		xmin=min(xmin,bodies[i].x);
		xmax=max(xmax,bodies[i].x);
		ymin=min(ymin,bodies[i].y);
		ymax=max(ymax,bodies[i].y);
	}
	
	rootnode = createnode(bodies+0,xmin,xmax,ymin,ymax);
	
	for(int i = 1; i < nbodies; i++)
	{
		insertbody(bodies+i, rootnode);
	}
	
  	omp_set_num_threads(THREADS);
	//#pragma omp parallel
	{	
		#pragma omp parallel for
		//#pragma omp for
		for(int i = rank; i < nbodies; i+=numtasks)  //sum accel, split up n bodies into p processes
		//# pragma omp critical
    	{
			sumforce(rootnode, bodies+i, G, treeratio);
		}
		//#pragma omp parallel for
    	for(int i = 0; i < nbodies; i++)
		{
			broadcastaccel(bodies+i,i%numtasks);
		}
		
	    #pragma omp parallel for
	    //#pragma omp for
		for(int i = 0; i < nbodies; i++)
		{
			// leap frog
      		bodies[i].x += (timestep/2)*bodies[i].vx; // pos += vel * dt/2
			bodies[i].y += (timestep/2)*bodies[i].vy;
      
      		bodies[i].vx += bodies[i].ax * timestep; // vel += accel * dt
			bodies[i].vy += bodies[i].ay * timestep;
			
			bodies[i].x += (timestep/2)*bodies[i].vx; // pos += vel * dt/2
			bodies[i].y += (timestep/2)*bodies[i].vy;
		} 		
	}
	destroytree(rootnode);
	return 0;
}

//Simulation loop,  increments the timer and runs timesteps.
int simulate(struct body * bodies, const int nbodies, const int stepnum, const double timestep, const double G, const double treeratio, const char * outfile)
{
	int run = 1;
  
  //# pragma omp parallel for
  for (int i=1;i<=stepnum;i++)
  {
		evolve(bodies, nbodies, timestep, G, treeratio);
		
		outputbodies(outfile, bodies, nbodies, timestep);
     
    if (i == stepnum)
      run = 0;
    
		MPI_Bcast(&run,1,MPI_INT,0,MPI_COMM_WORLD); 
		
	} 
	
	return 1;
	
}

int main(int argc, char * argv[])
{
	struct body * bodies;
	int nbodies;
	int stepnum;
	double timestep;
	double G = 6.67384e-11;
	double mass = 1000;
	double mzero = 1000000; //1e6 
	double treeratio = 0.25; // angle between 0 and 1
	//char * infile = NULL;
	char outfile[20];

	double start_time, end_time;
	// MPI stuff
	int provided;
	/* Initialize for MPI */
	if (MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided) != MPI_SUCCESS) {
		
		cerr << "MPI initialization error" << endl;
		exit(EXIT_FAILURE);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("nbody: MPI init on rank %d of %d\n",rank, numtasks); 
	
  
	if (argc < 5)
	{
	if(rank==0) 
	   
	{
		cerr << "Need to pass in number of bodies, steps, dt and outputfile" << endl 
		<< "  nbody <nbodies> <steps> <dt> <output>  " << endl
		<< "  for example './nbody 1000 400 0.1 output1.txt 2000 500 0.1 output2.txt'  " << endl;
	}
	MPI_Finalize();
	  
	return -1;
  }
	int num_regions = (argc-1)/4;
  
	//clock_t begin = clock(); // Timer before all
	//#pragma omp parallel for private(nbodies, stepnum, timestep, outfile) num_threads ( THREADS )
	for(int region=0;region<num_regions;region++){
	// scan the arguments
	sscanf(argv[region*4+1],"%i",&nbodies);
	sscanf(argv[region*4+2],"%i",&stepnum);
	sscanf(argv[region*4+3],"%lf",&timestep);
	sscanf(argv[region*4+4],"%[^\t]",&outfile);
	
	
	// initialization of the system
	if(rank==0)
	{
		
	bodies = initbodies(nbodies, mass, mzero);
		
    
	cout << "N bodies = " << nbodies << endl;
	cout << "Time step (dt) = " << timestep << endl;
	cout << "Theta for BH tree = " << treeratio << endl;
	cout << "Universal gravitational factor = " << G << endl;
	cout << "Simplified mass for particles = " << mass << endl;
	cout << "Center particle mass = " << mzero << endl;
	}
	
	
	broadcastbodies(&nbodies,&bodies);
	
	//************************SIMULATION***********************
	//// 1 - PRE-SIMULATION
  

	if(rank==0) 
	{
	start_time = MPI_Wtime();
	cout << "SIMULATION BEGIN\n";
	}

	// 2 - RUN SIMULATION
	simulate(bodies, nbodies, stepnum, timestep, G, treeratio, outfile);

	// 3 - POST-SIMULATION
	if(rank==0) 
	{
	cout << "SIMULATION END\n";
	end_time = MPI_Wtime();
	cout << "Simulation completed in " << end_time - start_time << " seconds.\n";  
	}
	
	free(bodies);
	}// End for regions
	MPI_Finalize();

	//clock_t end = clock(); // End up all
	//double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//cout << "The total execution time is" <<  time_spent << endl;
	return 0;
}




