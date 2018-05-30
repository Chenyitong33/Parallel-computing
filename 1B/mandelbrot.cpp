#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>

/* Constants for message-passing */
#define WORK_TAG    1           /* "work" message (master to worker) */
#define DATA_TAG    2           /* "data" message (worker to master) */
#define STOP_TAG    3           /* "stop" message (master to worker) */

// return 1 if in set, 0 otherwise
int inset(double real, double img, int maxiter){
	double z_real = real;
	double z_img = img;
	for(int iters = 0; iters < maxiter; iters++){
		double z2_real = z_real*z_real-z_img*z_img;
		double z2_img = 2.0*z_real*z_img;
		
		z_real = z2_real + real;
		z_img = z2_img + img;
		if(z_real*z_real + z_img*z_img > 4.0) 
			return 0;
	}
	return 1;
}

// count the number of points in the set, within the region
int mandelbrotSetCount(double real_lower, double real_upper, double img_lower, double img_upper, int num, int maxiter){
	int count=0;
	double real_step = (real_upper-real_lower)/num;
	double img_step = (img_upper-img_lower)/num;
	#pragma omp parallel for reduction(+:count) schedule(dynamic,1)
	for(int real=0; real<=num; real++){
		for(int img=0; img<=num; img++){
			count+=inset(real_lower+real*real_step,img_lower+img*img_step,maxiter);
		}
	}
	return count;
}

/*
 * ---- Program for worker process ---- 
 *
 */
void worker_pgm(int myID, int width, int height, double real_min, 
	double real_max, double imag_min, double imag_max, int maxiter) {

	MPI_Status status;
	int the_row;
	double scale_real, scale_imag;
	int *data_msg = (int *)malloc((width) * sizeof(int));
	
	
	scale_real = (real_max - real_min) / width;
	scale_imag = (imag_max - imag_min) / height;
	
	/* While master is still sending "work" (not "stop") messages .... */

	/* Receive message and check tag */
	while ( ((MPI_Recv(&the_row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
		&status)) == MPI_SUCCESS) && (status.MPI_TAG == WORK_TAG)) {

	/* Calculate points and return to master */
  
	for (int col = 0; col < width; ++col) {
      
		int inside = 1; //initialized, return 1 if in set, 0 otherwise
  
		/* Scale display coordinates to actual region */
		double real = real_min + ( col * scale_real);
		double img = imag_min + ( the_row * scale_imag);
      
      
		double z_real = real;
		double z_img = img;
		double real_tmp = 0;
		double img_tmp = 0;
		for(int iters = 0; iters < maxiter; iters++){
			//double z2_real = z_real*z_real-z_img*z_img;
			//double z2_img = 2.0*z_real*z_img;
      		
			real_tmp = z_real*z_real-z_img*z_img+real;
			img_tmp = 2.0*z_real*z_img+img;
			if (z_real == real_tmp && z_img == img_tmp) {
				inside = 1;
				break;
			}
			z_real = real_tmp;
			z_img = img_tmp;
			if(z_real*z_real + z_img*z_img > 4.0) { 
				inside = 0;
				break;
			}
         
		}                                  
  
	data_msg[col] = inside;//inset(real_step, img_step, maxiter);
	}
  
	MPI_Send(data_msg, width, MPI_INT, 0, DATA_TAG,
		MPI_COMM_WORLD);
	}

  /* END OF WORKER*/
  free(data_msg);
}

// main
int main(int argc, char *argv[]){
	double real_lower;
	double real_upper;
	double img_lower;
	double img_upper;
	int num;
	int maxiter;
	int num_regions = (argc-1)/6;
	
	int nprocs;
	int myid;
	
	/* Initialize for MPI */
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		fprintf(stderr, "MPI initialization error\n");
		exit(EXIT_FAILURE);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (nprocs < 2) {
		fprintf(stderr, "Number of processes must be at least 2\n");
		MPI_Finalize(); exit(EXIT_FAILURE);
	}
	
	#pragma omp parallel for private(real_lower, real_upper, img_lower, img_upper, num, maxiter)
	for(int region=0;region<num_regions;region++){
		// scan the arguments
	sscanf(argv[region*6+1],"%lf",&real_lower);
	sscanf(argv[region*6+2],"%lf",&real_upper);
	sscanf(argv[region*6+3],"%lf",&img_lower);
	sscanf(argv[region*6+4],"%lf",&img_upper);
	sscanf(argv[region*6+5],"%i",&num);
	sscanf(argv[region*6+6],"%i",&maxiter);		
				
	if (myid == 0) {
	int this_row, next_row;
	double start_time, end_time;
	int *data_msg = (int *)malloc((num) * sizeof(int));
	MPI_Status status;
	int tasks_not_done;
	int id;
	int count = 0;
    	
	/* Start timing  */
	/* could use MPI_Wtime but get_time is what sequential code uses */
	start_time = MPI_Wtime(); 
    	/* Set up for dynamic task assignment */
	next_row = 0;          /* next row to assign to a worker */
	tasks_not_done = 0;    /* count of workers still working */
	#pragma omp parallel for num_threads(nprocs-1)
	/* Send each worker first row to work on */
	for (int p = 0; p < (nprocs-1); ++p) {
		MPI_Send(&next_row, 1, MPI_INT, p+1, WORK_TAG, MPI_COMM_WORLD);
		++next_row;
		++tasks_not_done;
	}
    
	/* Receive results from workers and calculate points */

	while (tasks_not_done > 0) {
		#pragma omp critical
		/* Receive a row from a worker */
		{
		MPI_Recv(data_msg, num, MPI_INT, MPI_ANY_SOURCE,
			DATA_TAG, MPI_COMM_WORLD, &status);
		}  
		--tasks_not_done;
		id = status.MPI_SOURCE;
  
		/* More rows? */
		if (next_row < num) {
 			#pragma omp critical
			{
			/* If so, give this worker another row to work on */
			MPI_Send(&next_row, 1, MPI_INT, id, WORK_TAG, MPI_COMM_WORLD);
			++next_row;
			++tasks_not_done;
			}
		} else {
			#pragma omp critical
			{
			/* Otherwise shut this worker down */
			MPI_Send(&next_row, 0, MPI_INT, id, STOP_TAG, MPI_COMM_WORLD);
			}
		}
  
		/* Count received data */
		#pragma omp parallel for reduction(+:count) schedule(dynamic,1)
		for (int col = 0; col < num; ++col) {
			count += data_msg[col];
		}
	}

	/* End timing  */
	end_time = MPI_Wtime();

	/* Produce text output  */
	fprintf(stdout, "\n");
	fprintf(stdout, "A Mandelbrot set\n");
	fprintf(stdout, "number of worker processes = %d\n", nprocs-1);
	fprintf(stdout, "execution time in seconds = %g\n", end_time - start_time);
	fprintf(stdout, "count = %d\n", count);
	fprintf(stdout, "\n");

	free(data_msg);
	} else {
	worker_pgm(myid, num, num, real_lower, real_upper, img_lower, img_upper, maxiter);
	}// End if-else for master-slave 
	
	}// End for region
	/* Finish up */
	MPI_Finalize();
	return EXIT_SUCCESS;
}
