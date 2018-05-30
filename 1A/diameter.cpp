#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <unistd.h>
using namespace std;

#define MAX 10000
#define NOT_CONNECTED -1
#define THREADS 4

int distance[MAX][MAX];

int nodesCount;

void Initialize(){
    for (int i=0;i<MAX;++i){
        for (int j=0;j<MAX;++j){
            distance[i][j]=NOT_CONNECTED;

        }
        distance[i][i]=0;
    }
}

void floyd_warshall() {
    int i, j, k;
    //Floyd-Warshall
    for (k=1;k<=nodesCount;++k){
	   
        omp_set_num_threads(THREADS);
	
        #pragma omp for nowait
        for (i=1;i<=nodesCount;++i){
            if (distance[i][k]!=NOT_CONNECTED){
                for (j=1;j<=nodesCount;++j){
                    if (distance[k][j]!=NOT_CONNECTED && (distance[i][j]==NOT_CONNECTED || distance[i][k]+distance[k][j]<distance[i][j])){
                        distance[i][j]=distance[i][k]+distance[k][j];
			printf("Thread %d has completed iteration %d %d.\n", omp_get_thread_num(), i, j);
			printf("Thread number: %d.\n", omp_get_num_threads());
                    }
                }
            }
	//sleep(1);
        }

    }
 }

int main(){
    double start, stop;
    Initialize();
   
    scanf("%d", &nodesCount);
     
    int a, b, c; 
    while(scanf("%d %d %d", &a, &b, &c)!= EOF){
        if ( a > nodesCount || b > nodesCount){
            printf("Vertex index out of boundary.");
            return -1;
        }
        distance[a][b]=c;
    }
    //start parallelized algorithm
    start = omp_get_wtime();

    //Floyd-Warshall
    //floyd_warshall();
    omp_set_num_threads(THREADS);
    #pragma omp parallel
    {
    for (int k=1;k<=nodesCount;++k){
        
        #pragma omp for nowait
        for (int i=1;i<=nodesCount;++i){
            if (distance[i][k]!=NOT_CONNECTED){
                for (int j=1;j<=nodesCount;++j){
                    if (distance[k][j]!=NOT_CONNECTED && (distance[i][j]==NOT_CONNECTED || distance[i][k]+distance[k][j]<distance[i][j])){
                        distance[i][j]=distance[i][k]+distance[k][j];
			//printf("Thread %d has completed iteration %d %d.\n", omp_get_thread_num(), i, j);
			//printf("Thread number: %d.\n", omp_get_num_threads());
                    }
                }
            }
	//sleep(1);
        }
      
    }
    }

    int diameter=-1;

    //look for the most distant pair
    for (int i=1;i<=nodesCount;++i){
        for (int j=1;j<=nodesCount;++j){
            if (diameter<distance[i][j]){
                diameter=distance[i][j];
                //printf("%d-%d-%d\n", i, diameter, j);
            }
        }
    }

    printf("%d\n", diameter);

    /* all threads done */
    stop = omp_get_wtime();
    printf("Runtime is %f!\n",stop - start);
    return 0;

}
