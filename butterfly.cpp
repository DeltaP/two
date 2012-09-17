/*
 * Gregory Petropoulos
 *
 * This is my butterfly all reduce program for assignment 2
 *
 * To compile:  mpicxx -g -Wall -std=c99 -o butterfly butterfly.cpp
 * To run:  mpiexec -n 2 ./butterfly
 *          <> -> mandatory
 *          [] -> optional
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

// -----------------------------------------------------------------
// Ends the program on an error and prints message to node 0
void cleanup (int my_rank, const char *message) {
  if (my_rank==0) {printf("%s\n",message);}
  MPI_Finalize();                                       /* kills mpi                            */
  exit(0);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// my version of the all reduce
int BUTTERFLY_Reduce(double *vector, double *result, int count, string datatype, string operation, string comm, int my_rank, bool lth, bool verbose) {
  int i, condition, talk;
  double addtoresult[100];

  for(i = 0; i < 100; i++) {
    result[i] = vector[i];
  }

  if (lth == true) {                                    /* if the ordering is low to high       */
    condition = 1;
    while (condition < comm_sz ) {
      talk = condition ^ my_rank;
      if (talk < my_rank) {
        cout << "Condition " << condition << "|  Communication pair " << talk << " and " << my_rank << endl;
        //MPI_Send(vector, 100, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        //MPI_Recv(addtoresult, 100, MPI_DOUBLE, talk, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(i = 0; i < 100; i++) {
          result[i] += addtoresult[i];
        }
      }
      else {
        //MPI_Recv(addtoresult, 100, MPI_DOUBLE, talk, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Send(result, 100, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        for(i = 0; i < 100; i++) {
          result[i] = vector[i] + addtoresult[i];
        }
      }
      condition << 1;
    }
  }
  else {
    condition = comm_sz;
    while (condition < 1) {                             /* if the ordering is high to low       */
      talk = condition ^ my_rank;
      if (talk < my_rank) {
        cout << "Condition " << condition << "|  Communication pair " << talk << " and " << my_rank << endl;
        //MPI_Send(result, 1, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        //MPI_Recv(addtoresult, 1, MPI_DOUBLE, talk, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        result += addtoresult;
      }
      else {
        //MPI_Recv(addtoresult, 100, MPI_DOUBLE, talk, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Send(result, 100, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        for(i = 0; i < 100; i++) {
          result[i] += addtoresult[i];
        }
      }
      condition >> 1;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Checks if integer is a power of two
int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// the main program
int main(int argc, char *argv[]) {
  double vector[100], result[100];
  int i, my_rank, comm_sz;
  bool verbose = false, lth;
  string datatype = "MPI_DOUBLE", operation = "MPI_SUM", comm = "MPI_COMM_WORLD", flag, order;

  MPI_Init(&argc, &argv);                               /* start up MPI                         */
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);              /* get the number of processes          */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);              /* get my rank among all the processes  */

  if (argc < 1) {                                       /* too few arguments aborts the program */
    cleanup(my_rank, "Error:  Too few arguments");
  }
  else if (argc == 1) {                                      /* option to run with a b n as inputs   */
    order = argv[1];
    if (order.compare("lowtohigh") == 0) {lth = true;}
    else if (order.compare("hightolow") == 0) {lth == false;}
    else {cleanup(my_rank, "Error:  incorrect argument for order");}
  }
  else if (argc == 2) {                                      /* option to run -verbose               */
    flag = argv[1];
    order = argv[2];

    if (flag.compare("-verbose") != 0) {cleanup(my_rank, 
        "Error:  Wrong flag, only '-verbose' supported");}
    else if (flag.compare("-verbose") == 0) {verbose = true;}

    if (order.compare("lowtohigh") == 0) {lth = true;}
    else if (order.compare("hightolow") == 0) {lth == false;}
    else {cleanup(my_rank, "Error:  incorrect argument for order");}
  }
  else (argc > 3) {                                       /* too many arguments aborts the program*/
    cleanup(my_rank, "Error:  Too many arguments");
  }

  if( comm_sz != isPowerOfTwo(comm_sz) ) {
    cleanup(my_rank, 
        "You have not executed the program with a number of process that is a power of two");
  }



  for (i = 0, i < 100, i++) {
    vector[i] = i/100;
    result[i] = 0;
  }

  BUTTERFLY_Reduce(&vector, &result, 100, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, lth, verbose);
  

  cleanup(my_rank, "Program Complete");                 /* terminates the program               */
}
