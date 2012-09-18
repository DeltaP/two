/*
 * Gregory Petropoulos
 *
 * This is my butterfly all reduce program for assignment 2
 *
 * To compile:  mpicxx -g -Wall -std=c99 -o butterfly butterfly.cpp
 * To run:  mpiexec -n 8 ./butterfly [-verbose] <lth>
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

#define N 10
int My_rank, Comm_sz;


// -----------------------------------------------------------------
// prints a prefix
void printpre() {
    static int linenum = 0;
    printf ("[%d/%d:%4d] : ", My_rank, Comm_sz, linenum++); 
}
// -----------------------------------------------------------------



// Ends the program on an error and prints message to node 0
void cleanup (const char *message) {
  if (My_rank==0) {
    printpre();
    printf("%s\n",message);
  }
  MPI_Finalize();                                       /* kills mpi                            */
  exit(0);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// my version of the all reduce
void butterfly_reduce(double vector[], double result[], int count, string datatype, string operation, string comm, bool lth) {
  int i, condition, talk;
  double addtoresult[N];

  if (datatype.compare("MPI_DOUBLE") != 0) { cleanup("Error:  Only supports Double"); }
  if (operation.compare("MPI_SUM") != 0) { cleanup("Error:  Only supports Sum"); }
  if (comm.compare("MPI_COMM_WORLD") != 0) { cleanup("Error:  Wrong comm world"); }

  for (i = 0; i < count; i++) {
    result[i] = vector[i];
  }

  if (lth == true) {                                    /* if the ordering is low to high       */
    condition = 1;
    while (condition < Comm_sz ) {
      talk = condition ^ My_rank;
      if (talk < My_rank) {
        printpre();
        cout << "Process " << My_rank << ":  stage " << condition << ", sending to " << talk << endl;
        MPI_Send(&result[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        MPI_Recv(&addtoresult[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (i = 0; i < N; i++) {
          result[i] += addtoresult[i];
        }
      }
      else {
        printpre();
        cout << "Process " << My_rank << ":  stage " << condition << ", recieving from " << talk << endl;
        MPI_Recv(&addtoresult[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&result[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        for (i = 0; i < N; i++) {
          result[i] += addtoresult[i];
        }
      }
      condition <<= 1;
    }
  }
  else {
    condition = Comm_sz;
    while (condition > 1) {                             /* if the ordering is high to low       */
      condition >>= 1;
      talk = condition ^ My_rank;
      if (My_rank < talk) {
        printpre();
        cout << "Process " << My_rank << ":  stage " << condition << ", sending to " << talk << endl;
        MPI_Send(&result[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        MPI_Recv(&addtoresult[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (i = 0; i < N; i++) {
          result[i] += addtoresult[i];
        }
      }
      else {
        printpre();
        cout << "Process " << My_rank << ":  stage " << condition << ", recieving from " << talk << endl;
        MPI_Recv(&addtoresult[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&result[0], count, MPI_DOUBLE, talk, 0, MPI_COMM_WORLD);
        for (i = 0; i < N; i++) {
          result[i] += addtoresult[i];
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Checks if integer is a power of two
int isPowerOfTwo (int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// the main program
int main(int argc, char *argv[]) {
  double vector[N], result[N];
  int i;
  bool verbose = false;
  bool lth;
  string datatype = "MPI_DOUBLE";
  string operation = "MPI_SUM"; 
  string comm = "MPI_COMM_WORLD"; 
  string flag, order;

  MPI_Init(&argc, &argv);                               /* start up MPI                         */
  MPI_Comm_size(MPI_COMM_WORLD, &Comm_sz);              /* get the number of processes          */
  MPI_Comm_rank(MPI_COMM_WORLD, &My_rank);              /* get my rank among all the processes  */

  if (argc < 2) {                                       /* too few arguments aborts the program */
    cleanup("Error:  Too few arguments");
  }
  else if (argc == 2) {                                      /* option to run with a b n as inputs   */
    order = argv[1];
    if (order.compare("lowtohigh") == 0) {lth = true;}
    else if (order.compare("hightolow") == 0) {lth = false;}
    else {cleanup("Error:  incorrect argument for order");}
  }
  else if (argc == 3) {                                      /* option to run -verbose               */
    flag = argv[1];
    order = argv[2];

    if (flag.compare("-verbose") != 0) {cleanup( 
        "Error:  Wrong flag, only '-verbose' supported");}
    else if (flag.compare("-verbose") == 0) {verbose = true;}

    if (order.compare("lowtohigh") == 0) {lth = true;}
    else if (order.compare("hightolow") == 0) {lth = false;}
    else {cleanup("Error:  incorrect argument for order");}
  }
  else {                                       /* too many arguments aborts the program  */
    cleanup("Error:  Too many arguments");
  }

  if(isPowerOfTwo(Comm_sz) == 0) {
    cleanup( 
        "You have not executed the program with a number of process that is a power of two");
  }

  for (i = 0; i < N; i++) {
    vector[i] = i*0.001;
    result[i] = 0;
  }

  butterfly_reduce(vector, result, N, datatype, operation, comm, lth);

  if (My_rank==0) {
    printpre();
    cout << "The answer is [";
      for (i = 0; i < N-1; i++) {
        cout << result[i] << ", ";
      }
    cout << result[N-1] << "]\n";
  }

  cleanup("Program Complete");                 /* terminates the program               */
}
