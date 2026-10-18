#include <stdio.h>
#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "comm.h"
#include "vanilla.h"
#include <err.h>
#include <omp.h>




MPI_Comm Comm_world; // Comunicator (process inside the same comunicator can comunicate betwen them), GLobal comunicator, it keeps all the precesss
MPI_Comm Comm_world_x; // Comunicator for all process in the line x
MPI_Comm Comm_world_y; // Comunicator for all process in the line y
MPI_Status Comm_status; // Status of the last comunicator, it keeps messages ssuch how many of them ssend and received

void comm_init(int nx, int ny, int argc, char **argv) { // function to prepare all comunicators, it is called at the begining of the program
    int provided; // OpenMP threads exist, but only the master thread calls MPI: FUNNELED is the level we need
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided); // Initialize the MPI environment
    if (provided < MPI_THREAD_FUNNELED) fprintf(stderr, "Warning: MPI library does not support MPI_THREAD_FUNNELED\n");
    MPI_Comm_size(MPI_COMM_WORLD, &Ps); // turn on the comunication system, discover how many process is running
    MPI_Comm_rank(MPI_COMM_WORLD, &Rs); // Find out the identifier (rank) of THIS specific process, Rs (Rank Self) will receive a unique number between 0 and Ps-1
    vanilla_lattice_init(nx, ny, Rs, Ps); // Initialize the lattice, Each process is responsible for a part of the grid. Example: 100x100 grid with 4 processes → each one is responsible for 50x50.
    Nxp = (Rx + 1) % Px + Ry * Px; // Caclculate the rank of the right newighboor, 
    // Imagine the processes organized in a 2D grid:

    // [0] [1] [2]
    // [3] [4] [5]

    // If I am in process 1 (Rx=1, Ry=0):

    // (Rx + 1) % Px = (1 + 1) % 3 = 2 (next in line)
    // Ry * Px = 0 * 3 = 0
    // Nxp = 2 + 0 = 2 (process 2 is my neighbor to the right)
    Nxm = (Rx - 1 + Px) % Px + Ry * Px; // calculate the rank of the left neighbor
    // Continuing the previous example, if I am in process 1:

    // (Rx - 1 + Px) % Px = (1 - 1 + 3) % 3 = 3 % 3 = 0
    // Nxm = 0 + 0 = 0 (process 0 is my left neighbor)
    // The "+ Px" before the modulus ensures that we don't have negative numbers
    // If I am on the left edge (Rx=0), I go back to the right (Rx=Px-1)

    Nyp = Rx + ((Ry + 1) % Py) * Px; // calculate the rank of the top neightboor
    // [0] [1] [2] ← line 0 (Ry=0)
    // [3] [4] [5] ← line 1 (Ry=1)
    // If I am in process 1 (Rx=1, Ry=0):

    // (Ry + 1) % Py = (0 + 1) % 2 = 1
    // ((Ry + 1) % Py) * Px = 1 * 3 = 3
    // Nyp = 1 + 3 = 4 (process 4 is below me)


    Nym = Rx + ((Ry - 1 + Py) % Py) * Px; // calculate the rank of the bottom neighbor
    // Continuing the previous example, if I am in process 1 (Rx=1, Ry=0):

    // (Ry - 1 + Py) % Py = (0 - 1 + 2) % 2 = 1 % 2 = 1
    // ((Ry - 1 + Py) % Py) * Px = 1 * 3 = 3
    // Nym = 1 + 3 = 4 (process 4 is above me)

    Comm_world = MPI_COMM_WORLD; // Copy the global comunicator
    MPI_Comm_split(Comm_world, Ry, Rx, &Comm_world_x); // Split the global comunicator into comunicators for each row
    MPI_Comm_split(Comm_world, Rx, Ry, &Comm_world_y); // Split the global comunicator into comunicators for each column
}

void comm_finish() { // function to finish the comunication
    MPI_Comm_free(&Comm_world_x); // Free the comunicator for each row
    MPI_Comm_free(&Comm_world_y); // Free the comunicator for each column
    MPI_Finalize(); // Finish the comunication
}
void comm_barrier() { MPI_Barrier(Comm_world); } // function to synchronize all the processes, make all processors wait one for another

int comm_broadcast_i(int x) { // function to broadcast a integer to all processes (the number of processe, lattice size etc...))
    MPI_Bcast(&x, 1, MPI_INT, 0, Comm_world); // Broadcast the integer (number of process) to al ranks
    return x; // Return the integer
}

int comm_broadcast_s(char *buf, size_t bytes) { // function to broadcast a string to all processes (the flags when runs it)
    MPI_Bcast(buf, bytes, MPI_CHAR, 0, Comm_world); // Broadcast the string to all processes (the process 0 read and send to anothers)
    return 0;
}

int comm_broadcast_d(double *buf, size_t bytes) { // function to broadcast a double to all processes
    int count = bytes / sizeof(double); // Converts size in bytes to the number of elements of type double.
    MPI_Bcast(buf, count, MPI_DOUBLE, 0, Comm_world); // Broadcast the double to all processes (the process 0 read and send to anothers)
    return 0;
}

double comm_globsum(double x) { // function to sum all the values of the above functions
    double result;
    MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, MPI_SUM, Comm_world);
    return result;
}

int comm_global_or(int x) { // verify if at least one process has the true condition
    int result;
    MPI_Allreduce(&x, &result, 1, MPI_INT, MPI_LOR, Comm_world);
    return result;
}

int comm_master() { return (Rs == 0) ? 1 : 0; } // will print true for the zero process and false for anothers (process 0 is the master)
int comm_master_rank() { return 0; } // pass the master ID 
void comm_abort(int err) { MPI_Abort(Comm_world, err); } // abort if something goes wrong

int comm_mine(int coords[2]) { // verify for each process the point (x,y) belongs to it
    int x = coords[0], y = coords[1];
    return (x >= Ox && x < Ox + Lx && y >= Oy && y < Oy + Ly);
}

//#endif
