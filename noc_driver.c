#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "noc_header.h"

int main(int argc, char *argv[]) {

    int M_rows = 0, N_columns = 0;   // Mesh network dimensions
    int num_cores = 0;               // Total number of cores in the mesh network
    int num_io_pairs = 0;            // Number of io pairs available
    int num_test_cores = 0;          // Number of test cores
    NoC_node *noc_nodes;             // Pointer to NoC node structure array
    IO_pairs *io_pairs;              // Pointer to IO pairs structure array
    int *soln_vector;                // PSO solution vector
    FILE *fptr;                      // Input file pointer

    // Open and read input file
    fptr = fopen ("input.txt","r");
    if (fptr == NULL) {
        printf(" ERROR: Could not open the input file\n");
        return -1;
    }

    // All frequencies normalized wrt default test freq
    double freq[1] = {1.0};

    // Read network dimensions, number of i/o pairs
    fscanf (fptr,"%d\t%d", &M_rows, &N_columns);
    fscanf (fptr,"%d", &num_io_pairs);

    // Allocate memory for noc nodes array
    num_cores = M_rows * N_columns;
    noc_nodes = (NoC_node *) malloc (num_cores * sizeof (NoC_node));

    // Allocate memory for i/o pairs array
    io_pairs = (IO_pairs *) malloc (num_io_pairs * sizeof (IO_pairs));

    // Initialize all the nodes (NoC tiles) in the network by assigning cores and ids
    initialize_nodes (noc_nodes, num_cores, M_rows, N_columns);

    // Configure i/o pairs and store i/o pair info in io_pairs array for quick ref
    configure_io_pairs (fptr, io_pairs, noc_nodes, num_io_pairs);

    // Determine total test cores and read number of test patterns and scan chain lengths for each of them
    num_test_cores = num_cores - 2 * num_io_pairs;
    read_test_core_parameters (fptr, io_pairs, noc_nodes, num_cores, num_test_cores);

    // Close input file
    fclose (fptr);

    // Generate PSO particles
    srand(time(0));
    particle_swarm_optimization(noc_nodes, num_cores, M_rows, freq, 1/*num_freq*/, io_pairs, num_io_pairs);

    // Free allocated memory
    free (noc_nodes);
    free (io_pairs);

    return 0;
}
