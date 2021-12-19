// =================
// MACRO DEFINITIONS
// =================

#define MAX_NUM_CORES 30                           // Maximum number of cores allowed for a NoC -- just for array declaration convenience
#define MAX_IO_PAIRS 5                             // Maximum number of io pairs allowed for a NoC -- just for array declaration convenience
#define NUM_PSO_PARTICLES 1                        // Number of PSO particles considered for simulation
#define UNALLOCATED -1                             // To indicate UNALLOCATED field elements

// NoC node parameter values

#define TEST_CORE 1                                // core_type flag value 1 corresponds to test core
#define INPUT_CORE 2                               // core_type flag value 2 corresponds to input core
#define OUTPUT_CORE 3                              // core_type flag value 3 corresponds to output core       

// PSO evolution constants

#define ALPHA 0.5
#define BETA 0.5

// Router ports

// 16 Valid Router Statuses -- in accordance with XY routing
// For Output port --> Injection, Input ports allowed --> North, South, East, West, Unallocated (4)
// For Output port --> North, Input ports allowed --> Injection, South, East, West, Unallocated (4)
// For Output port --> East, Input ports allowed --> Injection, West, Unallocated (2) --> as Y-direction routing happens strictly after X
// For Output port --> South, Input ports allowed --> Injection, North, East, West, Unallocated (4)
// For Output port --> West, Input ports allowed --> Injection, East, Unallocated (2) --> as Y-direction routing happens strictly after X

// Input port numbers

#define INJECTION 0   
#define NORTH 1
#define EAST 2
#define SOUTH 3
#define WEST 4

// Output port number
// (Remaining output ports use the same indexing as above)

#define EJECTION 0

// Busyports struct 2D array second indices - indicating if the port and busytime information is for OUTPUT port corresponding to port index or INPUT port corresponding to port index

#define INPUT 0
#define OUTPUT 1

// To indicate test completion (completion of all preemption instances)

#define LAST_TEST_ADMINISTERED -73


#define FILLED 1
#define NOT_FILLED 0

// NoC node

typedef struct {
    int core_no;                                   // Core number corresponding to this node
    int router_no;                                 // Router number corresponding to this node
    int x_cord;                                    // x-coordinate of the core in the given mesh network
    int y_cord;                                    // y-coordinate of the core in the given mesh network
    int core_type;                                 // Flag to specify whether the core type is input/output/test
    int test_patterns;                             // Number of test patterns used to test the core (set to 0 for i/o)
    int scan_chain_length;                         // Scan chain length corresponding to this core
} NoC_node;

// Structs required to maintain schedule list per IO pair

struct _io_node {
    double starttime;                              // Test start time
    double endtime;                                // Test end time
    int test_core;                                 // Core under test
    struct _io_node *next;
};

typedef struct _io_node IO_node;

// IO schedule list head

typedef struct {
    int size;                                      // Total number of nodes in IO schedule list
    double max_busytime;                           // Maximum busytime of IO pair
    IO_node *head_node;                            // Pointer to head node of IO schedule list
} IO_head;

// IO pair

typedef struct {
    int io_pair_no;                                // IO pair index
    int input_core_no;                             // Input core index
    int output_core_no;                            // Output core index
    IO_head *io_head;                              // IO schedule list head
} IO_pairs;

// PSO particle

typedef struct {
    double mapping[4 * MAX_NUM_CORES];             // Mapping structure: | test core ids | io pairs assigned | test frequency | preemptions |
    double testtime;                               // Test time required for testing all the cores with the current mapping
    double SNR;                                    // Worst case SNR generated at the time of testing
    double communication_cost;                     // Communication cost for the given mapping (number of active MRs * number of test packets)
    double fitness;                                // Fitness function value calculated for the given mapping
    double lbest_mapping[4 * MAX_NUM_CORES];       // The mapping corresponding to the best fitness function value obtained this particle till now
    double lbest_fitness;                          // The best fitness function value obtained this particle till now
} PSO_particle;

// Global best particle

typedef struct {
    double gbest_mapping [4 * MAX_NUM_CORES];      // The mapping corresponding to the best fitness function value obtained this among all particles
    double gbest_fitness;                          // The best fitness function value obtained among all particles 
} Gbest_PSO_particle;

typedef struct {
    int port; 
    double busytime;
} Busyports;

// Resource matrix

typedef struct {
    double busytime;                               // If the resource status is BUSY, this field gives the time till which the resource remains occupied  
    Busyports busyports[5][2];                     // If the resource type is a router, this structure tracks all ACTIVE busy statuses and their busytimes
                                                   // First index gives port, second index gives the port type for which information is stored
                                                   // For instance, busyports[EAST][INPUT].port = WEST implies, the WEST is the INPUT port corresponding
                                                   // to OUTPUT port EAST, busytime for corresponding pair is stored in busyports[EAST][INPUT].busytimes
} Resource;

// Swap operators

typedef struct {
    int swap_idx1;                                 // Index of the first element to be swapped 
    int swap_idx2;                                 // Index of the second element to be swapped 
} Swap_operator;

// Sorted list of all starttimes and endtimes
// (Used to create CLAP input list)

struct _all_times {
    double time;
    struct _all_times *next;
};

typedef struct _all_times All_times;

// CLAP input list

typedef struct {
    int source_x;
    int source_y;
    int destination_x;
    int destination_y;
} Test_signals;

// CLAP input list

struct _clap_inputs {
    Test_signals test_signals[2 * MAX_IO_PAIRS];
    struct _clap_inputs *next;
};

typedef struct _clap_inputs Clap_inputs;

// Function declarations

// Initializes the node structs with ids, core type and coordinates info
void initialize_nodes (NoC_node *noc_nodes, int num_cores, int M_rows, int N_columns);

// Reads the input and output core indices from the file and stores info in IO pair struct array 
// Also configures the corresponding node structs as INPUT/OUTPUT CORES                      
void configure_io_pairs (FILE* in_file, IO_pairs *io_pairs, NoC_node *noc_nodes, int num_io_pairs);

// Reads the number of test patterns and scan chain lengths for each test core from the file and stores this info in node struct array
void read_test_core_parameters (FILE* in_file, IO_pairs *io_pairs, NoC_node *noc_nodes, int num_cores, int num_test_cores);

// Finds testtime for the given mapping assuming there are no resource conflicts involved (XY routing, circuit switching)
double find_individual_testtime (NoC_node *noc_nodes, int input_core, int output_core, int test_core, double frequency, double preemption_point);

// Finds the communication cost for a given PSO particle mapping --> consider hops (circuit switching scenario) --> use testtime (non-preemptive, single frequency)
void find_communication_cost (PSO_particle *pso_particle, NoC_node *noc_nodes, int num_cores, IO_pairs *io_pairs, int num_io_pairs);

// Maintains an ordered list of all starttimes and endtimes (used to generate CLAP input) 
void update_times_list (All_times **head, double time);

// Creates IO schedule lists
IO_head *create_IO_list_head ();

// Updates IO schedule lists
void update_IO_list (IO_head *head, double starttime, double endtime, int testcore);

// Creates an input list for the CLAP tool 
void create_clap_input_list (Clap_inputs *head, IO_pairs *io_pairs, int num_io_pairs);

// Populates the resource matrix with busytimes for all resources [links and router ports] in accordance with the XY-routing algorithm
void find_resource_busytimes (PSO_particle *pso_particle, NoC_node *noc_nodes, int N_columns, int num_cores, IO_pairs *io_pairs, int num_io_pairs);

// Initializes PSO particles by initializing the I/O core, frequencies and test core mapping; calculates the fitness value for each particle 
void init_pso_particles (PSO_particle *pso_particle, Gbest_PSO_particle *gbest_pso_particle, NoC_node *noc_nodes, int num_cores, double *freq, int num_freq, IO_pairs *io_pairs, int num_io_pairs, int N_columns);

// Generates random number between 0 and 1 (probability distribution: uniform)
double generate_random_number ();

// Swaps io pairs with given probability 
void swap_io_pair (int num_test_cores, double *a, double *b, double probability);

// Checks frequency validity for newly assigned test core, swaps frequencies
void swap_frequencies (int num_test_cores, double *a, double *b, double probability);

// Generates a sequence of swap operators for evolving a given particle's test core sequence
int generate_swap_operator_sequence (int num_test_cores, double *a, double *b, Swap_operator* swap_operator);

// Applies the sequence of swap operators on test core sequence a with give probability
void swap_test_core_sequence (int num_test_cores, double *a, Swap_operator* swap_operator, int num_swap_operators, double probability);

// Modifies the preemption points for test cores in a given particle (new position of a particle in continuous PSO)
// void modify_preemption_points (int num_test_cores, double *a, double *b, double *c);

// Simulates Particle Swarm Optimization algorithm to determine the mapping with minimum cost
void particle_swarm_optimization (NoC_node *noc_nodes, int num_cores, int M_rows, double *freq, int num_freq, IO_pairs *io_pairs, int num_io_pairs);

// Finds the maximum of two given numbers
double max (double a, double b);

// Prints the mapping and test schedule information for all PSO particles
void print_pso_particle_info (PSO_particle *pso_particle, int num_test_cores);

// Prints the mapping and test schedule information for the global best PSO particle
void print_global_best_info (Gbest_PSO_particle *gbest_pso_particle, int num_test_cores);

// Prints IO schedule lists
void print_IO_schedule_lists (IO_head* head);

