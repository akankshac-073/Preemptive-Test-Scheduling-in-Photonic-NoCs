#define MAX_NUM_CORES 20                           // Maximum number of cores allowed for a NoC -- just for array declaration convenience
#define NUM_PSO_PARTICLES 1                        // Number of PSO particles considered for simulation
#define UNALLOCATED -1                             // To indicate UNALLOCATED field elements

// NOC NODE PARAMETER VALUES

#define TEST_CORE 1                                // core_type flag value 1 corresponds to test core
#define INPUT_CORE 2                               // core_type flag value 2 corresponds to input core
#define OUTPUT_CORE 3                              // core_type flag value 3 corresponds to output core       

// PSO EVOLUTION CONSTANT VALUES

#define ALPHA 0.5
#define BETA 0.5

// ROUTER PORTS

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

// Busyports struct 2D array second indices 
// - indicating if the port and busytime information is for OUTPUT port corresponding to port index or INPUT port corresponding to port index

#define INPUT 0
#define OUTPUT 1

// Active MR counts for different routing paths MR_<input port>_<output port>
// Consider number of Mrs into test packets
#define MR_NORTH_INJECTION 1
#define MR_SOUTH_INJECTION 1 // check path
#define MR_EAST_INJECTION 1
#define MR_WEST_INJECTION 1
#define MR_INJECTION_NORTH 1 
#define MR_SOUTH_NORTH 0
#define MR_EAST_NORTH 1
#define MR_WEST_NORTH 1
#define MR_INJECTION_SOUTH 1
#define MR_NORTH_SOUTH 0
#define MR_EAST_SOUTH 1
#define MR_WEST_SOUTH 1
#define MR_INJECTION_EAST 1
#define MR_WEST_EAST 0
#define MR_INJECTION_WEST 1
#define MR_EAST_WEST 0

// Indicates test completion

#define LAST_TEST_ADMINISTERED -73

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

// IO pair

typedef struct {
    int io_pair_no;                                // IO pair index
    int input_core_no;                             // Input core index
    int output_core_no;                            // Output core index
} IO_pairs;

// Particle schedule linked list node

struct _schedule {
    double starttime;
    double endtime;
    struct _schedule* next;
};

typedef struct _schedule Schedule_node;

// Schedule HEAD structure

typedef struct {
    double min_time;
    double max_time; 
    Schedule_node *head_node;
} Schedule_head;

// PSO particle

typedef struct {
    double mapping [4 * MAX_NUM_CORES];            // Mapping structure: | test core ids | io pairs assigned | test frequency | preemptions |
    double testtime;                               // Test time required for testing all the cores with the current mapping
    double SNR;                                    // Worst case SNR generated at the time of testing
    double communication_cost;                     // Communication cost for the given mapping (number of active MRs * number of test packets)
    double fitness;                                // Fitness function value calculated for the given mapping
    double lbest_mapping [4 * MAX_NUM_CORES];      // The mapping corresponding to the best fitness function value obtained this particle till now
    double lbest_fitness;                          // The best fitness function value obtained this particle till now
    Schedule_head *schedule;                       // Linked list containing start and end times of testing for the given particle
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
    Busyports busyports[5][2];                     // If the resource type is a router, this structure tracks all ACTIVE busy statuses and their busytimes
                                                   // First index gives port, second index gives the port type for which information is stored
                                                   // For instance, busyports[EAST][INPUT].port = WEST implies, the WEST is the INPUT port corresponding
                                                   // to OUTPUT port EAST, busytime for corresponding pair is stored in busyports[EAST][INPUT].busytimes
    double busychannel;                            // If the resource type is a link, track the fraction of channel occupied for given resource
    double prev_busytime;                          // Required to calculate busytime when the resource is used simultaneously by multiple tests
                                                   // This field holds the previous busytime (max{busytime calculated by all prev occupying tests})
    double busytime;                               // If the resource status is BUSY, this field gives the time till which the resource remains occupied  
} Resource;

// Swap operators

typedef struct {
    int swap_idx1;                                 // Index of the first element to be swapped 
    int swap_idx2;                                 // Index of the second element to be swapped 
} Swap_operator;


// Function declarations

// Initializes the node structs with ids, core type and coordinates info
void initialize_nodes (NoC_node *noc_nodes, int num_cores, int M_rows, int N_columns);

// Reads the input and output core indices from the file and stores info in IO pair struct array 
// Also configures the corresponding node structs as INPUT/OUTPUT CORES                      
void configure_io_pairs (FILE* in_file, IO_pairs *io_pairs, NoC_node *noc_nodes, int num_io_pairs);

// Reads the number of test patterns and scan chain lengths for each test core from the file and stores this info in node struct array
void read_test_core_parameters (FILE* in_file, IO_pairs *io_pairs, NoC_node *noc_nodes, int num_cores, int num_test_cores);

// Initializes PSO particles by initializing the solution vector fields and calculating the resulting fitnessAlso initializes the local and global best particles for this stage 
void init_pso_particles (PSO_particle *pso_particle, Gbest_PSO_particle *gbest_pso_particle, NoC_node *noc_nodes, int num_cores, double *freq, int num_freq, IO_pairs *io_pairs, int num_io_pairs, int N_columns);

// Finds the particle fitness = w.testtime + (1-w).SNR -- w optimal, as testtime increases and SNR decreases
void find_particle_fitness ();

// Generates a random number between 0 and 1 (uniform distribution)
double generate_random_number ();

// Swaps io pairs in arrays a and b with given probability
void swap_io_pair (int num_test_cores, double *a, double *b, double probability);

// Swaps frequencies in arrays a and b with given probability
void swap_frequencies (int num_test_cores, double *a, double *b, double probability);

// Generates a sequence of swap operators for evolving a given particle's test core sequence 
int generate_swap_operator_sequence (int num_test_cores, double *a, double *b, Swap_operator* swap_operator);

// Modifies the preemption points for test cores in a given particle (new position of a particle in continuous PSO)
void modify_preemption_points (int num_test_cores, double *a, double *b, double *c);

// Simulates Particle Swarm Optimization algorithm to find the mapping with minimum cost
void particle_swarm_optimization(NoC_node *noc_nodes, int num_cores, int M_rows, double *freq, int num_freq, IO_pairs *io_pairs, int num_io_pairs);

// Utility function to print PSO particles' info
void print_pso_particle_info (PSO_particle *pso_particle, int num_test_cores);

// Utility function to print global best particle's info
void print_global_best_info (Gbest_PSO_particle *gbest_pso_particle, int num_test_cores);

