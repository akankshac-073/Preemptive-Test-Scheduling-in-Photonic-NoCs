#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "noc_header.h"

int times_list_size = 0;

// Assign core numbers and coordinates for each node struct

void initialize_nodes (NoC_node *noc_nodes, int num_cores, int M_rows, int N_columns) {		
    int i = 0;
    for (i = 0; i < num_cores; i++) {
        noc_nodes[i].core_no = i + 1;
        noc_nodes[i].router_no = i + 1;
        noc_nodes[i].x_cord = (i % N_columns);
        noc_nodes[i].y_cord = (i / N_columns) % M_rows;
        noc_nodes[i].core_type = TEST_CORE;
    }
}

// Reads the input and output core numbers from the file and configures corresponding nodes

void configure_io_pairs (FILE* in_file, IO_pairs *io_pairs, NoC_node *noc_nodes, int num_io_pairs) {

    int temp_in = 0;    // temporary variable to store input core no
    int temp_out = 0;   // temporary variable to store output core no

    for (int i = 0 ; i < num_io_pairs ; i++) {
        fscanf (in_file,"%d\t%d", &temp_in, &temp_out);

        // Update node structure array with core type
        noc_nodes[temp_in - 1].core_type = INPUT_CORE;
        noc_nodes[temp_out - 1].core_type = OUTPUT_CORE;

        // Populate the io pair structure array
        io_pairs[i].io_pair_no = i + 1;
        io_pairs[i].input_core_no = temp_in;
        io_pairs[i].output_core_no = temp_out;
        io_pairs[i].io_head = create_IO_list_head();
    }
}

// Reads the number of test patterns, scan chain length and max power budget for each test core from the file and configures corresponding nodes

void read_test_core_parameters (FILE* in_file, IO_pairs *io_pairs, NoC_node *noc_nodes, int num_cores, int num_test_cores) {

    for (int i = 0; i < num_cores; i++) {
        if (noc_nodes[i].core_type == TEST_CORE) {
            fscanf (in_file,"%d\t%d\n", &noc_nodes[i].test_patterns, &noc_nodes[i].scan_chain_length);
        }
    }
}

// Finds testtime for the given mapping assuming there are no resource conflicts involved (XY routing, circuit switching)

double find_individual_testtime (NoC_node *noc_nodes, int input_core, int output_core, int test_core, double frequency, double preemption_point) {

    int hop_length_ic = 0;        // Distance in terms of hops from input core to test core
    int hop_length_co = 0;        // Distance in terms of hops from output core to test core
    double testtime = 0;          // Testtime = [1 + Max{hop_length_ic, hop_length_co}] × no_test_patterns + [Min{hop_length_ic, hop_length_co}]

    // Input core to test core
    if ((noc_nodes[input_core - 1].x_cord == noc_nodes[test_core - 1].x_cord) || (noc_nodes[input_core - 1].y_cord == noc_nodes[test_core - 1].y_cord))
        hop_length_ic = 1;
    else
        hop_length_ic = 2;

    // Test core to output core
    if ((noc_nodes[test_core - 1].x_cord == noc_nodes[output_core - 1].x_cord) || (noc_nodes[test_core - 1].y_cord == noc_nodes[output_core - 1].y_cord))
        hop_length_co = 1;
    else
        hop_length_co = 2;

    // Testtime calculation - [ref Thermal-aware Test Scheduling Strategy for Network-on-Chip based Systems]
    if (hop_length_ic >= hop_length_co)
        testtime = (1 + hop_length_ic + (noc_nodes[test_core-1].scan_chain_length - 1)) * (int)(noc_nodes[test_core-1].test_patterns * preemption_point) + hop_length_co;
    else
        testtime = (1 + hop_length_co + (noc_nodes[test_core-1].scan_chain_length - 1)) * (int)(noc_nodes[test_core-1].test_patterns * preemption_point) + hop_length_ic;

    // Accounting for the frequency (frequency --> normalized wrt the base testing freq)
    testtime = testtime / frequency;

    return testtime;
}

// Finds the communication cost for a given PSO particle mapping --> consider hops (circuit switching scenario) --> use testtime (non-preemptive, single frequency)

void find_communication_cost (PSO_particle *pso_particle, NoC_node *noc_nodes, int num_cores, IO_pairs *io_pairs, int num_io_pairs) {

    int num_test_cores = num_cores - (2 * num_io_pairs);   // Number of test cores in the NoC mesh network
    double communication_cost = 0;                         // Communication cost = Number of active MRs * test packets
    int test_core = 0;                                     // Temporary variable to store test core number
    int input_core = 0;                                    // Temporary variable to store input core number
    int output_core = 0;                                   // Temporary variable to store output core number

    // For all test cores
    for (int i = 0; i < num_test_cores; i++) {

        test_core = (int)pso_particle->mapping[i];
        input_core = io_pairs[(int)(pso_particle->mapping[i + num_test_cores]) - 1].input_core_no;
        output_core = io_pairs[(int)(pso_particle->mapping[i + num_test_cores]) - 1].output_core_no;

        // Find communication costs

        // Input core to test core
        if ((noc_nodes[input_core - 1].x_cord == noc_nodes[test_core - 1].x_cord) || (noc_nodes[input_core - 1].y_cord == noc_nodes[test_core - 1].y_cord))
            pso_particle->communication_cost = 2 * noc_nodes[test_core-1].test_patterns;
        else
            pso_particle->communication_cost = 3 * noc_nodes[test_core-1].test_patterns;

        // Test core to output core
        if ((noc_nodes[test_core - 1].x_cord == noc_nodes[output_core - 1].x_cord) || (noc_nodes[test_core - 1].y_cord == noc_nodes[output_core - 1].y_cord))
            pso_particle->communication_cost += 2 * noc_nodes[test_core-1].test_patterns;
        else
            pso_particle->communication_cost += 3 * noc_nodes[test_core-1].test_patterns;
    }
}

// Maintains a sorted list of all starttimes and endtimes 
// (giving intervals for which CLAP inputs need to be generated)

void update_times_list (All_times **head, double time) {
    All_times* temp;
    All_times* new_node = (All_times*) malloc (sizeof (All_times));
    new_node->time = time;

    // If list is EMPTY or the new node is to be inserted at the head of the list
    if (*head == NULL || (*head)->time >= new_node->time) {
        new_node->next = *head;
        *head = new_node;
    }
    
    // Else, if the new node is to be inserted in the middle/ at the end of the list
    else {
        temp = *head;
        while (temp->next != NULL && temp->next->time < new_node->time)
            temp = temp->next;
        new_node->next = temp->next;
        temp->next = new_node;
    }

    // Increment the list size (global var)
    times_list_size++;
}

// Creates an EMPTY IO Schedule list

IO_head *create_IO_list_head () {

    // Allocate memory for the head struct
    IO_head *head;
    head = malloc (sizeof (IO_head));
    
    // Initialize queue parameters
    head->head_node = NULL;                       // Initialize head node pointer as NULL
    head->size = 0;                               // Initialize list size to 0
    head->max_busytime = 0.0;                     // Initialize max busytime as 0.0
    
    // Return pointer to run queue head
    return head;
}

// Updates IO Schedule lists

void update_IO_list (IO_head *head, double starttime, double endtime, int test_core) {

    IO_node* temp;

    // Allocate memory for the new node and set node parameters
    IO_node* new_node = (IO_node *) malloc (sizeof (IO_node));
    new_node->starttime = starttime;
    new_node->endtime = endtime;
    new_node->test_core = test_core;
    head->size++;
    head->max_busytime = new_node->endtime;

    // If the queue is empty, new node becomes the head node
    if (head->head_node == NULL) {
        new_node->next = head->head_node;
        head->head_node = new_node;
    }

    // Else, the new node is added to the end of the non-empty list
    else {
        temp = head->head_node;
        while (temp->next != NULL)
            temp = temp->next;
        new_node->next = temp->next;
        temp->next = new_node;
    }
}

// Populates the resource matrix with busytimes for all resources [links and router ports] in accordance with the XY-routing algorithm

void find_resource_busytimes (PSO_particle *pso_particle, NoC_node *noc_nodes, int N_columns, int num_cores, IO_pairs *io_pairs, int num_io_pairs) {

    int num_test_cores = num_cores - (2 * num_io_pairs);      // Number of test cores in the NoC mesh network
    Resource resource_matrix[num_cores][num_cores];           // Resource matrix contains staus and busytime info (FIXME: could cause segmentation fault)
    double individual_testtime = 0.0;                         // Individual testtime for a given test core assuming no resource conflicts
    int test_core = 0;                                        // Temporary variable to store test core number
    int input_core = 0;                                       // Temporary variable to store input core number
    int output_core = 0;                                      // Temporary variable to store output core number
    double frequency = 0.0;                                   // Temporary variable to store test frequency
    double preemption = 0.0;                                  // Temporary variable to store preemption point
    int ic_core = 0;                                          // Stores the intermediate core obtained after matching x coordinates of input core and test core
    int co_core = 0;                                          // Stores the intermediate core obtained after matching x coordinates of test core and output core 
    double max_busytime = 0.0;                                // Maximum busytime for all resources conflicting with the given input/output port
    int xrouting_port_ic = -1;                                // Temporary variable to store the (EAST/WEST) routing port used at the last core after X routing
    int xrouting_port_co = -1;                                // Temporary variable to store the (EAST/WEST) routing port used at the last core after X routing
    All_times *times_head = NULL; 
    int io_pair = 0;
    double starttime = 0.0;
    double endtime = 0.0;
    double io_busytime = 0.0;

    // Initializing the resource matrix
    for (int i = 0; i < num_cores; i++) {
        for (int j = 0; j < num_cores; j++) {
            resource_matrix[i][j].busytime = 0.0;

            // Initializing input output port assignments and busytimes
            for (int k = 0; k < 5; k++) {
                resource_matrix[i][j].busyports[k][INPUT].busytime = 0.0;
                resource_matrix[i][j].busyports[k][OUTPUT].busytime = 0.0;
                resource_matrix[i][j].busyports[k][INPUT].port = UNALLOCATED;
                resource_matrix[i][j].busyports[k][OUTPUT].port = UNALLOCATED;
            }
        }
    }

    // Populating resource matrix with busytime (latest time till which the resource is busy) for all test cores
    for (int i = 0; i < num_test_cores; i++) {

        // Use temporary variables to store core numbers and parameters -- just for convenience
        test_core = (int)pso_particle->mapping[i];
        io_pair = (int)(pso_particle->mapping[i + num_test_cores]);
        input_core = io_pairs[io_pair - 1].input_core_no;
        output_core = io_pairs[io_pair - 1].output_core_no;
        frequency = pso_particle->mapping[i + (2 * num_test_cores)];
        preemption = pso_particle->mapping[i + (3 * num_test_cores)];

        // Get max busytime corresponding to the IO pair used to test the given core
        io_busytime = io_pairs[io_pair - 1].io_head->max_busytime;

        // Find core individual testtime (assuming no resource conflicts)
        individual_testtime = find_individual_testtime (noc_nodes, input_core, output_core, test_core, frequency, preemption);
        printf(" %lf\n\n", individual_testtime);

        // ------------------------------------------------------------------- ROUTING LOGIC -------------------------------------------------------------------

        // ------------------------------------------------------------------ X ROUTING LOGIC ------------------------------------------------------------------

        // If input core x coordinate < test core x coordinate, update busytime for resources b/w cores
        if (noc_nodes[input_core - 1].x_cord < noc_nodes[test_core - 1].x_cord) {

            xrouting_port_ic = EAST;

            // Update LINK busytimes between the input core and its adjacent core (link direction +x)
            // resource_matrix[input_core - 1][input_core].busytime += individual_testtime;
            resource_matrix[input_core - 1][input_core].busytime += max(resource_matrix[input_core - 1][input_core].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[input_core - 1][input_core].busytime);

            // Update ROUTER busytimes and statuses (routing direction +x i.e. input port INJECTION to output port EAST)
            resource_matrix[input_core - 1][input_core - 1].busyports[EAST][INPUT].busytime = max(resource_matrix[input_core - 1][input_core - 1].busyports[EAST][INPUT].busytime, io_busytime);
            resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime = max(resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime, io_busytime);
            max_busytime = max(resource_matrix[input_core - 1][input_core - 1].busyports[EAST][INPUT].busytime, resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime);
            resource_matrix[input_core - 1][input_core - 1].busyports[EAST][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[input_core - 1][input_core - 1].busyports[EAST][INPUT].port = INJECTION;
            resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].port = EAST;
            if (resource_matrix[input_core - 1][input_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[input_core - 1][input_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[input_core - 1][input_core - 1].busytime);

            // From (input core + 1) to (input core + (difference in x-coordinates - 1))
            // * Index of input core + 1 is input_core
            for (int j = 0; j < (noc_nodes[test_core - 1].x_cord - noc_nodes[input_core - 1].x_cord - 1); j++) {
                
                // Update LINK busytimes between adjacent cores (link direction +x)
                // resource_matrix[input_core + j][input_core + j + 1].busytime += individual_testtime;
                resource_matrix[input_core + j][input_core + j + 1].busytime += max(resource_matrix[input_core + j][input_core + j + 1].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[input_core + j][input_core + j + 1].busytime);

                // Update ROUTER busytimes and statuses (routing direction +x i.e. input port WEST to output port EAST)
                resource_matrix[input_core + j][input_core + j].busyports[EAST][INPUT].busytime = max(resource_matrix[input_core + j][input_core + j].busyports[EAST][INPUT].busytime, io_busytime);
                resource_matrix[input_core + j][input_core + j].busyports[WEST][OUTPUT].busytime = max(resource_matrix[input_core + j][input_core + j].busyports[WEST][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[input_core + j][input_core + j].busyports[EAST][INPUT].busytime, resource_matrix[input_core + j][input_core + j].busyports[WEST][OUTPUT].busytime);
                resource_matrix[input_core + j][input_core + j].busyports[EAST][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[input_core + j][input_core + j].busyports[WEST][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[input_core + j][input_core + j].busyports[EAST][INPUT].port = WEST;
                resource_matrix[input_core + j][input_core + j].busyports[WEST][OUTPUT].port = EAST;
                if (resource_matrix[input_core + j][input_core + j].busytime < max_busytime + individual_testtime)
                    resource_matrix[input_core + j][input_core + j].busytime = max_busytime + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[input_core + j][input_core + j].busytime);
            }
        }

        // If input core x coordinate > test core x coordinate, update busytime for links b/w cores
        else if (noc_nodes[input_core - 1].x_cord > noc_nodes[test_core - 1].x_cord) {

            xrouting_port_ic = WEST;

            // Update LINK busytimes between the test core and its adjacent core (link direction -x)
            // resource_matrix[input_core - 1][input_core - 2].busytime += individual_testtime;
            resource_matrix[input_core - 1][input_core - 2].busytime += max(resource_matrix[input_core - 1][input_core - 2].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[input_core - 1][input_core - 2].busytime);

            // Update ROUTER busytimes and statuses (routing direction -x i.e. input port INJECTION to output port WEST)
            resource_matrix[input_core - 1][input_core - 1].busyports[WEST][INPUT].busytime = max(resource_matrix[input_core - 1][input_core - 1].busyports[WEST][INPUT].busytime, io_busytime);
            resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime = max(resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime, io_busytime);
            max_busytime = max (resource_matrix[input_core - 1][input_core - 1].busyports[WEST][INPUT].busytime, resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime);
            resource_matrix[input_core - 1][input_core - 1].busyports[WEST][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[input_core - 1][input_core - 1].busyports[WEST][INPUT].port = INJECTION;
            resource_matrix[input_core - 1][input_core - 1].busyports[INJECTION][OUTPUT].port = WEST;
            if (resource_matrix[input_core - 1][input_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[input_core - 1][input_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[input_core - 1][input_core - 1].busytime);

            // From (input core - 1) to (input core - (difference in x-coordinates - 1))
            // * Index of input core - 1 is input_core - 2
            for (int j = 0; j < (noc_nodes[input_core - 1].x_cord - noc_nodes[test_core - 1].x_cord - 1); j++) {

                // Update link busytimes between adjacent cores (link direction -x)
                // resource_matrix[input_core - j - 2][input_core - j - 3].busytime += individual_testtime;
                resource_matrix[input_core - j - 2][input_core - j - 3].busytime += max(resource_matrix[input_core - j - 2][input_core - j - 3].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[input_core - j - 2][input_core - j - 3].busytime);

                // Update ROUTER busytimes and statuses (routing direction -x i.e. input port EAST to output port WEST)
                resource_matrix[input_core - j - 2][input_core - j - 2].busyports[WEST][INPUT].busytime = max(resource_matrix[input_core - j - 2][input_core - j - 2].busyports[WEST][INPUT].busytime, io_busytime);
                resource_matrix[input_core - j - 2][input_core - j - 2].busyports[EAST][OUTPUT].busytime = max(resource_matrix[input_core - j - 2][input_core - j - 2].busyports[EAST][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[input_core - j - 2][input_core - j - 2].busyports[WEST][INPUT].busytime, resource_matrix[input_core - j - 2][input_core - j - 2].busyports[EAST][OUTPUT].busytime);
                resource_matrix[input_core - j - 2][input_core - j - 2].busyports[WEST][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[input_core - j - 2][input_core - j - 2].busyports[EAST][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[input_core - j - 2][input_core - j - 2].busyports[WEST][INPUT].port = EAST;
                resource_matrix[input_core - j - 2][input_core - j - 2].busyports[EAST][OUTPUT].port = WEST;
                if (resource_matrix[input_core - j - 2][input_core - j - 2].busytime < max_busytime + individual_testtime)
                    resource_matrix[input_core - j - 2][input_core - j - 2].busytime = max_busytime + individual_testtime;

               // Update test endtime
               endtime = max(endtime, resource_matrix[input_core - j - 2][input_core - j - 2].busytime);
            }
        }
                
        else
            xrouting_port_ic = EJECTION;

        // If test core x coordinate < output core x coordinate, update busytime for resources b/w cores
        if (noc_nodes[test_core - 1].x_cord < noc_nodes[output_core - 1].x_cord) {

            xrouting_port_co = EAST;

            // Update LINK busytimes between the test core and its adjacent core (link direction +x)
            // resource_matrix[test_core - 1][test_core].busytime += individual_testtime;
            resource_matrix[test_core - 1][test_core].busytime += max(resource_matrix[test_core - 1][test_core].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[test_core - 1][test_core].busytime);

            // Update ROUTER busytimes and statuses (routing direction +x i.e. input port INJECTION to output port EAST)
            resource_matrix[test_core - 1][test_core - 1].busyports[EAST][INPUT].busytime = max(resource_matrix[test_core - 1][test_core - 1].busyports[EAST][INPUT].busytime, io_busytime);
            resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime = max(resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime, io_busytime);
            max_busytime = max (resource_matrix[test_core - 1][test_core - 1].busyports[EAST][INPUT].busytime, resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime);
            resource_matrix[test_core - 1][test_core - 1].busyports[EAST][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[test_core - 1][test_core - 1].busyports[EAST][INPUT].port = INJECTION;
            resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].port = EAST;
            if (resource_matrix[test_core - 1][test_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[test_core - 1][test_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[test_core - 1][test_core - 1].busytime);

            // From (test core + 1) to (test core + (difference in x-coordinates - 1))
            // * Index of test core + 1 is input_core
            for (int j = 0; j < (noc_nodes[output_core - 1].x_cord - noc_nodes[test_core - 1].x_cord - 1); j++) {

                // Update LINK busytimes between adjacent cores (link direction +x)
                // resource_matrix[test_core + j][test_core + j + 1].busytime += individual_testtime;
                resource_matrix[test_core + j][test_core + j + 1].busytime += max(resource_matrix[test_core + j][test_core + j + 1].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[test_core + j][test_core + j + 1].busytime);

                // Update ROUTER busytimes and statuses (routing direction +x i.e. input port WEST to output port EAST)
                resource_matrix[test_core + j][test_core + j].busyports[EAST][INPUT].busytime = max(resource_matrix[test_core + j][test_core + j].busyports[EAST][INPUT].busytime, io_busytime);
                resource_matrix[test_core + j][test_core + j].busyports[WEST][OUTPUT].busytime = max(resource_matrix[test_core + j][test_core + j].busyports[WEST][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[test_core + j][test_core + j].busyports[EAST][INPUT].busytime, resource_matrix[test_core + j][test_core + j].busyports[WEST][OUTPUT].busytime);
                resource_matrix[test_core + j][test_core + j].busyports[EAST][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[test_core + j][test_core + j].busyports[WEST][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[test_core + j][test_core + j].busyports[EAST][INPUT].port = WEST;
                resource_matrix[test_core + j][test_core + j].busyports[WEST][OUTPUT].port = EAST;
                if (resource_matrix[test_core + j][test_core + j].busytime < max_busytime + individual_testtime)
                    resource_matrix[test_core + j][test_core + j].busytime = max_busytime + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[test_core + j][test_core + j].busytime);
            }
        }

        // If test core x coordinate > output core x coordinate, update busytime for links b/w cores
        else if (noc_nodes[test_core - 1].x_cord > noc_nodes[output_core - 1].x_cord) {
                
            xrouting_port_co = WEST;

            // Update LINK busytimes between the test core and its adjacent core (link direction -x)
            // resource_matrix[test_core - 1][test_core - 2].busytime += individual_testtime;
            resource_matrix[test_core - 1][test_core - 2].busytime += max(resource_matrix[test_core - 1][test_core - 2].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[test_core - 1][test_core - 2].busytime);

            // Update ROUTER busytimes and statuses (routing direction -x i.e. input port INJECTION to output port WEST)
            resource_matrix[test_core - 1][test_core - 1].busyports[WEST][INPUT].busytime = max(resource_matrix[test_core - 1][test_core - 1].busyports[WEST][INPUT].busytime, io_busytime);
            resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime = max(resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime, io_busytime);
            max_busytime = max (resource_matrix[test_core - 1][test_core - 1].busyports[WEST][INPUT].busytime, resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime);
            resource_matrix[test_core - 1][test_core - 1].busyports[WEST][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[test_core - 1][test_core - 1].busyports[WEST][INPUT].port = INJECTION;
            resource_matrix[test_core - 1][test_core - 1].busyports[INJECTION][OUTPUT].port = WEST;
            if (resource_matrix[test_core - 1][test_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[test_core - 1][test_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[test_core - 1][test_core - 1].busytime);

            // From (test core - 1) to (test core - (difference in x-coordinates - 1))
            // * Index of test core - 1 is input_core - 2
            for (int j = 2; j < (noc_nodes[test_core - 1].x_cord - noc_nodes[output_core - 1].x_cord - 1); j++) {

                // Update link busytimes between adjacent cores (link direction -x)
                // resource_matrix[test_core - j - 2][test_core - j - 3].busytime += individual_testtime;
                resource_matrix[test_core - j - 2][test_core - j - 3].busytime += max(resource_matrix[test_core - j - 2][test_core - j - 3].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[test_core - j - 2][test_core - j - 3].busytime);

                // Update ROUTER busytimes and statuses (routing direction -x i.e. input port EAST to output port WEST)
                resource_matrix[test_core - j - 2][test_core - j - 2].busyports[WEST][INPUT].busytime = max(resource_matrix[test_core - j - 2][test_core - j - 2].busyports[WEST][INPUT].busytime, io_busytime);
                resource_matrix[test_core - j - 2][test_core - j - 2].busyports[EAST][OUTPUT].busytime = max(resource_matrix[test_core - j - 2][test_core - j - 2].busyports[EAST][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[test_core - j - 2][test_core - j - 2].busyports[WEST][INPUT].busytime, resource_matrix[test_core - j - 2][test_core - j - 2].busyports[EAST][OUTPUT].busytime);
                resource_matrix[test_core - j - 2][test_core - j - 2].busyports[WEST][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[test_core - j - 2][test_core - j - 2].busyports[EAST][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[test_core - j - 2][test_core - j - 2].busyports[WEST][INPUT].port = EAST;
                resource_matrix[test_core - j - 2][test_core - j - 2].busyports[EAST][OUTPUT].port = WEST;
                if (resource_matrix[test_core - j - 2][test_core - j - 2].busytime < max_busytime + individual_testtime)
                    resource_matrix[test_core - j - 2][test_core - j - 2].busytime = max_busytime + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[test_core - j - 2][test_core - j - 2].busytime);
            }
        }
                
        else
            xrouting_port_co = EJECTION;

        // ------------------------------------------------------------------ Y ROUTING LOGIC ------------------------------------------------------------------
            
        // Core in XY routing intermediate step (i.e. after traversing X part)
        ic_core = input_core + (noc_nodes[test_core - 1].x_cord - noc_nodes[input_core - 1].x_cord);    // ic_core y coordinate = input core y coordinate
        co_core = test_core + (noc_nodes[output_core - 1].x_cord - noc_nodes[test_core - 1].x_cord);    // co_core y coordinate = test core y coordinate

        // If intermediate core y coordinate < test core y coordinate, update busytime for resources b/w cores 
        if (noc_nodes[input_core - 1].y_cord < noc_nodes[test_core - 1].y_cord) {

            // Update LINK busytimes between the input core and its adjacent core (link direction +y)
            // resource_matrix[ic_core - 1][ic_core - 1 + N_columns].busytime += individual_testtime;
            resource_matrix[ic_core - 1][ic_core - 1 + N_columns].busytime += max(resource_matrix[ic_core - 1][ic_core - 1 + N_columns].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[ic_core - 1][ic_core - 1 + N_columns].busytime);

            // Update ROUTER busytimes and statuses (routing direction +y i.e. input port EAST/WEST (xrouting_port_ic) to output port NORTH)
            resource_matrix[ic_core - 1][ic_core - 1].busyports[NORTH][INPUT].busytime = max(resource_matrix[ic_core - 1][ic_core - 1].busyports[NORTH][INPUT].busytime, io_busytime);
            resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime = max(resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime, io_busytime);
            max_busytime = max (resource_matrix[ic_core - 1][ic_core - 1].busyports[NORTH][INPUT].busytime, resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime);
            resource_matrix[ic_core - 1][ic_core - 1].busyports[NORTH][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[ic_core - 1][ic_core - 1].busyports[NORTH][INPUT].port = xrouting_port_ic;
            resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].port = NORTH;
            if (resource_matrix[ic_core - 1][ic_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[ic_core - 1][ic_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[ic_core - 1][ic_core - 1].busytime);

            // From (intermediate core + N columns) to (intermediate core + N columns * (difference in y-coordinates - 1))
            // * Index of intermediate core + N columns is intermediate_core - 1 + N columns
            for (int j = 0; j < (noc_nodes[test_core - 1].y_cord - noc_nodes[input_core - 1].y_cord - 1); j++) {

                // Update LINK busytimes between adjacent cores (link direction +y)
                // resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 2) * N_columns].busytime += individual_testtime;
                resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 2) * N_columns].busytime += max(resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 2) * N_columns].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 2) * N_columns].busytime);

                // Update ROUTER busytimes and statuses (routing direction +y i.e. input port SOUTH to output port NORTH))
                resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime = max(resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime, io_busytime);
                resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime = max(resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime, resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime);
                resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].port = SOUTH;
                resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].port = NORTH;
                if (resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busytime < max_busytime + individual_testtime)
                    resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busytime = max_busytime + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[ic_core - 1 + (j + 1) * N_columns][ic_core - 1 + (j + 1) * N_columns].busytime);
            }
        } 

        // If input core y coordinate > test core y coordinate, update busytime for resources b/w cores 
        else if (noc_nodes[input_core - 1].y_cord > noc_nodes[test_core - 1].y_cord) {
                
            // Update LINK busytimes between the input core and its adjacent core (link direction -y)
            // resource_matrix[ic_core - 1][ic_core - 1 - N_columns].busytime += individual_testtime;
            resource_matrix[ic_core - 1][ic_core - 1 - N_columns].busytime += max(resource_matrix[ic_core - 1][ic_core - 1 - N_columns].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[ic_core - 1][ic_core - 1 - N_columns].busytime);

            // Update ROUTER busytimes and statuses (routing direction -y i.e. input port EAST/WEST (xrouting_port_ic) to output port SOUTH)
            resource_matrix[ic_core - 1][ic_core - 1].busyports[SOUTH][INPUT].busytime = max(resource_matrix[ic_core - 1][ic_core - 1].busyports[SOUTH][INPUT].busytime, io_busytime);
            resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime = max(resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime, io_busytime);
            max_busytime = max (resource_matrix[ic_core - 1][ic_core - 1].busyports[SOUTH][INPUT].busytime, resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime);
            resource_matrix[ic_core - 1][ic_core - 1].busyports[SOUTH][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[ic_core - 1][ic_core - 1].busyports[SOUTH][INPUT].port = xrouting_port_ic;
            resource_matrix[ic_core - 1][ic_core - 1].busyports[xrouting_port_ic][OUTPUT].port = SOUTH;
            if (resource_matrix[ic_core - 1][ic_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[ic_core - 1][ic_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[ic_core - 1][ic_core - 1].busytime);

            // From (intermediate core - N columns) to (intermediate core - N columns * (difference in y-coordinates - 1))
            // * Index of intermediate core - N columns is intermediate_core - 1 - N columns
            for (int j = 0; j < (noc_nodes[input_core - 1].y_cord - noc_nodes[test_core - 1].y_cord - 1); j++) {

                // Update LINK busytimes between adjacent cores (link direction -y)
                // resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 2) * N_columns].busytime += individual_testtime;
                resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 2) * N_columns].busytime += max(resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 2) * N_columns].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 2) * N_columns].busytime);

                // Update ROUTER busytimes and statuses (routing direction -y i.e. input port NORTH to output port SOUTH))
                resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime = max(resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime, io_busytime);
                resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime = max(resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime, resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime);
                resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].port = NORTH;
                resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].port = SOUTH;
                if (resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busytime < max_busytime + individual_testtime)
                    resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busytime = max_busytime + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[ic_core - 1 - (j + 1) * N_columns][ic_core - 1 - (j + 1) * N_columns].busytime);
            }
        }

        // If test core y coordinate < output core y coordinate, update busytime for resources b/w cores 
        if (noc_nodes[test_core - 1].y_cord < noc_nodes[output_core - 1].y_cord) {

            // Update LINK busytimes between the input core and its adjacent core (link direction +y)
            // resource_matrix[co_core - 1][co_core - 1 + N_columns].busytime += individual_testtime;
            resource_matrix[co_core - 1][co_core - 1 + N_columns].busytime += max(resource_matrix[co_core - 1][co_core - 1 + N_columns].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[co_core - 1][co_core - 1 + N_columns].busytime);

            // Update ROUTER busytimes and statuses (routing direction +y i.e. input port EAST/WEST (xrouting_port_co) to output port NORTH)
            resource_matrix[co_core - 1][co_core - 1].busyports[NORTH][INPUT].busytime = max(resource_matrix[co_core - 1][co_core - 1].busyports[NORTH][INPUT].busytime, io_busytime);
            resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime = max(resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime, io_busytime);
            max_busytime = max (resource_matrix[co_core - 1][co_core - 1].busyports[NORTH][INPUT].busytime, resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime);
            resource_matrix[co_core - 1][co_core - 1].busyports[NORTH][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[co_core - 1][co_core - 1].busyports[NORTH][INPUT].port = xrouting_port_co;
            resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].port = NORTH;
            if (resource_matrix[co_core - 1][co_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[co_core - 1][co_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[co_core - 1][co_core - 1].busytime);

            // From (intermediate core + N columns) to (intermediate core + N columns * (difference in y-coordinates - 1))
            // * Index of intermediate core + N columns is intermediate_core - 1 + N columns
            for (int j = 0; j < (noc_nodes[output_core - 1].y_cord - noc_nodes[test_core - 1].y_cord - 1); j++) {

                // Update LINK busytimes between adjacent cores (link direction +y)
                // resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 2) * N_columns].busytime += individual_testtime;
                resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 2) * N_columns].busytime += max(resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 2) * N_columns].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 2) * N_columns].busytime);

                // Update ROUTER busytimes and statuses (routing direction +y i.e. input port SOUTH to output port NORTH))
                resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime = max(resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime, io_busytime);
                resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime = max(resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime, resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime);
                resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[NORTH][INPUT].port = SOUTH;
                resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busyports[SOUTH][OUTPUT].port = NORTH;
                if (resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busytime < max_busytime + individual_testtime)
                    resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busytime = max_busytime + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[co_core - 1 + (j + 1) * N_columns][co_core - 1 + (j + 1) * N_columns].busytime);
            }
        }

        // If test core y coordinate > output core y coordinate, update busytime for resources b/w cores 
        else if (noc_nodes[test_core - 1].y_cord > noc_nodes[output_core - 1].y_cord) {
            
            // Update LINK busytimes between the input core and its adjacent core (link direction -y)
            // resource_matrix[co_core - 1][co_core - 1 - N_columns].busytime += individual_testtime;
            resource_matrix[co_core - 1][co_core - 1 - N_columns].busytime += max(resource_matrix[co_core - 1][co_core - 1 - N_columns].busytime, io_busytime) + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[co_core - 1][co_core - 1 - N_columns].busytime);

            // Update ROUTER busytimes and statuses (routing direction -y i.e. input port EAST/WEST (xrouting_port_co) to output port SOUTH)
            resource_matrix[co_core - 1][co_core - 1].busyports[SOUTH][INPUT].busytime = max(resource_matrix[co_core - 1][co_core - 1].busyports[SOUTH][INPUT].busytime, io_busytime);
            resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime = max(resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime, io_busytime);
            max_busytime = max (resource_matrix[co_core - 1][co_core - 1].busyports[SOUTH][INPUT].busytime, resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime);
            resource_matrix[co_core - 1][co_core - 1].busyports[SOUTH][INPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].busytime = max_busytime + individual_testtime;
            resource_matrix[co_core - 1][co_core - 1].busyports[SOUTH][INPUT].port = xrouting_port_co;
            resource_matrix[co_core - 1][co_core - 1].busyports[xrouting_port_co][OUTPUT].port = SOUTH;
            if (resource_matrix[co_core - 1][co_core - 1].busytime < max_busytime + individual_testtime)
                resource_matrix[co_core - 1][co_core - 1].busytime = max_busytime + individual_testtime;

            // Update test endtime
            endtime = max(endtime, resource_matrix[co_core - 1][co_core - 1].busytime);

            // From (intermediate core - N columns) to (intermediate core - N columns * (difference in y-coordinates - 1))
            // * Index of intermediate core - N columns is intermediate_core - 1 - N columns
            for (int j = 0; j < (noc_nodes[test_core - 1].y_cord - noc_nodes[output_core - 1].y_cord - 1); j++) {

                // Update LINK busytimes between adjacent cores (link direction -y)
                // resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 2) * N_columns].busytime += individual_testtime;
                resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 2) * N_columns].busytime += max(resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 2) * N_columns].busytime, io_busytime) + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 2) * N_columns].busytime);

                // Update ROUTER busytimes and statuses (routing direction +y i.e. input port NORTH to output port SOUTH))
                resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime = max(resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime, io_busytime);
                resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime = max(resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime, io_busytime);
                max_busytime = max (resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime, resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime);
                resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].busytime = max_busytime + individual_testtime;
                resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[SOUTH][INPUT].port = NORTH;
                resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busyports[NORTH][OUTPUT].port = SOUTH;
                if (resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busytime < max_busytime + individual_testtime)
                    resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busytime = max_busytime + individual_testtime;

                // Update test endtime
                endtime = max(endtime, resource_matrix[co_core - 1 - (j + 1) * N_columns][co_core - 1 - (j + 1) * N_columns].busytime);
            }
        }

        // Update test starttime
        starttime = endtime - individual_testtime;

        // Update list of all starttimes and endtimes
        update_times_list(&times_head, starttime);
        update_times_list(&times_head, endtime);

        // Update IO list schedule
        update_IO_list(io_pairs[io_pair - 1].io_head, starttime, endtime, test_core);

        // Printing the resource matrix with calculated busytimes
        printf("\n\n");
        for (int i = 0; i < num_cores; i++) {
            for (int j = 0; j < num_cores; j++) {
                printf(" | %.2lf", resource_matrix[i][j].busytime);
                if (pso_particle->testtime < resource_matrix[i][j].busytime)
                    pso_particle->testtime = resource_matrix[i][j].busytime;
            } 
            printf(" |\n\n");
        }
    }

    // Printing the resource matrix with calculated busytimes
    printf("\n\n");
    for (int i = 0; i < num_cores; i++) {
        for (int j = 0; j < num_cores; j++) {
            printf(" | %.2lf", resource_matrix[i][j].busytime);
            if (pso_particle->testtime < resource_matrix[i][j].busytime)
                pso_particle->testtime = resource_matrix[i][j].busytime;
        } 
        printf(" |\n\n");
    }

    printf(" Total testtime for given mapping: %lf\n", pso_particle->testtime);
}

// void create_clap_input_list (Clap_inputs_head *head, IO_pairs *io_pairs, int num_io_pairs) {
    // Clap_list *temp = *head;

    // new_node->next = NULL;

    // if (*head == NULL) {
    //     Clap_list *new_node = (Clap_list*) malloc (sizeof (Clap_list));
    //     *head = new_node;
    //     new_node->clap_inputs[io_pair - 1].starttime = starttime;
    //     new_node->clap_inputs[io_pair - 1].endtime = endtime;
    //     new_node->clap_inputs[io_pair - 1].input_core = input_core;
    //     new_node->clap_inputs[io_pair - 1].test_core = test_core;
    //     new_node->clap_inputs[io_pair - 1].output_core = output_core;
    //     new_node->clap_inputs[io_pair - 1].filled = FILLED;
    // }
    // else {
    //     while (temp->next != NULL) {
    //         for (int i = 0; i < num_io_pairs; i++){
    //             if (new_node->clap_inputs[i].starttime < new_node->clap_inputs[io_pair - 1].starttime)
    //         }
    //         temp = temp->next;
    //     }
    // }
// }

// double calculate_worst_case_SNR (PSO_particle *pso_particle, NoC_node *noc_nodes, int N_columns, int num_cores, IO_pairs *io_pairs, int num_io_pairs) {
//     int num_test_cores = num_cores - (2 * num_io_pairs);    // Number of test cores in the NoC mesh network
//     int test_core = 0;                                      // Temporary variable to store test core number
//     int input_core = 0;                                     // Temporary variable to store input core number
//     int output_core = 0;                                    // Temporary variable to store output core number
//     double starttime = 0.0;                                 // Temporary variable to store test starttime
//     double endtime = 0.0;                                   // Temporary variable to store test endtime
//     int io_pair = 0;
//     Clap_list *head = NULL;                                 // First node of the clap inputs list

//     // Populating resource matrix with busytime (latest time till which the resource is busy) for all test cores
//     for (int i = 0; i < num_test_cores; i++) {

//         create_clap_input_list (&head, starttime, endtime, input_core, test_core, output_core, num_io_pairs, io_pair);     

        // CALL CLAP FUNCTIONS
// }

// Initializes PSO particles by initializing the I/O core, frequencies and test core mapping; calculates the fitness value for each particle

void init_pso_particles (PSO_particle *pso_particle, Gbest_PSO_particle *gbest_pso_particle, NoC_node *noc_nodes, int num_cores, double *freq, int num_freq, IO_pairs *io_pairs, int num_io_pairs, int N_columns) {

    int num_test_cores = num_cores - (2 * num_io_pairs);       // Number of test cores in the NoC mesh network
    int *temp_arr;                                             // Temporary array to store test core ids
    int ii = 0;                                                // Index to traverse through the temporary array
    int p = 0;                                                 // Index to traverse through PSO particle struct array
    int i = 0;                                                 // Index to traverse through NoC nodes struct array
    int j = 0;                                                 // Index to traverse through the mapping 
    double best_fitness = 128797218.0;                         // Temporary variable to store best fitness value
    int best_idx = 0;                                          // Temporary variable to store index of the particle with best fitness value
    int io_idx = 0;                                            // Temporary variable to store io pairs array index

    // Store all test core numbers temporarily in an array - for ease of access
    temp_arr = (int *) malloc (num_test_cores * sizeof (int));
    for (i = 0; i < num_cores; i++) {
        if (noc_nodes[i].core_type == TEST_CORE) { 
            temp_arr[ii] = noc_nodes[i].core_no;
            ii++;
        }
    }

    // For all particles
    for (p = 0; p < NUM_PSO_PARTICLES; p++) {

        // Create schedule lists for all particles
        // pso_particle[p].schedule = create_schedule_list();

        // Initialize mapping  
        // mapping  structure: | test core ids | io pairs assigned | test frequency | preemptions |

        // Initialize mapping  test core field to UNALLOCATED
        for (j = 0; j < num_test_cores; j++)
            pso_particle[p].mapping[j] = UNALLOCATED; 

        // Randomly assign test core nos. sequence to the first field -- gives order of testing
        for (ii = 0; ii < num_test_cores; ii++) {

            // Pick a random UNALLOCATED index for allocating a test core number
            do {
                j = rand() % num_test_cores;
            } while (pso_particle[p].mapping[j] != UNALLOCATED);

            pso_particle[p].mapping[j] = temp_arr[ii];
        }

        // Randomly assign io pairs to test cores
        for (j = num_test_cores; j < (2 * num_test_cores); j++)
            pso_particle[p].mapping[j] = (rand() % num_io_pairs) + 1;

        // Randomly assign a valid test frequency   
        for (j = (2 *num_test_cores); j < (3 * num_test_cores); j++) {
            pso_particle[p].mapping[j] = freq[(rand() % num_freq)];
        }

        // Assign preemption points - randomly generated value between 0 and 1
        for (j = 3 * num_test_cores; j < (4 * num_test_cores); j++)
            pso_particle[p].mapping[j] = drand48 () * 0.9 + 0.1;

        find_resource_busytimes (pso_particle, noc_nodes, N_columns, num_cores, io_pairs, num_io_pairs);
        // find SNR
        // fitness = w * (testtime) + (1 - w) * SNR;

        // Updating the best fitness and the corresponding particle index
        if (best_fitness > pso_particle[p].fitness) {
            best_fitness = pso_particle[p].fitness;
            best_idx = p;
        }
    }

    // Set the local best parameters
    for (p = 0; p < NUM_PSO_PARTICLES; p++) {

        // Local best mapping  - same as the initialized mapping 
        for (j = 0; j < 4 * num_test_cores; j++)
            pso_particle[p].lbest_mapping[j] = pso_particle[p].mapping[j];

        // Local best fitness value - same as the fitness value calculated for initialized mappings
        pso_particle[p].lbest_fitness = pso_particle[p].fitness;
    }

    // Set the global best particle parameters
    // Global best mapping 
    for (j = 0; j < 4 * num_test_cores; j++) {
        gbest_pso_particle->gbest_mapping[j] = pso_particle[best_idx].mapping[j];
        gbest_pso_particle->gbest_fitness = best_fitness;   
    }

    free (temp_arr);
}

// Generates random number between 0 and 1 (probability distribution: uniform)

double generate_random_number () {
    double x = 0.0; 
    int y = rand(); 
    x = (y % 1009) + 1;
    x = x / 1009;
    return x;
}

// Swaps io pairs with given probability 

void swap_io_pair (int num_test_cores, double *a, double *b, double probability) {
    double temp = 0.0;
    int i = 0;

    // For all test cores
    for (i = 0; i < num_test_cores; i++) {
        
        // Generate a random probability value between 0 and 1
        temp = generate_random_number();
        
        // If the given probability exceeds or equals this value -- SWAP
        if (temp <= probability) 
            a[num_test_cores + i] = b[num_test_cores + i];
    }
}

// Checks frequency validity for newly assigned test core, swaps frequencies

void swap_frequencies (int num_test_cores, double *a, double *b, double probability) {
    double temp = 0.0;
    int i = 0;

    // For all test cores
    for (i = 0; i < num_test_cores; i++) {

        // Generate a random probability value between 0 and 1
        temp = generate_random_number(); 

        // If the given probability exceeds or equals this value, and the frequency to be 
        // assigned lies within the valid range for this particle -- SWAP
        if (temp <= probability) {
            a[2 * num_test_cores + i] = b[2 * num_test_cores + i];
        }
    }
}

// Generates a sequence of swap operators for evolving a given particle's test core sequence 

int generate_swap_operator_sequence (int num_test_cores, double *a, double *b, Swap_operator* swap_operator) {
    int num_swap_operators = 0;          // Number of swap operators in the swap sequence
    double temp[num_test_cores];
    int i = 0;
    int j = 0;

    // Copy a (particle test core sequence) into temp
    for(i = 0; i < num_test_cores; i++)
        temp[i] = a[i]; 

    // Traverse the entire test core sequence of the particle and its local best
    for (i = 0; i < num_test_cores; i++) {

        // If the corresponding ith elements don't match
        if (temp[i] != b[i]) {
            
            // Find the jth element in temp that matches the ith element 
            // of local best particle's test core sequence 
            j = i;
            while (temp[j] != b[i]) 
                j++;

            // (Update temp for finding indices for modified sequence in the next iteration)
            temp[j] = temp[i];
            temp[i] = b[i];

            // We can now swap the ith element with jth element in the particle's
            // test core sequence to modify it to its local best sequence
            // Add these indices to the swap operator sequence
            swap_operator[num_swap_operators].swap_idx1 = i;
            swap_operator[num_swap_operators].swap_idx2 = j;

            num_swap_operators++;
        }
    }
    return num_swap_operators;
}

// Applies the sequence of swap operators on test core sequence a with give probability

void swap_test_core_sequence (int num_test_cores, double *a, Swap_operator* swap_operator, int num_swap_operators, double probability) {
    double temp;
    int swap_idx1 = 0;
    int swap_idx2 = 0;
    int temp1 = 0;
    int i = 0;

    for (i = 0; i < num_swap_operators; i++) {

        // Generate a random probability value between 0 and 1
        temp = generate_random_number();
        
        // If the given probability exceeds or equals this value -- SWAP
        if (temp <= probability) {
            swap_idx1 = swap_operator[i].swap_idx1;
            swap_idx2 = swap_operator[i].swap_idx2;
            temp1 = a[(swap_operator[i].swap_idx1)];
            a[(swap_operator[i].swap_idx1)] = a[(swap_operator[i].swap_idx2)];
            a[(swap_operator[i].swap_idx2)] = temp1;
        }
    }
}

// Modifies the preemption points for test cores in a given particle 
// (new position of a particle in continuous PSO)

// void modify_preemption_points (int num_test_cores, double *a, double *b, double *c) {
//     int c1 = 1;
//     int c2 = 1;
//     int i = 0;
//     double r1 = 0.0;
//     double r2 = 0.0;

//     for(i = (3 * num_test_cores); i < (4 * num_test_cores); i++) {	
//         r1 = drand48();
//         r2 = drand48();

//         if(r1 > ALPHA && r2 > BETA) {
//             a[i] = a[i] + (c1 * r1 * fabs(b[i] - a[i])) + (c2 * r2 * fabs (c[i] - a[i]));
//                 if (a[i] > 2.1f)
//                     a[i] = a[i] - 2;
//                 if (a[i] > 2.0f)
//                     a[i] = a[i] - 1.9;
//                 if (a[i] > 1.1f)
//                     a[i] = a[i] - 1;
//                 if (a[i]>1.0f)
//                     a[i] = a[i] - 0.9;
//         }

//         else if (r1 > ALPHA && r2 < BETA) {
//             a[i] = a[i] + (c1 * r1 * fabs(b[i] - a[i]));
//             if (a[i]>1.1f)
//                 a[i] = a[i] - 1;
//             if (a[i]>1.0f)
//                 a[i] = a[i] - 0.9;
//         }

//         else if (r1 < ALPHA && r2 > BETA) {
//             a[i] = a[i] + (c2 * r2 * fabs(c[i] - a[i]));
//             if (a[i] > 1.1f)
//                 a[i] = a[i] - 1;
//             if (a[i] > 1.0f)
//                 a[i] = a[i] - 0.9;
//         }
//     }
// }


// Particle swarm optimization

void particle_swarm_optimization (NoC_node *noc_nodes, int num_cores, int M_rows, double *freq, int num_freq, IO_pairs *io_pairs, int num_io_pairs) {
    int N_columns = num_cores / M_rows;                                 // Number of columns in NoC mesh network
    int num_test_cores = num_cores - (2 * num_io_pairs);                // Number of test cores in the NoC mesh network
    int num_pso_runs = 0;                                               // Number of PSO runs
    int max_generations = 0;                                            // Maximum number of generations
    PSO_particle *pso_particle;                                         // PSO particle struct array
    pso_particle = malloc (NUM_PSO_PARTICLES * sizeof (PSO_particle));  // Allocating memory for PSO particle struct array
    Gbest_PSO_particle gbest_pso_particle;                              // Global best PSO particle
    Swap_operator swap_operator[num_test_cores];                        // Sequence of swap operators for a given particle
                                                                        // The swap operator exchanges the values at positions
                                                                        // (swap_idx1, swap_idx2) generate a new particle
    int num_swap_operations = 0;                                        // Number of swap operators in the swap sequence


    // Initializes the PSO particles with randomized mapping, calculates respective costs and sets the initial local and global best
    init_pso_particles (pso_particle, (&gbest_pso_particle), noc_nodes, num_cores, freq, num_freq, io_pairs, num_io_pairs, N_columns);

    print_pso_particle_info (pso_particle, num_test_cores);
    print_global_best_info ((&gbest_pso_particle), num_test_cores);

    // while (1) {
        for (int p = 0; p < NUM_PSO_PARTICLES; p++) {

            swap_io_pair (num_test_cores, pso_particle[p].mapping, pso_particle[p].lbest_mapping, ALPHA);
            swap_io_pair (num_test_cores, pso_particle[p].mapping, gbest_pso_particle.gbest_mapping, BETA);
            
            swap_frequencies (num_test_cores, pso_particle[p].mapping, pso_particle[p].lbest_mapping, ALPHA);
            swap_frequencies (num_test_cores, pso_particle[p].mapping, gbest_pso_particle.gbest_mapping, BETA);

            num_swap_operations = generate_swap_operator_sequence (num_test_cores, pso_particle[p].mapping, pso_particle[p].lbest_mapping, swap_operator);
            swap_test_core_sequence (num_test_cores, pso_particle[p].mapping, swap_operator, num_swap_operations, ALPHA);

            num_swap_operations = generate_swap_operator_sequence (num_test_cores, pso_particle[p].mapping, gbest_pso_particle.gbest_mapping, swap_operator);
            swap_test_core_sequence (num_test_cores, pso_particle[p].mapping, swap_operator, num_swap_operations, BETA);
            
            // modify_preemption_points (num_test_cores, pso_particle[p].mapping, pso_particle[p].lbest_mapping, gbest_pso_particle.gbest_mapping);
        }

        for (int i = 0; i < num_io_pairs; i++)
            print_IO_schedule_lists (io_pairs[i].io_head);
    // }   
    print_pso_particle_info (pso_particle, num_test_cores);
    print_global_best_info((&gbest_pso_particle), num_test_cores);

}


// =================
// UTILITY FUNCTIONS
// =================

// Finds the maximum of two given numbers

double max (double a, double b) {
    if (a >= b)
        return a;
    else 
        return b;
}

// Prints the mapping and test schedule information for all PSO particles

void print_pso_particle_info (PSO_particle *pso_particle, int num_test_cores) {
    int p = 0;
    int i = 0;
    int j = 0;
    
    for (p = 0; p < NUM_PSO_PARTICLES; p++) {
        printf(" Particle %d\n", p + 1);

        printf(" Test core IDs: \n");
        for (i = 0; i < num_test_cores; i++)
            printf(" %d\t", (int)pso_particle[p].mapping[i]);

        printf("\n Corresponding IO pair IDs: \n");
        for (i = num_test_cores; i < (2 * num_test_cores); i++)
            printf(" %d\t", (int)pso_particle[p].mapping[i]);

        printf("\n Test frequencies (normalized): \n");
        for (i = (2 * num_test_cores); i < (3 * num_test_cores); i++)
            printf(" %.2lf\t", pso_particle[p].mapping[i]);

        printf("\n Preemption points: \n");
        for (i = (3 * num_test_cores); i < (4 * num_test_cores); i++)
            printf(" %.2lf\t", pso_particle[p].mapping[i]);

        printf("\n Particle fitness: %.2lf", pso_particle[p].fitness);
        printf("\n\n"); 
    }
}

// Prints the mapping and test schedule information for the global best PSO particle

void print_global_best_info (Gbest_PSO_particle *gbest_pso_particle, int num_test_cores) {
    int i = 0;
    printf(" Global best info\n");

    printf(" Test core IDs: ");
    for (i = 0; i < num_test_cores; i++)
        printf(" %d\t", (int)gbest_pso_particle->gbest_mapping[i]);

    printf("\n Corresponding IO pair IDs: ");
    for (i = num_test_cores; i < (2 * num_test_cores); i++)
        printf(" %d\t", (int)gbest_pso_particle->gbest_mapping[i]);

    printf("\n Test frequencies (normalized): ");
    for (i = (2 * num_test_cores); i < (3 * num_test_cores); i++)
        printf(" %.2lf\t", gbest_pso_particle->gbest_mapping[i]);

    printf("\n Preemption points: ");
    for (i = (3 * num_test_cores); i < (4 * num_test_cores); i++)
        printf(" %.2lf\t", gbest_pso_particle->gbest_mapping[i]);

    printf("\n Particle fitness: %.2lf", gbest_pso_particle->gbest_fitness);
    printf("\n\n");
}

// Print IO schedule lists

void print_IO_schedule_lists (IO_head* head) {

    IO_node *temp;
    temp = head->head_node;

    printf("\n");
    if (temp == NULL)
        printf (" The List is Empty\n");
    else {
        while (temp != NULL) {
            printf(" Test core: %d (%lf to %lf) -->",temp->test_core, temp->starttime, temp->endtime);
            temp=temp->next;
        }
        printf("\n");
    }
} 
