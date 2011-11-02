
int (*topology_fn)(signed int **, int, long);
int topology_id;

int topology0(signed int **adj, int no_neurons, long random_seed);
int topology1(signed int **adj, int no_neurons, long random_seed);
int topology2(signed int **adj, int no_neurons, long random_seed);
int topology10(signed int **adj, int no_neurons, long random_seed);
