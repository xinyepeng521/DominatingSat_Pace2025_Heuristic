# The Heuristic Solver for the PACE 2025 Challenge
This is a solver specifically designed for the Dominating Set problem in the heuristic track of the PACE 2025 Challenge. 
It fully supports the input and output formats as defined in the challenge for this particular problem.

## Solver Description
The solver operates through the following steps:
1. **Instance Reduction**: First, the instance is simplified using reduction rules.
2. **Initial Solution Construction**: A greedy principle is then applied to construct an initial solution.
3. **Iterative Local Search**: The solver attempts to reduce the number of vertices in the solution through iterative local search to find feasible solutions.
   
   In each iteration:
   - **Vertex Removal**: A vertex is removed from either the neighborhood of the last inserted vertex (if exists) or randomly selected vertices.
   - **Vertex Insertion**:
     - Randomly select an undominated vertex and insert either itself or its highest-scoring neighbor, OR
     - Globally insert the highest-scoring vertex
   - **Feasibility Check**: When the current solution size is one less than the best feasible solution size or redundant vertices exist, insert the highest-scoring vertex if it yields a feasible solution.
4. **Redundancy Handling**: When a feasible solution is found, redundant vertices (incrementally updated each iteration) are identified and removed.
5.**Cycle Escape**: A tabu mechanism using **solution hashing** is applied to escape cycles. Crucially, the hash of tabu solutions is incrementally updated **only for removed vertices**, not for inserted vertices.

## How to Build and Run

### Build

```shell
g++ -o solver -O3 -I . *.cpp
```

### Run
The solver reads input from standard input (stdin) and writes the output to standard output (stdout).

```shell
./solver < input_file > output_file
```

## Contact

For any questions or feedback, please contact:
778845121@qq.com

