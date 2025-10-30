To build the executable file, run the following command：g++ -o message_filt -O3 message_filt.cpp.
To use the kissat solver, you must replace the path to "kissat" in the system call command within message_filt.cpp with the installation path of kissat on your computer.
In SAT.cpp, we show how to construct our new SAT model with a prefix pair message.
In keccak_permutation.h, we implement deduce-and-sieve algorithm and keccak permutation.
In message_gen.cpp, we show how to generate prefix pairs fulfilling conditions with hash table.
In message_filt.cpp, we implement our whole framework with pre-filtering pairs by more conditions to find a collision.
In table.h we list our fixed conditions and 3-round differential characteristic.
The CNF file will be generated in test_sat.cnf.
The SAT solver result with running time will be generated in result.txt.
Other files are test files which we use to test correctness of code.
