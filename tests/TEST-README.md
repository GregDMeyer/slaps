SLAPS Test Suite
====

To run these tests, download `catch.hpp` from the Catch2 framework [here](https://github.com/catchorg/Catch2/releases/download/v2.2.2/catch.hpp), put it in this directory, and run `make`. This will generate three executables:

 - `matrix-tests`
 - `vector-tests`
 - `utils-tests`

Each of which can be run to do the tests.

They can be run equally well in parallel (with `upcxx-run` or sequentially. If run in parallel with more than one process, they automatically output the results to text files labeled with the process number, instead of stdout.

The Catch2 framework's repository is [here](https://github.com/catchorg/Catch2).
