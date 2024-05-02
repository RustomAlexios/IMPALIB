# **IMPALIB Test Framework**

=============

- Unit testing input and output files are saved in "test/ut_inputs" and "test/ut_results" folders, respectively.

- In order to perform unit-testing of a certain function $F$ with a function name $N$ in a class $C$ using the bash script files "unit_test_kc_mwm.sh" or "unit_test_tsp.sh", do the following:

  - Step $1$: 
    
    - Include the function's name $N$ in the variable "unit_tests" in "unit_test_kc_mwm.sh" or "unit_test_tsp.sh" files.

  - Step $2$: 
    
    - Navigate to "python_kc_mwm/test/" or "python_tsp/test/" and add the python unit test for function $F$ in class $C$:

    - If class $C$ file exists, open "ut_class.py". If class $C$ file does not exist, create a new file "ut_class.py", and make sure to import this class in "impalib_unit_tests.py".

    - Add a unit test function "ut_function()" in "ut_class.py". Note that "ut_function()" takes some specific arguments depending on the application, and will generate Python output files.

    - Navigate to "impalib_unit_tests.py" and call the function unit test "ut_class.ut_function()" based on its function name $N$.
    
    - Example: Function $F$ "degree_constraint_to_edge_ec_update()" with function name $N$ "DegreeConstraint2EdgeEcUpdate" is already implemented for class $C$ with file "ut_update_degree_constraint.py". Notice that $F$ is first called in "impalib_unit_tests.py" based on its function name. The unit test of $F$ in "ut_update_degree_constraint.py" will generate the Python input and output files.

    - Note that "impalib_unit_tests.py" file takes some inputs to process the python unit testing. These are user specific inputs (number of nodes in TSP, etc.)

  - Step $3$: 
    
    - Navigate to "test/include/":

    - If class $C$ file exists, open "ut_class.hpp". If class $C$ file does not exist, create a new file "ut_class.hpp", and make sure to add "#include "impalib/ut_class.hpp" in "impalib_unit_tests.hpp".

    - Add a unit test function "ut_function()" in "ut_class.hpp". Note that "ut_function()" takes specific arguments, and will generate C++ output files.

    - Next, navigate to "test/src/impalib_unit_tests_kc_mwm.cpp" or "test/src/impalib_unit_tests_tsp.cpp" and call the function unit test "ut_function()".

    - Example: Function $F$ "degree_constraint_to_edge_ec_update()" with function name $N$ "DegreeConstraint2EdgeEcUpdate" is already implemented for class $C$ with file "ut_update_degree_constraint.hpp". Notice that $F$ is called in "impalib_unit_tests_tsp.cpp" based on its function name. The unit test of $F$ in "ut_update_degree_constraint.hpp" will read the input files from "test/ut_inputs/" and generate the output files in "test/ut_results/". Note that some general variables are exported in "unit_test_kc_mwm.sh" or "unit_test_tsp.sh", and these will be read in the unit test functions. These inputs are common across various unit tests.

  - Step $4$:
    - "python3 ut_methods_utils.py" takes the function, and is called to read the outputs files of the Python and C++ codes, and check agreement or not.

    - Example: Based on the above examples, "ut_methods_utils.py" will read the Python and C++ output files for function name $N$ "DegreeConstraint2EdgeEcUpdate", and check whether there is agreement or not.

- To perform unit testing of Application $1$ or $2$, navigate to the root directory:
  - Run: ``cmake -B build ``
  - Run: ``cmake --build build ``
  - Run: ``cd build/test/src``
  - Run: ``./unit_test_kc_mwm.sh or ./unit_test_tsp.sh or ./unit_test_ksat.sh``

- Note that in unit testing, we loop over all unit tests names and Multiple sub-tests could be performed depending on the value of the "total_sub_tests" variable. 

### **License**

Distributed under the MIT License.
See accompanying file [`LICENSE`](https://github.com/RustomAlexios/IMPALIB/blob/main/LICENSE) or at
([https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT))
