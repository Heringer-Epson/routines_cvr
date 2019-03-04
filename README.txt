To compile the Cython routine, on this directory type:
>>python setup.py build_ext --inplace

To compile with a profiler flag, add the '@profile'
decorator to the target function and run master.py as
>>kernprof -l -v master.py


