### Package changes from previous WIPF version 0.1.0-1

Compared to version 0.1.0-1 of WIPF, the new version includes:

* The **WIPF()** function, which is generalizes the functions **WIPF2** and **WIPF3** 
for an indeterminate number N >= 2 of dimensions.

* The **df2array()** function, which is instrumental in organizing the information contained in a long-table data frame to create inputs for the WIPF functions. It transforms a data frame with N column factors and a column of values into a K-dimensional (K <= N) array, where the size of each dimension corresponds to the number of levels in the respective factor.

* The **array2df()** function, which is instrumental in organizing the outputs of the WIPF functions. It organizes the information in an array with N dimensions into a long-format data frame with N + 1 columns, where the last column contains the values of the array.
