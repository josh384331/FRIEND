# Fortran Recursive Interpolation Evaluation for N-dimensional Datasets (FRIEND)

This project contains a multi-dimensional recursive linear interpolator written in Fortran.

## How it Works

The code is organized into a module called `dataset_mod`, which defines a type called `dataset`. The `dataset` type represents a dataset with a table of values and provides various operations such as initialization, sorting, printing, and interpolation.  The main purpose of the code is to linearly interpolate multiple(currently 1) dependent variable over many independent variables.  It accomplishes this by recursivly calling the dataset_interp function until the code has calculated the weights for each independent variable, then solves the weighted average for all of the interpolations.

### Initialization

To initialize a dataset, use the `dataset_init` subroutine. It takes the following arguments:
- `data_table`: The table of data values.
- `n_IndepVars`: The number of independent variables.

### Sorting Data

To sort the data in the dataset, use the `dataset_sort_data` subroutine. It sorts the data in ascending order based on the first independent variable. 
    
this has not beed added yet, in the meantime, sort your tables assending by independant variable starting at the first independant variable. For example:

1,2,3\
1,2,4\
1,3,1\
1,3,2\
2,2,3\
2,2,4\
2,3,1\
2,3,2

where the first two columns are indpendenat and the third is dependant

### Printing Data

To print the data in the dataset, use the `dataset_print_data` subroutine. It prints each row of the table.

### Interpolation

To perform interpolation on the dataset, use the `dataset_interp` function. It takes the following arguments:
- `indep_Vars`: The independent variables.
- `i_indepVar`: The index of the independent variable to interpolate. Always choose 1 (I will make this automatic eventually)
- `rowi`: The starting row index for interpolation. This should probably be the start of your table (I will make this automatic eventually)
- `rowf`: The ending row index for interpolation. This should probably be the end of your table (I will make this automatic eventually)

The function returns the interpolated value.

## How to Operate

To use the multi-dimensional recursive linear interpolator, follow these steps:

1. Create a Fortran program and include the `dataset_mod` module.
2. Initialize a dataset object using the `dataset_init` subroutine, providing the data table and the number of independent variables.
3. Sort the data using the `dataset_sort_data` subroutine.
4. Print the data using the `dataset_print_data` subroutine (optional).
5. Perform interpolation using the `dataset_interp` function, providing the dataset object, the independent variables, the index of the independent variable to interpolate, and the row indices for interpolation.

Make sure to compile and run the Fortran program using a Fortran compiler.

## Example

Here's an example usage of the multi-dimensional recursive linear interpolator:

```fortran
