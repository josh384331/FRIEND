# MultiInterp

This project contains a multi-dimensional recursive linear interpolator written in Fortran.

## How it Works

The code is organized into a module called `dataset_mod`, which defines a type called `dataset`. The `dataset` type represents a dataset with a table of values and provides various operations such as initialization, sorting, printing, and interpolation.

### Initialization

To initialize a dataset, use the `dataset_init` subroutine. It takes the following arguments:
- `data_table`: The table of data values.
- `n_IndepVars`: The number of independent variables.

### Sorting Data

To sort the data in the dataset, use the `dataset_sort_data` subroutine. It sorts the data in ascending order based on the first independent variable.

### Printing Data

To print the data in the dataset, use the `dataset_print_data` subroutine. It prints each row of the table.

### Interpolation

To perform interpolation on the dataset, use the `dataset_interp` function. It takes the following arguments:
- `this`: The dataset object.
- `indep_Vars`: The independent variables.
- `i_indepVar`: The index of the independent variable to interpolate.
- `rowi`: The starting row index for interpolation.
- `rowf`: The ending row index for interpolation.

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
