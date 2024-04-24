# Fortran Recursive Interpolation Evaluation for N-dimensional Datasets (FRIEND)

This project contains a multi-dimensional recursive linear interpolator written in Fortran.

## How it Works

The code is organized into a module called `dataset_mod`, which defines a type called `dataset`. The `dataset` type represents a dataset with a table of values and provides various operations such as initialization, sorting, printing, and interpolation.  The main purpose of the code is to linearly interpolate multiple dependent variable over many independent variables.  It accomplishes this by recursivly calling the dataset_interp function until the code has calculated the weights for each independent variable, then solves the weighted average for all of the interpolations.

If the independent variables given are outside the range of the dataset, the code will select the closest value and continue the interpolation.  It will also update the `error` variable to reflect where the interpolation requires more data.  If `verbose` is enabled, the code will also print this information.

### Initialization

To initialize a dataset, use the `dataset_init` subroutine. It takes the following arguments:
- `data_table`: The table of data values as dimention(:,:) float
- `n_IndepVars`: The number of independent variables. as an integer
- `verbose`: An optional flag to disable warning messages from the program as a logical

### Sorting Data

To sort the data in the dataset, use the `dataset_sort_data` subroutine. It sorts the data in ascending order based on the first independent variable. 
    
this has not been added yet, in the meantime, sort your tables assending by independant variable starting at the first independant variable. For example:

1,2,3\
1,2,4\
1,3,1\
1,3,2\
2,2,3\
2,2,4\
2,3,1\
2,3,2

where the first two columns are independent and the third is dependent

### Printing Data

To print the data in the dataset, use the `dataset_print_data` subroutine. It prints each row of the table.

### Interpolation

To perform interpolation on the dataset, use the `dataset_interp` function. It takes the following arguments:
- `indep_Vars`: Float array of length `n_IndepVars` containing The independent variables. 
- `i_indepVar`: integer containing the index of the independent variable to interpolate. Always choose 1 (I will make this automatic eventually)
- `rowi`: Integer containing the starting row index for interpolation. This should probably be the start of your table (I will make this automatic eventually)
- `rowf`: Integer containing The ending row index for interpolation. This should probably be the end of your table (I will make this automatic eventually)
- `error`: Integer array of length `n_IndepVars` containing error codes. 0=within range 1=above range -1=below range

The function returns the interpolated values in an array of length equal to the number of dependant variables (stored as `n_depVars`)

## How to Operate

To use the multi-dimensional recursive linear interpolator, follow these steps:

1. Create a Fortran program and include the `dataset_mod` module.
2. Initialize a dataset object using the `dataset_init` subroutine, providing the data table and the number of independent variables.
3. Perform interpolation using the `dataset_interp` function, providing the independent values and the row indices for interpolation.

note: Step 3 can be performed as many times as you want after initializing the dataset with different independent values. 

## Example

The test cases used for development is coded in `main.f90` and is included inte the repo.  In addition, here's an example usage of the multi-dimensional recursive linear interpolator:

```fortran
program FRIEND_Example
    use dataset_mod
    implicit none
    type(dataset) :: ds
    real, dimension(4,4) :: data_table = reshape([1.0, 1.0, 2.0, 2.0 , 1.0, 2.0, 1.0, 2.0, 1.0, 3.0, 1.0,3.0,4.,5.,5.,6.], [4,4])
    integer :: n_IndepVars = 2
    real, dimension(2) :: indep_Vars = [1.5, 1.5]
    real, dimension(2) :: result
    
    ! init dataset
    call ds%init(data_table, n_IndepVars)
    ! interpolate
    result = ds%interp(indep_Vars,1, 1, 4)
end program FRIEND_Example
