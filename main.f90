program test_dataset_mod
    use dataset_mod
    implicit none
    type(dataset) :: ds
    real, dimension(4,3) :: data_table = reshape([1.0, 1.0, 2.0, 2.0 , 1.0, 3.0, 1.0, 3.0, 1.0, 3.0, 1.0,3.0], [4,3])
    integer :: n_IndepVars = 2
    real, dimension(2) :: indep_Vars = [1.0, 2.]
    real :: result

    ! Test dataset_init
    call ds%init(data_table, n_IndepVars)
    print*, "Expected n_rows: 3, Actual: ", ds%n_rows
    print*, "Expected n_depVars: 1, Actual: ", ds%n_depVars
    print*, "Expected n_IndepVars: 2, Actual: ", ds%n_IndepVars
    ! test dataset print
    print*, "Printing dataset"
    call ds%print_data()
    
    ! Test dataset sort
    print*, "Sorting dataset"
    call ds%sort_data()

    ! Test dataset_interp
    result = ds%interp(2, indep_Vars)
    print*, "Expected Interpolation result: 2, Interpolation result: ", result

    ! Check if the answer is correct
   
end program test_dataset_mod