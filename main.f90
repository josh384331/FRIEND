program test_dataset_mod
    use dataset_mod
    implicit none
    type(dataset) :: ds
    real, dimension(4,3) :: data_table = reshape([1.0, 1.0, 2.0, 2.0 , 1.0, 2.0, 1.0, 2.0, 1.0, 3.0, 1.0,3.0], [4,3])
    integer :: i
    integer :: n_IndepVars = 2
    real, dimension(:,:), allocatable :: indep_Vars_list
    real, dimension(2) :: indep_Vars = [1.5, 1.5]
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
    result = ds%interp(indep_Vars,1, 1, 4)
    print*, "Expected Interpolation result: 2, Interpolation result: ", result

    ! Test dataset_interp for many points
    print*, "Interpolating many points"
    allocate(indep_Vars_list(25,2))
    ! Check if the answer is correct
    indep_Vars_list = reshape([1.,1.,1.,1.,1.,1.25,1.25,1.25,1.25,1.25,1.5,1.5,1.5,1.5,1.5,1.75,1.75,1.75,1.75,1.75,2.,2.,2.,2.,2.,&
    1.,1.25,1.5,1.75,2.,1.,1.25,1.5,1.75,2.,1.,1.25,1.5,1.75,2.,1.,1.25,1.5,1.75,2.,1.,1.25,1.5,1.75,2.,1.,1.5,2.,2.5,3.,1.,1.5,2.,&
    2.5,3.,1.,1.5,2.,2.5,3.,1.,1.5,2.,2.5,3.,1.,1.5,2.,2.5,3.], [25,3])
    do i = 1, 25
        result = ds%interp(indep_Vars_list(i,:2),1, 1, 4)
        print*, "Expected Interpolation result: ", indep_Vars_list(i,3), "Interpolation result: ", result
    end do

end program test_dataset_mod
    

    
 