module dataset_mod
    ! use statements

    implicit none

    type :: dataset
        real, dimension(:,:), allocatable :: table
        integer :: n_IndepVars, n_rows, n_depVars
        integer, dimension(:), allocatable :: error
        logical :: verbose = .true.


        contains
            procedure :: init => dataset_init
            procedure :: sort_data => dataset_sort_data
            procedure :: print_data => dataset_print_data
            procedure :: interp => dataset_interp

    end type dataset




    contains 

    
    subroutine dataset_init(this, data_table, n_IndepVars, verbose)
        ! Initializes the dataset with the given data table.
            !! Inputs
            !! this - The instance of the dataset.
            !! data_table -  The data table to initialize the dataset with.
            !! n_IndepVars - The number of independent variables in the dataset.
            !! verbose -  Flag indicating whether to print verbose output.
        ! This subroutine initializes the dataset with the provided data table. It sets the number of independent variables and determines whether to print verbose output based on the value of the `verbose` flag.

        class(dataset) :: this
        real, dimension(:,:), intent(in) :: data_table
        integer, intent(in) :: n_IndepVars
        logical, optional, intent(in) :: verbose
        this%table = data_table
        this%n_depVars = size(data_table,2) - n_IndepVars
        this%n_indepVars = n_IndepVars
        this%n_rows = size(data_table,1)
        if (present(verbose)) then
            this%verbose = verbose
        end if

    end subroutine dataset_init


    subroutine dataset_sort_data(this)
        class(dataset) :: this


        call this%print_data()
        
    end subroutine dataset_sort_data

    subroutine dataset_print_data(this)
        class(dataset) :: this
        integer :: i, j

        do i = 1, this%n_rows
            write(*, '(100F10.2)') this%table(i,:)
        end do

    end subroutine dataset_print_data


    recursive function dataset_interp(this,indep_Vars,i_indepVar,rowi,rowf,error) result(ans)
    !> This function performs interpolation on a dataset.
        !! Inputs:
        !!   this - The dataset object.
        !!   indep_Vars - The independent variables.
        !!   i_indepVar - The index of the independent variable to interpolate.
        !!   rowi - The starting row index for interpolation.
        !!   rowf - The ending row index for interpolation.
        !!   error - An optional array of size n_indep_Vars  to store error codes. 0=within range 1=above range -1=below range
        !! Returns:
        !!   ans - The interpolated value.
        class(dataset) :: this
        integer, intent(in) :: i_indepVar, rowi, rowf
        real, dimension(:), intent(in) :: indep_Vars
        integer :: idx1, i, n, row1i, row1f, row2i, row2f
        real :: weight
        real, dimension(:), allocatable :: ans
        logical :: found = .false., high = .false., low = .false.
        integer, dimension(:), optional, intent(inout) :: error ! integer array of dimention i_indepVar to store error code
        
        ! perform some checks before interpolating
        !! check to make sure the number of independent variables is correct
        if (size(indep_Vars) /= this%n_indepVars) then
            print *, 'Error: Number of independent variables does not match number of layers'
            stop
        end if

        if (present(error)) then
            if (size(error) /= this%n_indepVars) then
                print *, 'Error: Size of error array is:', size(error), 'Expected:', this%n_indepVars
                stop
            end if
        end if

        !! check to make sure the independent variables are within the range of the table
        if (indep_Vars(i_indepVar) < this%table(rowi,i_indepVar) .or. indep_Vars(i_indepVar) > this%table(rowf,i_indepVar)) then
            if (this%verbose) then
                print *, 'Error: Independent variable out of range'
                print *, 'Requested Independent Variable', indep_Vars(i_indepVar)
                print *, 'Range of Independent Variable', this%table(rowi,i_indepVar), this%table(rowf,i_indepVar)
            end if
            
            if (indep_Vars(i_indepVar) < this%table(rowi,i_indepVar)) then
                high = .false.
                low = .true.
                if (present(error)) then
                    error(i_indepVar) = -1
                end if
            else 
                high = .true.
                low = .false.
                if (present(error)) then
                    error(i_indepVar) = 1
                end if
            end if
        else
            high = .false.
            low = .false.
            if (present(error)) then
                error(i_indepVar) = 0
            end if
        end if
        
        
        ! allocate ans array
        allocate(ans(this%n_depVars))
        ! interpolate table
        !! find weight and index
        idx1 = -1
        row1i = rowi
        row1f = -1
        row2i = rowi
        row2f = -1

        
        found = .false.
        ! loop through the rows of the table starting at rowi and ending at rowf
        do n = rowi, rowf
            ! print *, 'n, row1i, row1f, row2i, row2f,  ',n,  row1i, row1f, row2i, row2f
            ! print *, 'this%table(n,i_indepVar)', this%table(n,i_indepVar)
            ! if the independent variable is found, set the index and weight
            if (indep_Vars(i_indepVar) <= this%table(n,i_indepVar) .and. .not. found ) then
                ! print *, 'Found', n
                found = .true.
                idx1 = n-1
                row1f = n-1
                row2i = n
                weight = (indep_Vars(i_indepVar) - this%table(idx1,i_indepVar)) / (this%table(idx1+1,i_indepVar) &
                - this%table(idx1,i_indepVar))
            
            ! if the independent variable is not found and the value changes, set row1i to the next row
            else if (this%table(n,i_indepVar) /= this%table(row1i,i_indepVar) .and. .not. found )then
                row1i = n
            end if

            ! if the independent variable is found and the value changes, set row2f to the previous row
            if (found) then
                if (this%table(n,i_indepVar) /= this%table(row2i,i_indepVar)) then
                    row2f = n-1
                    exit
                end if
            end if
        end do

        ! if the last row is reached, set row2f to rowf
        if (row2f == -1) then
            row2f = rowf
        end if


        ! print *, 'idx1, weight, row1i, row1f, row2i, row2f', idx1, weight, row1i, row1f, row2i, row2f
        ! either recursively interpolate or linear interpolate
        if (i_indepVar == this%n_indepVars) then
            ! linear interpolate
            
            
            if (high .or. low) then
                if (high) then
                    do i = 1, this%n_depVars
                        ans(i) = this%table(row2f,i+this%n_indepVars)
                    end do
                else
                    do i = 1, this%n_depVars
                        ans(i) = this%table(row1i,i+this%n_indepVars)
                    end do
                end if
            else
                do i = 1, this%n_depVars
                    ans(i) = this%table(idx1,i+this%n_indepVars) * (1-weight) + this%table(idx1+1,i+this%n_indepVars) * weight
                end do
            end if
        else
            ! recursively interpolate
            if (high .or. low) then
                if (high) then
                    ans = dataset_interp(this,indep_Vars,i_indepVar+1,row2i,row2f)
                else
                    ans = dataset_interp(this,indep_Vars,i_indepVar+1,row1i,row2f)
                end if
            else
                ans = dataset_interp(this,indep_Vars,i_indepVar+1,row1i,row1f) * (1-weight)&
                + dataset_interp(this,indep_Vars,i_indepVar+1,row2i,row2f) * weight
            end if
        end if

    end function dataset_interp

end module dataset_mod