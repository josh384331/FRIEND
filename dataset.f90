module dataset_mod
    ! use statements

    implicit none

    type :: dataset
        real, dimension(:,:), allocatable :: table
        integer :: n_IndepVars, n_rows, n_depVars


        contains
            procedure :: init => dataset_init
            procedure :: sort_data => dataset_sort_data
            procedure :: print_data => dataset_print_data
            procedure :: interp => dataset_interp

    end type dataset




    contains 

    subroutine dataset_init(this,data_table,n_IndepVars)
        class(dataset) :: this
        real, dimension(:,:), intent(in) :: data_table
        integer, intent(in) :: n_IndepVars
        this%table = data_table
        this%n_depVars = size(data_table,2) - n_IndepVars
        this%n_IndepVars = n_IndepVars
        this%n_rows = size(data_table,1)

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


    recursive function dataset_interp(this,indep_Vars,i_indepVar,rowi,rowf) result(ans)
    !> This function performs interpolation on a dataset.
        !! Inputs:
        !!   this - The dataset object.
        !!   indep_Vars - The independent variables.
        !!   i_indepVar - The index of the independent variable to interpolate.
        !!   rowi - The starting row index for interpolation.
        !!   rowf - The ending row index for interpolation.
        !! Returns:
        !!   ans - The interpolated value.
        class(dataset) :: this
        integer, intent(in) :: i_indepVar, rowi, rowf
        real, dimension(:), intent(in) :: indep_Vars
        integer :: idx1, n, row1i, row1f, row2i, row2f
        real :: weight
        real :: ans
        logical :: found = .false.
        ! perform some checks before interpolating
        ! check to make sure the number of independent variables is correct
        if (size(indep_Vars) /= this%n_indepVars) then
            print *, 'Error: Number of independent variables does not match number of layers'
            stop
        end if
        ! check to make sure the independent variables are within the range of the table
        
        if (indep_Vars(i_indepVar) < this%table(rowi,i_indepVar) .or. indep_Vars(i_indepVar) > this%table(rowf,i_indepVar)) then
            print *, 'Error: Independent variable out of range'
            stop
        end if
        

        ! interpolate table
        !! find weight and index
        idx1 = -1
        row1i = rowi
        row1f = -1
        row2i = -1
        row2f = -1
        print*, 'Finding Weight'
        print *, 'i_indepVar (Column)', i_indepVar
        print *, 'indep_Vars(i_indepVar)', indep_Vars(i_indepVar)
        print *, 'this%table(:,i_indepVar)', this%table(:,i_indepVar)
        print *, 'rowi, rowf', rowi, rowf 
        print *, ""
        found = .false.
        ! loop through the rows of the table starting at rowi and ending at rowf
        do n = rowi, rowf
            print *, 'n, row1i, row1f, row2i, row2f,  ',n,  row1i, row1f, row2i, row2f
            print *, 'this%table(n,i_indepVar)', this%table(n,i_indepVar)
            ! if the independent variable is found, set the index and weight
            if (indep_Vars(i_indepVar) <= this%table(n,i_indepVar) .and. .not. found ) then
                print *, 'Found', n
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
            if (found .and. this%table(n,i_indepVar) /= this%table(row2i,i_indepVar)) then
                row2f = n-1
                exit
            end if
        end do

        ! if the last row is reached, set row2f to rowf
        if (row2f == -1) then
            row2f = rowf
        end if


        print *, 'idx1, weight, row1i, row1f, row2i, row2f', idx1, weight, row1i, row1f, row2i, row2f
        ! either recursively interpolate or linear interpolate
        if (i_indepVar == this%n_indepVars) then
            ! linear interpolate
            print *, 'Linear Interpolating'
            ans = this%table(idx1,this%n_indepVars+1) * (1-weight) + this%table(idx1+1,this%n_indepVars+1) * weight
        else
            ! recursively interpolate
            print *, 'Recursively Interpolating'
            ans = dataset_interp(this,indep_Vars,i_indepVar+1,row1i,row1f) * (1-weight)&
             + dataset_interp(this,indep_Vars,i_indepVar+1,row2i,row2f) * weight
        end if
    end function dataset_interp

end module dataset_mod