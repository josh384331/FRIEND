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


    recursive function dataset_interp(this,i_indepVar,indep_Vars) result(ans)
        class(dataset) :: this
        integer, intent(in) :: i_indepVar
        real, dimension(:), intent(in) :: indep_Vars
        integer :: idx1, i, n
        real :: weight
        real :: ans
        ! perform some checks before interpolating
        ! check to make sure the number of independent variables is correct
        if (size(indep_Vars) /= this%n_indepVars) then
            print *, 'Error: Number of independent variables does not match number of layers'
            stop
        end if
        ! check to make sure the independent variables are within the range of the table
        do i = 1, this%n_indepVars
            if (indep_Vars(i) < this%table(1,i) .or. indep_Vars(i) > this%table(this%n_rows,i)) then
                print *, 'Error: Independent variable out of range'
                stop
            end if
        end do

        ! interpolate table
        !! find weight and index
        print*, 'Finding Weight'
        print *, 'i_indepVar (Column)', i_indepVar
        print *, 'indep_Vars(i_indepVar)', indep_Vars(i_indepVar)
        print *, 'this%table(:,i_indepVar)', this%table(:,i_indepVar)
        do n = 1, this%n_rows
            if (this%table(n,i_indepVar) >= indep_Vars(i_indepVar)) then
                idx1 = n-1
                weight = (indep_Vars(i_indepVar) - this%table(idx1,i_indepVar)) / (this%table(idx1+1,i_indepVar) -&
                 this%table(idx1,i_indepVar))
                exit
            end if
        end do
        print *, 'idx1, weight', idx1, weight
        ! either recursively interpolate or linear interpolate
        if (i_indepVar == this%n_indepVars) then
            ! linear interpolate
            ans = this%table(idx1,this%n_indepVars+1) * (1-weight) + this%table(idx1+1,this%n_indepVars+1) * weight
        else
            ans = dataset_interp(this,i_indepVar+1,indep_Vars) * (1-weight) + dataset_interp(this,i_indepVar+1,indep_Vars) * weight
        end if
    end function dataset_interp

end module dataset_mod