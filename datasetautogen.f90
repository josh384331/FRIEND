
module dataset_mod
    implicit none

    type :: darray
        real, allocatable :: data(:)
    end type darray

    type :: dmatrix
        type(darray), allocatable :: rows(:)
    end type dmatrix

    type :: dataset_type
        character(len=256) :: mFilename
        integer :: mNumIndependentVars
        logical :: mHasSubset
        real, allocatable :: mX(:)
        type(dmatrix) :: mY
        character(len=256), allocatable :: mIndependentVarNames(:)
        character(len=256), allocatable :: mDependentVarNames(:)
        type(dataset_type), allocatable :: mSubset(:)
    contains
        procedure :: create_from_file
        procedure :: create_from_data
        procedure :: print_data
        procedure :: linear_interpolate
    end type dataset_type

contains

    subroutine create_from_file(this, filename)
        class(dataset_type), intent(inout) :: this
        character(len=*), intent(in) :: filename
        character(len=256), allocatable :: dataString(:,:)
        type(dmatrix) :: fileData
        type(darray) :: row
        integer :: totalColumns, i, j

        this%mFilename = filename

        call read_file(filename, dataString)

        this%mNumIndependentVars = int(dataString(1,3))
        totalColumns = size(dataString(3,:))

        write(*,'(a,i0)') '|  Number of independent variables = ', this%mNumIndependentVars
        write(*,'(a,i0)') '|  Number of dependent variables = ', totalColumns - this%mNumIndependentVars

        write(*,'(a)') '|  Independent Variable Names = '
        do i = 1, this%mNumIndependentVars
            write(*,'(a,a)') '|  ', dataString(2,i)
            this%mIndependentVarNames(i) = dataString(2,i)
        end do
        write(*,*)

        write(*,'(a)') '|  Dependent Variable Names = '
        do i = this%mNumIndependentVars+1, totalColumns
            write(*,'(a,a)') '|  ', dataString(2,i)
            this%mDependentVarNames(i-this%mNumIndependentVars) = dataString(2,i)
        end do
        write(*,*)

        write(*,'(a)') '|  Checking and Converting Data...'
        do i = 3, size(dataString,1)
           
            if (size(dataString(i,:)) /= totalColumns) then
                write(*,'(a,i0,a)') 'Row ', i, ' has an error:'
                error stop 'Invalid number of columns in row'
            end if
            do j = 1, totalColumns
                call row%data%append(real(dataString(i,j)))
            end do
            call fileData%rows%append(row)
        end do

        write(*,'(a)') '|  Sorting Data by the Independent Variables'
        call sort_data(fileData, this%mNumIndependentVars)

        call this%create_from_data(fileData, this%mNumIndependentVars, this%mIndependentVarNames, this%mDependentVarNames)
        write(*,'(a)') '|  Dataset Successfully Loaded.'
        write(*,*)
    end subroutine create_from_file

    subroutine create_from_data(this, inputData, numIndependentVars, independentVarNames, dependentVarNames)
        class(dataset_type), intent(inout) :: this
        type(dmatrix), intent(in) :: inputData
        integer, intent(in) :: numIndependentVars
        character(len=256), intent(in) :: independentVarNames(:), dependentVarNames(:)
        integer :: totalColumns, i, j, k
        type(dmatrix) :: subset
        type(dataset_type), allocatable :: newDataset

        this%mNumIndependentVars = numIndependentVars
        this%mIndependentVarNames = independentVarNames
        this%mDependentVarNames = dependentVarNames

        totalColumns = size(inputData%rows(1)%data)

        if (this%mNumIndependentVars == 1) then
            this%mHasSubset = .false.
            do i = 1, size(inputData%rows)
                call this%mX%append(inputData%rows(i)%data(1))
                type(darray) :: row
                do j = 2, totalColumns
                    call row%data%append(inputData%rows(i)%data(j))
                end do
                call this%mY%rows%append(row)
            end do
        else
            this%mHasSubset = .true.
            do i = 1, size(inputData%rows)
                if (all(this%mX /= inputData%rows(i)%data(1))) then
                    call this%mX%append(inputData%rows(i)%data(1))
                end do
            end do

            do k = 1, size(this%mX)
                type(dmatrix) :: subset
                do i = 1, size(inputData%rows)
                    if (inputData%rows(i)%data(1) == this%mX(k)) then
                        type(darray) :: row
                        do j = 2, totalColumns
                            call row%data%append(inputData%rows(i)%data(j))
                        end do
                        call subset%rows%append(row)
                    end if
                end do
                allocate(this%mSubset(k))
                call this%mSubset(k)%create_from_data(subset, this%mNumIndependentVars-1, this%mIndependentVarNames(2:), this%mDependentVarNames)
            end do
        end if
    end subroutine create_from_data

    subroutine print_data(this, spacer)
        class(dataset_type), intent(in) :: this
        integer, intent(in) :: spacer
        integer :: i, j

        if (.not. this%mHasSubset) then
            call print_hspace(spacer)
            write(*,'(a,1x)', advance='no') this%mIndependentVarNames(1)
            call sarray_print(this%mDependentVarNames)
        end if

        do i = 1, size(this%mX)
            call print_hspace(spacer)
            write(*,'(a,a,f0.2)') this%mIndependentVarNames(1), ' = ', this%mX(i)
            if (this%mHasSubset) then
                write(*,*)
                call this%mSubset(i)%print_data(spacer+3)
            else
                do j = 1, size(this%mY%rows(i)%data)
                    write(*,'(f0.2,1x)', advance='no') this%mY%rows(i)%data(j)
                end do
                write(*,*)
            end if
        end do
    end subroutine print_data

    function linear_interpolate(this, independentValues, verbose) result(ans)
        class(dataset_type), intent(in) :: this
        type(darray), intent(in) :: independentValues
        integer, intent(in) :: verbose
        type(darray) :: ans
        integer :: lower, i, subVerbose
        type(darray) :: vals1, vals2
        type(darray) :: independentSubset

        double precision :: value, ratio

        value = independentValues%data(1)
        if (verbose > 0) then
            call print_hspace(verbose)
            write(*,'(a,f0.2)') 'Interpolating for ', this%mIndependentVarNames(1), ' = ', value
            subVerbose = verbose + 3
        end if

        if (this%mHasSubset) then
            do i = 2, size(independentValues%data)
                call independentSubset%data%append(independentValues%data(i))
            end do
        end if

        lower = -1
        do i = 1, size(this%mX)
            if (value > this%mX(i)) then
                lower = i
            end if
        end do

        if ((lower <= -1) .or. (lower >= size(this%mX)-1)) then
            if (lower == -1) then
                lower = 1
            end if
            if (this%mHasSubset) then
                vals1 = this%mSubset(lower)%linear_interpolate(independentSubset, subVerbose)
            else
                vals1 = this%mY%rows(lower)
            end if
            ans = vals1
        else
            if (this%mHasSubset) then
                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,'(a,i0,a,f0.2)') this%mIndependentVarNames(1), ' lower index is ', lower, ' with value of ', this%mIndependentVarNames(1), ' = ', this%mX(lower)
                end if
                vals1 = this%mSubset(lower)%linear_interpolate(independentSubset, subVerbose)

                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,'(a,i0,a,f0.2)') this%mIndependentVarNames(1), ' upper index is ', lower+1, ' with value of ', this%mIndependentVarNames(1), ' = ', this%mX(lower+1)
                end if
                vals2 = this%mSubset(lower+1)%linear_interpolate(independentSubset, subVerbose)
            else
                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,'(a,i0,a,f0.2)') this%mIndependentVarNames(1), ' lower index is ', lower, ' with value of ', this%mIndependentVarNames(1), ' = ', this%mX(lower)
                end if
                vals1 = this%mY%rows(lower)

                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,'(a,i0,a,f0.2)') this%mIndependentVarNames(1), ' upper index is ', lower+1, ' with value of ', this%mIndependentVarNames(1), ' = ', this%mX(lower+1)
                end if
                vals2 = this%mY%rows(lower+1)
            end if

            ratio = (value - this%mX(lower)) / (this%mX(lower+1) - this%mX(lower))
            do i = 1, size(vals1%data)
                call ans%data%append(vals1%data(i) + ratio * (vals2%data(i) - vals1%data(i)))
            end do

            if (verbose > 0) then
                call print_hspace(verbose)
                write(*,'(a)') 'Solution for '//this%mIndependentVarNames(1)//' = '//trim(adjustl(to_string(value)))
                call print_hspace(verbose)
                write(*,'(a)') '   Dependent Variables = '
                call sarray_print(this%mDependentVarNames)
                call print_hspace(verbose)
                write(*,'(a)') 'Lower dependent values = '
                call darray_print(vals1)
                call print_hspace(verbose)
                write(*,'(a)') 'Upper dependent values = '
                call darray_print(vals2)
                call print_hspace(verbose)
                write(*,'(a,f0.2)') '                 Ratio = ', ratio
                call print_hspace(verbose)
                write(*,'(a)') '                Answer = '
                call darray_print(ans)
                write(*,*)
            end if
        end if
    end function linear_interpolate

end module dataset_mod

