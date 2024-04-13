! FILEPATH: /c:/Users/Josh/Documents/AeroLab_Projects/MultiInterp/test.f90

module DatasetModule
    implicit none
    private

    type :: Dataset
        private
        character(len=:), allocatable :: mFilename  ! Filename of the dataset
        integer :: mNumIndependentVars  ! Number of independent variables
        logical :: mHasSubset  ! Flag indicating if the dataset has subsets
        character(len=:), allocatable :: mIndependentVarNames(:)  ! Names of independent variables
        character(len=:), allocatable :: mDependentVarNames(:)  ! Names of dependent variables
        real(kind=8), allocatable :: mX(:)  ! Array storing the values of the first independent variable
        real(kind=8), allocatable :: mY(:,:)  ! Array storing the values of the dependent variables
        type(Dataset), allocatable :: mSubset(:)  ! Array of subsets of the dataset
    contains
        procedure :: create_from_data  ! Subroutine to create a dataset from input data
        procedure :: print_data  ! Subroutine to print the dataset
        procedure :: linear_interpolate  ! Function to perform linear interpolation
    end type Dataset

    interface Dataset
        module procedure :: init_Dataset  ! Interface for initializing a dataset
    end interface Dataset

    contains
            ! this subroutine will create the dataset from the data you need for input into the interpolator and should be changed to include the data you need
    subroutine init_Dataset(this)
        class(Dataset), intent(inout) :: this
        this%mFilename = ""  ! Initialize the filename to an empty string
        this%mNumIndependentVars = 0  ! Initialize the number of independent variables to 0
        this%mHasSubset = .false.  ! Initialize the subset flag to false
        this%mIndependentVarNames = []  ! Initialize the array of independent variable names to an empty array
        this%mDependentVarNames = []  ! Initialize the array of dependent variable names to an empty array
        this%mX = []  ! Initialize the array of values for the first independent variable to an empty array
        this%mY = []  ! Initialize the array of values for the dependent variables to an empty array
        this%mSubset = []  ! Initialize the array of subsets to an empty array
    end subroutine init_Dataset

    subroutine create_from_data(this, inputData, numIndependentVars, independentVarNames, dependentVarNames)
        class(Dataset), intent(inout) :: this
        real(kind=8), intent(in) :: inputData(:,:)  ! Input data array
        integer, intent(in) :: numIndependentVars  ! Number of independent variables
        character(len=*), intent(in) :: independentVarNames(:)  ! Names of independent variables
        character(len=*), intent(in) :: dependentVarNames(:)  ! Names of dependent variables

        integer :: i, j, k, totalColumns
        integer :: numDataPoints
        real(kind=8), allocatable :: row(:)
        real(kind=8), allocatable :: subset(:,:)
        character(len=:), allocatable :: tempNames(:)
        type(Dataset), allocatable :: newDataset

        this%mNumIndependentVars = numIndependentVars  ! Set the number of independent variables
        this%mIndependentVarNames = independentVarNames  ! Set the names of independent variables
        this%mDependentVarNames = dependentVarNames  ! Set the names of dependent variables

        totalColumns = size(inputData, 2)

        ! Create subsets if needed
        if (this%mNumIndependentVars == 1) then
            this%mHasSubset = .false.  ! No subsets for single independent variable
            numDataPoints = size(inputData, 1)
            allocate(this%mX(numDataPoints))  ! Allocate memory for the array of values for the first independent variable
            allocate(this%mY(numDataPoints, totalColumns - 1))  ! Allocate memory for the array of values for the dependent variables
            do i = 1, numDataPoints
                this%mX(i) = inputData(i, 1)  ! Store the values of the first independent variable
                this%mY(i, :) = inputData(i, 2:totalColumns)  ! Store the values of the dependent variables
            end do
        else
            this%mHasSubset = .true.  ! Dataset has subsets for multiple independent variables
            allocate(this%mX(0))  ! Allocate memory for the array of values for the first independent variable
            do i = 1, size(inputData, 1)
                if (.not. any(this%mX == inputData(i, 1))) then
                    this%mX = [this%mX, inputData(i, 1)]  ! Append unique values of the first independent variable
                end if
            end do

            do k = 1, size(this%mX)
                subset = reshape([inputData(i, 2:totalColumns) | i = 1, size(inputData, 1), inputData(:, 1) == this%mX(k)], [count(inputData(:, 1) == this%mX(k)), totalColumns - 1])  ! Create subset of data for each unique value of the first independent variable
                tempNames = this%mIndependentVarNames(2:)  ! Names of remaining independent variables
                allocate(newDataset)  ! Allocate memory for a new dataset
                call init_Dataset(newDataset)  ! Initialize the new dataset
                call create_from_data(newDataset, subset, this%mNumIndependentVars - 1, tempNames, dependentVarNames)  ! Recursively create subsets for remaining independent variables
                this%mSubset = [this%mSubset, newDataset]  ! Append the new dataset to the array of subsets
            end do
        end if
    end subroutine create_from_data

    subroutine print_data(this, spacer)
        class(Dataset), intent(in) :: this
        integer, intent(in) :: spacer

        integer :: i, j

        if (.not. this%mHasSubset) then
            call print_hspace(spacer)
            write(*,*) trim(this%mIndependentVarNames(1)), " ", this%mDependentVarNames  ! Print the names of the independent and dependent variables
        end if

        do i = 1, size(this%mX)
            call print_hspace(spacer)
            write(*,*) trim(this%mIndependentVarNames(1)), " = ", this%mX(i), " : "  ! Print the value of the first independent variable
            if (this%mHasSubset) then
                call print_data(this%mSubset(i), spacer + 3)  ! Recursively print the subsets
            else
                write(*,*) this%mY(i, :)  ! Print the values of the dependent variables
            end if
        end do
    end subroutine print_data

    recursive function linear_interpolate(this, independentValues, verbose) result(ans)
        class(Dataset), intent(in) :: this
        real(kind=8), intent(in) :: independentValues(:)  ! Array of independent variable values for interpolation
        integer, intent(in) :: verbose  ! Verbosity level for printing intermediate results
        integer :: lower, subVerbose
        real(kind=8), allocatable :: vals1(:), vals2(:), ans(:)
        real(kind=8), allocatable :: independentSubset(:)
        real(kind=8) :: value, ratio
        character(len=:), allocatable :: iVarName

        value = independentValues(1)  ! Value of the first independent variable for interpolation
        iVarName = this%mIndependentVarNames(1)  ! Name of the first independent variable
        lower = -1  ! Index of the lower bound for interpolation
        subVerbose = 0  ! Verbosity level for printing intermediate results within the recursive function
        vals1 = []  ! Array for storing the values of the dependent variables for the lower bound
        vals2 = []  ! Array for storing the values of the dependent variables for the upper bound
        ans = []  ! Array for storing the interpolated values
        allocate(independentSubset(0))  ! Allocate memory for the subset of independent variable values

        if (verbose > 0) then
            call print_hspace(verbose)
            write(*,*) "Interpolating for ", trim(iVarName), " = ", value  ! Print the current interpolation value
            subVerbose = verbose + 3
        end if

        ! Create subset of values to be interpolated if needed
        if (this%mHasSubset) then
            independentSubset = independentValues(2:)  ! Subset of independent variable values excluding the first variable
        end if

        ! Find index of lower bound
        do i = 1, size(this%mX)
            if (value > this%mX(i)) then
                lower = i
            end if
        end do

        if ((lower <= -1) .or. (lower >= size(this%mX) - 1)) then ! Case when value is out of bounds of data
            if (lower == -1) then ! set lower to 0 to represent lower bound value
                lower = 0
            end if
            if (this%mHasSubset) then
                vals1 = linear_interpolate(this%mSubset(lower), independentSubset, subVerbose)  ! Recursively interpolate for the lower bound subset
            else
                vals1 = this%mY(lower, :)  ! Use the values of the dependent variables for the lower bound
            end if
            ans = vals1  ! Set the interpolated values to the values for the lower bound
        else
            if (this%mHasSubset) then
                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,*) trim(iVarName), " lower index is ", lower, " with value of ", trim(iVarName), " = ", this%mX(lower)  ! Print the lower bound index and value
                end if
                vals1 = linear_interpolate(this%mSubset(lower), independentSubset, subVerbose)  ! Recursively interpolate for the lower bound subset

                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,*) trim(iVarName), " upper index is ", lower + 1, " with value of ", trim(iVarName), " = ", this%mX(lower + 1)  ! Print the upper bound index and value
                end if
                vals2 = linear_interpolate(this%mSubset(lower + 1), independentSubset, subVerbose)  ! Recursively interpolate for the upper bound subset
            else
                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,*) trim(iVarName), " lower index is ", lower, " with value of ", trim(iVarName), " = ", this%mX(lower)  ! Print the lower bound index and value
                end if
                vals1 = this%mY(lower, :)  ! Use the values of the dependent variables for the lower bound

                if (verbose > 0) then
                    call print_hspace(verbose)
                    write(*,*) trim(iVarName), " upper index is ", lower + 1, " with value of ", trim(iVarName), " = ", this%mX(lower + 1)  ! Print the upper bound index and value
                end if
                vals2 = this%mY(lower + 1, :)  ! Use the values of the dependent variables for the upper bound
            end if

            ! Perform the Interpolation
            ratio = (value - this%mX(lower)) / (this%mX(lower + 1) - this%mX(lower))  ! Calculate the interpolation ratio
            ans = vals1 + ratio * (vals2 - vals1)  ! Perform linear interpolation
            if (verbose > 0) then
                call print_hspace(verbose)
                write(*,*) "Solution for ", trim(iVarName), " = ", value  ! Print the interpolated value
                call print_hspace(verbose)
                write(*,*) "   Dependent Variables = ", this%mDependentVarNames  ! Print the names of the dependent variables
                call print_hspace(verbose)
                write(*,*) "Lower dependent values = ", vals1  ! Print the values of the dependent variables for the lower bound
                call print_hspace(verbose)
                write(*,*) "Upper dependent values = ", vals2  ! Print the values of the dependent variables for the upper bound
                call print_hspace(verbose)
                write(*,*) "                 Ratio = ", ratio  ! Print the interpolation ratio
                call print_hspace(verbose)
                write(*,*) "                Answer = ", ans  ! Print the interpolated values
                call print_hspace(verbose)
                write(*,*)
            end if
        end if

    end function linear_interpolate

end module DatasetModule

program Main
    use DatasetModule
    implicit none

    type(Dataset) :: myDataset

    ! Usage example
    call myDataset%create_from_file("/c:/Users/Josh/Documents/AeroLab_Projects/MultiInterp/test.txt")  ! Create dataset from file
    call myDataset%print_data(0)  ! Print the dataset
    write(*,*)
    write(*,*) "Interpolated values:"
    write(*,*) myDataset%linear_interpolate([1.0, 2.0, 3.0], 0)  ! Perform linear interpolation

end program Main
