module iomanip

contains

    ! convert number to string of length 10 with
    ! leading zeros
    function zfill(num) result(name)
        integer, intent(in) :: num
        character(len=10)   :: name

        write(name, fmt='(I10.10)') num
    end function zfill

end module iomanip
