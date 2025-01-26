module iomanip

contains

    ! convert number to string of length 10 with
    ! leading zeros
    function zfill(num) result(name)
        integer, intent(in) :: num
        character(len=10)   :: name

        write(name, fmt='(I10.10)') num
    end function zfill

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function int2string(num) result(name)
        integer, intent(in) :: num
        character(len=10)   :: name

        if (num < 10) then
            write(name, fmt='(I1.1)') num
        else if (num < 100) then
            write(name, fmt='(I2.2)') num
        else if (num < 1000) then
            write(name, fmt='(I3.3)') num
        else
            name = zfill(num)
        endif
    end function int2string

end module iomanip
