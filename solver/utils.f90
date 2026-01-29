!------------------------------------------------------------------------------
! Small string and filesystem utilities.
!------------------------------------------------------------------------------
module utils
  use iso_fortran_env, only: error_unit
  implicit none
contains
  subroutine lower_string(s)
    character(len=*), intent(inout) :: s
    integer :: i, c

    do i = 1, len(s)
      c = ichar(s(i:i))
      if (c >= 65 .and. c <= 90) then
        s(i:i) = char(c + 32)
      end if
    end do
  end subroutine lower_string

  subroutine read_next_line(unit, line, ios)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    integer, intent(out) :: ios
    integer :: p1, p2, cut

    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) return

      p1 = index(line,'!')
      p2 = index(line,'#')
      cut = 0
      if (p1 > 0) cut = p1
      if (p2 > 0 .and. (cut == 0 .or. p2 < cut)) cut = p2
      if (cut > 0) line = line(1:cut-1)

      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      return
    end do
  end subroutine read_next_line

  subroutine ensure_output_dir()
    call execute_command_line('mkdir -p output')
  end subroutine ensure_output_dir
end module utils
