!Aditya Barman, graduate student, Weizmann Institute of Science, July 11, 2025
! Email: atomicadi2023@gmail.com
program wigner_theory_clac
       implicit none
        real(kind=8), parameter :: h = 6.626D-34, k_B = 1.3806D-23, c= 2.99D10
        real(kind=8) :: imgf, T, nume, denom, val, k_wig 
      
     print *, "imagineary frequency of TS in cm-1 (eg. -672.34)"
     read(*,*) imgf
     print *, "Temperature i K (eg. 273.60)"
     read(*,*) T
     nume = h * (imgf*(-1)) * c
     denom = k_B * T
     val = (nume/denom) **2
     k_wig = 1.0D0 + (1.0D0/24.0D0)*val
     write(*,*) "The Wigner tunneling factor (k_wig) is", k_wig
end program wigner_theory_clac
