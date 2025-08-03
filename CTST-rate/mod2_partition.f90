module partition_clac_scratch
       implicit none
       real(8), parameter :: pi = 3.141592653589793, k_B = 1.380649e-23, h = 6.626e-34, c = 2.99e10, R = 8.314, kilo_fac = 1000.0
       real(8) :: trans_val, rot_val, vai_val, parti_val, rate_val
       contains

!===================Translational Partition function==========================
       subroutine trans_parti(mass_tot, vol_fac, Temp, trans_val)
          implicit none 
            real(8) :: mass_tot, vol_fac, Temp, trans_val
 
            trans_val = ((((2 * pi * mass_tot * k_B * Temp)**(3.0/2.0))/(h**3))/(vol_fac))
       end subroutine
!====================Rotational Partition function============================
       subroutine rot_parti(atomicity, mom_iner, sym_fac, Temp, rot_val)
          implicit none
           integer :: atomicity, sym_fac
            real(8) :: mom_iner, Temp, rot_val
 
            if(atomicity == 0) then
              rot_val = 1.0
            else
               rot_val = ((8 * (pi**2) * mom_iner * k_B * Temp))/(sym_fac*(h**2))
            end if
       end subroutine
!====================Vibrational Partition function============================
       subroutine vai_parti(atomicity, Temp, vai_val)
          implicit none
           integer :: i, atomicity, freq_num
            real(8) ::   fac_temp, freq(10), fac(10), Temp, vai_val
 
            if(atomicity == 0) then
              vai_val = 1.0
            else if(atomicity > 0) then
                read(10,*) freq_num
                if(freq_num == 1) then
                read(10,*) freq(1)
                fac(1) = (h * (freq(1) * c))/(k_B * Temp)
                vai_val = 1 / (1-(exp(-fac(1))))
                else if(freq_num == 2) then
                do i = 1, 2
                  read(10,*) freq(i)
                end do
                do i = 1, 2
                   fac(i) = (h * (freq(i) * c))/(k_B * Temp)
                end do
                vai_val = 1 / ((1-(exp(-fac(1)))) * ((1-(exp(-fac(2))))**2))
                end if
            end if
         end subroutine
!=========================Total Partition function=============================
         subroutine parti(trans_val, rot_val, vai_val, parti_val)
            implicit none
              real(8) :: trans_val, rot_val, vai_val, parti_val
              parti_val = (trans_val * rot_val * vai_val)
         end subroutine
!========================Rate Constant========================================
         subroutine rate(parti_fac, E_0, Temp, rate_val)
            implicit none
              real(8) :: parti_fac, E_0, Temp, rate_val
              rate_val = ((k_B * Temp)/h) * parti_fac * (exp(-((E_0*kilo_fac)/(R * Temp))))
         end subroutine

end module