module requirement_scratch
       implicit none
             integer :: i, num_atom
       real(8), parameter :: kg_fac = 1.67e-27, ang_met_fac = 1e-10
       real(8) :: M_A, M_B, mass_atom(20), mass_tot, dis_A_B, atom_dis(20)
     contains

!=============calc_mass==============================
       subroutine calc_mass(atomicity, mass_tot)
         implicit none
           integer :: i, atomicity
           real(8) :: mass_tot
      
           if (atomicity == 0) then 
           read(10,*) 
           read(10,*) M_A
           mass_tot = M_A
           mass_tot = mass_tot * kg_fac
           else if (atomicity == 1) then
           read(10,*)
           read(10,*) M_A
           read(10,*) M_B
           mass_tot = M_A + M_B
           mass_tot = mass_tot * kg_fac
           else if (atomicity == 2) then
           read(10,*) num_atom
           do i = 1, num_atom
             read(10,*) mass_atom(i)
           end do
           mass_tot = 0.0
           do i = 1, num_atom
             mass_tot = mass_tot + mass_atom(i)
           end do
           mass_tot = mass_tot * kg_fac
          else 
            print*, "Something wrong in your input"
          end if 
       end subroutine calc_mass

!===============calc_moment_iner=================================
       subroutine calc_moment_iner(atomicity, shape, mon_iner)
         implicit none
            integer :: atomicity, shape   
            real(8) :: red_mass, dis_A_B, mon_iner
          
         if(atomicity==0 .and. shape == 0) then
         mon_iner = 1.0
         else if(atomicity == 1 .and.  shape == 1) then
         red_mass = ((M_A * M_B)/(M_A + M_B)) * kg_fac
         read(10,*) dis_A_B
         mon_iner =   red_mass * (dis_A_B * ang_met_fac)**2
         else if(atomicity == 2 .and. shape == 1) then
         do i = 1, num_atom
           read(10,*) atom_dis(i)
         end do
         mon_iner = ((mass_atom(1)*mass_atom(2)*((kg_fac)**2))*((atom_dis(1) * ang_met_fac)**2) + (mass_atom(3)*mass_atom(1)*((kg_fac)**2))*((atom_dis(3) * ang_met_fac)**2) + (mass_atom(2)*mass_atom(3)*((kg_fac)**2))*((atom_dis(2) * ang_met_fac)**2))/mass_tot
         else
            print*, "Something wrong in your input"
         end if
      end subroutine calc_moment_iner
end module 
         

