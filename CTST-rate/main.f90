program rate_calc
       use requirement_scratch
       use partition_clac_scratch
       implicit none
       integer :: R_1_atomicity, R_2_atomicity, R_1_shape, R_2_shape, R_1_sym_fac, R_2_sym_fac, pos
       integer :: ac_atomicity, ac_shape, ac_sym_fac
       real(8) :: Temp, R_1_mass_tot, R_2_mass_tot,  mass_tot_f, mon_iner, R_1_mom_iner, R_2_mom_iner, red_mass, R_1_vol_fac, R_2_vol_fac
       real(8) :: ac_mass_tot, ac_mom_iner, ac_vol_fac, R_1_trans_val, R_1_rot_val, R_1_vai_val, R_1_parti_val
       real(8) :: R_2_trans_val, R_2_rot_val, R_2_vai_val, ac_trans_val, ac_rot_val, ac_vai_val, R_2_parti_val, ac_parti_val, E_0 
       real(8) :: parti_fac
       character (len=256) :: filename, output, line1
  
      
     call system("ls *.inp > reshima")
     open(unit=20, file='reshima')
     read(20,*) filename
     close(unit=20)
     open(unit=10, file=filename)

     read(10,'(A)') line1
     read(10,*) Temp
     read(10,*)
     read(10,*) R_1_atomicity
     call calc_mass(R_1_atomicity, mass_tot)
     R_1_mass_tot = mass_tot
      
     read(10,*) R_1_shape
     call calc_moment_iner(R_1_atomicity, R_1_shape, mon_iner)
     R_1_mom_iner = mon_iner

     read(10,*) R_1_vol_fac
     read(10,*) R_1_sym_fac

     call trans_parti(R_1_mass_tot, R_1_vol_fac, Temp, trans_val)
     R_1_trans_val = trans_val
     call rot_parti(R_1_atomicity, R_1_mom_iner, R_1_sym_fac, Temp, rot_val)
     R_1_rot_val = rot_val
     call vai_parti(R_1_atomicity, Temp, vai_val)
     R_1_vai_val = vai_val

     call parti(R_1_trans_val, R_1_rot_val, R_1_vai_val, parti_val)
     R_1_parti_val = parti_val
!--------------------------------------------------------------------
     read(10,*)
     read(10,*) R_2_atomicity
     call calc_mass(R_2_atomicity, mass_tot)
     R_2_mass_tot = mass_tot
      
     read(10,*) R_2_shape
     call calc_moment_iner(R_2_atomicity, R_2_shape, mon_iner)
     R_2_mom_iner = mon_iner

     read(10,*) R_2_vol_fac
     read(10,*) R_2_sym_fac

     call trans_parti(R_2_mass_tot, R_2_vol_fac, Temp, trans_val)
     R_2_trans_val = trans_val
     call rot_parti(R_2_atomicity, R_2_mom_iner, R_2_sym_fac, Temp, rot_val)
     R_2_rot_val = rot_val
     call vai_parti(R_2_atomicity, Temp, vai_val)
     R_2_vai_val = vai_val

     call parti(R_2_trans_val, R_2_rot_val, R_2_vai_val, parti_val)
     R_2_parti_val = parti_val
!-------------------------------------------------------------------
     read(10,*)
     read(10,*) ac_atomicity
     call calc_mass(ac_atomicity, mass_tot)
     ac_mass_tot = mass_tot
      
     read(10,*) ac_shape
     call calc_moment_iner(ac_atomicity, ac_shape, mon_iner)
     ac_mom_iner = mon_iner

     read(10,*) ac_vol_fac
     read(10,*) ac_sym_fac

     
     call trans_parti(ac_mass_tot, ac_vol_fac, Temp, trans_val)
     ac_trans_val = trans_val
     call rot_parti(ac_atomicity, ac_mom_iner, R_2_sym_fac, Temp, rot_val)
     ac_rot_val = rot_val
     call vai_parti(ac_atomicity, Temp, vai_val)
     ac_vai_val = vai_val

     call parti(ac_trans_val, ac_rot_val, ac_vai_val, parti_val)
     ac_parti_val = parti_val
!---------------------------------------------------------------------
     read(10,*)
     read(10,*) E_0

     parti_fac = (ac_parti_val)/(R_1_parti_val * R_2_parti_val)
     call rate(parti_fac, E_0, Temp, rate_val)

     close(unit=10)
     
!------------------------------------------------------------------ 
     filename = adjustl(trim(filename))
      pos = index(filename, '.inp')
      if(pos > 0) then
         output = filename(1:pos-1)// '.out'
      else
         output = filename // '.out'
      end if
      open(unit=30, file= output)
      write(30,*) "=================================================================================="
      write(30,*) "                      =        ========           =          ========             "
      write(30,*) "                    =              =            =                =                "
      write(30,*) "                   =               =              =              =                "
      write(30,*) "                   =               =               =             =                "
      write(30,*) "                    =              =               =             =                "
      write(30,*) "                      =            =             =               =                "
      write(30,*) "=================================================================================="
      write(30,*) "                                     Aditya Barman                                "
      write(30,*) "                            Email: atomicadi2023@gmai.com                         "
      write(30,*) "                  Current Institute: Weizmann Institute of Science, Israel        "
      write(30,*) "           Alma Mater: Malaviya National Institute of Technology Jaipur, India    "
      write(30,*) "                                    August 03, 2025                               "
      write(30,*) "=================================================================================="
      
      write(30, *) "The reaction is:", line1
      write(30,*)  "The temperature is", Temp, "K" 
      write(30,*) "The barrier hight of the reaction is", E_0, "kj mol-1"
      write(30,*)

      write(30,*) "__________________________Information about Reactant-1_____________________________"
      if(R_1_atomicity == 0) then
        write(30,*) "Reactant-1 is monatomic."
      else if(R_1_atomicity == 1) then
        write(30,*) "Reactant-1 is diatomic."
      else if(R_1_atomicity == 2) then
        write(30,*) "Reactant-1 is polyatomic."
      end if
     write(30,*) "Total mass of Reactant-1 is", R_1_mass_tot, "kg"
     write(30,*) "Moment of inertia of Reactant-1 is", R_1_mom_iner, "kg m-2"
     write(30,*) "Volume factor of Reactant-1 is", R_1_vol_fac, "m3"
     write(30,*) "Symmetry factor of Reactant-1 is", R_1_sym_fac
     write(30,*) "Translational partition function of Reactant-1 is", R_1_trans_val, "m-3"
     write(30,*) "Rotational partition function of Reactant-1 is", R_1_rot_val
     if(R_1_atomicity == 0) then
       write(30,*) "Reactant-1 has no any vibrational frequency."
     else
       write(30,*) "Vibrational partition function of Reactant-1 is", R_1_vai_val
     end if
     write(30,*) "Total partition function of Reactant-1 is:", R_1_parti_val, "m-3"
     write(30,*)

     write(30,*) "__________________________Information about Reactant-2_____________________________"
      if(R_2_atomicity == 0) then
        write(30,*) "Reactant-2 is monatomic."
      else if(R_2_atomicity == 1) then
        write(30,*) "Reactant-2 is diatomic."
      else if(R_2_atomicity == 2) then
        write(30,*) "Reactant-2 is polyatomic."
      end if
     write(30,*) "Total mass of Reactant-2 is", R_2_mass_tot, "kg"
     write(30,*) "Moment of inertia of Reactant-2 is", R_2_mom_iner, "kg m-2"
     write(30,*) "Volume factor of Reactant-2 is", R_2_vol_fac, "m-3"
     write(30,*) "Symmetry factor of Reactant-2 is", R_2_sym_fac
     write(30,*) "Translational partition function of Reactant-2 is", R_2_trans_val, "m3"
     write(30,*) "Rotational partition function of Reactant-2 is", R_2_rot_val
     if(R_2_atomicity == 0) then
       write(30,*) "Reactant-2 has no any vibrational frequency."
     else
       write(30,*) "Vibrational partition function of Reactant-2 is", R_2_vai_val
     end if
     write(30,*) "Total partition function of Reactant-2 is:", R_2_parti_val, "m-3"
     write(30,*)

     write(30,*) "__________________________Information about Activated Complex_____________________________"
      if(ac_atomicity == 0) then
        write(30,*) "Activated Complex is monatomic."
      else if(ac_atomicity == 1) then
        write(30,*) "Activated Complex is diatomic."
      else if(ac_atomicity == 2) then
        write(30,*) "Activated Complex is polyatomic."
      end if
     write(30,*) "Total mass of Activated Complex is", ac_mass_tot, "kg"
     write(30,*) "Moment of inertia of Activated Complex is", ac_mom_iner, "kg m-2"
     write(30,*) "Volume factor of Activated Complex is", ac_vol_fac, "m-3"
     write(30,*) "Symmetry factor of Activated Complex is", ac_sym_fac
     write(30,*) "Translational partition function of Activated Complex is", ac_trans_val, "m3"
     write(30,*) "Rotational partition function of Activated Complex is", ac_rot_val
     if(ac_atomicity == 0) then
       write(30,*) "Activated Complex has no any vibrational frequency."
     else
       write(30,*) "Vibrational partition function of Activated Complex is", ac_vai_val
     end if
     write(30,*) "Total partition function of Activated Complex is:", ac_parti_val, "m-3"
     write(30,*)

     write(30,*) "__________________________Rate Constant_____________________________"
     write(30,*) "Rate constant of the reaction is", rate_val, "m3 s-1"
     write(30,*) "____________________________________________________________________"
     write(30,*) "                  | Hurray! A successful calculation |              " 
     write(30,*) "---------------------------End of the file--------------------------"
      
     close(unit=30)
     call system("rm reshima")
end program rate_calc
     