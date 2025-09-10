! Aditya Barman, Graduate student, Weizmann Institute of Science, Israel; Aug 15, 2025
! Email: atomicadi2023@gmail.com

program rate_calc
       use mod_1_requirement_scratch
       use mod_2_partition_calc_scratch
        
       implicit none
       integer :: R_1_atomicity, R_2_atomicity, ac_atomicity, R_1_sym_fac, R_2_sym_fac, ac_sym_fac, pos
       integer :: R_1_freq_num,  R_2_freq_num, ac_freq_num 
       real(8) :: Temp, R_1_vol_fac, R_2_vol_fac, ac_vol_fac, R_1_freq(10), R_2_freq(10), ac_freq(10) 
       real(8) :: R_1_mass_atom(100), R_2_mass_atom(100), ac_mass_atom(100)
       real(8) :: R_1_mass_tot, R_2_mass_tot, R_1_mass_tot_kg, R_2_mass_tot_kg, R_1_atom_mat(100,3), R_2_atom_mat(100,3)
       real(8) :: R_1_COM(1,3), R_2_COM(1,3), R_1_I_A, R_1_I_B, R_1_I_C, R_2_I_A, R_2_I_B, R_2_I_C
       real(8) :: ac_mass_tot, ac_mass_tot_kg, ac_COM(1,3), ac_I_A, ac_I_B, ac_I_C
       real(8) :: R_1_I_A_unit, R_1_I_B_unit, R_1_I_C_unit, R_2_I_A_unit, R_2_I_B_unit, R_2_I_C_unit, ac_I_A_unit, ac_I_B_unit, ac_I_C_unit
       real(8) :: R_1_trans_val, R_2_trans_val, ac_trans_val, R_1_rot_val, R_2_rot_val, ac_rot_val, R_1_vib_val, R_2_vib_val, ac_vib_val
       real(8) :: R_1_parti_val, R_2_parti_val, ac_parti_val, E_0, parti_fac
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
     call calc_mass(R_1_atomicity, mass_tot, R_1_mass_atom, mass_tot_kg)
     R_1_mass_tot = mass_tot
     R_1_mass_tot_kg = mass_tot_kg
     call calc_iner(R_1_atomicity, R_1_mass_tot, R_1_mass_atom, R_1_COM, R_1_I_A, R_1_I_B, R_1_I_C)
     read(10,*) R_1_vol_fac
     read(10,*) R_1_sym_fac
 
     call trans_parti(R_1_mass_tot_kg, R_1_vol_fac, Temp, R_1_trans_val)
     call rot_parti(R_1_atomicity, R_1_sym_fac, Temp, R_1_I_A, R_1_I_B, R_1_I_C, R_1_I_A_unit, R_1_I_B_unit, R_1_I_C_unit, R_1_rot_val)
     call vai_parti(R_1_atomicity, Temp, R_1_freq_num, R_1_freq, R_1_vib_val)
     call parti(R_1_trans_val, R_1_rot_val, R_1_vib_val, R_1_parti_val)
     
!--------------------------------------------------------------------
     read(10,*)
     read(10,*) R_2_atomicity
     call calc_mass(R_2_atomicity, mass_tot, R_2_mass_atom, mass_tot_kg)
     R_2_mass_tot = mass_tot
     R_2_mass_tot_kg = mass_tot_kg
     call calc_iner(R_2_atomicity, R_2_mass_tot, R_2_mass_atom, R_2_COM, R_2_I_A, R_2_I_B, R_2_I_C)
     read(10,*) R_2_vol_fac
     read(10,*) R_2_sym_fac

     call trans_parti(R_2_mass_tot_kg, R_2_vol_fac, Temp, R_2_trans_val)
     call rot_parti(R_2_atomicity, R_2_sym_fac, Temp, R_2_I_A, R_2_I_B, R_2_I_C, R_2_I_A_unit, R_2_I_B_unit, R_2_I_C_unit, R_2_rot_val)
     call vai_parti(R_2_atomicity, Temp, R_2_freq_num, R_2_freq, R_2_vib_val)
     call parti(R_2_trans_val, R_2_rot_val, R_2_vib_val, R_2_parti_val)

!---------------------------------------------------------------------
     read(10,*)
     read(10,*) ac_atomicity
     call calc_mass(ac_atomicity, mass_tot, ac_mass_atom, mass_tot_kg)
     ac_mass_tot = mass_tot
     ac_mass_tot_kg = mass_tot_kg
     call calc_iner(ac_atomicity, ac_mass_tot, ac_mass_atom, ac_COM, ac_I_A, ac_I_B, ac_I_C)
     read(10,*) ac_vol_fac
     read(10,*) ac_sym_fac

     call trans_parti(ac_mass_tot_kg, ac_vol_fac, Temp, ac_trans_val)
     call rot_parti(ac_atomicity, ac_sym_fac, Temp, ac_I_A, ac_I_B, ac_I_C, ac_I_A_unit, ac_I_B_unit, ac_I_C_unit, ac_rot_val)
     call vai_parti(ac_atomicity, Temp, ac_freq_num, ac_freq, ac_vib_val)
     call parti(ac_trans_val, ac_rot_val, ac_vib_val, ac_parti_val)

!-----------------------------------------------------------------------
     read(10,*) 
     read(10,*) E_0

!-----------------------------------------------------------------------
     parti_fac = (ac_parti_val/(R_1_parti_val * R_2_parti_val))
     call rate(parti_fac, E_0, Temp, rate_val)

!==================================================================Print Scetion===============================================
    filename = adjustl(trim(filename))
      pos = index(filename, '.inp')
      if(pos > 0) then
         output = filename(1:pos-1)// '.out'
      else
         output = filename // '.out'
      end if

     open(unit=30, file= output)
     write(30,*) "================================================================================================================"
     write(30,*) "                  =====      =      ===    =====           |=                         "
     write(30,*) "                  =        =   =    =   =    =            / \                         "
     write(30,*) "                  =       =     =   =    =   =            |-|                         " 
     write(30,*) "                  =       =     =   =    =   =            |-|                         "      
     write(30,*) "                  =       =     =   =   =    =            | |  || || ||               "            
     write(30,*) "                  =====   =     =   =  =     =            | |--||-||-||--             "
     write(30,*) "                  =       =     =   =   =    =           -----------------            "
     write(30,*) "                  =       =     =   =   =    =            |=============|             "
     write(30,*) "                  =       =     =   =   =    =            |====| |==| |=|             "  
     write(30,*) "                  =       =     =   =   =    =            |=============|             "
     write(30,*) "                  =       =     =   =   =    =            |====|     |==|             "
     write(30,*) "                  =        =   =    =   =    =            |====|     |==|             "
     write(30,*) "                  =          =      =   =    =           =================            " 
     write(30,*) "================================================================================================================"
     write(30,*) "               FORT: Fortran Operator for Rate Theory, A CTST rate calculator                                   "
     write(30,*) "                                      Aditya Barman                                                             "
     write(30,*) "                              Email: atomicadi2023@gmail.com                                                    "
     write(30,*) "                  Graduate student @ Weizmann Institute of Science, Israel                                      "
     write(30,*) "             Alma mater: Malaviya National Institute of Technology Jaipur, India                                "
     write(30,*) "================================================================================================================"
     write(30,*) "The reaction is ", line1
     write(30,*) " "
     write(30,*) "________________________________Information about Reactant-1____________________________________________________"
     if(R_1_atomicity == 0) then
     write(30,*) "Reactant-1 is Monatomic."
     else if(R_1_atomicity == 1) then
     write(30,*) "Reactant-1 is Diatomic."
     else
         write(30,*) "Reactant-1 is Polyatomic."
     end if

     write(30,*) "Total mass of Reactant-1 is", R_1_mass_tot_kg, "kg."
     
     if(R_1_atomicity == 0) then
     write(30,*) "Again, Reactant-1 is monatomic."
     else
         write(30,*) "Moments of Inertia (kg m2) of Reactant-1 are: I_A=",R_1_I_A_unit, ", I_B=", R_1_I_B_unit, ", I_C=", R_1_I_C_unit
     end if

     write(30,*) "Volume factor of Reactant-1 is", R_1_vol_fac, "m3."
     write(30,*) "Symmetry factor of Reactant-1 is", R_1_sym_fac

     write(30,*) "Translational partition function of Reactant-1 is", R_1_trans_val, "m-3."

     if(R_1_atomicity == 0 .and. R_1_I_A == 0.0 .and. R_1_I_B == 0.0 .and. R_1_I_C == 0.0) then
     write(30,*) "Your Reactant-1 has only a single atom !!!"
     else if(R_1_atomicity > 0 .and. (abs(R_1_I_A)) < tolarance .and. (abs(R_1_I_B - R_1_I_C)) < tolarance) then
     write(30,*) "Reactant-1 is Linear."
     write(30,*) "Rotation partition function of Reactant-1 is", R_1_rot_val
     else if(R_1_atomicity > 0 .and. (abs(R_1_I_A - R_1_I_B)) < tolarance .and. (abs(R_1_I_B - R_1_I_C)) < tolarance) then
     write(30,*) "Reactant-1 belongs to Spherical top."
     write(30,*) "Rotation partition function of Reactant-1 is", R_1_rot_val
     else if(R_1_atomicity > 0 .and. R_1_I_A < R_1_I_B .and. (abs(R_1_I_B - R_1_I_C)) < tolarance) then
     write(30,*) "Reactant-1 belongs to Prolet Symmetric top."
     write(30,*) "Rotation partition function of Reactant-1 is", R_1_rot_val
     else if(R_1_atomicity > 0 .and. R_1_I_A < R_1_I_C .and. (abs(R_1_I_A - R_1_I_B)) < tolarance) then
     write(30,*) "Reactant-1 belongs to Oblet Symmetric top."
     write(30,*) "Rotation partition function of Reactant-1 is", R_1_rot_val
     else if(R_1_atomicity > 0 .and. R_1_I_A < R_1_I_B .and. R_1_I_B < R_1_I_C .and.  (R_1_I_B - R_1_I_A) > 0.4 .and. (R_1_I_C - R_1_I_A) > 0.4) then
     write(30,*) "Reactant-1 belongs to Asymmetric top."
     write(30,*) "Rotation partition function of Reactant-1 is", R_1_rot_val
     end if

     if(R_1_atomicity == 0) then
     write(30,*) "Again, your Reactant-1 has only a single atom !!!"
     else 
         if(R_1_freq_num == 1) then
         write(30,*) "Reactant-1 has", R_1_freq_num, "vibrational frequency."
         write(30,*) "The vibrational frequency is", R_1_freq(1), "cm-1."
         write(30,*) "Vibrational partition function of Reactant-1 is", R_1_vib_val
         else if(R_1_freq_num > 1) then
         write(30,*) "Reactant-1 has", R_1_freq_num, "vibrational frequencies."
         do i = 1, R_1_freq_num
           write(30,*) i, "frequency is", R_1_freq(i), "cm-1."
         end do
         write(30,*) "Vibrational partition function of Reactant-1 is", R_1_vib_val
         end if
     end if
   
     write(30,*) "Total partition function of Reactant-1 is",  R_1_parti_val, "m-3."
     write(30,*) "    "

     write(30,*) "________________________________Information about Reactant-2____________________________________________________"
     if(R_2_atomicity == 0) then
     write(30,*) "Reactant-2 is Monatomic."
     write(30,*) "Reactant-2 is Monatomic."
     else if(R_2_atomicity == 1) then
     write(30,*) "Reactant-2 is Diatomic."
     else
         write(30,*) "Reactant-2 is Polyatomic."
     end if

     write(30,*) "Total mass of Reactant-2 is", R_2_mass_tot_kg, "kg."
      
     if(R_2_atomicity == 0) then
     write(30,*) "Again, Reactant-2 is monatomic."
     else
         write(30,*) "Moments of Inertia (kg m2) of Reactant-2 are: I_A=",R_2_I_A_unit, ", I_B=", R_2_I_B_unit, ", I_C=", R_2_I_C_unit
     end if

     write(30,*) "Volume factor of Reactant-2 is", R_2_vol_fac, "m3."
     write(30,*) "Symmetry factor of Reactant-2 is", R_2_sym_fac

     write(30,*) "Translational partition function of Reactant-2 is", R_2_trans_val, "m-3."

     if(R_2_atomicity == 0 .and. R_2_I_A == 0.0 .and. R_2_I_B == 0.0 .and. R_2_I_C == 0.0) then
     write(30,*) "Your Reactant-2 has only a single atom !!!"
     else if(R_2_atomicity > 0 .and. (abs(R_2_I_A)) < tolarance .and. (abs(R_2_I_B - R_2_I_C)) < tolarance) then
     write(30,*) "Reactant-2 is Linear."
     write(30,*) "Rotation partition function of Reactant-2 is", R_2_rot_val
     else if(R_2_atomicity > 0 .and. (abs(R_2_I_A - R_2_I_B)) < tolarance .and. (abs(R_2_I_B - R_2_I_C)) < tolarance) then
     write(30,*) "Reactant-2 belongs to Spherical top."
     write(30,*) "Rotation partition function of Reactant-2 is", R_2_rot_val
     else if(R_2_atomicity > 0 .and. R_2_I_A < R_2_I_B .and. (abs(R_2_I_B - R_2_I_C)) < tolarance) then
     write(30,*) "Reactant-2 belongs to Prolet Symmetric top."
     write(30,*) "Reactant-2 partition function of Reactant-2 is", R_2_rot_val
     else if(R_2_atomicity > 0 .and. R_2_I_A < R_2_I_C .and. (abs(R_2_I_A - R_2_I_B)) < tolarance) then
     write(30,*) "Reactant-2 belongs to Oblet Symmetric top."
     write(30,*) "Rotation partition function of Reactant-2 is", R_2_rot_val
     else if(R_2_atomicity > 0 .and. R_2_I_A < R_2_I_B .and. R_2_I_B < R_2_I_C .and.  (R_2_I_B - R_2_I_A) > 0.4 .and. (R_2_I_C - R_2_I_A) > 0.4) then
     write(30,*) "Reactant-2 belongs to Asymmetric top."
     write(30,*) "Rotation partition function of Reactant-2 is", R_2_rot_val
     end if

     if(R_2_atomicity == 0) then
     write(30,*) "Again, your Reactant-2 has only a single atom !!!"
     else 
         if(R_2_freq_num == 1) then
         write(30,*) "Reactant-2 has", R_2_freq_num, "vibrational frequency."
         write(30,*) "The vibrational frequency is", R_2_freq(1), "cm-1."
         write(30,*) "Vibrational partition function of Reactant-2 is", R_2_vib_val
         else if(R_2_freq_num > 1) then
         write(30,*) "Reactant-2 has", R_2_freq_num, "vibrational frequencies."
         do i = 1, R_2_freq_num
           write(30,*) i, "frequency is", R_2_freq(i), "cm-1."
         end do
         write(30,*) "Vibrational partition function of Reactant-2 is", R_2_vib_val
         end if
     end if

     write(30,*) "Total partition function of Reactant-2 is",  R_2_parti_val, "m-3."
     write(30,*) "   "

     write(30,*) "_____________________________Information about Activated Complex________________________________________________"
     if(ac_atomicity == 0) then
     write(30,*) "Activated complex is Monatomic."
     write(30,*) "Activated complex is Monatomic."
     else if(ac_atomicity == 1) then
     write(30,*) "Activated complex is Diatomic."
     else
         write(30,*) "Activated complex is Polyatomic."
     end if

     write(30,*) "Total mass of activated complex is", ac_mass_tot_kg, "kg."

     if(ac_atomicity == 0) then
     write(30,*) "Again, activated complex is monatomic."
     else
         write(30,*) "Moments of Inertia (kg m2) of activated complex are: I_A=", ac_I_A_unit, ", I_B=", ac_I_B_unit, ", I_C=", ac_I_C_unit
     end if 
     
     write(30,*) "Volume factor of activated complex is", ac_vol_fac, "m3."
     write(30,*) "Symmetry factor of activated complex is", ac_sym_fac

     write(30,*) "Translational partition function of activated complex is", ac_trans_val, "m-3."

     if(ac_atomicity == 0 .and. ac_I_A == 0.0 .and. ac_I_B == 0.0 .and. ac_I_C == 0.0) then
     write(30,*) "Your activated complex has only a single atom !!!"
     else if(ac_atomicity > 0 .and. (abs(ac_I_A)) < tolarance .and. (abs(ac_I_B - ac_I_C)) < tolarance) then
     write(30,*) "Activated complex is Linear."
     write(30,*) "Rotation partition function of activated complex is", ac_rot_val
     else if(ac_atomicity > 0 .and. (abs(ac_I_A - ac_I_B)) < tolarance .and. (abs(ac_I_B - ac_I_C)) < tolarance) then
     write(30,*) "Activated complex belongs to Spherical top."
     write(30,*) "Rotation partition function of activated complex is", ac_rot_val
     else if(ac_atomicity > 0 .and. ac_I_A < ac_I_B .and. (abs(ac_I_B - ac_I_C)) < tolarance) then
     write(30,*) "Activated complex belongs to Prolet Symmetric top."
     write(30,*) "Rotation partition function of activated complex is", ac_rot_val
     else if(ac_atomicity > 0 .and. ac_I_A < ac_I_C .and. (abs(ac_I_A - ac_I_B)) < tolarance) then
     write(30,*) "Activated complex belongs to Oblet Symmetric top."
     write(30,*) "Rotation partition function of activated complex is", ac_rot_val
     else if(ac_atomicity > 0 .and. ac_I_A < ac_I_B .and. ac_I_B < ac_I_C .and.  (ac_I_B - ac_I_A) > 0.4 .and. (ac_I_C - ac_I_A) > 0.4) then
     write(30,*) "Activated complex belongs to Asymmetric top."
     write(30,*) "Rotation partition function of activated complex is", ac_rot_val
     end if

     if(ac_atomicity == 0) then
     write(30,*) "Again, your activated complex has only a single atom !!!"
     else 
         if(ac_freq_num == 1) then
         write(30,*) "Activated complex has", ac_freq_num, "vibrational frequency."
         write(30,*) "The vibrational frequincy is", ac_freq(1), "cm-1."
         write(30,*) "Vibrational partition function of activated complex is", ac_vib_val
         else if(ac_freq_num > 1) then
         write(30,*) "Activated complex has", ac_freq_num, "vibrational frequencies."
         do i = 1, ac_freq_num
           write(30,*) i, "frequency is", ac_freq(i), "cm-1."
         end do
         write(30,*) "Vibrational partition function of activated complex is", ac_vib_val
         end if
     end if
     
      write(30,*) "Total partition function of activated complex is",  ac_parti_val, "m-3."
     write(30,*) "    "
     
     write(30,*) "________________________Information about Temperature, Energy barrier and Rate constant_________________________"
     write(30,*) "Temperature is", Temp, "K."
     write(30,*) "The energy barrier of the reaction is", E_0, "kj mol-1."
     write(30,*) "The CTST rate constant of the reaction is", rate_val, "m3 s-1."
     write(30,*) "________________________________________________________________________________________________________________"
     write(30,*) "                                     | Hurray! A successful calculation! |                                       "
     write(30,*) "------------------------------------------------End of the file-------------------------------------------------"
     
     close(unit=30) 
     close(unit=10)
     call system("rm reshima")
end program rate_calc
     
