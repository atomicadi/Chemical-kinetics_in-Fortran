module mod_1_requirement_scratch
    implicit none
           real(8), parameter :: kg_fac = 1.67e-27           
           character(len=2) :: atom_symb(100)
           real(8) :: M_A, mass_atom(100), mass_tot, mass_tot_kg
           integer :: i, j, num_atom

           real(8) :: atom_mat(100,3), mass_atom_mat(100,3), mass_sum_atom_mat(1,3), COM(1,3), new_atom_mat(100,3)
           real(8) :: I_xx, I_yy, I_zz, I_xy, I_yz, I_xz, I_A, I_B, I_C
            
    contains

!==============Total Mass Calculation========================================
   subroutine calc_mass(atomicity, mass_tot, mass_atom, mass_tot_kg)
         implicit none
           integer :: i, atomicity
           real(8) :: mass_tot, mass_tot_kg, mass_atom(100)
      
           if (atomicity == 0) then 
           read(10,*) 
           read(10,*) M_A
           mass_tot = M_A
           mass_tot_kg = mass_tot * kg_fac
           else if(atomicity > 0) then   
               read(10,*) num_atom
               do i = 1, num_atom
                 read(10,*) mass_atom(i)
               end do
               mass_tot = 0.0
               do i = 1, num_atom
                  mass_tot = mass_tot + mass_atom(i)
               end do
               mass_tot_kg = mass_tot * kg_fac
          end if 
       end subroutine calc_mass

!================Moment of Inertia Calculation================================
   subroutine calc_iner(atomicity, mass_tot, mass_atom, COM, I_A, I_B, I_C)
         implicit none
           integer :: i, atomicity
           real(8) :: atom_mat(100,3), COM(1,3), mass_sum_atom_mat(1,3), I_A, I_B, I_C, mass_tot, mass_atom(100)
           integer, parameter :: dp = kind(1.0d0)
           integer :: n, lda, lwork, info
           double precision, allocatable :: A(:,:), w(:), work(:)

 
           if(atomicity == 0) then
            I_A = 0.0
            I_B = 0.0
            I_C = 0.0
           else if(atomicity > 0) then
             do i = 1, num_atom
              read(10,*) atom_symb(i), atom_mat(i,1), atom_mat(i,2), atom_mat(i,3)
             end do
             do i = 1, num_atom
               do j = 1,3
                 mass_atom_mat(i,j) = atom_mat(i,j) * mass_atom(i)
               end do
             end do
             mass_sum_atom_mat(1,:) = 0.0  
             do i = 1, num_atom
                do j = 1,3
                  mass_sum_atom_mat(1,j) = mass_sum_atom_mat(1,j) + mass_atom_mat(i,j)
                end do
             end do
             do j = 1, 3
                COM(1,j) = (mass_sum_atom_mat(1,j)/mass_tot)
             end do
             do i = 1, num_atom
                do j =1, 3
                  new_atom_mat(i,j) =  atom_mat(i,j) - COM(1,j) 
                end do
             end do
             
             I_xx = 0.0
               do i = 1, num_atom
                 I_xx = (mass_atom(i) * ((new_atom_mat(i,2)**2) + (new_atom_mat(i,3)**2))) + I_xx
               end do

             I_yy = 0.0
               do i = 1, num_atom
                  I_yy = (mass_atom(i) * ((new_atom_mat(i,1)**2) + (new_atom_mat(i,3)**2))) + I_yy
               end do

             I_zz = 0.0
                do i = 1, num_atom
                   I_zz = (mass_atom(i) * ((new_atom_mat(i,1)**2) + (new_atom_mat(i,2)**2))) + I_zz
                end do

             I_xy = 0.0
                do i = 1, num_atom
                  I_xy = (-(mass_atom(i) * ((new_atom_mat(i,1)) * (new_atom_mat(i,2))))) + I_xy
                end do

             I_xz = 0.0
                do i = 1, num_atom
                  I_xz = (-(mass_atom(i) * ((new_atom_mat(i,1)) * (new_atom_mat(i,3))))) + I_xz
                end do

             I_yz = 0.0
               do i = 1, num_atom
                 I_yz = (-(mass_atom(i) * ((new_atom_mat(i,2)) * (new_atom_mat(i,3))))) + I_yz
               end do
             
            if(I_xy == 0.0 .and. I_xz == 0.0 .and. I_yz == 0.0) then
            I_A = I_xx
            I_B = I_yy
            I_C = I_zz
            else
                n = 3
                lda = n
                allocate(A(n,n), w(n))
 
                A = reshape([ I_xx, I_xy, I_xz, &
                              I_xy, I_yy, I_yz, &
                              I_xz, I_yz, I_zz ], shape(A))

                lwork = -1
                allocate(work(1))
                call dsyev('V', 'U', n, A, lda, w, work, lwork, info)

               lwork = int(work(1))
               deallocate(work)
               allocate(work(lwork))

               call dsyev('V', 'U', n, A, lda, w, work, lwork, info)
               I_A = w(1)
               I_B = w(2)
               I_C = w(3)
            end if
        end if
   end subroutine calc_iner
end module