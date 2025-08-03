! Aditya Barman, Graduate Student, Weizmann Institute of Science, June 23, 2025
! Email: atomicadi2023@gmail.com
program RRKM_rate_calc
       implicit none
       integer :: i, j, k, l, m, n, a, b, sum_count, total_ro_count, count_ro_list(11)
          real, allocatable :: freq_list(:)
                       real :: ro_list(11) , E, E_0, val, temp_val, h, rrkm_rate
                   character(len=50) :: label

    h = 33.358e-12 !(in cm-1.s unit)
    open(unit=10,file='input.in')
    read(10,*) E
    read(10,*) E_0
    read(10,*) n
    allocate(freq_list(n))
    do i = 1, n
      read(10,*) freq_list(i)
    end do
    val = E - E_0
    do i= 0, 10
       ro_list(i+1) =  (E/10.0)*i 
    end do 
     
    b = 0
    do a = 1, 11
    if (n == 5) then
       do i = 0, 9
           do j = 0, 9
               do k = 0, 9
                 do l = 0, 9
                   do m = 0, 9
                   temp_val = i*freq_list(1) + j*freq_list(2) + k*freq_list(3) + l*freq_list(4) + m*freq_list(5)
                   if (temp_val <= ro_list(a)) then
                   b = b + 1
                   else
                       exit
                   end if 
                   end do
                 end do
               end do
            end do
          end do
          count_ro_list(a) = b
          b = 0
         

       else if (n == 4) then
         do i = 0, 9
           do j = 0, 9
               do k = 0, 9
                 do l = 0, 9
                  temp_val = i*freq_list(1) + j*freq_list(2) + k*freq_list(3) + l*freq_list(4)  
                   if (temp_val <= ro_list(a)) then
                   b = b + 1
                   else
                       exit
                   end if 
                   end do
                 end do
               end do
            end do
           count_ro_list(a) = b
           b = 0
          
  
       else if (n == 3) then
        do i = 0, 9
           do j = 0, 9
               do k = 0, 9
                 temp_val = i*freq_list(1) + j*freq_list(2) + k*freq_list(3) 
                 if (temp_val <= ro_list(a)) then
                 b = b + 1
                 else
                      exit
                  end if 
                end do
              end do
            end do
            count_ro_list(a) = b
            b = 0
          

        else if (n == 2) then
         do i = 0, 9
           do j = 0, 9
              temp_val = i*freq_list(1) + j*freq_list(2) 
               if (temp_val <= ro_list(a)) then
                b = b + 1
                else
                    exit
                end if 
              end do
           end do
          count_ro_list(a) = b
          b = 0

         end if
     end do

     total_ro_count = 0
     do i = 1, 11
       total_ro_count = total_ro_count + count_ro_list(i)
     end do

     sum_count = 0
     do i = 1, 11
       if (val >= ro_list(i)) then
          sum_count = sum_count + count_ro_list(i)
       end if
     end do

     rrkm_rate= sum_count/(h*total_ro_count)

     open(unit=20, file='output_rrkm.out')
     write(20, '(A)') "Authors: Philips Kumar Rai and Aditya Barman"
     write(20, '(A)') "Email: atomicadi2023@gmai.com"
     write(20,'(A)') "==================================Input data=================================="
     write(20, '(A, F12.6)') "Total energy (in cm-1) =", E
     write(20, '(A, f12.6)') "Threshold energy (in cm-1) =", E_0
     write(20,'(A)') "==================================Output data=================================="
     write(20, '(A)') "Density of states at various discrete energy levels upto total energy (E):" 
     do i = 1, 11
      write(label, '(F12.6)') ro_list(i) 
      write(20, '(A,A,A,I5)') 'ro_', trim(adjustl(label)),'=', count_ro_list(i)
     end do
     write(20, '(A, I5)') "The sum of states N(E-E_0) =", sum_count
     write(20, '(A, I5)') "Total density of states ro(E) =", total_ro_count
     write(20, *)    "RRKM rate constant k(E) in s-1 unit=", rrkm_rate
     write(20,*)
     write(20, *) "     ------------------End of the file--------------------         "
     close(unit=10)
     close(unit=20)
     deallocate(freq_list)
end program RRKM_rate_calc 
      
    
