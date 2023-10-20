! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use chem_def, only: chem_isos, ic12, ic13, in14, in15, io16, io17, io18
      
      implicit none
      logical :: active_thermalpulse
      integer :: TP_count, dredgeup
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (.not. restart) then
            TP_count = 0
            dredgeup = 0
         end if

      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 5
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         real(dp) :: age_at_TP, age_TP1, mcore_at_TP, mcore_TP1, surface_opacity, surface_opacity_TP1
         real(dp) :: c12_c13, c12_c13_TP1, n14_n15, n14_n15_TP1, c_o, c_o_TP1, surface_c12, prev_surface_c12
         integer, intent(out) :: ierr
         integer :: k, surf_type
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Counting TPs and retrieving stats
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         surface_c12 = s% xa(s% net_iso(ic12), 1)

         if (s% center_he4 .lt. 0.000001) then  

             if ((.not. active_thermalpulse) .and. (s% power_he_burn .gt. 1e4)) then 
                 active_thermalpulse = .true.
                 TP_count = TP_count + 1
                 write(*,*) 'starting thermal pulse', TP_count
                 do k = 1, s% nz
                    surface_opacity = s% opacity(1)
                 enddo
                 mcore_at_TP = s% he_core_mass
                 age_at_TP = s% star_age

                ! surface isotope ratios
                 c12_c13 = (s% xa(s% net_iso(ic12), 1)) / (s% xa(s% net_iso(ic13), 1))
                 n14_n15 = (s% xa(s% net_iso(in14), 1)) / (s% xa(s% net_iso(in15), 1))
                 c_o = (s% xa(s% net_iso(ic12), 1) + s% xa(s% net_iso(ic13), 1)) / (s% xa(s% net_iso(io16), 1) + s% xa(s% net_iso(io17), 1) + s% xa(s% net_iso(io18), 1))
               
                 if (c_o .le. 0.45) then
                    surf_type = 0 ! undef
                 else if ((c_o .gt. 0.45) .and. (c_o .le. 0.55)) then
                    surf_type = 1 ! M type                
                 else if ((c_o .gt. 0.55) .and. (c_o .le. 0.95)) then
                    surf_type = 2 ! S type        
                 else if ((c_o .gt. 0.95) .and. (c_o .le. 1.05)) then
                    surf_type = 3 ! SC type                           
                 else if (c_o .gt. 1.05) then
                    surf_type = 4 ! C type                
                 end if

                ! quantities at outset of TP-AGB phase
                 if (TP_count == 1) then
                 mcore_TP1 = s% he_core_mass
                 age_TP1 = s% star_age   
                 surface_opacity_TP1 = surface_opacity
                 c12_c13_TP1 = c12_c13
                 n14_n15_TP1 = n14_n15
                 c_o_TP1 = c_o
                 endif

                 open(1, file='agb_stats.dat', action='write', position='append')
                 write(1,*) 'TP_count=', TP_count, 'age_at_TP=', age_at_TP, 'mcore_at_TP=', mcore_at_TP, 'surface_opacity', surface_opacity, 'c12_c13', c12_c13, 'n14_n15', n14_n15, 'c_o', c_o, 'type', surf_type 
                 close(1)
             else
                 if (s% power_h_burn/s% power_he_burn .gt. 10) active_thermalpulse = .false. 
                  write(*,*) 'Thermal pulse over'        
             end if
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! measuring dredge up
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (active_thermalpulse .eqv. .true.) then
            dredgeup = 0
         else
            if (TP_count == 0) then
               dredgeup = 0
            else if (TP_count .gt. 0) then
               if (surface_c12 .gt. prev_surface_c12) then
                  dredgeup = TP_count
               else 
                  dredgeup = 0
               end if
            end if
         end if 

         prev_surface_c12 = surface_c12

   !      mag_conv_env = 0
   !      do k=1, s% nz
   !         if (mixing_type = 1) then
   !            mag_conv_env = mag_conv_env + 1
   !         end if
   !      enddo

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         names(1) = 'c12_c13'
         vals(1) = c12_c13

         names(2) = 'n14_n15'
         vals(2) = n14_n15

         names(3) = 'c_o'
         vals(3) = c_o

         names(4) = 'dredgeup'
         vals(4) = dredgeup

         names(5) = 'surf_type'
         vals(5) = surf_type        

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.

      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr, k
         real(dp) :: h_env_frac
         real(dp) :: tau_tp_agb, delta_mcore

         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Termination condition: third to last thermal pulse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         if (TP_count == 41) extras_finish_step = terminate ! leave? not?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Termination condition: Post-AGB EEP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         h_env_frac = (s% star_mass - s% he_core_mass) / s% star_mass
         write(*,*) 'fraction of H envelope left=', h_env_frac

         ! stop if h_env_frac is less than 0.2
         if (h_env_frac < 0.2) extras_finish_step = terminate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writing out final surface conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

         if (extras_finish_step == terminate) then
           ! tau_tp_agb = age_TP1 - s% star_age
           ! delta_mcore = mcore_TP1 - s% he_core_mass
           !  write(*,*) 'Total TP-AGB lifetime:', tau_tp_agb, 'years'
           !  write(*,*) 'H exhausted core mass grew by:', delta_mcore, 'Msun'     
         endif

      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      end module run_star_extras
      
