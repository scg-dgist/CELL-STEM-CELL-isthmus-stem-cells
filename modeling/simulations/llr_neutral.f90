! ---------------------------------------------------------------------------- !
!     Lateral neutral clonal dynamics                                          !
!     David Jorg                                                               !
!																																							 !
!     compile:																																 !
!				gfortran -cpp -Dsingle llr_neutral.f90 -ffree-line-length-none         !
!       -fbackslash -Ofast -o llr_neutral                                      !
! ---------------------------------------------------------------------------- !

program clone
  implicit none

	character(255),   parameter                           :: data_directory = './results/'

	integer,          dimension(:),         allocatable   :: mod2, rmod_numseg2

  real(kind=8),     parameter                           :: pi2 = 6.283185307179586

	integer,          parameter                           :: N_L = 1, N_R = 2

	real(kind=8)                                          :: induction_frequency, lambda, omega_in, omega_out, sim_time, eq_time, dt, potential_width
  integer                                               :: barrier_equivalent, num_segments, num_colors

  integer                                               :: num_runs, bins_clones, bins_angles

	! derived parameters
	integer                                               :: readout_points, num_segments_2, max_gland_size

  integer,          dimension(:,:),       allocatable   :: neighbors

  ! fields
  integer,          dimension(:),         allocatable   :: field, pres

  ! computation variables
	real(kind=8)                                          :: tp, die, total_t
	real(kind=8),     dimension(:),         allocatable   :: potential
	real(kind=8),     dimension(:),         allocatable   :: pp_field
	integer                                               :: update_field, expanding_field, gland_size, unlabelled_gland, num_clones, num_labelled_clones, barrier_number, dist
  integer,          dimension(:),         allocatable   :: relative_clone_sizes, global_unlabelled_glands, global_num_clones, global_num_labelled_clones, barrier_angles, barrier_list
	integer,          dimension(:,:),       allocatable   :: global_barrier_number_dist
	integer,          dimension(:,:,:),     allocatable   :: global_barrier_angle_dist

  integer,          dimension(:,:),       allocatable   :: global_relative_clone_sizes

	real(kind=8)                                          :: reaction, p_sum

  integer                                               :: rnd, num_exp_fields, candidate_field, n_steps
  integer,          dimension(0:2)                      :: potential_exp_fields
  integer,          dimension(:),         allocatable   :: barrier_segment_list
  logical                                               :: iterate, searching, equilibrated

  ! counters
  integer                                               :: i, ic, j, jc, k, r, c, d, rp, lrp
	real(kind=8)                                          :: time_step, t, t_readout

	character(255)                                        :: directory, prefix, filename

	! random number generator
	!type(fgsl_rng_type)                                   :: rnd_type
  !type(fgsl_rng)                                        :: rng
	!integer(8)                                            :: seed

  integer(4)                                            :: seed

	! command line arguments
	integer,          parameter                           :: req_num_args = 16
  integer                                               :: num_args
  character(255),   dimension(1:req_num_args)           :: cl_args

	! functions
	integer                                               :: rmod, ceiling2

  write(*,'(a)') achar(27)//'[30m[Clonal dynamics]' // achar(27)// '[0m'
  num_args=iargc()
	if (num_args.ne.req_num_args) then
		write(*,'(a)') 'Wrong number of arguments. Terminated.'
		call exit()
	end if
  do i=1,num_args
    call getarg(i, cl_args(i))
  end do

	read(cl_args(1),*) num_segments
	read(cl_args(2),*) induction_frequency
	read(cl_args(3),*) lambda
	read(cl_args(4),*) omega_in
	read(cl_args(5),*) omega_out
	read(cl_args(6),*) barrier_equivalent
	read(cl_args(7),*) potential_width
  read(cl_args(8),*) num_colors
	read(cl_args(9),*) sim_time
	read(cl_args(10),*) eq_time
	read(cl_args(11),*) dt
	read(cl_args(12),*) bins_clones
	read(cl_args(13),*) bins_angles
	read(cl_args(14),*) seed
	read(cl_args(15),*) num_runs
	prefix = trim(cl_args(16))

	! derived parameters
	readout_points = floor(sim_time / dt)
	num_segments_2 = 2 * num_segments

	max_gland_size = num_segments * (1 + barrier_equivalent)

	! allocate arrays
	allocate(mod2(0:2*num_segments_2), rmod_numseg2(0:2*num_segments_2))

	allocate(neighbors(1:num_segments_2,1:2))
	allocate(field(1:num_segments_2), pres(1:num_segments_2))
	allocate(potential(0:max_gland_size))

	allocate(pp_field(1:num_segments_2), barrier_list(1:num_segments))

	allocate(relative_clone_sizes(0:bins_clones), barrier_angles(1:num_segments**2))
	allocate(global_relative_clone_sizes(0:bins_clones,0:readout_points), global_unlabelled_glands(0:readout_points), global_num_clones(0:readout_points), global_num_labelled_clones(0:readout_points), global_barrier_number_dist(0:num_segments,0:readout_points), global_barrier_angle_dist(2:num_segments,0:bins_angles,0:readout_points))

  ! create output directory
	directory = trim(data_directory) // trim(prefix) // '/'
  call system('mkdir -p ' // trim(directory))

	! initialize random number generator
	call init_random_seed(seed)

	! static arrays
	do i=0,2*num_segments_2
		mod2(i) = mod(i, 2)
		rmod_numseg2(i) = rmod(i, num_segments_2)
	end do

	! neighborhood relations
	do i=1,num_segments_2
  	neighbors(i,N_L) = i - 1
  	neighbors(i,N_R) = i + 1
  end do
  neighbors(1,N_L) = num_segments_2
  neighbors(num_segments_2,N_R) = 1

	! initialize potential
	potential = 1.
	do i=0,floor(0.5 * real(num_segments,8))
		potential(i) = exp( - (real(i,8) - 0.5 * real(num_segments,8))**2. / (2. * potential_width**2.) )
	end do

	! --- start queue ---
	global_relative_clone_sizes = 0.
	global_unlabelled_glands = 0
	global_num_clones = 0
	global_num_labelled_clones = 0
	total_t = 0.
	n_steps = 0

	do r=1,num_runs
#ifdef single
		if (mod(r,25)==0) write (*,'(a,f5.1,a)',advance='no') '\b\b\b\b\b\b', (real(r,8)/real(num_runs,8))*100.0, '%'
#endif

		! initialize timers
		t = -eq_time
		equilibrated = .false.
		rp = -1
		lrp = -1

		! initialize fields
		field = 0
		pres = 0
		do i=1,num_segments_2,2
			pres(i) = 1
		end do

		! initial observables
		! -- clone sizes
		call get_clone_sizes(barrier_equivalent, num_segments, neighbors, field, pres, bins_clones, relative_clone_sizes, unlabelled_gland, num_clones, num_labelled_clones, gland_size)
		! -- barrier number
		call get_barrier_states(mod2, rmod_numseg2, gland_size, barrier_equivalent, num_segments_2, pres, bins_angles, barrier_number, barrier_angles, barrier_list)

		! time evolution
		do while (lrp<=readout_points)

			! label the system when equilibration of barrier segments is reached
			if (t>=0..and.(equilibrated.eqv..false.)) then
				t = 0.
				rp = 0
				equilibrated = .true.

				! initial conditions: label fields with a probability that matches the induction frequency
				field = 0
				do i=1,num_segments_2,2
					call random_number(die)
					if (die<induction_frequency) then
						call random_number(die)
						field(i) = ceiling(num_colors * die)
					end if
				end do

				! initial observables
				! -- clone sizes
				call get_clone_sizes(barrier_equivalent, num_segments, neighbors, field, pres, bins_clones, relative_clone_sizes, unlabelled_gland, num_clones, num_labelled_clones, gland_size)
				! -- barrier number
				call get_barrier_states(mod2, rmod_numseg2, gland_size, barrier_equivalent, num_segments_2, pres, bins_angles, barrier_number, barrier_angles, barrier_list)

				! create output files
#ifdef saveall
				! labeling trajectory
				write(filename, '(a12,i0.3,a4)') 'realization_', r, '.csv'
				open(unit=1,file=trim(directory) // trim(filename))

				write(1,'(f23.12,a)',advance='no') t, '    '
				do i=1,num_segments_2
					write(1,'(i3,a)',advance='no') field(i), ' '
				end do
				write(1,*)

				write(filename, '(a21,i0.3,a4)') 'relative_clone_sizes_', r, '.csv'
				open(unit=2,file=trim(directory) // trim(filename))

				write(2,'(f23.12,a)',advance='no') t, '    '
	  		do j=0,bins_clones
	  			write(2,'(i20,a)',advance='no') relative_clone_sizes(j), ' '
	  		end do
		  	write(2,*)

				! barrier trajectory
				write(filename, '(a9,i0.3,a4)') 'presence_', r, '.csv'
				open(unit=3,file=trim(directory) // trim(filename))

				write(3,'(f23.12,a)',advance='no') t, '    '
				do i=1,num_segments_2
					write(3,'(i3,a)',advance='no') pres(i), ' '
				end do
				write(3,*)
#endif
			end if

			! partial propensities
			! -- initialize partial propensities for proliferating segments
			pp_field = lambda

			! total propensity
			tp = sum(pp_field)

			! advance time
			call random_number(die)
			time_step = - (1./tp) * log(die)

			call random_number(die)
			reaction = tp * die

			t = t + time_step
			total_t = total_t + time_step
			n_steps = n_steps + 1

			! reactions
			p_sum     = 0.
			searching = .true.

			! -- loss and replacement in dynamic segments
			do i=1,num_segments_2,2
				p_sum = p_sum + pp_field(i)
				if (reaction<p_sum.and.searching) then

					! check status of left and right neighbors
					num_exp_fields       = 0
					potential_exp_fields = 0
					do d=1,2
						! if neighbor in d-direction is not blocked by barrier...
						if (pres(neighbors(i, d))==0) then
							! ... add to potential expansion sites
							num_exp_fields = num_exp_fields + 1
							potential_exp_fields(num_exp_fields) = neighbors(neighbors(i, d), d)
						end if
					end do

					! - if there is no potentially expanding field, nothing happens
					! - if there is only one potentially expanding field, this field will expand
					! - if there are two potentially expanding fields, one of these fields are picked randomly
					select case (num_exp_fields)
						case (1)
							field(i) = field(potential_exp_fields(1))
						case (2)
							call random_number(die)
							field(i) = field(potential_exp_fields(ceiling2(2. * die)))
					end select

					searching = .false.
					exit
				end if
			end do

			! data analysis
			rp = floor(t/dt)
			if (rp>lrp.and.equilibrated) then
				rp = min(rp, readout_points+1)

				do k=lrp+1,rp
					global_relative_clone_sizes(:,k) = global_relative_clone_sizes(:,k) + relative_clone_sizes(:)
					global_unlabelled_glands(k) = global_unlabelled_glands(k) + unlabelled_gland
					global_num_clones(k) = global_num_clones(k) + num_clones
					global_num_labelled_clones(k) = global_num_labelled_clones(k) + num_labelled_clones

					global_barrier_number_dist(barrier_number,k) = global_barrier_number_dist(barrier_number,k) + 1
					do i=1,num_segments_2
						if (barrier_angles(i)>0) then
							global_barrier_angle_dist(barrier_number,barrier_angles(i),k) = global_barrier_angle_dist(barrier_number,barrier_angles(i),k) + 1
						end if
					end do
				end do

				lrp = rp
			end if

			! compute observables
			! -- clone sizes
			call get_clone_sizes(barrier_equivalent, num_segments, neighbors, field, pres, bins_clones, relative_clone_sizes, unlabelled_gland, num_clones, num_labelled_clones, gland_size)
			! -- barrier number
			call get_barrier_states(mod2, rmod_numseg2, gland_size, barrier_equivalent, num_segments_2, pres, bins_angles, barrier_number, barrier_angles, barrier_list)

			! write results to files
#ifdef saveall
	 			write(1,'(f23.12,a)',advance='no') t, '    '
	 			do i=1,num_segments_2
	 				write(1,'(i3,a)',advance='no') field(i), ' '
	 			end do
	 			write(1,*)

				write(2,'(f23.12,a)',advance='no') t, '    '
	  		do j=0,bins_clones
	  			write(2,'(i20,a)',advance='no') relative_clone_sizes(j), ' '
	  		end do
		  	write(2,*)

	 			write(3,'(f23.12,a)',advance='no') t, '    '
	 			do i=1,num_segments_2
	 				write(3,'(i3,a)',advance='no') pres(i), ' '
	 			end do
	 			write(3,*)
#endif
		end do

#ifdef saveall
		close(1)
		close(2)
		close(3)
#endif

	end do
	! --- end of queue ---

	! write parameter file
  open(unit=1,file=trim(directory) // 'parameters.csv')
 	  call write_parameters(1, lambda, omega_in, omega_out, num_segments, num_colors, sim_time, eq_time, num_runs, induction_frequency, barrier_equivalent, potential_width, bins_clones, bins_angles, dt)
	close(1)

	! write clone size statistics
  open(unit=1,file=trim(directory) // 'global_relative_clone_sizes.csv')
  	do k=0,readout_points
  		do j=0,bins_clones
  			write(1,'(i20,a)',advance='no') global_relative_clone_sizes(j,k), ' '
  		end do
  		write(1,*)
  	end do
	close(1)

	open(unit=1,file=trim(directory) // 'global_unlabelled_glands.csv')
  	do k=0,readout_points
  		write(1,'(i20,a)') global_unlabelled_glands(k)
  	end do
	close(1)

	open(unit=1,file=trim(directory) // 'global_num_clones.csv')
  	do k=0,readout_points
  		write(1,'(i20,a)') global_num_clones(k)
  	end do
	close(1)

	open(unit=1,file=trim(directory) // 'global_num_labelled_clones.csv')
  	do k=0,readout_points
  		write(1,'(i20,a)') global_num_labelled_clones(k)
  	end do
	close(1)

	! write barrier statistics
	open(unit=1,file=trim(directory) // 'global_barrier_number_dist.csv')
  	do k=0,readout_points
			do j=0,num_segments
  			write(1,'(i20,a)',advance='no') global_barrier_number_dist(j,k), ' '
  		end do
  		write(1,*)
  	end do
	close(1)

	do i=2,num_segments/2 ! save only up to half of the possible number to save space
		write(filename, '(a26,i0.2,a4)') 'global_barrier_angle_dist_', i, '.csv'
		open(unit=1,file=trim(directory) // trim(filename))
	  	do k=0,readout_points
				do j=0,bins_angles
	  			write(1,'(i20,a)',advance='no') global_barrier_angle_dist(i,j,k), ' '
	  		end do
	  		write(1,*)
	  	end do
		close(1)
	end do

	write(*,*)
	write(*,*) trim(prefix)
	write(*,*) 'avg time step:', total_t/real(n_steps,8)
	write(*,*) ' - rel. to dt:', (total_t/real(n_steps,8))/dt
end program clone


subroutine get_clone_sizes(barrier_equivalent, num_segments, neighbors, field, pres, bins_clones, relative_clone_sizes, unlabelled_gland, num_clones, num_labelled_clones, gland_size)
	implicit none

	integer,    parameter                                     :: N_L = 1, N_R = 2

	integer,															       intent(in)   :: barrier_equivalent, num_segments, bins_clones
	integer,    dimension(1:2*num_segments,1:2), intent(in)   :: neighbors
	integer,    dimension(1:2*num_segments),     intent(in)   :: field, pres

	! second index of clone_boundaries indicates left/right boundary
	integer,    dimension(1:2*num_segments,1:2)               :: clone_boundaries
	integer,    dimension(1:2*num_segments)                   :: clone_label
	logical                                                   :: barrier
	integer                                                   :: num_segments_2, cur_clone_size, first_clone_size, i, j, k, d, left_boundary, right_boundary, bin, n

	integer, dimension(0:bins_clones),            intent(out) :: relative_clone_sizes
	integer,                                      intent(out) :: unlabelled_gland, num_clones, num_labelled_clones, gland_size

	! functions
	integer                                                 :: rmod

	num_segments_2 = 2 * num_segments

	clone_boundaries     = -1
	relative_clone_sizes = 0

	num_clones           = 1
	num_labelled_clones  = 0
	cur_clone_size       = 1

	unlabelled_gland     = 0

	! determine the current size of the gland
	gland_size = num_segments
	do i=2,num_segments_2,2
		gland_size = gland_size + pres(i) * barrier_equivalent
	end do

	! obtain clone boundaries
	do i=1,num_segments_2
		if (pres(i)==1) then
			! left clone boundary is detected ...
			! if left neighbor is present and has different label
			! if left neighbor is absent (only possible if barrier) and the next-nearest left neighbor has different label
			n = neighbors(i, N_L)
			if ( ((pres(n)==1).and.(field(n).ne.field(i))).or.((pres(n)==0).and.(field(neighbors(n, N_L)).ne.field(i))) ) then
				clone_boundaries(num_clones, N_L) = i
				clone_label(num_clones) = field(i)
			end if

			! symmetric for right boundary
			n = neighbors(i, N_R)
			if ( ((pres(n)==1).and.(field(n).ne.field(i))).or.((pres(n)==0).and.(field(neighbors(n, N_R)).ne.field(i))) ) then
				clone_boundaries(num_clones, N_R) = i
				num_clones = num_clones + 1
			end if
		end if
	end do


	! catch monoclonality
	if (num_clones==1) then
		! unlabelled clones will not be counted
		! field(1) is as good as any field sind gland is monoclonal
		if (field(1).ne.0) then
			relative_clone_sizes(bins_clones) = relative_clone_sizes(bins_clones) + 1
			num_labelled_clones = 1
		else
			unlabelled_gland = 1
		end if
	else
		! check whether first left clone boundary has been detected
		if (clone_boundaries(1, N_L)==-1) then
			clone_boundaries(1, N_L) = clone_boundaries(num_clones, N_L)
			clone_label(1) = clone_label(num_clones)
		end if

		! last right boundary leads to one superficial increment in num_clones
		num_clones = num_clones - 1

		! determine clone sizes
		do i=1,num_clones

			! unlabelled clones will not be counted
			if (clone_label(i).ne.0) then
				cur_clone_size = 0

				left_boundary  = clone_boundaries(i, N_L)
				right_boundary = clone_boundaries(i, N_R)

				! check if index of left bdry is larger than right bdry
				if (left_boundary>right_boundary) right_boundary = right_boundary + num_segments_2

				do j=left_boundary,right_boundary
					! real index of segment
					k = rmod(j, num_segments_2)

					if (mod(k, 2) == 1) then
						! if normal segment, add unit size
						cur_clone_size = cur_clone_size + 1
					else
						! if barrier segment and present, add barrier size
						cur_clone_size = cur_clone_size + pres(k) * barrier_equivalent
					end if
				end do

				bin = nint(real(cur_clone_size,8) / real(gland_size,8) * real(bins_clones,8))

				relative_clone_sizes(bin) = relative_clone_sizes(bin) + 1

				num_labelled_clones = num_labelled_clones + 1
			end if

		end do

	end if
end subroutine get_clone_sizes


subroutine get_barrier_states(mod2, rmod_numseg2, gland_size, barrier_equivalent, num_segments_2, pres, bins_angles, barrier_number, barrier_angles, barrier_list)
	implicit none

	integer,															       intent(in)   :: gland_size, barrier_equivalent, num_segments_2, bins_angles
	integer,    dimension(1:num_segments_2),     intent(in)   :: pres
	integer,    dimension(0:2*num_segments_2),   intent(in)   :: mod2, rmod_numseg2

	! second index of clone_boundaries indicates left/right boundary
	integer                                                   :: i, ic, j, jc, k, kc, m
	integer                                                   :: dist, barrier_angle

	integer,                                      intent(out) :: barrier_number
	integer, dimension(1:(num_segments_2/2)**2),  intent(out) :: barrier_angles
	integer, dimension(1:num_segments_2/2),       intent(out) :: barrier_list

	! functions
	integer                                                   :: rmod


	! barrier number
	m = 0
	barrier_list = 0
	do i=2,num_segments_2,2
		if (pres(i)==1) then
			m = m + 1
			barrier_list(m) = i
		end if
	end do
	barrier_number = m

	! barrier angles
	barrier_angles = 0
	m = 0
	if (barrier_number>1) then
		do i=1,barrier_number
			ic = barrier_list(i)

			do j=i+1,barrier_number
				jc = barrier_list(j)

				call get_barrier_distance(mod2, rmod_numseg2, gland_size, barrier_equivalent, num_segments_2, pres, ic, jc, dist)

				m = m + 1
				barrier_angles(m) = ceiling(real(dist,8) / real(gland_size,8) * real(bins_angles,8))
			end do

		end do
	end if
end subroutine get_barrier_states


subroutine get_barrier_distance(mod2, rmod_numseg2, gland_size, barrier_equivalent, num_segments_2, pres, ic, jc, dist)
	implicit none

	integer,															       intent(in)   :: barrier_equivalent, num_segments_2, gland_size
	integer,    dimension(1:num_segments_2),     intent(in)   :: pres
	integer,    dimension(0:2*num_segments_2),   intent(in)   :: mod2, rmod_numseg2

	integer                                                   :: ic, jc, k, kc, md, d1, d2

	integer,	   													        intent(out) :: dist

	! functionsg
	integer                                                   :: rmod

	if (ic>jc) jc = jc + num_segments_2

	dist = 0
	do k=ic+1,jc
		kc = rmod_numseg2(k)
		dist = dist + mod2(kc) + (1 - mod2(kc)) * pres(kc) * barrier_equivalent
	end do

	dist = min(dist, gland_size - dist)
end subroutine get_barrier_distance


! note: modulo which starts from 1 instead of 0 is given by rmod(n,m) = mod(n-1,m)+1
function rmod(n, m)
	implicit none

	integer                                               :: rmod
	integer,															   intent(in)   :: n, m

	rmod = modulo(n - 1, m) + 1
	return
end function


! modified ceiling function which satisfies ceiling2(0.) = 1 instead of ceiling2(0.) = 0
! to avoid singular out-of-bound errors when random numbers return exactly 0.
function ceiling2(q)
	implicit none

	integer                                               :: ceiling2
	real(kind=8),												     intent(in)   :: q

	if (q==0.) then
		ceiling2 = 0 ! temporary: check for occurrences of 0
	else
		ceiling2 = ceiling(q)
	end if

	return
end function


subroutine write_parameters(output_unit, lambda, omega_in, omega_out, num_segments, num_colors, sim_time, eq_time, num_runs, induction_frequency, barrier_equivalent, potential_width, bins_clones, bins_angles, dt)
	implicit none

	integer,															   intent(in)   :: output_unit
	integer,															   intent(in)   :: num_segments, num_colors, num_runs, bins_clones, bins_angles, barrier_equivalent
	real(kind=8),                            intent(in)   :: lambda, omega_in, omega_out, sim_time, eq_time, induction_frequency, potential_width, dt

  write(output_unit,'(a20,f23.10)')  'lambda',               lambda
  write(output_unit,'(a20,f23.10)')  'omega_in',             omega_in
  write(output_unit,'(a20,f23.10)')  'omega_out',            omega_out

  write(output_unit,'(a20,i7)')      'num_segments',	  		 num_segments
  write(output_unit,'(a20,f23.10)')  'potential_width',	     potential_width
  write(output_unit,'(a20,i7)')      'barrier_equivalent',	 barrier_equivalent
  write(output_unit,'(a20,i7)')      'num_colors',					 num_colors
  write(output_unit,'(a20,f23.10)')  'induction_frequency',  induction_frequency
  write(output_unit,'(a20,f23.12)')  'sim_time',						 sim_time
	write(output_unit,'(a20,f23.12)')  'eq_time', 						 eq_time
  write(output_unit,'(a20,i7)')      'num_runs',	  				 num_runs
  write(output_unit,'(a20,i7)')      'bins_clones',	 				 bins_clones
	write(output_unit,'(a20,i7)')      'bins_angles',	 				 bins_angles
	write(output_unit,'(a20,f23.12)')  'dt',	        				 dt
end subroutine write_parameters


SUBROUTINE init_random_seed(m)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: m
  INTEGER :: i, n
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  seed = m + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE


subroutine exit_call(msg)
	implicit none
	character(255), intent(in) :: msg

	write(*,*) trim(msg)
	call exit()
end subroutine
