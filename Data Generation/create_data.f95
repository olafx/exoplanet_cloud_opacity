
! ==================================================================================
! Create a trainings set for a neural netwrok using a number of 'length' randomised
! points with the dust size limits 'size_min' and 'size_max'. The dust sepecies are
! given by the ones loaded into mie_data and the wavelengths are limited by the\
! input wavelength file.
program net_data_set
    ! import parameters from modules
    use mie, only: mie_calc, mie_getter, mie_init, d_wl
    implicit none
    ! input parameter
    integer      :: net_data_len  ! number of training points
    real(8)      :: dsi_min       ! minimum dust size
    real(8)      :: dsi_max       ! maximum dust size
    real(8)      :: std_max       ! maximum dust size dist
    real(8)      :: wav_min       ! minimum wavelength
    real(8)      :: wav_max       ! maximum wavelength
    real(8)      :: rng_param_a   ! param for mixture ratio RNG
    real(8)      :: rng_param_b   ! param for mixture ratio RNG
    real(8), dimension(16) :: vs_max ! species mask
    character(len=8) :: typ ! type of processing
    ! internal parameter
    integer                            :: i, j, k, l, nr_wl, spec1, spec2
    real(8)                            :: r_dsi, dsi_net, r_std, std, r_wav, wav_net
    real(8), dimension(:), allocatable :: q_sca_net, q_ext_net, wave
    real(8), dimension(16)             :: r_vsx
    
    ! set input parameters
    open(unit=85, file='data.txt', action='read', form='formatted')
    read(85, *); read(85, *) rng_param_a
    read(85, *); read(85, *) rng_param_b
    read(85, *); read(85, *) typ
    read(85, *); read(85, *) net_data_len
    read(85, *); read(85, *) dsi_min
    read(85, *); read(85, *) dsi_max
    read(85, *); read(85, *) std_max
    read(85, *); read(85, *) wav_min
    read(85, *); read(85, *) wav_max
    read(85, *); read(85, *); read(85, *) vs_max
    close(85)
    open(unit=85, file='wavelengths.txt', action='read', form='formatted')
    read(85, *) nr_wl
    close(85)
    
    ! initalise mie and read out important value
    call mie_init (mute         = .False., &
                   wave_file    = 'wavelengths.txt', &
                   dnf          = 'nk_origin/aa_dust_names.txt', &
                   nkp          = 'nk_origin/', &
                   do_holo      = .True., &
                   max_hol_frac = 0.85d0, &
                   n_vol        = 16, &
                   do_size_dist = .True., &
                   func         = 'logn', &
                   std_res      = 10, &
                   completness  = 0.95d0 )
                   
    call mie_getter(nwl = nr_wl, wavelengths = wave)
    
    ! initalise
    dsi_net = -1
    allocate(q_sca_net(nr_wl))
    allocate(q_ext_net(nr_wl))
	
    
    ! open file to store trainings set
    open(unit=81, file='input_set.txt', action='write', &
    status ='replace', form='formatted')
    open(unit=82, file='training_set.txt', action='write', &
    status ='replace', form='formatted')
    
    
! ==================================================================================
! calculate random mie values within the given parameters
    if (typ == 'mie_rand') then
		! calculate Mie values for a number of 'length' random points
		do i = 1, net_data_len
            print *, i
			! --------------------- generate random data ---------------
			! calculate a rondom dust size
			call random_seed()
			call random_number(r_dsi)
			dsi_net = (log10(dsi_max) - log10(dsi_min)) * r_dsi + log10(dsi_min) 
			dsi_net = 10**dsi_net ! transform from log to normal space
			
			! calcualte a random standard deviation
			call random_seed()
			call random_number(r_std)
			std = (log10(std_max) - (-7.d0)) * r_dsi + (-7.d0)
			std = 10**std ! transform from log to normal space
			
			! randomize wavelength
			call random_seed()
			call random_number(r_wav)
			wav_net = (log10(wav_max) - log10(wav_min)) * r_wav + log10(wav_min)
			wav_net = 10**wav_net ! transform from log to normal space
			d_wl(1) = wav_net
			
			! calcuate a random volume mixing ratio
			do j = 1, 16
				call random_seed()
				call random_number(r_vsx(j))
				r_vsx(j) = r_vsx(j) * vs_max(j) ! apply maximum fraction mask
			end do
            
            ! uniform [0, 1] to uniform [a, b]
            r_vsx = (rng_param_b - rng_param_a) * (r_vsx + rng_param_a) 
            ! uniform [0, 1] to non-uniform [0, inf)
            r_vsx = tan(r_vsx * (2.0d0 * atan(1.0d0)))

			r_vsx = r_vsx/sum(r_vsx) ! normalise random volumefraction to 1
			
			
			! calculate training data
			call mie_calc(r_vsx, dsi_net, Q_sca_net, Q_ext_net, r_std)
			
			! safe training data
			do j = 1, nr_wl
				! input data
				write(81, *) d_wl(j), dsi_net, r_std, r_vsx
				! output data
				write(82, *) Q_sca_net(j), Q_ext_net(j)
			end do
		end do
	end if
    
! ==================================================================================
! calculate a grid of mie values within the given parameters
    if (typ == 'mie_grid') then
		! set only two volume species
		i = 0
		do l = 1, 16
			if ( vs_max(l) > 1.d-6) then
				if (i == 0) then
					spec1 = l
					i = 1
				else
					spec2 = l
					exit
				end if
			end if
		end do
		
		r_vsx = 0.d0
		
		! calculate Mie values for a number of 'length' random points
		do i = 1, net_data_len
		do j = 1, net_data_len
		do k = 1, net_data_len
			! --------------------- generate random data ---------------
			! calculate a rondom dust size
			dsi_net = (log10(dsi_max) - log10(dsi_min)) * (i-1)/(net_data_len-1) &
			          + log10(dsi_min) 
			dsi_net = 10**dsi_net ! transform from log to normal space
			
			! calcualte a random standard deviation
			std = (log10(std_max) - (-7.d0)) * (j-1)/(net_data_len-1) + (-7.d0)
			r_std = 10**std ! transform from log to normal space
			
			! calcuate a random volume mixing ratio
			r_vsx(spec1) = (k-1.d0) / (net_data_len-1.d0) * vs_max(spec1)
			r_vsx(spec2) = 1.d0 - r_vsx(spec1)
			r_vsx = r_vsx/sum(r_vsx) ! normalise random volumefraction to 1
			
			! calculate training data
			call mie_calc(r_vsx, dsi_net, Q_sca_net, Q_ext_net, r_std)
			
			! safe training data
			do l = 1, nr_wl
				! input data
				write(81, *) wave(l), dsi_net, r_std, r_vsx
				! output data
				write(82, *) Q_sca_net(l), Q_ext_net(l)
			end do
			
		end do
		end do
		end do
	end if
	
	
    ! close file
    close(81); close(82)
    
end program net_data_set
