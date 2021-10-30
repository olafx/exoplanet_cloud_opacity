

module rosland
	
	private
	! internal variables to remmebre array sizes
	integer :: nr_wl = -1  ! number of wavlengths
	integer :: nr_vol = -1 ! number of dust species
	real(8), dimension(:), allocatable :: wave ! wavelengths
	
	! parameters
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0) ! pi, you know it
    real(8), parameter :: c  = 2.99792458d8        ! light speed [m/s]
    real(8), parameter :: h  = 6.6260755d-34       ! plank constant [Js]
    real(8), parameter :: kb = 1.380658d-23        ! boltzman constant [J/K]
	
	! public functions to be called for rosland mean calculaiton
	public :: rosland_mean ! calculated rosland mean based on mie calc
	public :: rosland_init ! set up rosland mean claculation
	
	! internal functions
	private :: plank_derivative_function ! what the name says
	
	contains
	
    ! ==================================================================================
    ! Initalises Rosland mean calculation. Needs to be run after mie_init and before
    ! rosland_mean.
	subroutine rosland_init()
		use mie, only: mie_getter
		implicit none
		
		! initalise values from mie calc
		call mie_getter(nwl = nr_wl, nvol = nr_vol, )
	
	end subroutine rosland_init
	
	

	subroutine rosland_mean(vs, r_dust, std_dust, n_dust, temp, &
	                        lam_low, lam_up, kc)
	    use mie, only: mie_calc
		implicit none
		
		! input variables
		real*8 :: r_dust ! dust particle radius
		real*8 :: std_dust ! dust deviation
		real*8 :: n_dust ! particle number dust density
		real*8 :: temp ! dust temperature
		real*8 :: lam_low ! lower wavelength limit
		real*8 :: lam_up  ! upper wavelength limit
		real*8, dimension(nr_vol) :: vs ! volume fraction of particles
		
		! output variables
		real*8 :: kc ! output opacity, averaged from lam_low to lam_up
		
		! Internal variables
		integer :: j
        real(8) :: ros_mean_plank, ros_mean
		real*8, dimension(nr_wl) :: Q_sca, Q_ext, plank_der, kext
		
		! mie calculation for each grid point
		call mie_calc(vs, r_dust Q_sca, Q_ext, std_dust)
		
		! calculate opacities from extinction coeficinets
		kext = n_dust* pi * r_dust**2 * Q_ext
		
		! loop initalisation for Rosseland mean
		ros_mean_plank = 0.d0
		ros_mean       = 0.d0
		
		! calculate plank derivatives for integration
		call plank_derivative_function(temp, plank_der)
		
		! calculate Rosseland mean with trapezioid integration
		do j = 2, nr_wl
			if (nr_wl(j-1) > lam_low .and. nr_wl(j) < lam_up) then
				ros_mean       = ros_mean +                            &
								 ( kext*plank_der(j)     &
								  -kext*plank_der(j-1) ) &
								 / (m_wl(j) - m_wl(j-1))
				ros_mean_plank = ros_mean_plank +        &
								 ( plank_der(j)          &
								  -plank_der(j-1) )      &
								 / (m_wl(j) - m_wl(j-1))
			end do
		end do
		
		! assign Rossland mean opacity
		kc(i1,i2,i3) = ros_mean_plank / ros_mean
	
	
	end subroutine rosland_mean
	
	

    ! ==================================================================================
    ! Derivative of the Plank distribution evluated on the wavelength grid.
    subroutine plank_derivative_function(T, u)
        implicit none
        ! input parameter
        real(8), intent(in)       :: T ! Temperature
        ! output parameter
        real(8), dimension(m_nwl) :: u ! output function values
        ! internal parameter
        integer :: i
        
        do i = 1, m_nwl
            u(i) = exp(h*c/m_wl(i)/kb/T) / ((exp(h*c/m_wl(i)/kb/T) - 1)**2 * T**2)
            u(i) = (2*h**2*c**3 /m_wl(i)**6/kb) * u(i)
        end do
        
    end subroutine plank_derivative_function


end module rosland
