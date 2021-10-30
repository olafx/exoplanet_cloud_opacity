! ======================================================================================
! Module to sample size distributions. It defines a linear grid with dust sizes around
! the mean vlaue of the distribution. Each point will then be weighted according to a 
! given distribution. The set up is modular and allows to add your own distributions 
! at the bottom. Distributions currently implemented: Gaussian (dist_gaus), Lognormal
! (dist_logn).
!
! Date:    18.01.21
! Author:  Sven Kiefer
! Contact: kiefersv.mail@gmail.com
! ======================================================================================

module size_dist
    implicit none
    
    private
    real(8), parameter :: pi = 4.d0 * atan(1.d0)
    
    ! internal variables
    real(8) :: m_r_min   = 1.0d-5  ! minimal radii, needed to prevent zero devision
    real(8) :: m_stepper = 1.0d-5  ! step size in standard deviations
    real(8) :: m_acc     = 1.0d-25 ! stopping accuracy for delta integration values
    logical :: m_mute    = .True. ! mute all outputs of this  
    
    public :: distributer ! Calcuates dust size spacing and weighting
    public :: dist_gaus   ! Gaussian distribution
    public :: dist_logn   ! lognormal distribution
    public :: dist_setter ! sett internal variable (optional)
    
    contains
    
    ! ==================================================================================
    ! Setter for internal variables. This function does not have to be used as in most
    ! cases the default values work well.
    subroutine dist_setter(r_min, stepper, acc, mute)
        implicit none
        real(8), optional :: r_min, stepper, acc
        logical, optional :: mute
        
        !set optional parameters
        if (present(r_min))   m_r_min   = r_min   ! minimal dust radii
        if (present(stepper)) m_stepper = stepper ! step size in std
        if (present(acc))     m_acc     = acc     ! stopping accuracy
        if (present(mute))    m_mute    = mute    ! mute all prints
        
    end subroutine dist_setter
    
    
    ! ==================================================================================
    ! A function which claculates the necessary dust radii to consider while analysing
    ! a dust size distribution (incl. their weights). The function f needs to be
    ! normalised to 1 and takes the mean a as first input parameter and its 'standard
    ! deviation' b as second parameter.
    subroutine distributer(f, a, b, completness, resolution, r_dist, weights)
        implicit none
        ! input variables
        real(8), intent(in) :: a           ! mean of the function
        real(8), intent(in) :: b           ! standard deviation of the function (or eqv)
        real(8), intent(in) :: completness ! value between 0 and 1
        integer, intent(in) :: resolution  ! number of radii considered
        interface
            real(8) function f(x, a, b)            ! Function which describes the
                real(8), intent(in) :: x, a, b     !   size distribution of the dust
            end function f                         !   partricles.
        end interface
        ! output variables
        real(8), dimension(resolution) :: r_dist  ! dust size grid
        real(8), dimension(resolution) :: weights ! weights of the dust sizes
        !internal variables
        integer :: i
        real(8) :: r_low, r_hig, ratio, integ_value, delta_int, &
                   integ_low, integ_hig, sb
        
        ! input check
        if (completness > 1 .or. completness < 0) then
            print *, ''
            print *, '--------------------- ERROR ---------------------'
            print *, 'The variable completness has to be set to a value'
            print *, 'between 0 and 1. Please set it accordingly.'
            print *, 'Current value: ', completness
            stop     '---------------------------------------------'
        end if
        
        ! set up integration loop
        i = 0             ! itteration variable
        integ_value = 0   ! Integration value
        sb = m_stepper*b  ! step size of integration
        weights = 0.0d0   ! initalisation of weights
        
        ! Default case (if std is to low, return only mean with weight 1)
        if (b < m_acc) then
            r_dist = a
            weights(1) = 1.0d0
            
            ! print warning message if module is not muted
            if (.not. m_mute) &
                print *, 'WARNING: The standard deviation given is smaller', &
                         'than the allowed accuracy. Default to mean value', &
                         'only.'
                          
            ! stop subroutine
            return
        end if
        
        do while (integ_value <= completness)
            ! starting preparation
            i = i + 1
            
            ! calcualte integration values (Simpson rule)
            integ_hig = (f(a + sb*i,             a, b)   +  &
                         f(a + sb*(2*i-1)/2.0d0, a, b)*4 +  &
                         f(a + sb*(i-1),         a, b)   )* &
                         sb/6.0d0
            integ_low = (f(a - sb*i,             a, b)   +  &
                         f(a - sb*(2*i-1)/2.0d0, a, b)*4 +  &
                         f(a - sb*(i-1),         a, b)   )* &
                         sb/6.0d0
            
            ! addition to current integration value
            delta_int = integ_hig + integ_low
            
            ! breaking if stuck
            if (delta_int < m_acc) then
                print *, ''
                print *, '--------------------- ERROR ---------------------'
                print *, 'The integration steps were to low. Please set a'
                print *, 'larger stepper variable to prevent this.'
                print *, ''
                print *, 'stepper:             ', m_stepper
                print *, 'Inetgration step:    ', delta_int
                print *, 'Completness reached: ', integ_value
                stop     '---------------------------------------------'
            end if
            
            ! add newly calculated integrtation values to the total
            integ_value = integ_value + delta_int
        
        end do
        
        ! calculate the lower to upper boundry of the dust size space
        r_hig = a + sb*i
        r_low = a - sb*i
        if (r_low < m_r_min) r_low = m_r_min ! if to low, defualt to r_min
        
        ! calculate the dust size grid between the boundaries with given resolution
        do i = 1, resolution
            ratio = real((i-1),8)/real((resolution-1),8) ! r step size
            r_dist(i) = ratio*(r_hig - r_low) + r_low    ! r value grid
            weights(i) = f(r_dist(i), a,b)               ! r weight
        end do
        
        ! normalise weights if possible. If not, return only the default case (mean 
        ! value with weight 1).
        if (sum(weights) < m_acc) then
            r_dist     = a
            weights    = 0.0d0
            weights(1) = 1.0d0
            
            ! print warning message if module is not muted
            if (.not. m_mute) &
                print *, 'WARNING: The weights could not be properly &
                      &   assignet. Default to mean value only.'
        else
            ! normalisation
            weights =  weights / sum(weights)
        end if
        
    end subroutine distributer
    
    
    ! ==================================================================================
    ! Gaussain distribution
    real(8) function dist_gaus(x, a, b)
        implicit none
        real(8), intent(in) :: x, a, b
        
        dist_gaus = exp(-0.5d0*((x-a)/b)**2) / (b*sqrt(2.0d0*pi))
        
    end function dist_gaus
    
    
    ! ==================================================================================
    ! Log-normal distribution
    real(8) function dist_logn(x, a, b)
        implicit none
        real(8), intent(in) :: x, a, b
        real(8) :: mu
        
        if (x <  m_acc) then
			dist_logn = 0.d0
			return
		end if
        
        mu = log(a) - (b**2)/2.0d0
        dist_logn = exp(-0.5d0*((log(x)-mu)/b)**2)/(x*b*sqrt(2.0d0*pi))
        
    end function dist_logn
    
    
    ! ==================================================================================
    ! Here you can add further distributions. All distributions must take as input
    ! exactly the mean value a and the standard deviation b. You can also import and use
    ! numerical distributions. Nethertheless, make sure to give it a approxiamte value
    ! of the mean and an estimated value of a deviation b.

end module size_dist
