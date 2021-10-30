! ======================================================================================
! Module to calculate Mie-Theory. Supports Size distribution calculations and
! Holospheres. The Mie-Calcualtion is importet from Miex (Wolf et al. 2004). The
! Holosphere calculation is done using the work of 
!
! Requires: Miex (Wolf et al. 2004)
!           DMiLay (with the argument mute added)
!
! Date:    18.01.21
! Author:  Sven Kiefer
! Contact: kiefersv.mail@gmail.com
! ======================================================================================


module mie
    implicit none
    
    private
    
    ! --------------------------------- Default values ---------------------------------
    !           The Values below are reasonable default values for the mie 
    !           calculation. All can (and should) be set using mie_init.
    ! ----------------------------------------------------------------------------------
    
    ! default values for various parts of the mie calculation
    real(8) :: m_max_hol_frac  = 0.85d0 ! maximum fraction of vacuum in holo spheres
    real(8) :: m_completness   = 0.99d0 ! default completness for dust size dist
    integer :: m_std_res       = 20     ! default resolution for dust size dist
    character(len=4) :: m_func = 'logn' ! distribution function ('gaus' or 'logn')
    
    ! default values for qu-polate and net-polate
    real(8) :: m_std_grid_min  = -7.d0   ! minimum dust size dist to consider (log)
    
    ! defualt values for internal switches to determine the type of mie calculation
    logical :: m_do_holo       = .False. ! Holosphere switch
    logical :: m_do_size_dist  = .False. ! Size distribution switch
    logical :: m_mute          = .False. ! mute all non error prints
    
    
    ! ------------------------------- Internal settings --------------------------------
    !           Do not change vlaues below as they are necessary default values
    !           for the functions called during the mie calculation.
    ! ----------------------------------------------------------------------------------

    ! trigger to remember current state of the mie initalisation
    logical :: m_init         = .False. ! Initalisation trigger
    logical :: m_emt_init     = .False. ! initalisation trigger of eff. medium calc
    integer :: m_n_vol        = 0       ! numer of species (do not change this value)

    ! fixed parameters
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    integer, parameter :: n_dust = 16 ! number of dust species from drift    

    ! storage variables for use in various functions
    integer                                 :: d_nwl = -1 ! Number of wavelengths
    real(8), dimension(:), allocatable, public      :: d_wl ! Wavelengths for mie calculation
    character(len=7), dimension(n_dust)     :: d_dust_names ! name of the dust species
    real(8), dimension(:,:), allocatable    :: d_no, d_ko ! refractive indices
    complex(8), allocatable, dimension(:,:) :: d_minc ! m refrective indicies
    complex(8), allocatable, dimension(:,:) :: d_einc ! epsilon ref indicies

    
    
    ! ----------------------------------- Functions ------------------------------------
    !           Use these functions to calcualte mie Theory. 
    ! ----------------------------------------------------------------------------------
    
    ! Public functions to be used for mie calculation
    public :: mie_init   ! Function to intalise all necessary parameter for mie calc
    public :: mie_calc   ! Mie calculation
    public :: mie_getter ! Returns the initalisation state of this module
    
    
    ! Private functions for internal handling of data
    private :: init_wavelengths ! Read in wavelength data
    private :: init_nk          ! Read in refractive index data
    private :: calc_mie_holo    ! calculate mie including holospheres (uses DMiLay)
    private :: calc_emt         ! calculate effective medium theory
    private :: e2m, m2e         ! convertion between m, epsilon (m = sqrt(epsilon))
    private :: NR               ! Newton Raphson minimization
    private :: Bruggeman        ! Bruggeman formula used in Newton-Raphson minimization
    private :: gauss            ! Linear equation solver
    private :: eff_pullback     ! calculate next effective medium value in itteration
    
    
    contains
    
    ! ==================================================================================
    ! Mie calculation set up. Needs to be called at least onece before the mie calc
    ! can be run. Can be run without any input variables in which case mie calcuation
    ! will be initalised without Holospheres and without dust size distributions.
    subroutine mie_init (mute, wave_file, dnf, nkp, &
                         do_holo, max_hol_frac, n_vol, &
                         do_size_dist, func, std_res, std_grid_min, completness)
        implicit none
        ! input variables
        logical, optional, intent(in) :: do_holo        ! Switch to activate holospheres
        logical, optional, intent(in) :: do_size_dist   ! Switch for dust size dist
        logical, optional, intent(in) :: mute           ! mute all prints
        real(8), optional, intent(in) :: completness    ! Completness of size dist calc
        real(8), optional, intent(in) :: max_hol_frac   ! max hollosphere vacuum fraction
        real(8), optional, intent(in) :: std_grid_min   ! lowest dust size dist
        integer, optional, intent(in) :: std_res        ! dust dist grid resolution
        integer, optional, intent(in) :: n_vol          ! number of species intended
        character(len=4), optional, intent(in) :: func  ! name of the dust dist function
        character(*), optional, intent(in) :: wave_file ! Wavelength file
        character(*), optional, intent(in) :: dnf       ! dustname file
        character(*), optional, intent(in) :: nkp       ! nk path
        
        ! Turn prints on or off
        if (present(mute)) m_mute = mute
        
        ! set number of species if desired
        if (present(n_vol)) m_n_vol = n_vol
        
        ! set smallest dust size dist
        if (present(std_grid_min)) m_std_grid_min = std_grid_min
        
        ! Module initalisation (with trigger to only do it once)
        if (.not. m_init) then
            ! read in wavelength array
            if (present(wave_file)) then
                call init_wavelengths(wave_file) ! read in from wave file
            else
                call init_wavelengths()          ! read in from default
            end if
            
            ! read in refractive index data
            if(     present(dnf).and.     present(nkp)) call init_nk(dnf, nkp)
            if(.not.present(dnf).and.     present(nkp)) call init_nk(nk_path=nkp)
            if(     present(dnf).and..not.present(nkp)) call init_nk(dust_name_file=dnf)
            if(.not.present(dnf).and..not.present(nkp)) call init_nk()
            
            ! Set initalisation swith to True
            m_init = .True.
        end if
        
        ! check if valid distribution functions is selcted
        if (.not.any(func==(/ 'gaus', 'logn' /))) then
            print *, ''
            print *, '--------------------- ERROR ---------------------'
            print *, 'Selected Distribution function not one of gaus or'
            print *, 'logn. Please select one of them.'
            print *, 'Currently selected: ', func
            stop     '---------------------------------------------'
        end if
        
        ! Activate holospheres
        if (present(do_holo))      m_do_holo      = do_holo
        if (present(max_hol_frac)) m_max_hol_frac = max_hol_frac
        
        ! Warning if holo sphere has to use defualt parameters
        if (present(do_holo).and.(.not.present(max_hol_frac))) then
            if (.not.m_mute) then
                print *, 'WARNING: Not all hollowsphere parameters were initalised. ', &
                         'The following default values are used:'
                if (.not.present(max_hol_frac)) print *, 'Max hol frac:  ', max_hol_frac
            end if
        end if
        
        ! Activate size distributions
        if (present(do_size_dist)) m_do_size_dist = do_size_dist
        if (present(completness))  m_completness  = completness
        if (present(std_res))      m_std_res      = std_res
        if (present(func))         m_func         = func
        
        ! Warning if dust size has to use defualt parameters
        if (present(do_size_dist).and.(.not.present(completness).or. &
                                       .not.present(std_res)    .or. &
                                       .not.present(func)            )) then
            if (.not.m_mute) then
                print *, 'WARNING: Not all dust size parameters were initalised. The', &
                         'following default values are used:'
                if (.not.present(completness)) print *, 'Completness:  ', completness
                if (.not.present(std_res))     print *, 'Resolution:   ', std_res
                if (.not.present(func))        print *, 'Distribution: ', func
            end if
        end if
        
        ! check if all values with bounderies are within them
        If (m_completness < 0.d0 .or. m_completness > 1.d0) then
            print *, ''
            print *, '--------------------- ERROR ---------------------'
            print *, 'Completness needs to be between 0 and 1. Instead'
            print *, 'the current value is:'
            print *, m_completness
            stop     '---------------------------------------------'
        end if
        If (m_max_hol_frac < 0.d0 .or. m_max_hol_frac > 1.d0) then
            print *, ''
            print *, '--------------------- ERROR ---------------------'
            print *, 'The maximum hollowsphere fraction needs to be'
            print *, 'between 0 and 1. Instead the current value is:'
            print *, m_max_hol_frac
            stop     '---------------------------------------------'
        end if
        
    end subroutine mie_init
    
    
    ! ==================================================================================
    ! Calculate Mie Theory according to the description on the top of this file. To be
    ! run mie_init should be run first to ensure correct set up.
    subroutine mie_calc (Vs, r_dust, Q_sca, Q_ext, r_std)
        use mie_routines, only: shexqnn2             ! Import Mie calculation from miex
        use size_dist,    only: distributer          ! Import size distribution calc
        use size_dist,    only: dist_gaus, dist_logn ! size distributions
        implicit none
        
        ! ----------------------------------------------------------------------- Set up
        ! variables input
        real(8), dimension(:), intent(in)    :: Vs     ! Volume mixing ratios
        real(8), intent(in)                  :: r_dust ! dust size
        real(8), optional, intent(in)        :: r_std  ! deviation of dust size
        
        ! variables output
        real(8), dimension(d_nwl), intent(out) :: Q_sca ! scattering coeficient
        real(8), dimension(d_nwl), intent(out) :: Q_ext ! extinction coeficient
        
        ! variables internal
        integer :: ier, nang, i, j
        logical :: doSA
        real(8) :: iQext, iQsca, iQabs, iQbk, iQpr, ialbedo, ig, &
                   rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg, &
                   r_qsca, r_qext, Qabs, Qbk, Qpr, albedo, g, x_size
        real(8), dimension(m_std_res)         :: r_arr, weights
        real(8), dimension(d_nwl)             :: q_sca_sum, q_ext_sum
        complex(8), dimension(d_nwl)          :: eeff
        complex(8), dimension(:), allocatable :: SA1, SA2
        
        ! Initalisation
        doSa = .false.
        nang = 1
        q_sca_sum = 0
        q_ext_sum = 0
        
        ! check initalisation and variable assignment
        if (m_init .eqv. .False.) then
            if (.not. m_mute) then
                print *, 'WARNING: The mie calculation was not initalised. Now ', &
                         'calling mie_init with default values.'
            end if
            call mie_init()
        end if
        
        if ((m_do_size_dist .eqv. .True.) .and. (.not. present(r_std))) then
            print *, ''
            print *, '--------------------- ERROR ---------------------'
            print *, 'Mie calculation is set to include size dist. This'
            print *, 'requires to input the r_std parameter, which was'
            print *, 'not given in mie_calc.'
            stop     '---------------------------------------------'
        
        end if
        
        
        
        ! --------------------------------------------------------- Refractiv Index calc
        ! Use effective mediume theory according to Bruggemann (1934) to calculate the
        ! effective reflactive indices of a dust particle.
        call calc_emt(vs, size(vs), eeff)
        
        
        ! -------------------------------------------------------------- Mie calculation
        ! Calculate extinction and scattering efficiencies using Mie theory (1908)
        do i = 1, d_nwl
            
            ! -------------------- no holospheres / no dust size distribution
            if (.not. m_do_holo .and. .not. m_do_size_dist) then
                ! calculate mie thory
                x_size = 2 * pi * r_dust / d_wl(i) ! size parameter
                call shexqnn2(eeff(i), x_size, q_ext_sum(i), q_sca_sum(i), Qabs, &
                              Qbk, Qpr, albedo, g, ier, SA1, SA2, doSA, nang)
            end if
            
            
            ! -------------------- with holospheres / no dust size distribution
            if (m_do_holo .and. .not. m_do_size_dist) then
                ! precalculate mie theory
                x_size = 2 * pi * r_dust / d_wl(i) ! size parameter
                call shexqnn2(eeff(i), x_size, q_ext_sum(i), q_sca_sum(i), Qabs, &
                              Qbk, Qpr, albedo, g, ier, SA1, SA2, doSA, nang)
                
                ! prepare input varaibales
                iQext = q_ext_sum(i); iQsca = q_sca_sum(i); iQabs = Qabs; iQbk = Qbk
                iQpr = Qpr; ialbedo = albedo; ig = g
                
                ! call holosphere function (wirtten by dominic samra)
                call calc_mie_holo(eeff(i), d_wl(i), r_dust/10000.0d0, &
                                   iQext, iQsca, iQabs, iQbk, iQpr, ialbedo, ig, &
                                   rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg)
                
                ! assigne output variables
                q_ext_sum(i) = rQext; q_sca_sum(i) = rQsca
            end if
            
            
            ! -------------------- no Holosphers / with dust size distribution
            if (.not. m_do_holo .and. m_do_size_dist) then
                ! calculate correct distribution factors
                if (m_func == 'gaus') then
                    call distributer(dist_gaus, r_dust, r_std, m_completness, &
                                     m_std_res, r_arr, weights)
                else if (m_func == 'logn') then
                    call distributer(dist_logn, r_dust, r_std, m_completness, &
                                     m_std_res, r_arr, weights)
                end if
                    
                ! calculate mie scattering for all r_dust grid points
                do j = 1, m_std_res
                    x_size = 2 * pi * r_arr(j) / d_wl(i) ! size parameter
                    call shexqnn2(eeff(i), x_size, r_qext, r_qsca, Qabs, &
                                  Qbk, Qpr, albedo, g, ier, SA1, SA2, doSA, nang)
                    
                    ! weighted sum over all considered dust sizes
                    q_sca_sum(i) = q_sca_sum(i) + r_qsca * weights(j)
                    q_ext_sum(i) = q_ext_sum(i) + r_qext * weights(j)
                end do
            end if
            
            
            ! -------------------- with Holosphers / with dust size distribution
            if (m_do_holo .and. m_do_size_dist) then
                ! calculate correct distribution factors
                if (m_func == 'gaus') then
                    call distributer(dist_gaus, r_dust, r_std, m_completness, &
                                     m_std_res, r_arr, weights)
                else if (m_func == 'logn') then
                    call distributer(dist_logn, r_dust, r_std, m_completness, &
                                     m_std_res, r_arr, weights)
                else
                    stop 'No valid distribution was selceted. Check variable func.'
                end if
                    
                ! calculate mie scattering for all r_dust grid points
                do j = 1, m_std_res
                    x_size = 2 * pi * r_arr(j) / d_wl(i) ! size parameter
                    call shexqnn2(eeff(i), x_size, r_qext, r_qsca, Qabs, &
                                  Qbk, Qpr, albedo, g, ier, SA1, SA2, doSA, nang)
                          
                    ! prepare input varaibales
                    iQext = r_qext; iQsca = r_qsca; iQabs = Qabs; iQbk = Qbk
                    iQpr = Qpr; ialbedo = albedo; ig = g
                    
                    ! call holosphere function (wirtten by dominic samra)
                    call calc_mie_holo(eeff(i), d_wl(i), r_dust/10000.0d0, &
                                       iQext, iQsca, iQabs, iQbk, iQpr, ialbedo, ig, &
                                       rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg)
                    
                    ! assigne output variables
                    r_qext = rQext; r_qsca = rQsca
                    q_sca_sum(i) = q_sca_sum(i) + r_qsca * weights(j)
                    q_ext_sum(i) = q_ext_sum(i) + r_qext * weights(j)
                end do
            end if
        end do
        
        ! remnant to prevent data overwrite, no idea why but dont delet it
        q_sca = q_sca_sum
        q_ext = q_ext_sum

    end subroutine mie_calc
     
     
    ! ==================================================================================
    ! Returns the current initalisation state of the Mie calculation
    subroutine mie_getter (do_holo, do_size_dist, std_res, completness, func, nwl, &
                           nvol, wavelengths)
        implicit none
        ! variables
        logical, optional :: do_holo       ! True if holospheres are considered
        logical, optional :: do_size_dist  ! True if dust size dist are conisdered
        integer, optional :: std_res       ! Grid resolution of dust dist 
        integer, optional :: nwl           ! number of wavelengths
        integer, optional :: nvol          ! number of dust species
        real(8), optional :: completness   ! Completness factor for dust dist
        real(8), optional, allocatable, dimension(:) :: wavelengths ! wavelength array
        character(len=4), optional :: func ! Dust distribution function
        
        ! Return called variables
        if (present(do_holo))      do_holo      = m_do_holo
        if (present(do_size_dist)) do_size_dist = m_do_size_dist
        if (present(std_res))      std_res      = m_std_res
        if (present(nwl))          nwl          = d_nwl
        if (present(nvol))         nvol         = m_n_vol
        if (present(completness))  completness  = m_completness
        if (present(func))         func         = m_func
        if (present(wavelengths))  wavelengths  = d_wl
    end subroutine
   


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This code breaks because of illegal
    ! value assignment. TODO: Check it 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ==================================================================================
    ! Calculate mie theory including holospheres. This subroutine was writen by Dominic
    ! Samra and rewritten by Sven Kiefer. 
    subroutine calc_mie_holo(N_efft, wave, asize, &
                             iQext, iQsca, iQabs, iQbk, iQpr, ialbedo, ig, &
                             rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg)
        ! import libraries
        ! use dmilay_f95, only: DMiLay
        implicit none
        
        ! input variables
        complex(8), intent(in) :: N_efft ! effective refractive index of the dust
        real(8), intent(in)    :: wave   ! wavelength [cm]
        real(8), intent(in)    :: asize  ! dust particle size [??]
        real(8), intent(in)    :: iQext, iQsca, iQabs, &  ! Input mie calculation for
                                  iQbk, iQpr, ialbedo, ig !  backup in case holo fails
        ! output variables
        real(8),intent(out)    :: rQext, rQsca, rQabs, &  ! Output of mie holo calc or
                                  rQpr, rQbk, ralbedo, rg !  if fail, return input
        ! internal variables
        logical                :: ok
        integer                :: nFrac, ifrac, numang, maxang
        real(8)                :: wave_cm, Qext1, Qsca1, Qbk1, g, rcore, rshell, fnew, &
                                  fold, sext0, ssca0, sext1, ssca1,  gtot, stotext, &
                                  stotsca, g0, g1, wvno
        complex                :: rindsh, rindco
        real(8), dimension(3)  :: mu
        real(8),dimension(3,2) :: M1, M2, S21, D21
        real,dimension(320)    :: sext_loop,ssca_loop,g_loop


        ! Preparation step: Initalisation and conversion
        numang = 3               ! [TODO: what is this?]
        maxang = 10              ! [TODO: what is this?]
        mu = (/0.d0,0.5d0,1.d0/) ! [TODO: what is this?]
        
        wave_cm = (wave/10000.0d0) ! wavelength to centimeter for Dmilay
        wvno    = 2.d0*pi/wave_cm  ! convert to wave number as required by Dmilay
        rindsh  = cmplx(real(N_efft),-1.d0*dimag(N_efft)) ! effective mantle (img sign)
        rindco  = cmplx(1.d0,-0.d0)   ! vacuum, img sign flipped as dmilay is defined
        
        ! number of holosphere calculation accodring to fmax to reache suf. precision
        nFrac = 15
        if (m_max_hol_frac >= 0.9d0)   nFrac = 30 
        if (m_max_hol_frac >= 0.99d0)  nFrac = 60
        if (m_max_hol_frac >= 0.999d0) nFrac = 120
        
        ! initalisation of zero values for loop      
        sext0   = pi*asize**2*iQext 
        ssca0   = pi*asize**2*iQsca
        g0      = pi*asize**2*iQsca*ig
        stotext = 0.d0
        stotsca = 0.d0
        gtot    = 0.d0
        fold    = 0.d0 
        
        ! loop over different hollow sphere volume fractions
        do ifrac = 1, nFrac
            ! calculate current f value spaced on sqrt distribution
            fnew = m_max_hol_frac*(real(ifrac,8)/real(nFrac,8))**0.5d0
            
            ! calculate inner and outer sphere radius according to vacuum fraction
            rcore  = (asize**3/(1.d0/fnew-1.d0))**(1.d0/3.d0)
            rshell = (asize**3/(1.d0-fnew)     )**(1.d0/3.d0) 
            
            ! mie calculation including volume fraction
            call dmilay(rcore, rshell, wvno, rindsh, RINDCO, &
            & MU, NUMANG, Qext1, Qsca1, Qbk1, g, M1, M2, S21, D21, &
            & MAXANG, ok, m_mute)
                        
            if (.not. ok) then
                if (.not. m_mute) then
                    print *, 'WARNING: Hollowsphere calculation failed.', &
                             'Defaulting to input values.'
                end if
                
                ! set output values equal to input values
                rg = ig; rQext = iQext; rQsca = iQsca; rQabs = iQabs
                rQbk = iQbk; rQpr = iQpr; ralbedo = ialbedo
                
                ! stop subroutine
                return
            end if
                        
            ! check and convert g as dmilay uses differnt definition
            if (Qsca1 >  0.d0) g = g/Qsca1  
            if (Qsca1 <= 0.d0) g = 0.d0
            
            ! calculate current f value spaced on sqrt distribution
            fnew = m_max_hol_frac*(real(ifrac,8)/real(nFrac,8))**0.5d0
            
            ! assign current loop integration values
            sext1 = pi*rshell**2*Qext1
            ssca1 = pi*rshell**2*Qsca1
            g1    = pi*rshell**2*Qsca1*g
            
            ! trapezium integration
            stotext = stotext + 0.5d0*(sext0 + sext1)*(fnew-fold)
            stotsca = stotsca + 0.5d0*(ssca0 + ssca1)*(fnew-fold)
            gtot    = gtot    + 0.5d0*(g0    + g1   )*(fnew-fold)
            
            ! store old values for next integration step
            fold  = fnew
            sext0 = sext1
            ssca0 = ssca1
            g0    = g1
        enddo
                    
        
        ! prepare output if Mie-Calculation worked for all steps
        rQext = stotext / (pi*asize**2) / m_max_hol_frac
        rQsca = stotsca / (pi*asize**2) / m_max_hol_frac
        rQabs = rQext-rQsca  
        rQbk  = iQbk ! always equal to input, not calculated in this subroutine
        rQpr  = iQpr ! always equal to input, not calculated in this subroutine
        ralbedo = rQsca/rQext 
        rg    = gtot / stotsca 

    end subroutine calc_mie_holo


    ! ==================================================================================
    ! Calculates (n,k) constants of mixed material cloud particle using effective medium
    ! theory of Bruggemann or (if that fails) the LLL method. References: Helling et al.
    ! (2008), Lee et al. (2015b), Lee et al. (2016)
    subroutine calc_emt (vol, n_vol, m_eff)
        implicit none
        
        ! input variables
        integer, intent(in)                       :: n_vol ! number of volume species
        real(8), dimension(n_Vol), intent(in)     :: vol   ! volume fractions
        ! output variables
        complex(8), dimension(d_nwl), intent(out) :: m_eff ! effective ref. index
        ! internal variables
        logical                     :: errflag
        integer                     :: n, l
        complex(8)                  :: m_eff0, e_eff0
        complex(8), dimension(d_nwl) :: e_eff

        ! Inatlise arrays if not already correctly inatalised
        if (.not. m_emt_init .or. m_n_vol .ne. n_vol) then
            ! allocate arrays
            if (allocated(d_minc)) deallocate(d_minc)
            allocate(d_minc(n_vol, d_nwl))
            if (allocated(d_einc)) deallocate(d_einc)
            allocate(d_einc(n_vol, d_nwl))
            
            ! save current number of dust species
            m_n_vol = n_vol
            
            ! loop over species and wavelengths to fill arrays
            do l = 1, d_nwl
                do n = 1, m_n_vol
                    ! fill in non vacuum entries
                    if (n .le. n_dust) then
                        d_minc(n,l) = cmplx(d_no(n,l), d_ko(n,l), 8)
                        d_einc(n,l) = m2e(d_minc(n,l)) ! convert m to epsilon
                    
                    ! fill in vacuum entries
                    else       
                        d_minc(n,l) = cmplx(1, 0, 8)
                        d_einc(n,l) = m2e(d_minc(n,l)) ! convert m to epsilon
                    end if
                end do
            end do
            
            ! switch initalisation off
            m_emt_init = .True.
        end if
        
        ! loop over all wavelengths
        do l = 1, d_nwl
            ! loop initalistion
            m_eff0 = (0.d0, 0.d0) ! starting guess
            errflag = .False.     ! reset error flag trigger
            
            ! calculate starting condition with simple linear dependence
            do n = 1, n_vol
             m_eff0 = m_eff0 + vol(n) * d_minc(n,l)
            end do
        
            ! Call Newton-Raphson minimization
            call NR(d_minc(:,l), vol(:), m_eff0, errflag, m_eff(l))
            
            ! if fails (errflag) use LLL method
            if (errflag .eqv. .True.) then
                if (.not.m_mute) then
                    print *, 'WARNING: Bruggeman effective medium theory failed.', &
                             'Falling back to LLL method.'
                end if
                
                ! sumation initalisatoin
                e_eff0 = (0.d0,0.d0)
                
                !sum over all volume fractions
                do n = 1, n_vol
                    e_eff0 = e_eff0 + vol(n)*d_einc(n,l)**(1.d0/3.d0)
                end do
                
                ! set output values
                e_eff(l) = e_eff0**3
                m_eff(l) = e2m(e_eff(l))
            end if
            
        ! end wavelength loop
        end do

    end subroutine calc_emt


    ! ==================================================================================
    ! transform epsilon refrective indices to m refrective indices (m = sqrt(epsilon))
    pure complex(8) function e2m(e)
        implicit none
        ! input variables
        complex(8), intent(in) :: e
        ! internal variables
        real(8) :: ereal, eimag, n, k, sqrte2

        ! calcualte transformation
        ereal = real(e, 8)
        eimag = aimag(e)
        sqrte2 = sqrt(ereal*ereal + eimag*eimag)
        n = sqrt(0.5d0 * ( ereal + sqrte2))
        k = sqrt(0.5d0 * (-ereal + sqrte2))
        
        ! set output
        e2m = cmplx(n, k, 8)
    end function e2m


    ! ==================================================================================
    ! transform epsilon refrective indices to m refrective indices (m = sqrt(epsilon))
    pure complex(8) function m2e(m)
        implicit none
        ! input variables
        complex(8), intent(in) :: m
        ! internal variables
        real(8) :: ereal, eimag, n, k
        
        ! calcualte transformation
        n = real(m, 8)
        k = aimag(m)
        ereal = n*n - k*k
        eimag = 2.d0 * n * k
        
        ! set output
        m2e = cmplx(ereal, eimag, 8)
    end function m2e


    ! ==================================================================================
    ! Newton raphson minimizer for Bruggemann effective medium theory
    subroutine NR (M_inc, V_inc, M_eff0, unphysical, M_eff)
        implicit none
        
        ! input variables
        complex(8), intent(in)               :: M_eff0     ! inital ref. guess
        complex(8), dimension(:), intent(in) :: M_inc      ! refractie indces
        real(8),    dimension(:), intent(in) :: V_inc      ! volume fractions
        ! output variables
        logical, intent(out)                 :: unphysical ! errorflag
        complex(8), intent(out)              :: M_eff      ! effective ref. index
        ! internal variables
        integer                 :: itmax, it
        real(8)                 :: qual, acc, de1, de2
        real(8), dimension(2)   :: deff, eff_it, eff_new, FF, FF1, FF2, FF3, FF4
        real(8), dimension(2,2) :: DF
        
        ! Initalisation: set internal parameters
        unphysical = .False. ! start by assuming nothing broke yet
        itmax      = 30      ! maximum number of itterations
        acc        = 1.d-13  ! accuracy to be reached
        
        ! set starting value to estimation
        M_eff = M_eff0
        
        ! start Newton raphson itteration loop
        do it = 1, itmax
            ! real and img part of effective refractive index (eri)
            eff_it(1) = real (M_eff,8)
            eff_it(2) = aimag(M_eff)
            
            ! call bruggeman to calc quality of current eri (FF should be 0+i0)
            call Bruggeman(M_eff, M_inc, V_inc, FF)
            
            ! calculate absolute value of the complex FF
            qual = FF(1)*FF(1) + FF(2)*FF(2)
            
            ! breaking condition if minimisation succesfull
            if (abs(qual) < acc) return ! stop subroutine
            
            ! If not finished yet, find next point. Start by setting step size.
            de1 = eff_it(1)*1.0d-5
            de2 = eff_it(2)*1.0d-5
            
            ! Calcualate adjoint Bruggeman values to calculate deviation
            call Bruggeman(M_eff + cmplx(de1,  0.d0, 8), M_inc, V_inc, FF1)
            call Bruggeman(M_eff - cmplx(de1,  0.d0, 8), M_inc, V_inc, FF2)
            call Bruggeman(M_eff + cmplx(0.d0, de2,  8), M_inc, V_inc, FF3)
            call Bruggeman(M_eff - cmplx(0.d0, de2,  8), M_inc, V_inc, FF4)
            
            ! calcualte deviations in real and imag part
            DF(1,1) = (FF1(1)-FF2(1)) / (2.d0*de1)
            DF(1,2) = (FF3(1)-FF4(1)) / (2.d0*de2)
            DF(2,1) = (FF1(2)-FF2(2)) / (2.d0*de1)
            DF(2,2) = (FF3(2)-FF4(2)) / (2.d0*de2)
            
            ! calculate deviations
            call gauss(2, 2, DF, deff, FF)
            deff = -deff ! sign corrections
            
            ! callculate new eff itteration values using deff as offset
            call eff_pullback(eff_it, deff, FF, eff_new, unphysical, M_inc, V_inc)
            
            ! check breaking condition
            if (unphysical) then
                if (.not. m_mute) then
                    print *, 'WARNING: Bruggeman mixing faild due to unphysical', &
                             'effective refrective indice values. Using LLL instead'
                end if
                
                ! exit Bruggeman function
                return
            end if
            
            
            ! set new eri values for next itteration step
            M_eff = cmplx(eff_new(1), eff_new(2), 8)
        end do
        
        ! Bruggeman exceeded number of itteration steps, thus failed
        if (.not. m_mute) then
            print *, 'WARNING: Bruggeman mixing faild due to too many itteration', &
                     'steps needed. Using LLL instead.'
        end if
        
        ! set results as unphysical as minimisation failed
        unphysical = .True.
        
    end subroutine NR


    ! ==================================================================================
    ! Bruggeman formula to calcualte brugeman value and its derivative needed for the
    ! Newton raphson minimization.
    pure subroutine Bruggeman(M_eff,M_inc,V_inc,FF)
        implicit none
        
        ! input variables
        real(8),    dimension(:), intent(in) :: V_inc ! Volume fractions of dust species
        complex(8), dimension(:), intent(in) :: M_inc ! Refractive indices of species
        complex(8), intent(in)               :: M_eff ! eff. ref. index guess
        ! output variables
        real(8), dimension(2), intent(out)   :: FF    ! derivation value 
        ! internal variables
        integer    :: i
        complex(8) :: fun, mm2, mmi2
        
        ! [TODO: why do you square here]
        mm2 = M_eff**2
        
        ! set up loop
        fun = cmplx(0.d0, 0.d0, 8)
        
        ! sum over all dust species
        do i = 1, size(V_inc)
            mmi2 = M_inc(i)**2
            fun = fun + V_inc(i)*(mmi2 - mm2)/(mmi2 + 2.d0*mm2)
        end do
        
        ! set output parameters
        FF(1) = real(fun, 8)
        FF(2) = aimag(fun)

    end subroutine Bruggeman


    ! ==================================================================================
    ! Solve a linear system of equations ((a))*(x)=(b) for x with dimensions
    ! a=(n)-vector, x=(N)-vector, a=(NxN)-matrix.
    pure subroutine gauss (Nd, N, a_in, x, b_in)
        implicit none
        
        ! input variables
        integer, intent(in)                   :: Nd    ! Dimensions of the Matrix
        integer, intent(in)                   :: N     ! Dimensions of the vectors
        real(8), dimension(Nd,Nd), intent(in) :: a_in  ! input matrix
        real(8), dimension(Nd), intent(in)    :: b_in  ! input vector
        ! output variables
        real(8), dimension(Nd), intent(out)   :: x  ! output vector
        ! internal variables
        real(8)                   :: c, amax
        integer                   :: i, j, k, kmax
        real(8), dimension(Nd,Nd) :: a  ! input matrix
        real(8), dimension(Nd)    :: b  ! input vector
        
        ! define internal variables to prevent value change of input variabels
        a = a_in
        b = b_in
        
        ! switch rows to start with the largest entry
        do i = 1, N-1
            kmax = i
            amax = abs(a(i,i))
            
            ! find largest entry in row of matrix on top right side
            do k = i+1, N
                if (abs(a(k,i)) > amax) then
                    amax = abs(a(k,i))
                    kmax = k
                endif
            end do

            ! switch entries to have row start with largest entry
            if (kmax /= i) then
                do j = 1, N
                    c = a(i,j)
                    a(i,j) = a(kmax,j)
                    a(kmax,j) = c
                end do
                c = b(i)
                b(i) = b(kmax)
                b(kmax) = c
            end if

            ! bring matrix to triangle-shape
            do k = i+1, N
                c = a(k,i) / a(i,i)
                a(k,i) = 0.d0
                do j = i+1, N
                    a(k,j) = a(k,j) - c * a(i,j)
                end do
                b(k) = b(k) - c * b(i)
            end do
        end do
        
        ! solve for x
        do i = N, 1, -1
            c = 0.d0
            if (i < N) then
                do j = i+1, N
                    c = c + a(i,j) * x(j)
                end do
            end if
            x(i) = (b(i) - c) / a(i,i)
        end do

    end subroutine gauss


    ! ==================================================================================
    ! Calculate the next eff iteration value. If no valid value could be found, the
    ! flag unphysical is set to flase. (eri = effective refractiv index)
    pure subroutine eff_pullback (eff, deff, Fold, eff_new, unphysical, M_inc, V_inc)
        implicit none
        
        ! input parameters
        real(8), dimension(:),intent(in)    :: V_inc ! volume fraction
        complex(8), dimension(:),intent(in) :: M_inc ! refractive indicies of species  
        real(8), dimension(2), intent(in)   :: eff   ! old eri value
        real(8), dimension(2), intent(in)   :: deff  ! setp to next eff value
        real(8), dimension(2), intent(in)   :: Fold  ! current Brug. min. value
        ! output parameters
        logical, intent(out)                :: unphysical ! False, if stepper failed
        real(8), dimension(2), intent(out)  :: eff_new    ! new eri value
        ! internal parameters
        integer               :: itmax, it
        real(8)               :: fac, qnew, qold
        complex(8)            :: m_eff_new
        real(8), dimension(2) :: Fnew
        
        ! setting internal variables
        itmax      = 20     ! maximum number of iterations
        fac        = 1.d0   ! fraction of itteration steps (lowerd if failed)
        unphysical = .False. ! start by assuming it did not fail
        
        ! calculate qulity of old value, new value must be better than it
        qold = Fold(1)*Fold(1) + Fold(2)*Fold(2)
        
        ! loop over diverent fraction of deff step to find good new eff value
        do it = 1, itmax
            ! calculate new eff value
            eff_new = eff + fac*deff
            
            ! already prepare next itteration step as 70% of the last step size
            fac = fac*0.7d0
            
            ! check if real and img part are physical (larger than 0)
            ! If true, skip to next itteration step
            if ((eff_new(1) < 0.d0).and.(eff_new(2) < 0.d0)) then
                unphysical = .True.
                cycle
            end if
            
            ! calculate new Brug minimisation value 
            m_eff_new = cmplx(eff_new(1), eff_new(2), 8)
            call Bruggeman(m_eff_new, M_inc, V_inc, Fnew)
            qnew = Fnew(1)*Fnew(1) + Fnew(2)*Fnew(2)
            
            ! confirm that new eff value is physical
            unphysical = .False.
            
            ! break if new value is better than old one
            if (qnew < qold) return
            
        end do

    end subroutine eff_pullback


    ! ==================================================================================
    ! Read in wavelenghts from file (from Dominic Samra, rewritten by Sven Kiefer)
    subroutine init_wavelengths(wave_file)
        implicit none
        ! input variables
        character(*), optional, intent(in) :: wave_file ! wave info
        ! internal variables
        integer :: i
        character(:), allocatable :: wave_file_wk ! wave info
        
        ! Initalisation: set default arguments
        if (.not. present(wave_file)) then
            if (.not. m_mute) then
                print *, 'WARNING: No wavelength file path has been given. searching', &
                         'for wavelength information at wavelengths.txt'
            end if
            wave_file_wk = 'wavelengths.txt'
        else
            wave_file_wk = wave_file
        end if
        
        ! open wavelength file
        open(14, file=wave_file_wk, form='formatted', action='read')
        
        ! read in number of wavelngths (not yet implemented)
        read(14,*) d_nwl
        
        ! allocated wavelength array (if necessary, delet old allocation)
        if (allocated(d_wl)) deallocate(d_wl)
        allocate(d_wl(d_nwl))
        
        ! wavelength read in loop
        do i = 1,d_nwl
            read(14,*) d_wl(i)
        end do
        
        ! close file and finish read in
        close(14)
        
    end subroutine init_wavelengths


    ! ==================================================================================
    ! Read in nk data.
    subroutine init_nk(dust_name_file, nk_path)
        implicit none
        
        ! input variables
        character(*), optional, intent(in) :: dust_name_file ! dust names
        character(*), optional, intent(in) :: nk_path        ! nk data
        ! internal variables
        logical            :: exists, conducting
        integer            :: i, k, l, n, l1, n_lines 
        real(8)            :: fac
        character(len=200) :: zeile, fname, info, nk_path_wk, dust_name_file_wk
        character(len=10)  :: cond
        real(8), allocatable, dimension(:) :: wl_work, n_work, k_work
        
        ! Initalisation: set default arguments
        if (.not. present(dust_name_file)) then
            if (.not. m_mute) then
                print *, 'WARNING: No dust name file path has been given. Searching', &
                         'for it at nk_origin/aa_dust_names.txt'
            end if
            dust_name_file_wk = 'nk_origin/aa_dust_names.txt'
        else
            dust_name_file_wk = dust_name_file
        end if
        
        if (.not. present(nk_path)) then
            if (.not. m_mute) then
                print *, 'WARNING: No path to nk files has been given. Searching', &
                         'for it in the folder nk_origin.'
            end if
            nk_path_wk = 'nk_origin/'
        else
            nk_path_wk = nk_path
        end if
        
        ! allocated arrays
        if (allocated(d_no)) deallocate(d_no)
        allocate(d_no(n_dust, d_nwl))
        if (allocated(d_ko)) deallocate(d_ko)
        allocate(d_ko(n_dust, d_nwl))
        
        ! open dust names file
        open(12, file=trim(dust_name_file_wk), status='old')
        
        ! check if correct file was loaded
        read(12,*) i
        if (n_dust .ne. i) then
            print *, ''
            print *, '--------------------- ERROR ---------------------'
            print *, 'Dustchem file does not contain the right number'
            print *, 'of dust species. Please use the correct file.'
            print *, 'Expected dust species: ', n_dust
            print *, 'Given dust species:    ', i
            stop     '---------------------------------------------'
        end if
        
        ! loop over dust names and read them in
        do i = 1, n_dust
            read(12,*) cond
            k = index(cond,'[s]')
            d_dust_names(i) = trim(cond(:k-1)) 
        end do
        
        ! close file
        close(12)
        
        ! loop over all dust species to read in nk data
        do n = 1, n_dust
            ! define path to nk files and check if the file exists
            fname = trim(nk_path_wk)//trim(d_dust_names(n))//'[s]_ori.dat'
            inquire(file=fname,exist=exists)
            
            ! default case if nk data was not found
            if (.not. exists) then
                if (.not. m_mute) then
                    print *, 'WARNING: Missing (n,k)-data for ', &
                              trim(d_dust_names(n)), ' using vaccum data instead.'
                end if
                d_no(n,:) = 1.0
                d_ko(n,:) = 0.0 !eHERE vacuum values
                
                ! skip to the next loop step
                cycle
            end if
            
            ! open nk data file
            open(13, file=fname, form='formatted', action='read')
            
            ! read in nk data from file in whole lines
            read(13,*) n_lines, conducting
            read(13,*) ; read(13,'(A200)') info; read(13,*) ; read(13,*)
            
            ! allocated working arrays to store temporary data
            allocate(wl_work(n_lines), n_work(n_lines), k_work(n_lines))
            
            ! split lines into individual elements
            do l = 1, n_lines
                read(13,*) wl_work(l),n_work(l),k_work(l)
                ! enforce 0 boundary
                n_work(l) = max(0.d0, n_work(l))
                k_work(l) = max(0.d0, k_work(l))
            enddo
            
            ! close file
            close(13)
            
            ! loop over wavelengths to read and/or extrapolate nk data
            do l = 1, d_nwl
                ! If required wavelength is less than availible data - keep constant
                if (d_wl(l) < wl_work(1)) then
                    d_no(n,l) = n_work(1)
                    d_ko(n,l) = k_work(1)
                
                ! If required wavelength is greater than availible data - extrapolate
                else if (d_wl(l) > wl_work(n_lines)) then
                    ! Non conducting: n is constant - k is linear decreasing
                    if (conducting .eqv. .False.) then
                        d_no(n,l) = n_work(n_lines)
                        d_ko(n,l) = k_work(n_lines)*wl_work(n_lines)/d_wl(l)
                        
                    ! Conducting: n and k are log-log extrapolated
                    else if (conducting .eqv. .True.) then
                        ! find index were wl is 0.7 of max value
                        do l1 =  n_lines, 1, -1
                            if (wl_work(l1) < 0.7d0*wl_work(n_lines)) then
                                exit                                         
                            endif                                            
                        enddo
                        
                        ! log log interpolation
                        fac = log(d_wl(l)    /wl_work(n_lines)) / &
                              log(wl_work(l1)/wl_work(n_lines))
                        d_no(n,l) = exp(log(n_work(n_lines)) + &
                                        fac*log(n_work(l1)/n_work(n_lines)))
                        d_ko(n,l) = exp(log(k_work(n_lines)) + &
                                        fac*log(k_work(l1)/k_work(n_lines)))
                    endif
                
                ! Data is availible in the wavelength range - log-log interpolation
                else
                    ! Loop across work arrays untill straddle point, then interpolate
                    do l1 = 1, n_lines - 1
                        ! serach for straddle point
                        if (d_wl(l) >= wl_work(l1) .and. d_wl(l) <= wl_work(l1+1)) then
                            ! log log interpolation
                            fac = log(d_wl(l)      /wl_work(l1)) / &
                                  log(wl_work(l1+1)/wl_work(l1))
                                  
                            d_no(n,l) = exp(log(n_work(l1)) + &
                                            fac*log(n_work(l1+1)/n_work(l1)))
                            
                            ! default case, if k values are out of bounds
                            if (k_work(l1) <= 0.0d0 .or. k_work(l1+1) <= 0.d0) then
                                d_ko(n,l) = 0.d0
                            else
                                d_ko(n,l) = exp(log(k_work(l1)) + &
                                                fac*log(k_work(l1+1)/k_work(l1)))
                            endif
                            
                            ! break after interpoloation found
                            exit
                        endif
                    enddo
                    
                ! end interploation selection
                end if
            ! end wavelength loop
            end do
            
            ! print warinings if data was extrapolated
            if (.not. m_mute) then
                if (d_wl(1) < wl_work(1) .or. d_wl(d_nwl) > wl_work(n_lines)) then
                    print *, 'WARNING: extrapolated missing (n, k) data for ', &
                              d_dust_names(n), ' (conducting=', conducting, ')'
                end if
            end if
            
            ! deallocate working arrays
            deallocate(wl_work, n_work, k_work)
            
        ! end dust loop
        end do

    end subroutine init_nk

    
end module mie
