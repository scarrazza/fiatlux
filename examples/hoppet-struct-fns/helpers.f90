module helpers
  use convolution
  use pdf_tabulate
  use pdf_representation
  use types; use consts_dp
  use sub_defs_io
  implicit none


  integer, parameter, private :: maxlen = 1000

  interface write_lhapdf_with_photon
     module procedure write_lhapdf_with_photon_filename, write_lhapdf_with_photon_iunit
  end interface write_lhapdf_with_photon

  real(dp), parameter :: pdf4lhc_yvalues(100) = (/&
     & 2.07232659e+01_dp,  2.03511309e+01_dp,  1.99789958e+01_dp,  1.96068609e+01_dp,&
       1.92347259e+01_dp,  1.88625910e+01_dp,  1.84904561e+01_dp,  1.81183211e+01_dp,&
       1.77461861e+01_dp,  1.73740511e+01_dp,  1.70019162e+01_dp,  1.66297812e+01_dp,&
       1.62576463e+01_dp,  1.58855114e+01_dp,  1.55133764e+01_dp,  1.51412414e+01_dp,&
       1.47691064e+01_dp,  1.43969715e+01_dp,  1.40248365e+01_dp,  1.36527016e+01_dp,&
       1.32805665e+01_dp,  1.29084316e+01_dp,  1.25362966e+01_dp,  1.21641617e+01_dp,&
       1.17920267e+01_dp,  1.14198917e+01_dp,  1.10477568e+01_dp,  1.06756218e+01_dp,&
       1.03034868e+01_dp,  9.93135190e+00_dp,  9.55921695e+00_dp,  9.18708196e+00_dp,&
       8.81494701e+00_dp,  8.44281205e+00_dp,  8.07067703e+00_dp,  7.69854208e+00_dp,&
       7.32640709e+00_dp,  6.95427213e+00_dp,  6.58213724e+00_dp,  6.21000222e+00_dp,&
       5.83786721e+00_dp,  5.46573229e+00_dp,  5.09359729e+00_dp,  4.72146233e+00_dp,&
       4.34932737e+00_dp,  3.97719239e+00_dp,  3.60505748e+00_dp,  3.23292247e+00_dp,&
       2.86078752e+00_dp,  2.48865254e+00_dp,  2.30258509e+00_dp,  2.13396244e+00_dp,&
       1.98971280e+00_dp,  1.86367206e+00_dp,  1.75175412e+00_dp,  1.65111063e+00_dp,&
       1.55967641e+00_dp,  1.47590651e+00_dp,  1.39861483e+00_dp,  1.32687095e+00_dp,&
       1.25993146e+00_dp,  1.19719311e+00_dp,  1.13815956e+00_dp,  1.08241757e+00_dp,&
       1.02961938e+00_dp,  9.79469591e-01_dp,  9.31715154e-01_dp,  8.86137713e-01_dp,&
       8.42547269e-01_dp,  8.00777849e-01_dp,  7.60683385e-01_dp,  7.22134709e-01_dp,&
       6.85017099e-01_dp,  6.49227984e-01_dp,  6.14675596e-01_dp,  5.81277309e-01_dp,&
       5.48958511e-01_dp,  5.17651601e-01_dp,  4.87295133e-01_dp,  4.57833096e-01_dp,&
       4.29214285e-01_dp,  4.01391775e-01_dp,  3.74322449e-01_dp,  3.47966600e-01_dp,&
       3.22287582e-01_dp,  2.97251500e-01_dp,  2.72826945e-01_dp,  2.48984753e-01_dp,&
       2.25697798e-01_dp,  2.02940808e-01_dp,  1.80690268e-01_dp,  1.58923987e-01_dp,&
       1.37621402e-01_dp,  1.16763168e-01_dp,  9.63311258e-02_dp,  7.63082063e-02_dp,&
       5.66783467e-02_dp,  3.74264119e-02_dp,  1.85381241e-02_dp,  0.00000000e+00_dp /)
  
contains


  !----------------------------------------------------------------------
  function lhapdf_qmin(imem) result(qmin)
    integer, intent(in) :: imem
    real(dp) :: qmin, q2min
    call getq2min(imem,q2min)
    write(0,*) q2min
    qmin = sqrt(q2min)
  end function lhapdf_qmin

  !----------------------------------------------------------------------
  function lhapdf_qmax(imem) result(qmax)
    integer, intent(in) :: imem
    real(dp) :: qmax, q2max
    call getq2max(imem,q2max)
    write(0,*) q2max
    qmax = sqrt(q2max)
  end function lhapdf_qmax
  
  
  !----------------------------------------------------------------------
  function lhapdf_qmass(iflv) result(res)
    integer, intent(in) :: iflv
    real(dp)            :: res
    call getqmass(iflv,res)
  end function lhapdf_qmass
  
  !----------------------------------------------------------------------
  ! this assumes a 6-column format where column 1 is x and column 6 is
  ! xf_{gamma/p}(x)
  function read_photon_from_file(filename,grid) result(photon)
    character(len=*), intent(in) :: filename
    type(grid_def),   intent(in) :: grid
    !---------------------------------
    integer             :: iy, iread
    character(len=maxlen) :: line
    real(dp)            :: photon(0:grid%ny), xVals(0:grid%ny), numbers(6)
    real(dp), parameter :: tolerance = 1e-7_dp
    integer, save :: warn_n = 5
    real(dp) :: renorm1, renorm2
    
    iread = get_new_device()
    open(unit=iread,file=trim(filename),status='old')
    !iread = idev_open_opt("-read-photon",status='OLD')
    xVals = xValues(grid)
    write(0,*) 'reading photon from file ',trim(filename)
    do iy = 0, grid%ny
       ! skip comment lines
       do while (.true.)
          read(iread,'(a)') line
          if (line(1:1) /= '#') exit
       end do
       ! then get & check the numbers
       read(line,*) numbers
       if (abs(log(numbers(1)/xVals(iy))) > tolerance) then
          write(0,*) 'grid  x value = ', xVals(iy)
          write(0,*) 'input x value = ', numbers(1)
          call wae_error('read_photon', 'x value is inconsistent with grid')
       end if
       photon(iy) = numbers(6)
    end do
    close(iread)

    ! next we check for NaNs
    if (photon(0) /= photon (0)) photon(0) = zero
    do iy = 1, grid%ny
       if (isNaN(photon(iy))) then
          call wae_warn(warn_n,'read_photon_from_file: trying to correct NaN at iy,x = ',&
               &intval=iy, dbleval=xVals(iy))
          if (xVals(iy) < xVals(iy-1) .and. xVals(iy+1) < xVals(iy) &
               &.and. .not. (isNan(photon(iy+1)) .or. isNan(photon(iy-1)))) then
             ! interpolate between the two values (after deweighting by (1-x)^4)
             if (xVals(iy-1) == one) then
                renorm1 = zero
             else
                renorm1 = photon(iy-1)/(one-xVals(iy-1))**4
             end if
             renorm2 = photon(iy+1)/(one-xVals(iy+1))**4
             photon(iy) = half*(renorm1+renorm2)*(one-xVals(iy))**4
          else
             call wae_error('read_photon_from_file: failed to correct NaN at iy,x = ',&
                  &intval=iy,dbleval=xVals(iy))
          end if
       end if
       
    end do
    !write(0,*) photon(:)
    
  end function read_photon_from_file

  !----------------------------------------------------------------------
  function isNaN(val)
    real(dp), intent(in) :: val
    logical isNaN
    
    isNaN = (val /= val)
  end function isNan
  
  
  !----------------------------------------------------------------------
  ! looks through a file's header for the tag param//' =' and then
  ! reads in the value that follows it
  function string_param_from_header(filename,param) result(string)
    use sub_defs_io
    character(len=*), intent(in) :: filename, param
    character(len=maxlen)          :: string
    !--------------------------------------------------------------------
    character(len=maxlen) :: tag,line
    integer             :: len_tag, ios, ix, iunit

    tag     = trim(param)//' ='
    len_tag = len_trim(tag)
    
    iunit = get_new_device()
    open(unit=iunit,file=filename,status='old')
    do
       ! check we can read the line
       read(iunit,'(a)',iostat=ios) line
       if (ios /= 0) call wae_error('string_param_from_file', &
            &'problem reading line from'//trim(filename)//' (EOF?)')

       if (len_trim(line) == 0) then
          ! skip blank lines
          cycle
       else if (line(1:1) == '#') then
          ! process lines starting with a hash
          ix = index(line,trim(tag))
          if (ix == 0) cycle

          ! now get just the piece we want, removing anything after ", "
          string = adjustl(line(ix+len_tag:))
          ix = index(string,', ')
          if (ix /= 0) string = string(1:ix-1)
          exit
       else
          ! we have reached the end of the header without finding the info...
          call wae_error('string_param_from_file', 'could not find tag "'//trim(tag)//'"')
       end if
       
    end do

    close(iunit)
  end function string_param_from_header

  !----------------------------------------------------------------------
  function int_param_from_header(filename,param) result(res)
    character(len=*), intent(in) :: filename, param
    integer                      :: res
    !-----------------------------------
    character(len=maxlen) :: string
    integer :: ios
    string = string_param_from_header(filename,param)
    read(string,*,iostat=ios) res
    if (ios /= 0) call wae_error('int_param_from_header','problem converting "'//trim(string)//'" to int')
  end function int_param_from_header
  
  !----------------------------------------------------------------------
  function dble_param_from_header(filename,param) result(res)
    character(len=*), intent(in) :: filename, param
    real(dp)                     :: res
    !-----------------------------------
    character(len=maxlen) :: string
    integer :: ios
    string = string_param_from_header(filename,param)
    read(string,*,iostat=ios) res
    if (ios /= 0) call wae_error('dble_param_from_header','problem converting "'//trim(string)//'" to dble')
  end function dble_param_from_header
  
  !----------------------------------------------------------------------
  ! Fills a PDF from LHAPDF -- if the PDF has enough components it tries
  ! also to fill the photon component.
  subroutine fill_from_lhapdf(grid, this_pdf, Q)
    type(grid_def), intent(in) :: grid
    real(dp), intent(out) :: this_pdf(0:,-6:)
    real(dp), intent(in)  :: Q
    !---------
    real(dp) :: xv(0:grid%ny)
    integer  :: i

    xv = xValues(grid)
    this_pdf = zero
    do i = 0, grid%ny
       if (ubound(this_pdf,dim=2) >= 8) then
          call EvolvePDFPhoton(xv(i), Q, this_pdf(i,-6:6), this_pdf(i,8))
       else
          call EvolvePDF(xv(i), Q, this_pdf(i,-6:6))
       end if
       !write(iunit,*) xv(i), this_pdf(i,8)
    end do
  end subroutine fill_from_lhapdf


  !-----------------------------------------------------------------
  !! return zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function zeta_of_y(y, a) result(zeta)
    real(dp), intent(in) :: y, a
    real(dp)             :: zeta
    zeta = y + a*(one - exp(-y))
  end function zeta_of_y

  
  !-----------------------------------------------------------------
  !! return inversion of zeta = ln 1/x + a*(1-x) +b*(1-x^4) (x = exp(-y))
  function y_of_zetaext(zeta, a, b) result(y)
    real(dp), intent(in) :: zeta, a, b
    real(dp)             :: y, x, diff_from_zero, deriv
    integer             :: iter
    real(dp), parameter :: eps = 1e-12_dp
    integer,  parameter :: maxiter = 100

    ! starting condition (and soln if a = 0 and b=0)
    y = zeta 
    if (a /= zero .or. b /= zero) then
       do iter = 0, maxiter
          x = exp(-y);
          diff_from_zero = zeta - (y + a*(one-x) + b*(one-x**4));
          ! we have found good solution
          if (abs(diff_from_zero) < eps) exit
          deriv = -one  - a*x - four*b*x**4;
          y = y - diff_from_zero / deriv;
       end do
    end if
    
    if (iter > maxiter) write(0,*) "y_of_zeta reached maxiter"

  end function y_of_zetaext

  !-----------------------------------------------------------------
  !! return zeta = ln 1/x + a*(1-x) + b * (1-x)^4 (x = exp(-y))
  function zetaext_of_y(y, a, b) result(zeta)
    real(dp), intent(in) :: y, a, b
    real(dp)             :: zeta, x
    x = exp(-y)
    zeta = y + a*(one - x) + b * (one - x**4)
  end function zetaext_of_y
  
  
  !-----------------------------------------------------------------
  !! return inversion of zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function y_of_zeta(zeta, a) result(y)
    real(dp), intent(in) :: zeta, a
    real(dp)             :: y, x, diff_from_zero, deriv
    integer             :: iter
    real(dp), parameter :: eps = 1e-12_dp
    integer,  parameter :: maxiter = 100

    ! starting condition (and soln if a = 0)
    y = zeta 
    if (a /= zero) then
       do iter = 0, maxiter
          x = exp(-y);
          diff_from_zero = zeta - y - a*(one-x);
          ! we have found good solution
          if (abs(diff_from_zero) < eps) exit
          deriv = -one  - a*x;
          y = y - diff_from_zero / deriv;
       end do
    end if
    
    if (iter > maxiter) write(0,*) "y_of_zeta reached maxiter"

  end function y_of_zeta


  !----------------------------------------------------------------------
  subroutine write_moments(unit, grid, pdf,label)
    integer,          intent(in) :: unit
    type(grid_def),   intent(in) :: grid
    real(dp),         intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !-----------------
    real(dp) :: moments(ncompmin:ubound(pdf,dim=2))
    
    write(unit,"(a)", advance="no"), "# total momentum "//trim(label)//" (& components)"
    moments = TruncatedMoment(grid, pdf, one)
    write(unit,"(40f10.7)") sum(moments), moments

    write(unit,"(a)", advance="no"), "# total number "//trim(label)//" (& components)"
    moments(1:6) = TruncatedMoment(grid, pdf(:,1:6)-pdf(:,-1:-6:-1), zero)
    write(unit,"(40f11.7)") sum(moments(1:6)), moments(1:6)

  end subroutine write_moments
  

  !----------------------------------------------------------------------
  subroutine write_lhapdf_with_photon_filename(filename, table, pdf_type, iy_increment)
    character(len=*),  intent(in) :: filename
    type(pdf_table),   intent(in) :: table
    character(len=*),  intent(in) :: pdf_type
    integer, optional, intent(in) :: iy_increment
    integer :: iunit

    iunit = get_new_device()
    open(unit=iunit,file=trim(filename))
    call write_lhapdf_with_photon_iunit(iunit, table, pdf_type, iy_increment)
    close(iunit)
  end subroutine write_lhapdf_with_photon_filename
  
  !----------------------------------------------------------------------
  ! 
  subroutine write_lhapdf_with_photon_iunit(iunit, table, pdf_type, iy_increment)
    integer,           intent(in) :: iunit
    type(pdf_table),   intent(in) :: table
    character(len=*),  intent(in) :: pdf_type
    integer, optional, intent(in) :: iy_increment
    
    integer   :: flav_indices(14), flav_pdg_ids(14)
    real(dp)  :: flav_rescale(14)
    if (.not. table%nf_info_associated) then
       call wae_error('write_lhapdf_with_photon_iunit', 'currently need nf_info_associated')
    end if
    if (table%nfhi == 5) then
       flav_indices(1:12) = (/-5, -4 , -3, -2, -1,  0, 1, 2, 3, 4, 5,  8/)
       flav_pdg_ids(1:12) = (/-5, -4 , -3, -2, -1, 21, 1, 2, 3, 4, 5, 22/)
       flav_rescale(1:12) = one
       call write_lhapdf_dat_file(iunit, table, pdf_type, &
            &             flav_indices(1:12), flav_pdg_ids(1:12), flav_rescale(1:12), iy_increment)
    else if (table%nfhi == 6) then
       flav_indices(1:14) = (/-6, -5, -4 , -3, -2, -1,  0, 1, 2, 3, 4, 5, 6,  8/)
       flav_pdg_ids(1:14) = (/-6, -5, -4 , -3, -2, -1, 21, 1, 2, 3, 4, 5, 6, 22/)
       flav_rescale(1:14) = one
       call write_lhapdf_dat_file(iunit, table, pdf_type, &
            &             flav_indices(1:14), flav_pdg_ids(1:14), flav_rescale(1:14), iy_increment)
    end if
    
  end subroutine write_lhapdf_with_photon_iunit
  
  !----------------------------------------------------------------------
  ! A first attempt at outputting a table to an lhapdf dat file
  !
  ! - pdf_type can be one of: central, error
  !
  ! - flav_indices indicates which indices to use from our PDF array
  !
  ! - flav_pdf_ids indicates the corresponding PDG IDs
  !
  ! - flav_rescale specifies by how much each flavour should be
  !   rescaled (usually by one, unless you want to use an entry that
  !   has a sum of flavour and anti-flavour and so rescale by 50% to get
  !   just the flavour)
  !
  ! - iy_increment (default 1) indicates the increment to use in going
  !   up in y. The rule is that every time y goes up and iy has gone
  !   up by at least iy_increment, you output an x point.
  subroutine write_lhapdf_dat_file(iunit, table, pdf_type, &
                                   & flav_indices, flav_pdg_ids, flav_rescale,&
                                   & iy_increment)
    use sub_defs_io
    use assertions
    integer,           intent(in) :: iunit
    character(len=*),  intent(in) :: pdf_type
    type(pdf_table),   intent(in) :: table
    integer,           intent(in) :: flav_indices(:), flav_pdg_ids(:)
    real(dp),          intent(in) :: flav_rescale(:)
    integer, optional, intent(in) :: iy_increment
    !------------------------------------
    real(dp) :: flavs(lbound(table%tab,2):ubound(table%tab,2))
    real(dp) :: xVals(0:table%grid%ny), xVals_orig(0:table%grid%ny), last_x, val, yVal
    real(dp) :: zeta, zetamax, dzeta
    ! these parameters appear to give OK accuracy, <=10^{-4} for x<0.5
    ! and <=10^{-3} for x<0.9, at least for Q=100
    real(dp), parameter :: dzeta_def = 0.30_dp, zeta_a = 10.0_dp, zeta_b = 5.0_dp
    integer  :: iyVals(0:table%grid%ny), nx_max, iy, ix, iQ, iseg, iQlo, iQhi, ipdg, iflv
    type(pdfseginfo), pointer :: seginfos(:)
    integer  :: iy_inc, last_iy, n_since_last_spacing_change

    iy_inc = default_or_opt(-1, iy_increment)
    
    ! First work out what x values to use.
    ! In this simple version, we simply take the grid points (eliminating duplications)
    xVals_orig = xValues(table%grid)
    nx_max = -1
    last_x = two
    last_iy = -iy_inc
    n_since_last_spacing_change = 0
    write(iunit,'(a,i4)') '# iy_increment = ', iy_inc
    if (iy_inc > 0) then
       do iy = 0, table%grid%ny
          if (xVals_orig(iy) > last_x * 0.999999999_dp) then
             n_since_last_spacing_change = 0
          else
             n_since_last_spacing_change = n_since_last_spacing_change + 1
             ! skip some points, except just after a transition to a looser
             ! spacing
             if ( iy >= (last_iy + iy_inc) .or. &
                  & n_since_last_spacing_change < iy_inc**2) then
                last_x  = xVals_orig(iy)
                last_iy = iy
                nx_max = nx_max + 1
                xVals (nx_max) = xVals_orig(iy)
                iyVals(nx_max) = iy
             end if
          end if
       end do
    else if (iy_inc == 0) then
       ! use the PDF4LHC15 points; note that we need to reverse their order
       do iy = size(pdf4lhc_yvalues), 1, -1
          yVal = pdf4lhc_yvalues(iy)
          if (yVal < table%grid%ymax) then
             nx_max = nx_max + 1
             xVals(nx_max) = exp(-yVal)
          end if
       end do
       ! finish off with the last grid point
       nx_max = nx_max + 1
       xVals(nx_max) = exp(-table%grid%ymax)
    else if (iy_inc == -1) then
       write(iunit,'(a,f10.5,a,f10.5,a,f10.5)') &
            & '# dzeta_def = ', dzeta_def,&
            & ', zeta_a = ', zeta_a,&
            & ', zeta_b = ', zeta_b

       ! use our own smoothly changing spacing
       zetamax = zetaext_of_y(table%grid%ymax, zeta_a, zeta_b)
       nx_max = ceiling(zetamax/dzeta_def)
       dzeta = zetamax / nx_max
       do ix = 0, nx_max
          xVals(ix) = exp(-y_of_zetaext(dzeta*ix, zeta_a, zeta_b))
       end do
    else
       call wae_error('write_lhapdf_dat_file','value of iy_inc could not be handled; it was',&
            &intval = iy_inc)
    end if
    
    write(0,*) 'nx = ', nx_max+1, ', nQ = ', size(table%Q_vals)


    ! the official header
    write(iunit,'(a,a)') 'PdfType: ',trim(pdf_type)
    write(iunit,'(a)'  ) 'Format: lhagrid1'
    write(iunit,'(a)'  ) '---'

    ! handle different seginfo scenarios
    if (table%nf_info_associated) then
       ! if the table has seginfo, then just point to it
       seginfos => table%seginfo
    else
       ! otherwise create a fictitious seginfo with the info we need...
       allocate(seginfos(1:1))
       seginfos(1)%ilnlnQ_lo = lbound(table%Q_vals,1)
       seginfos(1)%ilnlnQ_hi = lbound(table%Q_vals,1)
    end if

    ! now loop over the actual or "invented" seginfos
    do iseg = lbound(seginfos,1), ubound(seginfos,1)
       ! first we output info about x structure, Q structure and flavours
       write(iunit,'(4000es14.7)') xVals(nx_max:0:-1)
       iQlo = seginfos(iseg)%ilnlnQ_lo
       iQhi = seginfos(iseg)%ilnlnQ_hi
       write(iunit,'(4000es14.7)') table%Q_vals(iQlo:iQhi)
       write(iunit,'(100i4)')      flav_pdg_ids

       ! then write out the PDF itself
       ! Note that LHAPDF wants _increasing_ x values
       do ix = nx_max, 0, -1
          iy = iyVals(ix)
          do iQ = iQlo, iQhi
             flavs = table%tab(:,:,iQ) .atx. (xVals(ix).with.table%grid)
             do ipdg = 1, ubound(flav_pdg_ids,1)
                val = flavs(flav_indices(ipdg)) * flav_rescale(ipdg)
                if (val == zero) then
                   write(iunit,'(a)',advance='no') '  0'
                else
                   write(iunit,'(es15.7)',advance='no') val
                end if
             end do
             write(iunit,'(a)') ''
          end do
       end do

       ! and finish with a yaml subdocument separator
       write(iunit,'(a)') '---'
    end do


    ! cleaning
    if (.not. table%nf_info_associated) deallocate(seginfos)
    close(iunit)
  end subroutine write_lhapdf_dat_file
  
  
end module helpers
