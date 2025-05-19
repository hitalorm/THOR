
!
!
!      *********************
!      ** penEasy Imaging **
!      *********************
!
!
!             ** DISCLAIMER **
!
!    This software and documentation (the "Software") were developed at the Food and
!    Drug Administration (FDA) by employees of the Federal Government in the course
!    of their official duties. Pursuant to Title 17, Section 105 of the United States
!    Code, this work is not subject to copyright protection and is in the public
!    domain. Permission is hereby granted, free of charge, to any person obtaining a
!    copy of the Software, to deal in the Software without restriction, including
!    without limitation the rights to use, copy, modify, merge, publish, distribute,
!    sublicense, or sell copies of the Software or derivatives, and to permit persons
!    to whom the Software is furnished to do so. FDA assumes no responsibility
!    whatsoever for use by other parties of the Software, its source code,
!    documentation or compiled executables, and makes no guarantees, expressed or
!    implied, about its quality, reliability, or any other characteristic. Further,
!    use of this code in no way implies endorsement by the FDA or confers any
!    advantage in regulatory decisions.  Although this software can be redistributed
!    and/or modified freely, we ask that any derivative works bear some notice that
!    they are derived from it, and any modified versions bear some notice that they
!    have been modified.
!
!
!                        @file    sourceRectangularBeam_2010-09-02.f
!                        @author  Andreu Badal (Andreu.Badal-Soler@fda.hhs.gov)
!                        @date    2010/09/02
!
!


!*******************************************************************
!*                   SOURCE RECTANGULAR BEAM                       *
!*                                                                 *
!* Short description:                                              *
!*                                                                 *
!*   Generation of primary particle states for radiation transport *
!*   calculations with PENELOPE.                                   *
!*                                                                 *
!*   A rectangular beam (also called pyramidal or fan beam) is     *
!*   emitted from a point focal spot. The particle directions are  *
!*   uniformly distribution on the unit sphere, and collimated to  *
!*   a rectangular field. The particle energy is sampled following *
!*   the input energy spectrum.                                    *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routines STORES,RAND                                       *
!*   from PENVOX.F                                                 *
!*   -> common /PARTVOX/                                           *
!*   -> routine LOCATEX                                            *
!*   -> routine STEPX                                              *
!*   from other penEasy libraries:                                 *
!*   -> routines GETLINE,TALLY                                     *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005,2006                                                     *
!*                                                                 *
!* -> Source created by A. Badal adapting the penEasy source       *
!*    "sourceBoxIsotropicGaussSpectrum.f" by J. Sempau and the     *
!*    pyramidal source model from PENMAIN.F (from penelope 2008)   *
!*    by F. Salvat.                                                *
!*******************************************************************


      subroutine RECTBEAMsource
!*******************************************************************
!*    Input:                                                       *
!*                                                                 *
!*    Output:                                                      *
!*      through /track/ and sec stack                              *
!*******************************************************************
	  use TRACK_mod
	  implicit none
	  
	  real*8 uold,vold,ivx,ivy,ivz,sdx,sdy,sdz
	  
      
      !integer*4 kpar,ibody,mat,ilb, ipol
      !real*8 e,x,y,z,u,v,w,wght
      !common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 xvox,yvox,zvox,absvox
      
      
      logical warned,srcpoint,active
      integer parsrc,matsrc,nspc, MAX_ANGLE_BINS, MAX_ENERGY_BINS,
     &        alias, max_angle_bin_used
      parameter (MAX_ANGLE_BINS=128, MAX_ENERGY_BINS=512)
      real*8 shots,usrc,vsrc,wsrc,xsrc,ysrc,zsrc,espc,cutoff,
     &       rot_beam, cos_theta_low, D_cos_theta, phi_low, D_phi,    
     &       max_height_at_y1cm     
      
      common /partvox/ uold,vold,ivx,ivy,ivz,sdx,sdy,sdz,
     &                 xvox,yvox,zvox,absvox
	  common /srcrect/ rot_beam(9),espc(MAX_ENERGY_BINS), 
     &           cutoff(MAX_ANGLE_BINS*MAX_ENERGY_BINS),
     &           shots,cos_theta_low, D_cos_theta, phi_low, D_phi,
     &           max_height_at_y1cm, usrc,vsrc,wsrc,xsrc,ysrc,zsrc,
     &           alias(MAX_ANGLE_BINS*MAX_ENERGY_BINS), parsrc,
     &           matsrc,nspc,max_angle_bin_used, warned,srcpoint,active
      integer seeki_walker, bin_angle_energy, bin_ang, bin_e
      integer*4 ncross
      real*8 PI,DOSPI,HALFPI, randno,rand,INFTY,dsef, u1,v1,
     &       phi_sampled, sin_theta_sampled
      parameter (PI=3.1415926535897932d0, DOSPI=2.0d0*PI, 
     &           HALFPI=PI*0.5d0, INFTY=1.0d30)
      external rand

      if (.not.active) return

      kpar = parsrc  ! Particle type


      ! -- Sample the fan beam angle and the particle energy from the input angle/energy spectra:
   55 continue   ! ["do" -> "while()" loop:]   

        ! -- Sample energy and angle bin at the same time from the input spectra, using Walker's aliasing methode:   !!Walker!!
        bin_angle_energy = seeki_walker(cutoff,alias,rand(2.0d1),
     &                                 (nspc*max_angle_bin_used))
        
        bin_ang = int( (bin_angle_energy-1)/nspc )       ! Angle bin, starting at 0.
        bin_e   = bin_angle_energy - bin_ang*nspc
            
        
        
        ! -- Fan beam sampling: PENMAIN ALGORITHM, corrected to emit a non-deformed square field.
        w = cos_theta_low + rand(1.0d0)*D_cos_theta      !!FANBEAM!!
        sin_theta_sampled = sqrt(1.0d0-w*w)


        !!DeBuG!! Original uniform fan beam:  phi_sampled = phi_low + rand(2.0d0)*D_phi
        
        ! - Sample the phi angle assuming a symmetric bow-tie filter.
        !   The first digit of the sampled number will select if we are on the "right" (<90deg) or the "left" (>90deg) of the fan beam.
        !   The remaining decimal part (0,1) will sample the actual angle with a uniform dist in the bin.
        !   Reject the sampled angle if it is lower than the minimum angle from the x axis (ie, larger than the maximum aperture from the y axis). 
        !   This may happen if the input aperture falls inside an angular bin defined in the spectrum.
        randno = 10.0d0*rand(2.2d1)
        if (randno<5.0d0) then   ! Emit to the "right"
          phi_sampled=(HALFPI-bin_ang*D_phi)-(randno-int(randno))*D_phi     ! Scale the lowest phi to the minimum angle at the sampled angular bin
          if (phi_sampled.lt.phi_low) goto 55          ! Reject the sampled angle if it is lower than the minimum angle from the x axis (ie, larger than the maximum aperture from the y axis). This may happen if the input aperture falls within an spectrum angular bin.
        else                     ! Emit to the "left"
          phi_sampled=(HALFPI+bin_ang*D_phi)+(randno-int(randno))*D_phi     ! The phi will be above 90 deg: assuming symmetric bow-tie filter!     !!DeBuG!!
          if (phi_sampled.gt.(PI-phi_low)) goto 55     ! Reject the sampled angle, larger than the input aperture
		  !goto 55
        end if
              
        u = sin_theta_sampled * cos(phi_sampled)
        v = sin_theta_sampled * sin(phi_sampled)   
      
	  
      if (abs(w/v).gt.max_height_at_y1cm) goto 55      !!DeBuG!!  Force square field for any phi rejecting directions that would go above a maximum precomputed height!!   !!DeBuG!! 
        
        
      ! Sample the energy uniformly within the sampled energy bin:
      e = espc(bin_e) + rand(3.0d1)*(espc(bin_e+1)-espc(bin_e))           !!DeBuG!! Sampling uniformly within the energy interval using a new rand bc we don't have the original cumul prob.
        
    
      if (rot_beam(1).gt.-1000.0d0) then
        ! - Initial beam not pointing to (0,1,0), apply rotation:
        u1 = u
        v1 = v
        u = rot_beam(1)*u1+rot_beam(2)*v1+rot_beam(3)*w
        v = rot_beam(4)*u1+rot_beam(5)*v1+rot_beam(6)*w
        w = rot_beam(7)*u1+rot_beam(8)*v1+rot_beam(9)*w
      endif  
	  
!	  if (v.lt.0) goto 55
    


      ! -Assign position (point focal spot):
      x = xsrc
      y = ysrc
      z = zsrc
	  


      
	  call locate
	  call locatevox
      if (mat.eq.0) then              ! The particle was born in vacuum
        call stepx(infty,dsef,ncross) ! Advance up to the object or infinity
        call tally(7,dsef)            ! Inform tallies about the maiden flight
      endif

      ! Kinetic energy (same code as box isotropic gaussian source):

      wght = 1.0d0                              ! Init statistical weight
      ilb(1) = 1                                ! Tag as source particle (i.e. 1st generation)
      ilb(2) = 0                                ! Clear other labels
      ilb(3) = 0
      ilb(4) = 0
      ilb(5) = 0                                ! Optional label (transferred to descendants)
	  
	  ipol = 0

      call stores(e,x,y,z,u,v,w,wght,kpar,ilb,ipol)  ! Push particle to stack
      call tally(0,e)                           ! Deposit its kinetic E
      end



      subroutine RECTBEAMinisrc(activated,emax)
!*******************************************************************
!*    Initializes the source.                                      *
!*                                                                 *
!*    Output:                                                      *
!*      activated -> TRUE if the source is active.                 *
!*      emax -> max source energy (eV)                             *
!*******************************************************************
      implicit none
      logical activated
      real*8 emax

      logical warned,srcpoint,active
      integer parsrc,matsrc,nspc, MAX_ANGLE_BINS, MAX_ENERGY_BINS,
     &        alias, max_angle_bin_used
      parameter (MAX_ANGLE_BINS=128, MAX_ENERGY_BINS=512)
      real*8 shots,usrc,vsrc,wsrc,xsrc,ysrc,zsrc,espc,cutoff,
     &       rot_beam, cos_theta_low, D_cos_theta, phi_low, D_phi,    
     &       max_height_at_y1cm
      common /srcrect/ rot_beam(9),espc(MAX_ENERGY_BINS), 
     &           cutoff(MAX_ANGLE_BINS*MAX_ENERGY_BINS),
     &           shots,cos_theta_low, D_cos_theta, phi_low, D_phi,
     &           max_height_at_y1cm, usrc,vsrc,wsrc,xsrc,ysrc,zsrc,
     &           alias(MAX_ANGLE_BINS*MAX_ENERGY_BINS), parsrc,
     &           matsrc,nspc,max_angle_bin_used, warned,srcpoint,active
     
      !! Common from sourceBoxIsotropicGaussSpectrum.f (required for image_model_3):
      logical warned0,srcpoint0,active_srcbig
      integer parsrc0,matsrc0,nspc0, dim
      parameter (dim=1000)      
      real*8 shots0,usrc0,vsrc0,wsrc0,cossrc0,espc0,pspc0,despc0,rot0
      real*8 xsrc0,ysrc0,zsrc0,dxsrc0,dysrc0,dzsrc0
      common /srcbig/rot0(3,3),espc0(dim),pspc0(dim),despc0(dim),shots0
     &           ,cossrc0,usrc0,vsrc0,wsrc0,xsrc0,ysrc0,zsrc0,dxsrc0,
     &           dysrc0,dzsrc0,parsrc0,matsrc0,nspc0,warned0,srcpoint0,
     &           active_srcbig
     
      character*80 buffer, fname
      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION SOURCE RECTANGULAR BEAM v.2010-09-02]')
      parameter (eos='[END OF RECTANGULAR BEAM SECTION]')
      
      integer i,j,error, num_angle_bins, finduf, in, position_ptr
      real*8 prob_angle_energy(MAX_ANGLE_BINS*MAX_ENERGY_BINS),
     &       prob_angles_tmp(MAX_ANGLE_BINS), total_prob_all_ang,
     &       total_prob_ang(MAX_ANGLE_BINS), phi_rot,alpha_rot,
     &       theta_aperture,phi_aperture,mc2,twomc2,PI,norm,ZERO,
     &       RAD2DEG,DEG2RAD, roll_rot, cp, sp, ca, sa, cr, sr
     
      parameter (PI=3.1415926535897932d0,DEG2RAD=PI/180.0d0)
      parameter (RAD2DEG=180.0d0/PI)      
      parameter (mc2=5.10998918d5,twomc2=2.0d0*mc2,ZERO=1.0d-9)
      

      write(*,'(a)') ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'RECTBEAMinisrc:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(1x,a3)') buffer
      if (adjustl(buffer(1:3)).eq.'ON') then
        active = .true.
        activated = active
      else if (buffer(1:3).eq.'OFF') then
        active = .false.
        activated = active
        write(*, '(a)')
     &    '>>>> Source rectangular Beam is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'RECTBEAMinisrc:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'RECTBEAMinisrc:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Type of particle:
      write(*,'(a)') 'Particle type:'
      read(*,*) parsrc
      write(*,'(1x,i1)') parsrc
      if (parsrc.ne.1.and.parsrc.ne.2.and.parsrc.ne.3) then
        write(*,'(a)') 'RECTBEAMinisrc:ERROR: invalid particle type'
        stop
      endif

      ! Energy:

      ! -- Read the source energy spectrum form an external file      !!BOW-TIE!!
      write(*,'(a)') 
      read(*,'(a80)') fname
      fname = adjustl(fname)
      fname = fname(1:scan(fname,' ')) ! Clip at 1st blank
      write(*,'(a,1x,a)')'Reading the angle/energy spectra'//
     &      ' input file: ', fname

      in = finduf()
      open(in,file=fname,status='old',iostat=error)
      if (error.ne.0) then
        write(*,'(a)')
     &      'Energy spectra source ERROR: cannot open the input file?'
        stop
      endif
          
      !! The spectrum data is read below, after reading the maximum fan angles.

            
      !!DeBuG!! Gaussian source disabled to speed up the execution: too many useless if's!

      if (parsrc.eq.3) emax = emax+twomc2  ! Allowance for e+ annihilation

      ! Position:
      write(*,'(a)') 'Focal spot coordinates (cm):'
      read(*,*) xsrc,ysrc,zsrc
      write(*,'(3(1x,1pe12.5))') xsrc,ysrc,zsrc

      ! Direction:
      write(*,'(a)')
     &'Direction vector for the rectangular beam central axis (u,v,w):'
      read(*,*) usrc,vsrc,wsrc
      write(*,'(3(1x,1pe12.5))') usrc,vsrc,wsrc

      norm = sqrt(usrc**2+vsrc**2+wsrc**2)
      if (norm.lt.ZERO) then
        write(*,'(a)') 'RECTBEAMinisrc:ERROR: null direction not valid.'
        stop
      else
        if (abs(norm - 1.0d0).gt.ZERO) then
          usrc = usrc/norm
          vsrc = vsrc/norm
          wsrc = wsrc/norm
          write(*,'(a)')
     &       'Normalized direction vector of the center of the beam:'
          write(*,'(3(1x,1pe12.5))') usrc,vsrc,wsrc
        endif
      endif

      write(*,'(a)') 'Polar and azimuthal aperture angles'//
     &               ' for the rectangular beam source (deg):'
      read(*,*) theta_aperture, phi_aperture, roll_rot            !  Read the fan beam aperture and the rectangular beam tilt (roll)
      write(*,'(2(1x,1pe12.5))') theta_aperture, phi_aperture
      if (theta_aperture.lt.0.0d0.or.theta_aperture.gt.180.0d0) then
        write(*,'(a)')
     &    'RECTBEAMinisrc:ERROR: polar aperture must be in [0,180].'
        stop
      endif
      if (phi_aperture.lt.0.0d0.or.phi_aperture.gt.360.0d0) then
        write(*,'(a)') 'RECTBEAMinisrc:ERROR: azimuthal aperture'//
     &             ' must be in [0,360].'
        stop
      endif
      write(*,'(a)') 'Roll angle (rectangular beam tilt):'      
      write(*,'(1x,1pe12.5)') roll_rot
      

! *** RECTANGULAR BEAM INITIALIZATION: aperture initially centered at (0,+1,0), ie, THETA_0=90, PHI_0=90
      cos_theta_low = cos((90.0d0 - 0.5d0*theta_aperture)*DEG2RAD)
      D_cos_theta   = -2.0d0*cos_theta_low     ! Theta aperture is symetric above and below 90 deg
      phi_low = (90.0d0 - 0.5d0*phi_aperture)*DEG2RAD
      
        !!BOWTIE!! The angle interval will be set by the angle/energy spectrum (see below).
        !!DeBuG!! Original:  D_phi   = phi_aperture*DEG2RAD
      
      max_height_at_y1cm = tan(0.5d0*theta_aperture*DEG2RAD)         !!DeBuG!! Force square field forcing a maximum Z at y=1cm;   !!DeBuG!! Is this correct???
  
  
            
      !! -- Rotation around X axis (alpha degrees, with alpha=0 on +Y) and Z axis (phi degrees, with phi=0 on +X):
      !!    The rotations move the vector (0,1,0) [==> alpha=0, phi=90] to the input direction vector (u,v,z).
      
      if (vsrc.gt.(1.0d0-ZERO) .and. abs(roll_rot).lt.1.0d-6) then
        ! Beam pointing to (0,1,0) and not tilted (roll=0): mark that rotation is not needed with a silly value
        rot_beam(1) = -1001.0d0
      else
        ! - Beam not pointing to (0,1,0), prepare rotation:
        alpha_rot = 0.5d0*PI - acos(wsrc)       ! Rotation about X: acos(wsrc)==theta, theta=90 for alpha=0, ie, +Y.
        if (abs(usrc).lt.1.0d-7 .and. abs(vsrc).lt.1.0d-7) then
          phi_rot = 0.0d0    ! Use when u and v are zero (atan(0/0) is not valid).
        else
          phi_rot = atan2(vsrc,usrc) - 0.5d0*PI   ! Rotation about Z:  initial phi = 90 (+Y).  [ATAN2(v,u) = TAN(v/u), with the angle in the correct quadrant.
        endif


CC !!DeBuG!! Original rotation for the fan beam (Rz*Rx) without roll rotation:

!  !! ** Rotation around X (alpha) and then around Z (phi): Rz*Rx
!         salpha = sin(alpha_rot)
!         calpha = cos(alpha_rot)
!         sphi   = sin(phi_rot)
!         cphi   = cos(phi_rot)!
!         rot_beam(1) = cphi
!         rot_beam(2) = -calpha*sphi
!         rot_beam(3) = salpha*sphi
!         rot_beam(4) = sphi
!         rot_beam(5) = calpha*cphi
!         rot_beam(6) = -salpha*cphi
!         rot_beam(7) = 0.0d0
!         rot_beam(8) = salpha
!         rot_beam(9) = calpha


CC --- Rotate the sampled vector 3 times:        !!DeBuG!! roll rotation
C      1) Roll rotation around the +Y axis (rho angle) to tilt the rectangular beam.
C      2) Rotation around X (alpha) to set the polar angle to the input direction.
C      3) Rotation around Z (phi) to set the azimuthal angle to the input direction.

        sa = sin(alpha_rot)
        ca = cos(alpha_rot)
        sp = sin(phi_rot)
        cp = cos(phi_rot)
        sr = sin(roll_rot*DEG2RAD)
        cr = cos(roll_rot*DEG2RAD)
        
        !! Rotation around Y (rho), around X (alpha) and then around Z (phi): Rz*Rx*Ry
        rot_beam(1) = cp*cr - sa*sp*sr
        rot_beam(2) =-ca*sp
        rot_beam(3) = cp*sr + sa*sp*cr
        rot_beam(4) = sp*cr + sa*cp*sr
        rot_beam(5) = ca*cp
        rot_beam(6) = sp*sr - sa*cp*cr
        rot_beam(7) =-ca*sr
        rot_beam(8) = sa
        rot_beam(9) = ca*cr
        
        
        write(*,'(a)')
     & 'Rotation from +Y to the input direction (around X and Z axis):'
        write(*,'(2(1x,1pe12.5))') alpha_rot*RAD2DEG, phi_rot*RAD2DEG
      endif
      
      
      
      ! -- Read the spectrum file:            
      do
        read(in,'(a80)') buffer          
        buffer = adjustl(buffer)
        if (buffer(1:1).ne.'#' .and. buffer.ne."") exit     ! skip all comments and blanks                    
      end do
      
            
      read(buffer,*) num_angle_bins, D_phi                       !!BOW-TIE!!  
      write(*,'(a)')"No. angular bins and bin width [degrees]:"
      write(*,'(1x,i4,1x,1pe12.5)') num_angle_bins, D_phi     

      if (num_angle_bins.gt.MAX_ANGLE_BINS) then
        write(*,'(a,i3)') 'RECTBEAMinisrc:ERROR: too many angular'//
     &   ' bins in the input spectrum. Increase MAX_ANGLE_BINS=',
     &   MAX_ANGLE_BINS
        stop       
      endif
      
      if ((num_angle_bins*D_phi).lt.0.0d0 .or.
     &    (num_angle_bins*D_phi).gt.360.0d0) then
        write(*,'(a,i4)') 'RECTBEAMinisrc:ERROR: the maximum '//
     &     'semi-aperture for the last angular bin in the spectrum '//
     &     'is outside [0,180]. Max angle = ', (num_angle_bins*D_phi)
        stop
      endif
      
      
        !!DeBuG!! OLD, using fix width for angles now       read(in,*) (angle_list(i), i=1,num_angle_bins)    ! Read the angle values with an implicit loop 
            

      max_angle_bin_used = 1+int((phi_aperture*0.5d0-0.00001d0)/D_phi)     !!BOW-TIE!!  Find the angular bin that contains the requested maximum semiaperture phi angle.
      
      if (max_angle_bin_used.gt.num_angle_bins) then        
        write(*,'(a)') 'RECTBEAMinisrc:ERROR: the input phi aperture'//
     &   ' of the fan beam is larger than the maximum angular span of'//
     &   ' angle/energy spectrum.'
        write(*,'(a)') ' Reduce the input fan beam or increse the'//
     &   ' angular width of the spectrum.'
        stop       
      endif

     
      D_phi = D_phi*DEG2RAD     ! Store the angle bin width in radians (used in the angle sampling).
         

      !! --Read the multiple spectra from the external file. 
      !!   The energy-angle bin will be sampled using Walker from the array of probabilities for each energy and angle bin.
      !!   The final energy and angle will be then sampled using a uniform distrib inside the bin. The sampling will be rejected if the angle is larger than the inoput fan beam.  !!BOW-TIE!!

      do i = 1, (MAX_ANGLE_BINS*MAX_ENERGY_BINS)         !!BOW-TIE!!  Initialize the probabilities to 0.
        prob_angle_energy(i) = 0.0d0
        cutoff(i)            = 0.0d0
        alias(i)             = 0
      end do


      do j = 1, MAX_ANGLE_BINS    ! Prepare auxiliar array to report the prob of each angle bin:
        total_prob_ang(j) = 0.0d0
      enddo
       

      nspc = 0

      do
        nspc = nspc+1
        read(in,*) espc(nspc), (prob_angles_tmp(i), i=1,num_angle_bins)    ! Read the lowest energy in the bin and the prob. for the different angles   !!BOW-TIE!!                
        
        if (espc(nspc).lt.0.0d0) then
          write(*,'(a)') 'RECTBEAMinisrc:ERROR: negative energy'
          stop
        else if (espc(nspc).lt.espc(max(nspc-1,1))) then
          write(*,'(a)') 'RECTBEAMinisrc:ERROR: decreasing energy'
          stop
        endif

        if (prob_angles_tmp(1).lt.0.0d0) exit    ! End of spectrum: negative prob input.
          
        if (nspc.ge.MAX_ENERGY_BINS) then
          write(*,'(a,i4)')'RECTBEAMinisrc:ERROR: too many bins'//
     &   ' in spectrum; enlarge MAX_ENERGY_BINS=',MAX_ENERGY_BINS
          stop
        endif


        ! - Copy the angular probabilities for the current energy bin to the angle/energy prob array.
        !   The array is organized giving first the energies for the first angle bin, then the energies for the second angle, etc.
        !   Skip the angular bins above the input aperture (prob will be zero for them). 
        do i = 1, max_angle_bin_used
          prob_angle_energy( (i-1)*MAX_ENERGY_BINS + nspc )
     &           = prob_angles_tmp(i)                                     !!BOW-TIE!!  !!Walker!!
        end do  
        
        do j = 1, max_angle_bin_used    ! Count total prob of each angle bin:
          total_prob_ang(j) = total_prob_ang(j) + prob_angles_tmp(j)
        enddo
        
      enddo    ! [End reading spectra]
      
      
      emax = espc(nspc)  ! Set max energy
      
      nspc = nspc - 1    ! Discount the last bin, which had negative prob just to mark the maximum energy.
      

      ! -- Now that we know how many energy and angle bins will be used, re-locate the input 
      !    probabilities to the beginning of the array ("nspc*max_angle_bin_used" bins used only):
      position_ptr = nspc + 1
      do i = 2, max_angle_bin_used    ! (The first angle is already in position)
        do j= 1, nspc
          prob_angle_energy(position_ptr) = 
     &               prob_angle_energy((i-1)*MAX_ENERGY_BINS + j)
          position_ptr = position_ptr + 1            
        end do  
      end do      
      do i = position_ptr, (MAX_ENERGY_BINS*MAX_ANGLE_BINS)
        prob_angle_energy(i) = 0.0d0     ! Fill the rest of the array with zeros (just to be safe).
      end do
        
                
      
      write(*,'(a)') 'No. of energy bins read:'
      write(*,'(2(1x,i0))') nspc
      if (nspc.lt.1) then
        write(*,'(a)')
     &    'RECTBEAMinisrc:ERROR: at least 1 bin must be defined'
        stop
      endif
      write(*,'(a)')'No. of angle bins read and maximum angle bin'//
     &              ' used for the input azimuthal aperture:'
      write(*,'(2(1x,i0))') num_angle_bins, max_angle_bin_used
      
      total_prob_all_ang = 0.0d0    ! Total prob for all angles:
      do j = 1, max_angle_bin_used
        total_prob_all_ang = total_prob_all_ang + total_prob_ang(j)
      enddo
      write(*,'(a)')'Relative emission probabability for each'//
     &                                         ' angular bin:'
      do j = 1, max_angle_bin_used
        write(*,'(a,i3,a,1pe10.3,a,1pe10.3,a,1pe14.8)')
     &                '     Bin ',j,':  ', D_phi*(j-1)*RAD2DEG,
     &                ' --', D_phi*j*RAD2DEG,' deg. --> Prob = ', 
     &                total_prob_ang(j)/total_prob_all_ang
      enddo
      write(*,'(a)')' '
        

          !   SUBROUTINE IRND0(W,F,K,N)     (subroutine from penelope.f)
          !
          !  Initialisation of Walker's aliasing algorithm for random sampling
          !  from discrete probability distributions:
          !
          !  Input arguments:
          !    N ........ number of **different** values of the random variable.
          !    W(1:N) ... corresponding point probabilities (not necessarily
          !               normalised to unity).
          !  Output arguments:
          !    F(1:N) ... cutoff values.
          !    K(1:N) ... alias values.

      ! - Init. Walker sampling of the angular/energy bin:
      call IRND0( prob_angle_energy, cutoff, alias, 
     &            (nspc*max_angle_bin_used) )         !!Walker!! Calling PENELOPE's function to init the Walker method

      write(*,'(a)')"Walker aliasing sampling algorithm"//
     &                " correctly initialized!"
      write(*,'(a)')' '
        
        

      ! Init performance vars:
      shots = 0.0d0
      warned = .false.

      !! -- Copy the rectangular beam source data in the standard source common block to allow image_tally (model 3) to use this data:
      xsrc0 = xsrc
      ysrc0 = ysrc
      zsrc0 = zsrc
      nspc0 = nspc
      do j = 1, nspc0
        espc0(j) = espc(j)   ! Copy the input energy spectrum
!         pspc0(j) = pspc(j)
!         despc0(j)= despc(j)
      enddo


      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'RECTBEAMinisrc:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ', eos
        write(*,'(a,a)') '  found instead:     ', buffer
        stop
      endif
      write(*,'(a)')
     &   '>>>> RECTANGULAR BEAM source initialization finished >>>>'
      end





      integer function seeki_walker(cutoff, alias, randno, n)
!********************************************************************
!*    Finds the interval (x(i),x(i+1)] containing the input value   *
!*    using Walker's aliasing method.                               *
!*                                                                  *
!*    Input:                                                        *
!*      cutoff(1..n) -> interval cutoff values for the Walker method*
!*      cutoff(1..n) -> alias for the upper part of each interval   *
!*      randno       -> point to be located                         *
!*      n            -> no. of data points                          *
!*    Output:                                                       *
!*      index i of the semiopen interval where randno lies          *
!*    Comments:                                                     *
!*      -> The cutoff and alias values have to be previously        *
!*         initialised calling the penelope subroutine IRND0.       *
!*                                                                  *
!*                                      AB, 2009-11-19              *
!********************************************************************
      implicit none

      integer n, int_part, alias(n)
      real*8 randno, cutoff(n), fraction_part, RN

      RN = randno*dble(n)+1.0D0       ! Find initial interval:
      int_part = int(RN)              !   -- Integer part
      fraction_part = RN-int_part     !   -- Fractional part

      if (fraction_part.lt.cutoff(int_part)) then    ! Check if we are in the aliased part
        seeki_walker = int_part
      else
        seeki_walker = alias(int_part)
      endif

      end
      

!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

