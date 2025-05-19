!*******************************************************************
!*******************************************************************
!*                          TALLY                                  *
!*               PIXELATED IMAGING DETECTOR                        *
!*                                                                 *
!* Short description:                                              *
!*   Tallies the pixelated image corresponding to the radiation    *
!*   flux arriving at a detector.                                  *
!* Rev:                                                            *
!*   2015-02-06                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module IMGmod
!*******************************************************************
!*    Vars for the IMG tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active,isLab,isEinteg,isPCM,isTrap   !!! logical isBiased enables or disables the Charge Sharing simulation option
      integer detmode,repformat,dispersionmode  !! dispersion Mode: with or without charge sharing
      integer*8 bodydet,matdet,nx,ny,nebin,ijmax,ijkmax,trap,nlast
	  integer*8 ccc, isaSe,isBiased
      real*8 lab2det(2,3),r0(3),dx,dy,idx,idy, EField   !!! EField: Bias Electric Field (V/cm)
      real*8 fano,eh,peso,e_tau,h_tau,thick,z_top,z_low
      real*8 unclimit,ethr,emin,emax,ebin,iebin, totalenergy
      real*8, parameter ::
     &       fwhm2sig2=1.0/(8.0*log(2.0d0)),
     &       npixmax=1.024e9    ! Max num pixels=32k squared (do not exceed int*4)
      real*8 :: memIMG=0.0      ! Memory used by arrays (bytes)
      real*8, allocatable ::
     &  wghtlast(:),e_counts(:),e_counts2(:),std(:)
	  real*8, allocatable :: image(:),image2(:),edeptmp(:),nlast_1(:)
      real*8, allocatable :: edeptmp_f(:)
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine IMGtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode: Identifies the state of the calling routine          *
!*      arg: energy loss (mode<0) or history no. (mode=1)          *
!*******************************************************************
      use IMGmod
      use TRACK_mod
      use ctrsimMod, only : nhist  ! Current history number
      use SECSTcommonMod, only : nsec,ilbs
      use PENELOPE_mod, only: eabs
      implicit none
      integer mode
      real*8 arg

      integer*8 maxpix
      parameter(maxpix=4500000) 
      real*8 edep,edep2,Dy_2,ipixarea,
     &       omega,theta,phi,rot_inv,xshift,yshift,zshift,
     &       rho, z_d, dist
      integer*8 imgmat, imgmodel, binout
      common /scoimg/ edep(maxpix),edep2(maxpix),
     &                Dy_2,ipixarea,
     &                omega,theta,phi,rho,rot_inv(9),xshift,yshift,
     &                zshift,imgmat,imgmodel,binout

     
      real*8 EPS, ONEPLUS, xt,yt,zt
      parameter(EPS=0.0000000001d0, ONEPLUS=1.0d0+EPS)

      integer*8 pix,i,j,ij,k,ijk
      real*8 rprima(3),xdet,ydet,zdet

      !integer out2


      
      !out2 = 300
      !open(out2,file='debug_z_value.dat')

      if (.not.active) return



      ! Score:
	  
	  if (isBiased.eq.1) then
		if (mat.eq.matdet)then
			if(mode.le.0.or.mode.eq.4) then   ! Energy has been deposited in detection material
				if (arg.eq.0.0) return                ! Nothing to deposit (e.g., elastic scattering)
      
			if (kpar.eq.2.and.isAse.eq.1.and.mat.eq.matdet) then
				if (mode.eq.4.or.ilb(1).eq.1.and.e.gt.0) then
					eh=(6+300/((abs(EField)/1e4)**0.9))*(0.38+4.8/((e/1e3)**0.5))
					!write(out2,'(2(1x,es12.5))') e, eh
				endif
			endif



		!Energy Integrator Detector
			if (isEinteg) then
				
				if (kpar.eq.1) then !selection of secondary electrons depositing energy
			    
				! the energy deposited on mode equal -1 
					if(totalenergy.lt.eh.and.mode.eq.-1)then
						totalenergy = totalenergy + ABS(arg)
					endif
					if(totalenergy.gt.eh)then
						z_d = dist(Efield, ccc, z, z_top, z_low)
			        !z_d = abs(thick - z)
						!write(out2,'(2(1x,es12.5))') z,z_d
						peso = 1
						call dispersion_model(totalenergy, x, y, z_d)
						totalenergy = 0.0
					endif

					if(mode.eq.-3.or.mode.eq.-4.or.mode.eq.-5) then
						z_d = dist(Efield, ccc, z, z_top, z_low)
						!write(out2,'(2(1x,es12.5))') z,z_d
						
						!z_d = abs(thick - z)
						peso = 1
						call dispersion_model(arg, x, y, z_d)
									
					! When a new particle is created its energy is subtracted from the
					! EHP counter
					else if(ilb(2).eq.1.and.mode.eq.-99) then
						!z_d = abs(thick - z)
						z_d = dist(Efield, ccc, z, z_top, z_low)
						peso = -1
						call dispersion_model(ABS(arg), x, y, z_d)
					endif
					
				endif 

			!Photon Counting Detector
			else if (isPCM) then
				if(mat.eq.matdet.and.kpar.eq.2) then
					!write(out2,'(1(1x,es12.5),i3,i3,i3)') e, mode, ilb(1),nlast	
					if(mode.eq.4.or.ilb(1).eq.0) then
						nlast = ilb(1)+1 !mode=4 went from one material to another
						!write(out2,'(1(1x,es12.5),i3,i3,i3)')e,mode, ilb(1),nlast
					endif

					if(ilb(2).eq.2.and.ilb(1).eq.nlast.and.ilb(4).eq.0) then
						if(mode.eq.-1.or.mode.eq.-2) then 
							nlast = ilb(1)+1
							!write(out2,'(1(1x,es12.5),i3,i3,i3)') e, mode, ilb(1),nlast	
						endif
					endif
				endif 


				! if(ilb(2).eq.2.and.mode.eq.-99)then
						! if (kpar.eq.2) then
						! ilb(5) = 1
						! do pix = 1, ijmax
							!write(out2,'(1x,es12.5)') edeptmp_f(pix)
		  
							! if (edeptmp_f(pix).ge.ethr) then
								! image(pix)  = image(pix) + 1           ! Add one count to this pixel   !Forcei
								! image2(pix) = image2(pix) + 1**2
			 
							! endif
					  
							! edeptmp_f(pix) = 0.0
						! enddo
			
					! else if(kpar.eq.1.and.ilb(1).eq.nlast) then
						! ilb(5) = 0

						! do pix = 1, ijmax
							! if (edeptmp_f(pix).ge.ethr) then
								! image(pix)  = image(pix) + 1           ! Add one count to this pixel   !Forcei
								! image2(pix) = image2(pix) + 1**2
							! endif
						! edeptmp_f(pix) = 0.0
						! enddo
					! endif
				! endif
			
				if (kpar.eq.1.and.ilb(1).gt.nlast) then !selection of secondary electrons depositing energy
					
					if(totalenergy.lt.eh.and.mode.eq.-1)then !stack electrons whose Eloss < W+-
						totalenergy = totalenergy + ABS(arg)
					endif

					if(totalenergy.gt.eh)then
						z_d = dist(Efield, ccc, z, z_top, z_low)
						write(out2,'(2(1x,es12.5))') z,z_d
						!z_d = abs(z_top - z)
						peso = 1
						call dispersion_model(totalenergy, x, y, z_d)
						totalenergy = 0.0
					endif

					if(mode.eq.-3.or.mode.eq.-4.or.mode.eq.-5) then
						z_d = dist(Efield, ccc, z, z_top, z_low)
						write(out2,'(2(1x,es12.5))') z,z_d
						!z_d = abs(z_top - z)
						peso = 1
						call dispersion_model(ABS(arg), x, y, z_d)
									
					else if(ilb(2).eq.1.and.mode.eq.-99) then
						z_d = dist(Efield, ccc, z, z_top, z_low)
						!z_d = abs(z_top - z)
						peso = -1
						call dispersion_model(ABS(arg), x, y, z_d)
					endif
					
				endif  
			   
			endif
      endif

      !if(mode.eq.6.and.isEinteg)then
      !  write(out2,'(1x,es12.5)') e_counts(1)
      !  e_counts(1) = 0

      !endif

      if(mode.eq.6.and.isPCM)then
		 !write(out2,'(2(1x,es12.5))') eabs(1,1),eabs(2,1)
            do pix = 1, ijmax
            
			if (edeptmp(pix).ge.ethr) then
                image(pix)  = image(pix) + 1           ! Add one count to this pixel   !Forcei
    		    image2(pix) = image2(pix) + 1**2
            endif
            edeptmp(pix) = 0.0
			edeptmp_f(pix) = 0.0
			
		 enddo
		 nlast = 0
		 ilb(5) = 0
		 totalenergy = 0
      endif
	  endif
	  
	  else
		  if (mat.eq.matdet.and.mode.le.0) then   ! Energy has been deposited in detection material
			if (arg.eq.0.0) return                ! Nothing to deposit (e.g., elastic scattering)

			if (isLab) then                       ! Slab and Sdet basis vectors coincide, no rotation needed
			  xdet = x-r0(1)                      ! Position vector w.r.t. pixels origin, only x&y coords are needed
			  ydet = y-r0(2)
			else                                  ! Both a translation and a rotatation are needed
			  rprima = (/x,y,z/)-r0               ! Translation
			  xdet = dot_product(lab2det(1,:),rprima) ! Rotation to obtain coords in Sdet frame
			  ydet = dot_product(lab2det(2,:),rprima)
			endif

			! Check if particle is inside pixels grid:
			i = xdet*idx+1.0                      ! Pixel index
			if (i.lt.1.or.i.gt.nx) return         ! Out of pixels grid
			j = ydet*idy+1.0                      ! Pixel index
			if (j.lt.1.or.j.gt.ny) return         ! Out of pixels grid
			ij = i+(j-1)*nx                       ! Mapping of i,j into a single index

			if (isEinteg) then                    ! E-integrating mode
			  if (nhist.gt.nlast_1(ij)) then               ! Apply "deferred bookkeeping": transfer temporary counter when new history
				!if (isEspread) edeptmp(ij)=espread(edeptmp(ij)) ! Gaussian energy spread
				image(ij)   = image(ij) +edeptmp(ij)     ! Transfer temp counter to permanent
				image2(ij)  = image2(ij)+edeptmp(ij)**2  ! ...and the variance too
				edeptmp(ij) = arg*wght                   ! Reset temp counter
				nlast_1(ij)   = nhist+0.5                  ! Reset last visit--add 1/2 to prevent round-off problems
			  else                                       ! Same history as previous event, continue temp scoring
				edeptmp(ij) = edeptmp(ij)+arg*wght       ! Score in temp counter
			  endif

			else if (nhist.gt.nlast_1(ij)) then     ! Photon modes; visit of a new history occurred
			  if (isPCM) then                     ! Photon counting mode (PCM)
				!if (isEspread) edeptmp(ij)=espread(edeptmp(ij))  ! Gaussian energy spread
				if (edeptmp(ij).gt.ethr) then                    ! Deposited energy exceeds threshold
				  ! AutoNote: Gaussian spread could be reduced to prob above the ETHR (Erf) for this case.
				  image(ij)  = image(ij) +wghtlast(ij)           ! Add one count to this pixel
				  image2(ij) = image2(ij)+wghtlast(ij)**2
				endif
			  endif
			  edeptmp(ij)  = arg                       ! Reset temp counter
			  nlast_1(ij)    = nhist+0.5                 ! Reset last visit--add 1/2 to prevent round-off problems
			  wghtlast(ij) = wght                      ! Reset weight to current value; shouldn't change for same hist

			else                                  ! Photon modes; same history as previous event, continue temp scoring
			  edeptmp(ij)  = edeptmp(ij)+arg           ! In photon modes WGHT is applied to the added count, not to deposited energy

			endif

		  endif
	  endif
      end


      subroutine IMGreport(mode,n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      mode: -1 if called at end-of-simulation, 0 for dump only,  *
!*             1 for tally report only, 2 for both                 *
!*      n: no. of histories simulated                              *
!*      cputim: elapsed CPU time                                   *
!*    Output:                                                      *
!*      uncdone: 2 if uncert reached, 1 if not defined, 0 else     *
!*               Not computed if mode=0                            *
!*******************************************************************
      use IMGmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      character*80 buffer
      integer out,finduf,errorcode
      integer*8 nchan,i,j,k,ij,ijk
      real*8 q,q2,sigma,eff,avesig,maxq,uncert,x,y,xmiddle,ymiddle,invn
      real*8 e,emiddle
      uncdone = 1
      if (.not.active) return
    

      ! Dump temp counters:
      do ij=1,ijmax
       if (nlast_1(ij).lt.0.5) cycle                     ! Skip empty pixels
        select case (detmode)
        case (1)                                      ! E-integrating mode
          image(ij)   = image(ij) +edeptmp(ij)          ! Transfer temp counter to permanent
          image2(ij)  = image2(ij)+edeptmp(ij)**2
        case (2)                                      ! Photon counting mode (PCM)
          if (edeptmp(ij).gt.ethr) then                 ! Deposited energy exceeds threshold
            image(ij)  = image(ij) +wghtlast(ij)        ! Add one count
            image2(ij) = image2(ij)+wghtlast(ij)**2     ! Variance
          endif
          wghtlast(ij) = 0.0                            ! Reset last weight
        end select

		! Reset all:
                                      ! Reset last visited to avoid recounting in next report
        edeptmp(ij)  = 0.0                              ! Reset counter
        nlast_1(ij)    = 0.0                              ! Reset last visited to avoid recounting in next report
		
      enddo
	  totalenergy = 0.0
        
		
      if (mode.ne.1) call IMGdump(1)  ! Sim dump file, notice that temp counter and nlast have been cleared before
      if (mode.eq.0) return           ! Dump only, no report

      ! Filename:
      select case (detmode)
        case (1)
          buffer = 'tallyPixelImageDetectEI'
        case (2)
          buffer = 'tallyPixelImageDetectPC'
      end select
      
      select case (repformat)
        case (1)
          buffer = trim(buffer)//'.dat'
        case (2)
          buffer = trim(buffer)//'-matrix.dat'
      end select

      out = finduf()
      open(out,file=buffer,iostat=errorcode)

      if (errorcode.ne.0) then
        write(*,*) ''
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'IMGreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        
        return
      endif



        write(out,'(a)')
     &  '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)')
     &  '# [SECTION REPORT PIXELATED IMAGING DETECTOR]'

        

        select case (detmode)
        case (1)
          write(out,'(a)')
     &     '# Energy integrating mode, units are electron-hole pairs.'
          write(out,'(a)') '# Electric Field (V/cm):'
          write(out,'(a,es12.5)') '# ',EField
        case (2)
          write(out,'(a)')
     &     '# Photon counting mode, units are electron-hole pairs.'
          write(out,'(a)') '#'
          write(out,'(a)') '# Energy deposition threshold (eV):'
          write(out,'(a,es12.5)') '# ',ethr
          write(out,'(a)') '# Electric Field (V/cm):'
          write(out,'(a,es12.5)') '# ',EField

        end select
        write(out,*) '#'
        write(out,'(a)')
     &    '# No. of pixels and pixel size (cm) in x & y directions:'
        write(out,'(a,1x,i0,1x,es12.5)') '# ',nx,dx
        write(out,'(a,1x,i0,1x,es12.5)') '# ',ny,dy
        write(out,'(a)')
     &    '# (x & y with respect to detector reference frame.)'
        write(out,'(a)') '#'

      select case (repformat)
      case (1)  ! Columnar
        write(out,'(a)')
     &    '# For plotting purposes, two values per '//
     &    'coordinate are given, namely, '//
     &    'the low end and the middle point of each interval.'
        write(out,'(a)') '#'
        write(out,'(a,$)')
     &    '# xPixelIndex : xLow(cm) : xMiddle(cm) : '//
     &      'yPixelIndex : yLow(cm) : yMiddle(cm) : '
        select case (detmode)
        case (1)
        write(out,'(a)') 'e-hCounts/hist : +-2sigma'
        
        case (2)
        write(out,'(a)') 'phton_Counts/hist) : +-2sigma'

        end select
      case (2)  ! Matrix format
        write(out,'(a)') '# Pixel x index:'
        write(out,'(a,$)') '#'
        do i=1,nx
          write(out,'(i12,a,$)') i,' '  ! Extra padding space to match data columns
        enddo
        write(out,*) ''                 ! Terminates non-advancing output
      end select                        ! Nothing to do for case3 binary

      ! Write data:
      invn   = 1.0/n
      nchan  = 0
      avesig = 0.0
      maxq   = 0.5*maxval(e_counts)                ! 1/2 of the max score

      do j=1,ny
        if (repformat.eq.1) then                ! Compute coords
          y = dy*(j-1)
          ymiddle = y+dy*0.5
        endif
        do i=1,nx
          if (repformat.eq.1) then
            x = dx*(i-1)
            xmiddle = x+dx*0.5
          endif
          ij = i+(j-1)*nx                       ! Mapping of i,j into a single index
            select case (detmode)
            case(1)
				if (isBiased.eq.1) then
					!q = e_counts(ij)*invn
					q = (edeptmp(ij) + edeptmp_f(ij))*invn
	!                sigma = (e_counts2(ij)*invn)*invn
					sigma = sqrt(max(e_counts(ij),0.0))
					sigma = sigma*invn
				else
					q = image(ij)*invn
					sigma = sqrt(max((image2(ij)-q**2)*invn,0.0))
				endif
            case(2)
                q = image(ij)*invn
                sigma = (image2(ij)*invn - q**2)*invn
                sigma = sqrt(max(sigma,0.0))
            end select
      
            
            select case (repformat)
            case (1)                            ! Columnar
              write(out,'(2(1x,i5,2(1x,es12.5)),1x,es12.5,1x,es12.5)')
     &          i,x,xmiddle,j,y,ymiddle,q,2.0*sigma
            case (2)                            ! Matrix
              write(out,'(1x,es12.5,$)') q
            case (3)                            ! Binary
              write(out) q,2.0*sigma
            end select
            if (e_counts(ij).gt.maxq) then         ! Evaluate average uncertainty for scores above 1/2 max score
              avesig = avesig+(sigma/q)**2
              nchan = nchan+1
            endif
        enddo
        write(out,*) ''
        
      enddo

      uncdone = 0
      if (nchan.gt.0) then
        uncert = 200.0*sqrt(avesig/nchan)
        if (uncert.lt.unclimit) uncdone = 2  ! Uncertainty reached
      else
        uncert = 0.0                    ! Uncertainty assumed not reached when score is nil
      endif

     
        write(out,'(a)') ' '
        write(out,'(a)') '# Performance report'
        write(out,'(a)') '#   Memory used (MB):'
        write(out,'(a,es12.5)') '#   ',memIMG*1.0e-6
        write(out,'(a)') '#   Random seeds:'
        write(out,'(a,i0)') '#   ',seed1
        write(out,'(a,i0)') '#   ',seed2
        write(out,'(a)') '#   No. of histories simulated [N]:'
        write(out,'(a,f18.0)') '#   ',n
        write(out,'(a)') '#   CPU time [t] (s):'
        write(out,'(a,es12.5)') '#   ',cputim
        if (cputim.gt.0.0) then
          write(out,'(a)') '#   Speed (histories/s):'
          write(out,'(a,es12.5)') '#   ',n/cputim
        endif
        write(out,'(a)')
     &  '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
        write(out,'(a,es12.5)') '#   ',uncert
        eff = n*uncert**2
        if (eff.gt.0.0) then
          write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
          write(out,'(a,es12.5)') '#   ',1.0/eff
        endif
        eff = cputim*uncert**2
        if (eff.gt.0.0) then
          write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
          write(out,'(a,es12.5)') '#   ',1.0/eff
        endif
        
      close(out)
      end


      subroutine IMGinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays.                          *
!*******************************************************************
      use IMGmod
      use dataTypesMod
      use PENELOPE_mod, only : eabs,maxmat
      use CSIMPHcommonMod, only : eabsh       ! Protons
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PIXELATED IMAGING DETECTOR v.2017-02-17]')
      parameter (eos='[END OF IMG SECTION]')
      character*80 buffer,upcase
      integer errorcode
      integer*8 ijwght
      real*8 pixlen
      real*8, parameter :: einf=1.1e9  ! Above max E in PENELOPE

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'SDDinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))      ! Adjust and capitalize
      buffer = buffer(1:scan(buffer,' '))   ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Pixelated Imaging Detector is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'IMGinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'IMGinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Detection material set to:'
      read(*,*) matdet
      write(*,'(1x,i0)') matdet
      if (matdet.le.0.or.matdet.gt.maxmat) then
        write(*,'(a)') 'IMGinitally:ERROR: invalid entry.'
        stop
      endif

      


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Added block for Charge Sharing Simulation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,'(a)') 'Include dispersion effects:'
      read(*,*) isBiased
      select case (isBiased)
		case (1)
			write(*,'(a)') '  Yes'
		case (0)
			write(*,'(a)') '  No'
	  end select
      
	  write(*,'(a)') 'Bias Electric Field (V/cm):'
      read(*,*) EField
      write(*,'(1x,es12.5)') EField
      
	  write(*,'(a)') 'Charge carrier collected:'
      read(*,*) ccc
      select case (ccc)
		case (0)
			write(*,'(a)') '  Collecting electrons'
		case (1)
			write(*,'(a)') '  Collecting holes'
		case (2)
			write(*,'(a)') '  Collecting electrons and holes'
	  end select
      
      write(*,'(a)') 'Fano Factor'
      read(*,*) fano
      write(*,'(1x,es12.5)') fano
      
      write(*,'(a)') 'Is W variable?'
      read(*,*) isaSe
      select case (isaSe)
		case (0)
			write(*,'(a)') '  Fixed Energy to create EHP'
		case (1)
			write(*,'(a)') '  Variable Energy to create EHP'
	  end select
      
	  
      write(*,'(a)') 'generated e-h pairs'
      read(*,*) eh
      write(*,'(1x,es12.5)') eh

      write(*,'(a)') 'Trapping:'
      read(*,*) trap
      select case (trap)
      case (0)
        write(*,'(a)') ' OFF'
        isTrap = .false.       ! Turn on charge trapping
      case (1)
        write(*,'(a)') ' ON'
        isTrap = .true.
      end select
      
	  write(*,'(a)') 'Electron mobility-lifetime (cm^2/V):'
      read(*,*) e_tau
      write(*,'(1x,es12.5)') e_tau
      write(*,'(a)') 'Hole mobility-lifetime (cm^2/V):'
      read(*,*) h_tau
      write(*,'(1x,es12.5)') h_tau
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      write(*,'(a)') 'Pixel size (cm) & No. pixels in x direction:'
      read(*,*) dx,nx
      write(*,'(1x,es12.5,1x,i0)') dx,nx
      if (dx.lt.0.0.or.nx.lt.0) then
        write(*,'(a)') 'IMGinitally:ERROR: invalid value.'
        stop
      endif
      if (dx.eq.0.0.and.nx.eq.0) then
        write(*,'(a)') 'IMGinitally:ERROR: '//
     &    'pixel size and no. pixels are both zero.'
        stop
      endif
      pixlen = dx*nx  ! Length of the pixels grid along x
      if (pixlen.ne.0.0) then
        write(*,'(a)') 'Pixels grid x-size is (cm):'
        write(*,'(1x,es12.5)') pixlen
      endif

      write(*,'(a)') 'Pixel size (cm) & No. pixels in y direction:'
      read(*,*) dy,ny
      write(*,'(1x,es12.5,1x,i0)') dy,ny
      if (dy.lt.0.0.or.ny.lt.0) then
        write(*,'(a)') 'IMGinitally:ERROR: invalid value.'
        stop
      endif
      if (dy.eq.0.0.and.ny.eq.0) then
        write(*,'(a)') 'IMGinitally:ERROR: '//
     &    'pixel size and no. pixels are both zero.'
        stop
      endif
      pixlen = dy*ny  ! Length of the pixels grid along y
      if (pixlen.ne.0.0) then
        write(*,'(a)') 'Pixels grid y-size is (cm):'
        write(*,'(1x,es12.5)') pixlen
      endif

      call setDetectFrame              ! Set detector reference frame (Sdet) and pixels size&num
      write(*,'(a)') 'Number of megapixels total:'
      write(*,'(1x,es12.5)') real(nx)*real(ny)*1.0e-6
      if (real(nx)*real(ny).gt.npixmax) then  ! Set a sensible max to ensure there will not be an int*4 overflow
        write(*,'(a,es12.5,a)') 'IMGinitally:ERROR: ',
     &    npixmax*1.0e-6,' max no. megapixels exceeded.'
        stop
      endif

      ! Detector mode:
      write(*,'(a)') 'Detector mode:'
      read(*,*) detmode
      select case (detmode)
      case (1)
        write(*,'(a)') ' Energy integrating'
        isEinteg = .true.       ! Switches help faster execution in tallying routine
        write(*,'(a)') 'E threshold ignored, not used in this mode.'
        read(*,'(a80)') buffer
        nebin = 1
      case (2)
        write(*,'(a)') ' Photon counting'
        isEinteg = .false.
        isPCM = .true.
        write(*,'(a)') 'Energy deposition threshold (eV):'
        read(*,*) ethr
        write(*,'(1x,es12.5)') ethr
        nebin = 1
      case default
        write(*,'(1x,i0)') detmode
        write(*,'(a)') 'IMGinitally:ERROR: invalid mode.'
        stop
      end select
      if (real(nx*ny)*nebin.gt.npixmax) then  ! Ensure there will not be an int*4 overflow
        write(*,'(a,es12.5,a)') 'IMGinitally:ERROR: ',
     &    npixmax*1.0e-6,' max no. (Mpixels x E bins) exceeded.'
        stop
      endif


      ! Format
      write(*,'(a)') 'Report format:'
      read(*,*) repformat
      select case (repformat)
      case (1)
        write(*,'(a)') ' Columnar'
      case (2)
        write(*,'(a)') ' Matrix'
        if (detmode.eq.3) then
          write(*,'(a)') 'IMGinitally:ERROR: '//
     &      'Matrix format is incompatible with '//
     &      'energy discriminating mode.'
          stop
        endif
      case (3)
        write(*,'(a)') ' Binary'
      case default
        write(*,'(1x,i0)') repformat
        write(*,'(a)') 'IMGinitally:ERROR: invalid format.'
        stop
      end select

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'IMGinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      ! Allocate scoring arrays:
      ijmax = nx*ny                          ! Size of 2D arrays
      select case (detmode)
      case (1)                               ! Energy integrating
        ijkmax = ijmax                       ! Note: DUMP uses this var even for 2D arrays
        ijwght = 1                           ! WGHTLAST defined for simplicity, unused for this mode
      case (2)                               ! Photon counting
        ijkmax = ijmax                       ! All arrays are 2D
        ijwght = ijmax                       ! WGHTLAST used
      case (3)                               ! Photon energy discriminating
        ijkmax = ijmax*nebin                 ! Size of 3D arrays
        ijwght = ijmax
      end select
      ! Memory usage:
      mem = sizeofReal8*real(ijkmax*2.0+ijmax*2.0+ijwght)
      write(*,'(a)') 'Memory used to allocate arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6
      allocate(image(ijkmax),image2(ijkmax),
     &         e_counts(ijmax),e_counts2(ijmax),wghtlast(ijwght),
     &         std(ijmax),edeptmp(ijmax),nlast_1(ijmax),
     &         edeptmp_f(ijmax),stat=errorcode)
	 
      if (errorcode.ne.0) then
        write(*,'(a)') 'IMGinitally:ERROR: not enough memory.'
        stop
      endif
      memIMG = mem  ! Store for later use by report

      ! Init arrays:
      image    = 0.0
      image2   = 0.0
      edeptmp  = 0.0
      edeptmp_f  = 0.0
      totalenergy = 0.0
      wghtlast = 0.0
      e_counts = 0.0
      e_counts2 = 0.0
      std = 0.0
	  nlast = 0
	  nlast_1 = 0
      peso = 1
      ! Sim restart file:
      call IMGdump(0)

      write(*,'(a)') '>>>> IMG tally initialization finished >>>>'



      end


      subroutine setDetectFrame
!*******************************************************************
!*    Identifies detector reference frame and prepares rot matrix. *
!*******************************************************************
      use IMGmod
      use PENGEOMcommonsMod
      use PENGEOM_mod, only : nxg,mater,nbody
      use dataTypesMod
      implicit none

      integer*8 i,isurf,sidepoint,surf,surfmax
      real*8 normx(3),scalarx,normy(3),scalary,normtop(3),scalartop
      real*8 imodu,detbox,pixgrid, normz(3),scalarz
      real*8, parameter :: eps=1.0d-13,oneplus=1.0d0+eps

	  
      ! Search last body containing the detector material
      bodydet = 0                     ! Init void
      do i=nbody,1,-1                 ! Reverse search, start from top
        if (mater(i).eq.matdet) then  ! This body is detector
          write(*,'(a)')
     &      'Detection material found in body/module no. '//
     &      '(PENGEOM internal numbering): '
          write(*,'(1x,i0)') i
          bodydet = i
          exit                        ! Only 1st ocurrence matters
        endif
      enddo
      if (bodydet.eq.0) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detection material not found.'
        stop
      endif
      surfmax = ksurf(bodydet,nxg)    ! Num. surfaces limiting this body

      ! Search for detector box faces:
      write(*,'(a)') 'Searching for detector box faces...'
      isurf = 0                     ! Init
      do i=1,surfmax                ! For all surfaces limiting BODYDET
        sidepoint = 2*kflag(bodydet,i)-3  ! Side pointer -1,+1,+3,etc
        if (sidepoint.gt.+1) cycle  ! Skip surfaces limiting inserted bodies
        isurf = i                   ! Found x-face
        exit
      enddo
      if (isurf.eq.0) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detector box faces not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)   ! Surface number
      write(*,'(a)') 'Found x-normal face; it is surface no. '//
     &  '(PENGEOM internal numbering):'
      write(*,'(1x,i0)') surf
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'x-normal face is not a plane.'
        stop
      endif
      normx = (/ax(surf),ay(surf),az(surf)/)  ! Surface normal
      imodu = sidepoint/sqrt(sum(normx**2))   ! Signed inverse modulus
      normx = normx*imodu                     ! Inwards-oriented unit normal
      scalarx = a0(surf)*imodu                ! Scalar term, normalized
      write(*,'(a)') 'Searching for the other 3 faces...'

      ! Seach for x-normal top face:
      isurf = isurf+1                         ! Next surface ID should be x-normal top face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x-normal top face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x-normal top face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'x-normal top face is not a plane.'
        stop
      endif
      normtop = (/ax(surf),ay(surf),az(surf)/)
      imodu = -sidepoint/sqrt(sum(normtop**2))
      normtop = normtop*imodu                 ! Outwards-oriented unit normal, should coincide with NORMX
      scalartop = a0(surf)*imodu              ! Scalar term, normalized

      ! Check if parallel and correctly oriented:
      if (abs(dot_product(normx,normtop)-1.0d0).gt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x-normal and x-normal top faces are not parallel or '//
     &    'their side pointers are incorrectly defined.'
        stop
      endif

      ! Set grid size along x:
      detbox = abs(scalarx-scalartop)        ! Detector box size = Distance between planes
      write(*,'(a)') 'Detector box x-size is (cm):'
      write(*,'(1x,es12.5)') detbox
      if (detbox.lt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detector box size is too small.'
        stop
      endif
      pixgrid = nx*dx                         ! Pixels grid size, as defined by user
      if (pixgrid.ne.0.0) then                ! Pixels grid size explicitly given by user
        if (abs(pixgrid-detbox).gt.eps*detbox) then
          write(*,'(a)') '(pixels grid and detector box sizes '//
     &      'do not coincide; detector box x-size is ignored.)'
        endif
      else if (nx.eq.0) then                  ! Infer num. pixels from geometry
        write(*,'(a)') 'Number of pixels along x set to:'
        if (detbox*oneplus/dx.gt.npixmax) then
          write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'overflow in no. pixels; redefine detector box parameters.'
          stop
        endif
        nx = detbox*oneplus/dx
        write(*,'(1x,i0)') nx
      else                                    ! user entered dx=0 => infer pixel size
        dx = detbox/nx
        write(*,'(a)') 'Pixel size along x set to (cm):'
        write(*,'(1x,es12.5)') dx
      endif
      idx = 1.0/dx

      ! Seach for y-normal face:
      isurf = isurf+1                         ! Next surface ID should be y-normal face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'y-normal face is not a plane.'
        stop
      endif
      normy = (/ax(surf),ay(surf),az(surf)/)
      imodu = sidepoint/sqrt(sum(normy**2))
      normy = normy*imodu                     ! Inwards-oriented unit normal
      scalary = a0(surf)*imodu                ! Scalar term, normalized

      ! Seach for y-normal top face:
      isurf = isurf+1                         ! Next surface ID should be y-normal top face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal top face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal top face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'y-normal top face is not a plane.'
        stop
      endif
      normtop = (/ax(surf),ay(surf),az(surf)/)
      imodu = -sidepoint/sqrt(sum(normtop**2))
      normtop = normtop*imodu                 ! Outwards-oriented unit normal, should coincide with NORMY
      scalartop = a0(surf)*imodu              ! Scalar term, normalized

      ! Check if parallel and correctly oriented:
      if (abs(dot_product(normy,normtop)-1.0d0).gt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal and y-normal top faces are not parallel or '//
     &    'their side pointers are incorrectly defined.'
        stop
      endif

      ! Set grid size along y:
      detbox = abs(scalary-scalartop)         ! Detector box size = Distance between planes
      write(*,'(a)') 'Detector box y-size is (cm):'
      write(*,'(1x,es12.5)') detbox
      if (detbox.lt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detector box size is too small.'
        stop
      endif
      pixgrid = ny*dy
      if (pixgrid.ne.0.0) then                ! Pixels grid size given by user
        if (abs(pixgrid-detbox).gt.eps*detbox) then
          write(*,'(a)') '(pixels grid and detector box sizes '//
     &      'do not coincide; detector box y-size is ignored.)'
        endif
      else if (ny.eq.0) then                  ! Infer num. pixels
        write(*,'(a)') 'Number of pixels along y set to:'
        if (detbox*oneplus/dy.gt.npixmax) then
          write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'overflow in no. pixels; redefine detector box parameters.'
          stop
        endif
        ny = detbox*oneplus/dy
        write(*,'(1x,i0)') ny
      else                                    ! Infer pixel size
        dy = detbox/ny
        write(*,'(a)') 'Pixels size along y set to (cm):'
        write(*,'(1x,es12.5)') dy
      endif
      idy = 1.0/dy

      ! Check perpendicularity of x&y faces:
      if (abs(dot_product(normx,normy)).gt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x & y faces are not perpendicular.'
        stop
      endif


      !!!!!!!!!!!!!!!! added script for thickness calculation!!!!!!!!!!

      ! Seach for z-normal face:
      isurf = isurf+1                         ! Next surface ID should be z-normal face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'z-normal face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'z-normal face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'z-normal face is not a plane.'
        stop
      endif
      normz = (/ax(surf),ay(surf),az(surf)/)
	  
      imodu = sidepoint/sqrt(sum(normz**2))
      normz = normz*imodu                     ! Inwards-oriented unit normal
      scalarz = a0(surf)*imodu                ! Scalar term, normalized

      ! Seach for y-normal top face:
      isurf = isurf+1                         ! Next surface ID should be z-normal top face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'z-normal top face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'z-normal top face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'z-normal top face is not a plane.'
        stop
      endif
      normtop = (/ax(surf),ay(surf),az(surf)/)      
      imodu = -sidepoint/sqrt(sum(normtop**2))
      normtop = normtop*imodu                 ! Outwards-oriented unit normal, should coincide with NORMY
      scalartop = a0(surf)*imodu              ! Scalar term, normalized
	  
      
      ! Check if parallel and correctly oriented:
      if (abs(dot_product(normz,normtop)-1.0d0).gt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'z-normal and z-normal top faces are not parallel or '//
     &    'their side pointers are incorrectly defined.'
        stop
      endif

      ! Set grid size along z:
      detbox = abs(scalarz-scalartop)         ! Detector box size = Distance between planes
      write(*,'(a)') 'Detector box z-size is (cm):'
      write(*,'(1x,es12.5)') detbox
      z_top = -scalartop
	  z_low = -scalarz
	  thick = abs(z_low - z_top)

      
      if (detbox.lt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detector box size is too small.'
        stop
      endif
      !z_top = detbox

      ! Lab->Pixels change of bases (only x&y coords will be needed):
      lab2det(1,:) = normx
      lab2det(2,:) = normy

      isLab = .false.
      if (abs(normx(1)-1.0d0).lt.eps.and.
     &    abs(normy(2)-1.0d0).lt.eps) then
        isLab = .true.
        write(*,'(a)') 'Detector reference frame coincides with lab '//
     &    'frame; no rotation applied for scoring.'
      endif

      ! Find an arbitrary point r on the line L defined by the two orthogonal planes
      !   defining the detector box. Point r is selected such that its distance to the
      !   origin in minimal, that is, it lies on the plane perpendicular to L
      !   that passes through the origin. Thus, it satisfies the eqs.
      !      normx.r+scalarx  = 0  (x-face plane)
      !      normy.r+scalary  = 0  (y-face plane)
      !     (normx X normy).r = 0  (orthogonal plane through origin)
      !
      !   ('.' represents scalar product, 'X' represents cross product)
      !   The solution of the above system for r is:
      r0 = -scalarx*normx-scalary*normy
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!  ------------ DETECTOR MODELS --------------


      subroutine dispersion_model(edeptmp_in, x0, y0, zd)
!*******************************************************************
!*  imgmodel == 4  : "Depth-dependent model, with PSF for optical  *
!*                    photons".                                    *
!*                                                                 *
!*   The detector model '4' creates the image with optical photons *
!*   created with a depth-of-interaction-dependent model developed *
!*   by Iacovos Kyprianou based on MANTIS simulations [Kyprianou   *
!*   et al., Proc. SPIE, Vol. 6913, p. 69130 (2008)].              *
!*   The optical photons are generated sampling and combining 3    *
!*   gaussian functions with a weight and sigma that depends on    *
!*   depth of the interaction.                                     *
!*                                                                 *
!*   This function is only called when energy is deposited inside  *
!*   the detector (mode<0):                                        *
!*         arg  ==> Deposited energy in this interaction.          *
!*     xo,yo,zo ==> interaction point in the reference system      *
!*                  where the detector is located on the XZ plane, *
!*                  and centered at the origin.                    *
!*                                                                 *
!*******************************************************************
      use IMGmod
      use TRACK_mod
      use ctrsimMod, only : nhist  ! Current history number
      implicit none
      

      integer*8 maxpix
      parameter(maxpix=4500000)
      real*8 edep,edep2,Dy_2,ipixarea,
     &       omega,theta,phi,rot_inv,xshift,yshift,zshift,
     &       rho
      integer*8 imgmat, imgmodel, binout

      common /scoimg/ edep(maxpix),edep2(maxpix),
     &                Dy_2,ipixarea,
     &                omega,theta,phi,rho,rot_inv(9),xshift,yshift,
     &                zshift,imgmat,imgmodel,binout
           
      integer*8 pix, k, l, i_disp, j_disp, ij_disp, ijk_disp, npairs_k,
     &          i,j,ij
      real*8 edeptmp_in, sample_gauss, efficiency, n,
     &       x_eh,y_eh,num_ehpairs_r, x0, y0, zd, tr
      real*8 rprima(3), sigSpread, weightnorm, numPairsEn, rand
      external rand
      !integer out3



      !out3 = 200
      !open(out3,file='x_eh-y_eh.dat')

      
	  
      numPairsEn = edeptmp_in/eh      ! E(eV)/eh(eV) --> generated e-h pairs
	  !write(out3,'(3(1x,es12.5))') eh,edeptmp_in,numPairsEn
      num_ehpairs_r = sample_gauss(numPairsEn, sqrt(fano*numPairsEn))   ! Fano distribution of the e-h pairs
	  
	  if (ccc.eq.0)then
		tr = (1 - EXP(-zd/(e_tau*abs(EField))))*(e_tau*abs(EField))/thick
	  else
 		tr = (1 - EXP(-zd/(h_tau*abs(EField))))*(h_tau*abs(EField))/thick
	  endif


      if (num_ehpairs_r.lt.0.5) return   ! Do not create a negative (or 0) number of optical photons.
            
      ! Spread as a function of depth (z):
      sigSpread=sqrt(2d0*(1.3806d-23)*298d0*zd/((1.6d-19)*abs(EField)))  !! Dispersion as a function of the electric field
      !write(out3,'(5(1x,es12.5))') zd,sigSpread,tr,e_tau,h_tau
	  
      npairs_k = idnint(num_ehpairs_r)      !rounds to an integer

	 do l = 1, npairs_k   !! (For each e-h pair)

		  !call sample_gauss2D(x0,y0,sigSpread,sigSpread,x_eh,y_eh)
		  x_eh = sample_gauss(x0,sigSpread)
          y_eh = sample_gauss(y0,sigSpread)
		  ! Count the optical photon in the pixel:


          !!!!!!! Same as the original tally to calculate the pixel index  !!!!!!!

          if (isLab) then                       ! Slab and Sdet basis vectors coincide, no rotation needed
            x_eh = x_eh-r0(1)                      ! Position vector w.r.t. pixels origin, only x&y coords are needed
            y_eh = y_eh-r0(2)
          else                                  ! Both a translation and a rotation are needed
            rprima = (/x_eh,y_eh,zd/)-r0               ! Translation
            x_eh = dot_product(lab2det(1,:),rprima) ! Rotation to obtain coords in Sdet frame
            y_eh = dot_product(lab2det(2,:),rprima)
          endif

          ! Check if particle is inside pixels grid:
          i_disp = x_eh*idx+1.0                      ! Pixel indexDetector box z-size is (cm):
          if (i_disp.lt.1.or.i_disp.gt.nx) cycle          ! Out of pixels grid
          j_disp = y_eh*idy+1.0                      ! Pixel index
          if (j_disp.lt.1.or.j_disp.gt.ny) cycle          ! Out of pixels grid
          ij_disp = i_disp+(j_disp-1)*nx                       ! Mapping of i,j into a single index

		  !write(out3,'((1x,es12.5),i3)') edeptmp_in,l
			
          
          !write(out3,'(2(1x, i9))') l, npairs_k
          !write(out3,'(1x,es12.5)') edeptmp(ij_disp)

          !wghtlast(ij_disp) = wght

          if(isTrap.and.peso.gt.0)then
            !write(out3,'(2(1x,es12.5))') tr, rand(1.d0)
            if(rand(1.d0).le.tr)then
                e_counts(ij_disp) = e_counts(ij_disp)+peso
                std(ij_disp) = std(ij_disp)+ 1
                !e_counts2(ij_disp) = e_counts2(ij_disp)+peso**2
                if(ilb(5).eq.0) then
                    edeptmp(ij_disp)  = edeptmp(ij_disp)+eh*peso
                else
                    edeptmp_f(ij_disp)  = edeptmp_f(ij_disp)+eh*peso
                endif
            endif
          else
            e_counts(ij_disp) = e_counts(ij_disp)+peso
            std(ij_disp) = std(ij_disp)+ 1
            !e_counts2(ij_disp) = e_counts2(ij_disp)+peso**2
            if(ilb(5).eq.0) then
                edeptmp(ij_disp)  = edeptmp(ij_disp)+eh*peso
            else
                edeptmp_f(ij_disp)  = edeptmp_f(ij_disp)+eh*peso
            endif
          endif


          
		  
      enddo

      end


      subroutine IMGdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode: 1 to write dump file, else to read from it.          *
!*******************************************************************
      use IMGmod
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) e_counts(1:ijmax),e_counts2(1:ijmax)
      else
        if (restartfile.lt.0) return  ! No restart file open
        write(dumpfile) e_counts(1:ijmax),e_counts2(1:ijmax)
      endif
      end


      function dist(Efield,ccc,z,z_top, z_low)
!********************************************************************
!*  Sample a Gaussian distribution centered in the input mean value *
!*  and with the input sigma, using the Box-Muller algorithm.       *
!*     -- Return value: sampled gaussian value                      *
!********************************************************************
      real*8 dist,Efield,z,z_low,z_top
	  integer*8 ccc
	  
	  !colecting electrons
	  if (Efield.lt.0.and.ccc.eq.0) dist = abs(z-z_low)
	  if (Efield.gt.0.and.ccc.eq.0) dist = abs(z-z_top)
	  if (Efield.lt.0.and.ccc.eq.1) dist = abs(z-z_top)
 	  if (Efield.gt.0.and.ccc.eq.1) dist = abs(z-z_low)
     
     
      return
      end  
	  
      function sample_gauss(mean,sigma)
!********************************************************************
!*  Sample a Gaussian distribution centered in the input mean value *
!*  and with the input sigma, using the Box-Muller algorithm.       *
!*     -- Return value: sampled gaussian value                      *
!********************************************************************
      implicit none
      real*8 sample_gauss, mean, sigma, PIx2, rand
      external rand
      parameter (PIx2=(2.0d0*3.14159265358979d0))
      
      sample_gauss = mean + 
     &   sigma * sqrt(-2.0d0*log(rand(1.d0)))*sin(PIx2*rand(2.d0))
     
      return
      end  
        
      subroutine sample_gauss2D(mean1,mean2,sigma1,sigma2,gauss1,gauss2)
!********************************************************************
!*  Sample two Gaussian distributions with the input mean and sigma *
!*  using the Box-Muller algorithm.                                 *
!*  The sampled values as returned in arguments gauss1 and gauss2   *
!********************************************************************
      implicit none
      real*8 mean1,mean2,sigma1,sigma2,PIx2,rand,
     &          gauss1,gauss2,radious,theta
      external rand
      parameter (PIx2=(2.0d0*3.14159265358979d0))
      radious=sqrt(-2.0d0*log(rand(1.d0)))
      theta=PIx2*rand(2.d0)
      gauss1 = mean1 + sigma1*radious*cos(theta)
      gauss2 = mean2 + sigma2*radious*sin(theta)
      return
      end

!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
