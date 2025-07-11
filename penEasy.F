!*******************************************************************
!*                         penEasy                                 *
!*                                                                 *
!* Short description:                                              *
!*   General-purpose main program for the PENELOPE system.         *
!*   Please refer to the README.txt file for detailed instructions.*
!*                                                                 *
!*   For a list of dependencies from PENELOPE, see                 *
!*   ~/documentation/dependencies.txt                              *
!*                                                                 *
!* Josep Sempau                                                    *
!* email: josep.sempau@upc.es                                      *
!* Universitat Politecnica de Catalunya, Barcelona, Spain          *
!* SEE COPYRIGHT NOTICE IN FILE README.txt                         *
!*******************************************************************

! PENELOPE library:
#include "penelope.f"
#include "pengeom.f"
#include "rita.f"
#include "penvared.f"
! Protons:
#ifdef PROTONS
#include "penh.f"
#endif

! penEasy library:
#include "penaux.F"
#include "penvox.F"
#include "penvr.F"
#include "timing.F"

! Source models:
#include "sourceBoxIsotropicGaussSpectrum.F"
#include "sourcePhaseSpaceFile.F"
#include "sourceRectangularBeam_2010-09-02.f"

! Tallies:
#include "tallyCylindricalDoseDistrib.F"
#include "tallyEnergyDeposition.F"
#include "tallyFluenceTrackLength.F"
#include "tallyParticleCurrentSpectrum.F"
#include "tallyParticleTrackStructure.F"
#include "tallyPhaseSpaceFile.F"
#include "tallyPixelImageDetect_THOR.F"
#include "tallyPulseHeightSpectrum.F"
#include "tallySpatialDoseDistrib.F"
#include "tallySphericalDoseDistrib.F"
#include "tallyVoxelDoseDistrib.F"


!*******************************************************************
!*    MAIN                                                         *
!*******************************************************************
      program main
      use TRACK_mod, only : mat,e,kpar
      use ctrsimMod, only : nhist
      use dsmaxMod
      use forcingMod, only : isforcing
      implicit none
      logical endsim,absorb
      integer*4 ncross,icol,left
      real*8 ds,dsef,de

      write(*,'(a)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') '>>>> This is penEasy v.2015-05-30 >>>>'
      write(*,'(a)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      call init      ! Initialize the PENELOPE/penEasy system and timers

      history: do              ! Each iteration simulates a new history
        nhist = nhist+1.0d0    ! Update history counter
        call cleans            ! Empty the stack, just in case
        call source            ! Put primary particles (from same history) in stack
        call tally(1,nhist)    ! The simulation of this history begins

        particle: do                       ! Each iteration simulates a new particle
          call secpar(left)                ! Retrieve a particle from the stack
          if (left.eq.0) exit particle     ! Stack was empty
          call locatevox                   ! Locate in voxels geometry (info not in PENELOPE's stack)
          call tally(-99,-e)               ! The simulation of this particle begins
          if (absorb()) cycle particle     ! Check particle absorption
          call start                       ! Reset transport mechanics
          call forcing                     ! Set interaction forcing status

          interact: do                     ! Each iteration simulates an interaction
            if (absorb()) cycle particle   ! Check particle absorption

! Following lines will be processed only if compiled with option -DPROTONS:
#ifdef PROTONS
            if (kpar.eq.4) then            ! It's a proton
              call jumph(dsmax(mat),ds)    ! Get distance to next interaction (note: forcing not implemented)
            else                           ! Photon-electron transport
#endif
            if (isforcing) then            ! Interaction forcing is active
              call jumpf(dsmax(mat),ds)    ! Get distance DS until next interaction
            else
              call jump(dsmax(mat),ds)     ! Same without interaction forcing
            endif
#ifdef PROTONS
            endif
#endif

            call stepx(ds,dsef,ncross)     ! Advance up to interaction point or interface
            if (ncross.eq.0) then          ! Moved a distance DS, no interface crossed
              call tally(3,ds)             ! Notify tallies
            else                           ! Moved a distance DSEF, interface found
              call tally(4,dsef)           ! Notify tallies
              if (mat.eq.0) cycle particle ! New material is vacuum => gone
              call start                   ! New material => reset transport mechanics
              call forcing                 ! Set interaction forcing status
              call russian                 ! Apply Russian roulette
              call splitting               ! Apply particle splitting
              cycle interact
            endif

#ifdef PROTONS
            if (kpar.eq.4) then            ! Proton transport
              call knockh(de,icol)         ! Simulate an interaction
            else                           ! Photon-electron transport
#endif
            if (isforcing) then            ! Interaction forcing (IF) is active
              call knockf(de,icol)         ! Simulate an interaction with IF
            else
              call knock(de,icol)          ! Same with no IF
            endif
#ifdef PROTONS
            endif
#endif

            call tally(-int(icol),de)      ! Tally kinetic energy released
          enddo interact
        enddo particle

        call tally(6,nhist)                ! End-of-history bookkeeping
        if (endsim()) exit history         ! Simulation is finished
      enddo history

      call report(nhist)                   ! Write final report
      end


      subroutine init
!*******************************************************************
!*    Initializes the simulation system.                           *
!*                                                                 *
!*    Comments:                                                    *
!*      -> The history counter is set to >0 in case of a restart   *
!*      -> The order in which init routines are called must not be *
!*         changed, since some of them provide data for others and *
!*         write data to a dump file. This file must be written in *
!*         the same order in which is will be read afterwards.     *
!*******************************************************************
      use ctrsimMod, only : nhist
      use RSEEDcommonMod
      implicit none

      integer restartf
      integer*4 nmatgeo
      real*8 emax,realtime,cputime,rtime,utime,dmem
      real*8 :: mem=0.0         ! Memory usage

      call initime              ! Write date on the screen
      call treset               ! Reset simulation timer to compute init() timing

      call iniconfig(restartf)  ! Simulation config
      nhist = 0.0d0             ! Clear number of simulated histories
      if (restartf.ge.0) read(restartf) nhist,seed1,seed2  ! Simulation restart
      call inisource(emax,dmem) ! Source models
      mem = mem+dmem            ! Used memory
      call inigeo(nmatgeo,dmem) ! Geometry: PENGEOM & penVox
      mem = mem+dmem
      call inipen(emax,nmatgeo) ! PENELOPE init
      call inimassvox(dmem)     ! Init voxels masses (cannot do before PENELOPE inits)
      mem = mem+dmem
      call initally(emax,dmem)  ! Tallies; tallies read from restart file if needed
      mem = mem+dmem
      call iniforce(emax)       ! Interaction forcing
      call inisplit             ! Particle splitting
      call inirussia            ! Russian roulette

      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'init: INITIALIZATION ENDED'
      write(*,'(a,f9.2,a)') 'Elapsed real time:',realtime(),' s'
      write(*,'(a,f9.2,a)') 'Elapsed CPU time :',cputime(),' s'
      write(*,'(a,es7.1,a)') "Memory used by penEasy arrays: ",
     &  mem*1.0e-6,' MB'

      call treset               ! Reset simulation timers

      if (restartf.ge.0) then   ! This is a sim restart
        read(restartf) rtime,utime
        close(restartf)
        call trestart(rtime,utime) ! Re-reset CPU and realtime timers in case of a restart
        write(*,*) ''
        write(*,'(a)') 'Simulation restarted with the following '//
     &                 'data from the previous run:'
        write(*,'(a,f18.0)') '  No. of histories: ',nhist
        write(*,'(a,2(1x,i0))') '  Seed1,Seed2:',seed1,seed2
        write(*,'(a,2(1x,es12.5))') '  RealTime(s),CPUtime(s):',
     &    rtime,utime
      endif

      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      flush(6)                  ! Flush stdout
      end


      subroutine report(n)
!*******************************************************************
!*    Reports final results.                                       *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated.                           *
!*******************************************************************
      use RSEEDcommonMod
      implicit none
      real*8 n

      integer unc
      real*8 cputime,realtime,nowcpu

      nowcpu = cputime()                 ! Set a reference to measure report timing
      call tallyreport(-1,n,nowcpu,unc)  ! EOS report, each tally will report its data

      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'report: SIMULATION ENDED'
      write(*,'(a)')
     & 'Results have been written to the corresponding DAT files.'
      select case(unc)
      case(0)
        write(*,'(a)')
     &   'The requested uncertainty has NOT been reached.'
      case(1)
        continue  ! Uncertainty limit not defined, no message written
      case default
        write(*,'(a)')
     &   'The requested uncertainty has been reached.'
      end select
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      write(*,*) ''
      write(*,'(a)') 'Last random seeds:'
      write(*,'(2(1x,i0))') seed1,seed2

      write(*,'(a)') 'Elapsed real time (s), excluding init:'
      write(*,'(1x,es12.5)') realtime()

      write(*,'(a)') 'Elapsed CPU time (s), excluding init:'
      write(*,'(1x,es12.5)') nowcpu

      write(*,'(a)') 'Each report update took (in CPU s):'
      write(*,'(1x,es12.5)') cputime()-nowcpu

      write(*,'(a)') 'No. of histories simulated:'
      write(*,'(1x,f18.0)') n

      if (nowcpu.gt.0.0) then
        write(*,'(a)') 'CPU Speed (histories/s):'
        write(*,'(1x,es12.5)') n/nowcpu
      endif

      call endtime  ! Report date and say goodbye
      end


!*******************************************************************
!*******************************************************************
!*    Source routines start here.                                  *
!*    Source models require:                                       *
!*     i) an initialization routine that must be called by         *
!*        INISOURCE                                                *
!*    ii) a particle generation routine that must be called        *
!*        by SOURCE                                                *
!*******************************************************************
!*******************************************************************

      subroutine inisource(emax,mem)
!*******************************************************************
!*    Init routines for source models.                             *
!*                                                                 *
!*    Output:                                                      *
!*      emax: max source energy (eV).                              *
!*      mem: memory used by source arrays.                         *
!*    Comments:                                                    *
!*      - Note that, if more than one source is defined and        *
!*        active, then EMAX must be the maximum of the value       *
!*        reported by all of them.                                 *
!*******************************************************************
      implicit none
      real*8 emax,mem
      integer nsrc
      logical active,activepsf
      real*8 dmem

      nsrc = 0
      mem = 0.0
      call BIGSinisrc(active,emax,dmem)
	  if (active) nsrc = nsrc+1
      !write(*,'(a)') nsrc
	  mem = mem+dmem                   ! Used memory
      call PSFSinisrc(active,emax)
	  if (active) nsrc = nsrc+1
	!  write(*,'(a)') nsrc
	  
	  call RECTBEAMinisrc(active,emax)
	  if (active) nsrc = nsrc+1
	 ! write(*,'(a)') nsrc
	  
      if (nsrc.gt.1) then
        write(*,'(a)') 'inisource:ERROR: PSF source ON is '//
     &    'incompatible with other sources ON.'
        stop
      endif
      end


      subroutine source
!*******************************************************************
!*    Source models.                                               *
!*******************************************************************
      implicit none
      call BIGSsource
      call PSFSource
  	  call RECTBEAMsource 
      end


!*******************************************************************
!*******************************************************************
!*    Tally routines start here.                                   *
!*    Tallies require:                                             *
!*    i) an initialization routine that must be called by INITALLY *
!*    ii) a tally routine that must be called by TALLY             *
!*    iii) a reporting routine that must be called by TALLYREPORT  *
!*    iv) a dump routine to allow restarting of the simulation.    *
!*                                                                 *
!*    Notice that the ordering of the tally initialization routines*
!*    must coincide with the ordering of the corresponding sections*
!*    in the input file.                                           *
!*******************************************************************
!*******************************************************************

      subroutine initally(emax,mem)
!*******************************************************************
!*    Init tallying routines.                                      *
!*                                                                 *
!*    Input:                                                       *
!*      emax: max source energy (eV).                              *
!*    Output:                                                      *
!*      mem: memory used by arrays in tallies.                     *
!*******************************************************************
      implicit none
      real*8 mem,emax,dmem

      mem = 0.0
      call VDDinitally(dmem)
      mem = mem+dmem
      call SDDinitally(dmem)
      mem = mem+dmem
      call CDDinitally(dmem)
      mem = mem+dmem
      call SPDinitally(dmem)
      mem = mem+dmem
      call EDPinitally(dmem)
      mem = mem+dmem
      call PHSinitally(dmem)
      mem = mem+dmem
      call IMGinitally(dmem)
      mem = mem+dmem
      call FTLinitally(dmem)
      mem = mem+dmem
      ! call FPDinitally(dmem,emax)  !! Under construction
      ! mem = mem+dmem
      call PSFinitally(emax)
      call PCSinitally(dmem)
      mem = mem+dmem
      call PTSinitally
      end


      subroutine tally(mode,arg)
!*******************************************************************
!*    Tallying routines.                                           *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      call VDDtally(mode,arg)  ! Must be 1st tally to be called
      call SDDtally(mode,arg)
      call CDDtally(mode,arg)
      call SPDtally(mode,arg)
      call EDPtally(mode,arg)
      call PHStally(mode,arg)
      call IMGtally(mode,arg)
      call FTLtally(mode,arg)
      ! call FPDtally(mode,arg)  !! Under construction
      call PSFtally(mode,arg)
      call PCStally(mode)
      call PTStally(mode,arg)
      end


      subroutine tallyreport(mode,n,cputim,unc)
!*******************************************************************
!*    Calls report routines for all tallies and PSF sources.       *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      unc -> 0: uncert not reached, 1: undefined, >1: reached    *
!*    Comments:                                                    *
!*      -> The sequence of calls to the report routines MUST be the*
!*         same used for the corresponding initialization routines *
!*         in INISOURCE and INITALLY (see above). This is to ensure*
!*         proper manipulation of the simulation dump file.        *
!*      -> PSF sources need to report the total num. of particles  *
!*         read to the dump file, which may be needed in case of a *
!*         restart.                                                *
!*******************************************************************
      use dumpMod
      use RSEEDcommonMod
      implicit none
      integer mode,unc,uncdone
      real*8 n,cputim

      integer finduf,error
      real*8 realtime,cputime

      ! Write dump file if needed:
      if (mode.ne.1.and.dumpfile.ne.-1) then  ! Create sim dump file
        dumpfile = finduf()              ! Find a valid unit for the file
        open(dumpfile,file=dumpfilen,status='replace',access='stream',
     &       iostat=error)
        if (error.ne.0) then
          write(*,'(a)')
     &      'tallyreport:ERROR: unable to open dump file.'
          dumpfile = -2                  ! Tells others that file coundn't be opened
        else
          write(dumpfile) n,seed1,seed2  ! Write history state to dump file
        endif
      endif

      ! Write source reports:
      call PSFSreport(mode)              ! PSF source dumps the num. of particles read

      ! Write tally reports:
      unc = 1                            ! Unknown uncertainty status at this point
      call VDDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call SDDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call CDDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call SPDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call EDPreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call PHSreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call IMGreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call FTLreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      ! call FPDreport(mode,n,cputim,uncdone)  !! Under construction
      ! unc = unc*uncdone
      call PSFreport(mode,n,cputim)      ! No uncertainty for this tally
      call PCSreport(mode,n,cputim,uncdone)
      unc = unc*uncdone

      if (mode.ne.1.and.dumpfile.ge.0) then
        write(dumpfile) realtime(),cputime() ! Write timings to dump file
        close(dumpfile)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
