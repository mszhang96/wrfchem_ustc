  !...................................................................................
  ! reorganized by Kai Zhang @ PNNL 
  !...................................................................................
  !
  ! ICE FREEZING PARAMETERIZATION FILES START HERE
  !
  ! ..................................................................................
  !
  ! Code Developer
  ! Donifan Barahona, GA TECH
  ! donifan@gatech.edu
  !
  ! ....................................
  !    DESCRIPTION
  ! ....................................
  !
  ! Parameterization of ICE crystal number concentration 
  ! for large scale models. 
  ! Donifan Barahona, Athanasios Nenes
  !    
  ! J. Geophys. Res., 111, D11211, 2008
  ! Atmos. Chem. Phys., 9, 369-381, 2009 
  ! Atmos. Chem. Phys., 9, 5933-5948, 2009  
  ! 
  ! SI units unless otherwise specified. 
  ! 
  ! INPUT:  
  ! 
  !      T      - Temperature (K)
  !      P      - Pressure  (Pa)
  !     nSO4    - Deliquesced aerosol number concentration (m-3)
  !     ndust   - Total dust number concentration (m-3)  
  !                 If an array is passed then information on average sizes must be included in the INspec function 
  !!    nbc     - Total black carbon number concentration (m-3)
  !
  ! OUTPUT:  
  !
  !     smax  - Maximum supersaturation wrt ice (absolute) 
  !     nice  - Total ice crystal concentration (cm-3)
  !     nhet  - Ice crystal concentration from heteroegeneous freezing (cm-3)         
  !     nlim  - Maximum ice nuclei concentration that allows homogeneous freezing (cm-3)
  !
  ! OPTIONS:
  !  
  !   Calculation options are set using several global variables. 
  !   Freezing mechanism can be set to 
  !  
  !     1) homogeneous only (purehom_ice = .TRUE.)
  !     2) heterogeneous only (purehet_ice=.TRUE.)
  !     3) both mechanisms included (purehet_ice = purehom_ice = .FALSE.) 
  !   
  !   Heterogeneous freezing is described using the function INSPEC_ice. Several correlations are included 
  !   and selected using the variable  typeofspec_ice (integer) as follows (notice that additional parameters 
  !   may need to be set depending on the spectrum):
  ! 
  !    *  -1;   Assume all IN as monodisperse with a single freezing threeeshold
  !    *   1;   Use the correlation of Meyers et al. 1992.  
  !    *   2;   Use the correlation of Phillips et al. 2007.
  !    *   3;   Use the CNT spectra of Barahona and Nenes, 2009.
  !    *   4;   Use the correlation of Phillips et al. 2008.
  !    *   5+;  Additional formulations can be included by modifying INSPEC_ice.
  !
  !   Three options are used to specify updraft velocity: 
  ! 
  !   (1) Using a single value; set use_av_v= .TRUE. and pass the updraf using updraft_ice   
  !   (2) Integrating over a normal PDF of updrafts: set use_av_v= .FALSE. and linear_V= .FALSE. and 
  !       pass the variables vmin_ice, vmax_ice, sigmav_ice  and miuv_ice
  !   (3) Integrating over a normal PDF of updrafts with sigmav_ice linearly 
  !       dependent on temperature: set linear_V= .TRUE. 
  !
  !  WRITTEN BY DONIFAN BARAHONA
  !
  !  Xiangjun Shi at PNNL made some changes (noted ++sxj) for bugs and coding
  !  June 2010
  !
  !=======================================================================

   module module_icenucleation_bn

     use shr_kind_mod,  only:  r8 => shr_kind_r8

     implicit none
     private
     save

     public IceParam_BN         ! ! Public interfaces

     real(r8),parameter :: amw_ice     = 0.029         !Air molecular Weight 
     real(r8),parameter :: wmw_ice     = 0.018         !Water molecular weight
     real(r8),parameter :: rgas_ice    = 8.314         !Universal gas constant
     real(r8),parameter :: grav_ice    = 9.81          !gravity
     real(r8),parameter :: cpa_ice     = 1005.1        !Air thermal gas capacity
     real(r8),parameter :: pi_ice      = 3.1415927   
     real(r8),parameter :: depcoef_ice = 0.1           !Default deposition coefficient
     real(r8),parameter :: thaccom_ice = 0.7           !Default thermal accommodation coefficient
     real(r8),parameter :: ddry_ice    = 4.0e-8        !Droplet mean geometric diameter, m
     logical,parameter  :: linear_V    = .FALSE.       !TRUE to use linear relation for sigma V from Penner & Wang ACPD 2009   

   contains
   

     subroutine IceParam_BN (T, P, nso4, ndust, nbc,  nhet, nice, smax, nlim, &
                          miuv_ice,vmin_ice,vmax_ice,sigmav_ice, &
                          updraft_ice,shdust_ice,effdust_ice,shbc_ice,effbc_ice, &
                          nin_ice,sh_ice,use_av_v, &
                          purehom_ice,purehet_ice,typeofspec_ice )
       implicit none

       real(r8),intent(in) :: T, P, nso4, ndust, nbc  !SCI unit k Pa #/m3
       real(r8),intent(in) :: miuv_ice             !Mean velocity for the pdf of updraft (m/s)
       real(r8),intent(in) :: vmin_ice             !Lower limit of integration when using a pdf of updraft (m/s)
       real(r8),intent(inout) :: vmax_ice          !Upper limit of integration when using a pdf of updraft (m/s)
       real(r8),intent(inout) :: sigmav_ice        !Standard deviation of the pdf of updraft distribution (m/s)
       real(r8),intent(in) :: updraft_ice          !When using fixed updraft velocity (m/s)
       real(r8),intent(in) :: shdust_ice           !Maximum freezing threeshold dust 
       real(r8),intent(in) :: effdust_ice          !Maximum freezing efficiency dust
       real(r8),intent(in) :: shbc_ice             !Maximum freezing threeshold bc only used in CNT
       real(r8),intent(in) :: effbc_ice            !Maximum freezing efficiency bc Only used in the Monodisperse approximation:       
       real(r8),intent(in) :: nin_ice              !Assumed Nin concentration (m-3)
       real(r8),intent(in) :: sh_ice               !Assumed Nin freeizng threeshold- monodisperse
       logical, intent(in) :: use_av_v             !FALSE to integrate over a pdf of updraft
       logical, intent(in) :: purehom_ice          !True supresses homog nucleation
       logical, intent(in) :: purehet_ice          !True supresses heter nucleation 
       integer, intent(in) :: typeofspec_ice       !Integer 1 to 4 see function NcSpec. -1 is monodisperse
       real(r8), intent(out) :: nice,  nhet, smax, nlim    ! #/m3

       real(r8) :: T_ice, P_ice
       real(r8) :: dhs_ice ,smax_ice,np_ice
       real(r8) :: rv_ice, denwat_ice, denice_ice, vpresw_ice, vpresi_ice
       real(r8) :: denair_ice, diff_ice, aircond_ice
       real(r8) :: alfa_ice, beta_ice, shom_ice, koft_ice 
       real(r8) :: dliq_ice, g1_ice, g2_ice, gdoin_ice, z_ice, doin_ice        
       real(r8) :: ndust_ice, lambda_ice, kdust_ice, kbc_ice 
       real(r8) :: del0dust_ice, si0dust_ice 
       real(r8) :: del0bc_ice, si0bc_ice, nbc_ice

       ! Conditions 

       P_ice=P
       T_ice=T
       ndust_ice=ndust
       np_ice=nso4
       nbc_ice=nbc

       ! Ff using a sigma_V(T) correlation 

       if (linear_V) then
       ! linear_V TRUE to usea linear relation for sigma V from Penner & Wang ACPD 2009
          if (T .lt. 235.) then
             sigmav_ice=0.25-0.006*(235.-T)  
          end if
          if (T .lt. 195.) then
             sigmav_ice=0.01
          end if
          !correct upper limit for integration
          vmax_ice=min(vmin_ice+(4.0*sigmav_ice), 0.5_r8)  
       end if 

       ! Set ice properties and physical constants 
       call prop_ice(T_ice, P_ice, &
                     rv_ice,dhs_ice,vpresw_ice,vpresi_ice, &
                     denice_ice,denwat_ice,denair_ice,diff_ice, &
                     aircond_ice,alfa_ice,beta_ice,shom_ice, &
                     koft_ice,dliq_ice,g1_ice,g2_ice, &
                     doin_ice,gdoin_ice,z_ice,lambda_ice, &
                     kdust_ice,si0dust_ice,del0dust_ice, &
                     Kbc_ice,  si0bc_ice,  del0bc_ice  )

       ! Calculation of smax, nice
       ! use_av_v          !FALSE to integrate over a pdf of updraft

       if (use_av_v) then

         call nice_param(updraft_ice, T_ice, P_ice, nice, smax, nhet, nlim, &
                           np_ice,dliq_ice,alfa_ice,beta_ice,z_ice,shom_ice,g1_ice,g2_ice, &
                           shdust_ice,ndust_ice,effdust_ice,kdust_ice,si0dust_ice,del0dust_ice, &
                           shbc_ice,nbc_ice,effbc_ice,kbc_ice,si0bc_ice,del0bc_ice, &
                           typeofspec_ice, &
                           vpresi_ice,vpresw_ice,sh_ice,nin_ice, &
                           denair_ice,denice_ice,koft_ice,gdoin_ice, &
                           purehom_ice,purehet_ice, &
                           lambda_ice)

       else

         call nice_Vdist(T_ice, P_ice, nice, smax, nhet, nlim, &
                           np_ice,dliq_ice,alfa_ice,beta_ice,z_ice,shom_ice,g1_ice,g2_ice, &
                           shdust_ice,ndust_ice,effdust_ice,kdust_ice,si0dust_ice,del0dust_ice, &
                           shbc_ice,nbc_ice,effbc_ice,kbc_ice,si0bc_ice,del0bc_ice, &
                           typeofspec_ice, &
                           vpresi_ice,vpresw_ice,sh_ice,nin_ice, &
                           denair_ice,denice_ice,koft_ice,gdoin_ice, &
                           purehom_ice,purehet_ice, &
                           lambda_ice, &
                           vmax_ice,vmin_ice,miuv_ice,sigmav_ice)

       end if

       ! write(*,*) "use_av_v",use_av_v
       ! convert to cm-3
       !nice=nice*1.0e-6 
       !nhet=nhet*1.0e-6  
       !nlim=nlim*1.0e-6 

     end subroutine IceParam_BN
  
!------------------------------------------------------------------------------------------
!    Subroutine nice_Vdist. Calculates the ice crystal number concentration
!    at the maximum supersaturation using a PDF of updraft integrating with a 
!    sixth order Gauss-Legendre quadrature  
!     Inputs:  T_ice, and P_ice all SI units)
!------------------------------------------------------------------------------------------

     subroutine nice_Vdist(T_ice, P_ice, nice, smax, nhet, nlim, &
                           np,dliq,alfa,beta,z,shom,g1,g2, &
                           shdust,ndust,effdust,kdust,si0dust,del0dust, &
                           shbc,nbc,effbc,kbc,si0bc,del0bc, &
                           typeofspec, &
                           vpresi,vpresw,sh_mono,nin_mono, &
                           denair,denice,koft,gdoin, &
                           purehom,purehet, &
                           lambda, &
                           vmax,vmin,miuv,sigmav )
       implicit none

       real(r8), intent(in) :: T_ice, P_ice
       real(r8), intent(in) :: np,dliq,alfa,beta,z,shom,g1,g2
       real(r8), intent(in) :: shdust,ndust,effdust,kdust,si0dust,del0dust
       real(r8), intent(in) :: shbc,nbc,effbc,kbc,si0bc,del0bc
       integer , intent(in) :: typeofspec
       real(r8), intent(in) :: vpresi,vpresw,sh_mono,nin_mono
       real(r8), intent(in) :: denair,denice,koft,gdoin
       logical,  intent(in) :: purehom,purehet
       real(r8), intent(in) :: lambda
       real(r8), intent(in) :: vmax,vmin,miuv,sigmav
       real(r8), intent(out) :: nice, smax, nhet, nlim 

       integer :: INDEX
       real(r8)  normv_ice
       real(r8) :: quadx(6), wpar, sum1, quadw(6), dp, sum2, sum3, sum4, x1, x2 
       data  quadx/0.23861918_r8,-0.23861918_r8, 0.66120939_r8,-0.66120939_r8, 0.93246951_r8,-0.93246951_r8/
       data  quadw/0.46791393_r8, 0.46791393_r8, 0.36076157_r8, 0.36076157_r8, 0.17132449_r8, 0.17132449_r8/

       ! calculate the integral in the denominator
       x1=(vmin-miuv)/(sqrt(2.0)*sigmav)
       x2=(vmax-miuv)/(sqrt(2.0)*sigmav)

       ! cummulative width of the distribution of velocities 
       normv_ice=(ERFAPP(x2)-ERFAPP(x1))*0.5
       sum1=0.
       sum2=0.
       sum3=0.
       sum4=0.

       ! use a Gauss-Legendre Quadrature
       do INDEX =1, 6
          wpar=0.5*(((vmax-vmin)*quadx(INDEX))+(vmax+vmin))                  
          call nice_param(wpar, T_ice, P_ice, nice, smax, nhet, nlim, &
                           np,dliq,alfa,beta,z,shom,g1,g2, &
                           shdust,ndust,effdust,kdust,si0dust,del0dust, &
                           shbc,nbc,effbc,kbc,si0bc,del0bc, &
                           typeofspec, &
                           vpresi,vpresw,sh_mono,nin_mono, &
                           denair,denice,koft,gdoin, &
                           purehom,purehet, &
                           lambda)
          call gausspdf(wpar, dp, &
                        miuv,sigmav,normv_ice)
          sum1=sum1+(nice*dp*quadw(INDEX))
          sum2=sum2+(smax*dp*quadw(INDEX))
          sum3=sum3+(nhet*dp*quadw(INDEX))
          sum4=sum4+(nlim*dp*quadw(INDEX))
       enddo

       nice=sum1*(vmax-vmin)*0.5
       smax=sum2*(vmax-vmin)*0.5
       nhet=sum3*(vmax-vmin)*0.5
       nlim=sum4*(vmax-vmin)*0.5

     contains

       ! Approximation to the error function
       real(r8) function ERFAPP(x)
          real(r8),  intent(in) :: x 
          real(r8) :: ax  
          ax=x*x*(1.27324+(0.147*x*x))/(1.0+(0.147*x*x))
          ERFAPP=SQRT(1.0-exp(-ax))
       end function ERFAPP

     end subroutine nice_Vdist

!------------------------------------------------------------------------------------------
!    Subroutine nice_param. Calculates the ice crystal number concentration
!    at the maximum supersaturation assuming that all crystals 
!    freeze at the same size. Inputs: wpar_ice, T_ice, and P_ice all SI units)
!    Output nice, nhet, nlim_ (m-3), and smax.
!------------------------------------------------------------------------------------------

     subroutine nice_param(wpar_ice,T_ice, P_ice, nice, smax, nhet, nlim_, &
                           np,dliq,alfa,beta,z,shom,g1,g2, &
                           shdust,ndust,effdust,kdust,si0dust,del0dust, &
                           shbc,nbc,effbc,kbc,si0bc,del0bc, &
                           typeofspec, &
                           vpresi,vpresw,sh_mono,nin_mono, &
                           denair,denice,koft,gdoin, &
                           purehom,purehet, &
                           lambda)
       implicit none
       real(r8), intent(in) :: T_ice, P_ice, wpar_ice     ! wpar  updraft velocity
       real(r8), intent(in) :: np,dliq,alfa,beta,z,shom,g1,g2
       real(r8), intent(in) :: shdust,ndust,effdust,kdust,si0dust,del0dust
       real(r8), intent(in) :: shbc,nbc,effbc,kbc,si0bc,del0bc
       integer , intent(in) :: typeofspec
       real(r8), intent(in) :: vpresi,vpresw,sh_mono,nin_mono
       real(r8), intent(in) :: denair,denice,koft,gdoin
       logical,  intent(in) :: purehom,purehet
       real(r8), intent(in) :: lambda
       real(r8), intent(out):: nice, smax, nhet, nlim_

       !++sxj
       real(r8) lambda_local
       real(r8) AUX1, AUX2, G
       real(r8) DPMAX, MIU, MONVOL, FDS, NLIM, DLIM,  DSTAR, DS, NSTAR 
       real(r8) NHOM, FC, PHIDO, AUXNC, SIZECORR, DSH, NHET_, F1, F2 
       real(r8) SAUX, SUP, SLOW, SHALF, FHALF, DPMIN, GAM
       integer :: INDEX 

       !czhao initialize
       nice=0.0
       smax=0.0
       nhet=0.0
       nlim_=0.0

       MONVOL=np*1.0e-6*ddry_ice*ddry_ice*ddry_ice
       ! ddry_ice    = 4.0e-8                !Droplet mean geometric diameter, m
       ! np_ice=nso4
       AUX1=1.6397e-14*T_ice-3.1769e-12      ! !from Barahona and Nenes JGR, 2008  (eq 32) !
       DPMAX=AUX1*(MONVOL**(-0.373))*(wpar_ice**(-0.05))  ! !from Barahona and Nenes JGR, 2008  (eq 32) !
       if (DPMAX.GT.1.0e-4) then
          DPMAX=1.0e-4  !from Barahona and Nenes JGR, 2008
       endif      
       DPMIN=dliq+(0.02/sqrt(alfa*wpar_ice*g1))  !MInimum size for DPmax (added on 09/11/09)
       DPMAX=max(DPMAX,DPMIN)   

       ! Barahona and Nenes JGR, 2008  (eq 25)
       AUX1=DPMAX-dliq
       AUX2=dlog((g2+(g1*DPMAX))/(g2+(g1*dliq)))     
       G=1.3346*((g1*AUX1)-(g2*AUX2))/(AUX1*g1*g1) ! Barahona and Nenes JGR, 2008  (eq 25)
       lambda_local=lambda/sqrt(wpar_ice)  ! (polydisperse eq 29)  !++sxj  bug 
       !write(*,*) "lambda_local",lambda_local,wpar_ice
       NSTAR=((g1*alfa*wpar_ice)**1.5)/beta/z/sqrt(2.0) ! (polydisperse eq 26)
       GAM=g2/g1 ! ???

       ! IS HOMOGENEOUS FREEZING HAPPENING?
       FDS=1.0 !CORRECTION TO Nc FROM HET FREEZING
       NHOM=0.
       if (T_ice .GE. 235.0) then  !no homogeneous above 235 K
          NHOM=0.0   
          goto 685
       endif
       if (purehet) then
          NHOM=0.0   
          goto 685
       endif 
       ! calculate limiting NIN for combined hom_het    
       if (typeofspec .ge. 0.0) then !polydisperse expressions
          call INSPEC_ice(shom, T_ice, NHET_, DSH,& 
                             shdust,ndust,effdust,kdust,si0dust,del0dust, &
                             shbc,nbc,effbc,kbc,si0bc,del0bc, & 
                             typeofspec,vpresi,vpresw)
          SIZECORR=EXP(-2.0/lambda_local/shom)
          DSTAR=((4.0*DSH*DSH/3.0)+(2.0*DSH*(shom-DSH))) /(shom-DSH+1.0)
!         write(0,*) 'czhao check nlim in icepara poly',nstar,shom,dstar,sizecorr
!         print*,'czhao check nlim in icepara poly',nstar,shom,dstar,sizecorr
          NLIM=NSTAR*(shom+1.0)/shom/sqrt(DSTAR)/SIZECORR     
       else !monodisperse approximation (eq 18) 
          DSH=shom-sh_mono
          DSTAR=((4.0*DSH*DSH/3.0)+(2.0*DSH*(shom-DSH)))/(shom-DSH+1.0)      
          DLIM=-GAM+sqrt((GAM*GAM)+(2.0*DSTAR/g1/alfa/wpar_ice))
!         write(0,*) 'czhao check nlim in icepara mono',alfa,wpar_ice,shom,z,beta,shom
!         print*, 'czhao check nlim in icepara mono',alfa,wpar_ice,shom,z,beta,shom
          NLIM=alfa*wpar_ice*(shom+1.0)/z/beta/shom
          NLIM=NLIM*((g1*DLIM)+g2)/DLIM/DLIM
          NHET_=nin_mono
       end if
       nlim_=NLIM
       FDS=1.0-((NHET_/NLIM)**1.5)  ! (monodisperse eq 28)
       if (purehom) then
          NHET_=0.0
          FDS=1.0
       end if
       if (FDS .LE. 0.) then  !Homogeneous nucleation completely inhibited by IN
          NHOM=0.   
          goto 685
       end if

       ! ********FRACTION OF FROZEN DROPLETS********************
       MIU=FDS*alfa*(shom+1.)*wpar_ice*koft/shom  !!!
       PHIDO=((pi_ice*G/MIU/2.)**0.5)*(G/MIU)
       AUXNC=2.*denair/beta/koft/denice/pi_ice/np
       FC=AUXNC/PHIDO
       ! calculating hom Nc
       if (np .gt.0.) then
          if (FC .le. 0.32) then
             NHOM=np*EXP(-FC)*(1.0-EXP(-FC)) ! From Barahona and Nenes, JGR, 2008
          else
             NHOM=np/(1.0+EXP(2.0*(1.0-FC)))   !correction needed for convective clouds (very high smax) (09/11/09)
          end if
       else  
          NHOM=0.0
       end if
       smax=shom
       nhet=NHET_
       goto 686  !finish        
       !PURE HETEROEGENEOUS FREEZING
       !find interval for bisection          
685    smax=0.0
       nhet=0.0
       SAUX=0.01 
       if (typeofspec .lt. 0.0) then
       ! typeofspec_ice    !Integer 1 to 4 see function NcSpec. -1 is monodisperse
          SAUX=sh_mono+0.00000000001 !minimun smax in monodisperse case
       end if
       F1=FINDSMAX(SAUX)   
       do INDEX =1, 20
          SAUX=SAUX+0.1
          F2=FINDSMAX(SAUX)      
          if (F2*F1 .lt.0.0) goto 677
          F2=F1
       end do 
       if (F2*F1 .gt.0.0) then
          nhet=0.0
          smax=SAUX    !No NIN present in pure heterogeneous mode smax>200%
          goto 686
       end if
       !Perform bisection
677    SUP=SAUX
       SLOW=SAUX-0.1
       do INDEX=1,100
          SHALF=0.5*(SUP+SLOW) 
          FHALF=FINDSMAX(SHALF)
          !++sxj
          !if (SIGN(1.0,F1)*SIGN(1.0,FHALF) .LE. 0.0) then 
          !if (  ( (F1.ge.0).and.(FHALF.le.0) ).or.( (F1.le.0).and.(FHALF.ge.0) )  ) then
          if (SIGN(1.0_r8,F1)*SIGN(1.0_r8,FHALF) .LE. 0.0) then
          !--sxj 
            F2    = FHALF
            SUP   = SHALF
          else
            F1    = FHALF
            SLOW   = SHALF
          endif
          if (ABS(SLOW-SUP) .LE. 1e-3) goto 678 
       enddo
678    smax=SHALF       
       if (typeofspec .ge. 0.0) then
          call INSPEC_ice(smax, T_ice, nhet, DSH,& 
                             shdust,ndust,effdust,kdust,si0dust,del0dust, &
                             shbc,nbc,effbc,kbc,si0bc,del0bc, & 
                             typeofspec,vpresi,vpresw)
       else
          nhet=nin_mono ! monodisperse approximation
       end if    
686    nice=NHOM+nhet             ! Adding Heterogeneous 
       !if (isnan(nice)) then     ! avoiding errors
       !   nice=0.0
       !end if
       !++sxj 
       if ((nice<0.).or.(nice>10e30)) nice=0. 
       !--sxj
     contains  
       real(r8) function FINDSMAX(SX)
          real(r8), intent(in) :: SX
          real(r8) :: tao
          if (typeofspec .ge. 0.0) then !polydisperse expressions
             call INSPEC_ice(SX, T_ice, NHET_, DSH,&
                             shdust,ndust,effdust,kdust,si0dust,del0dust, &
                             shbc,nbc,effbc,kbc,si0bc,del0bc, & 
                             typeofspec,vpresi,vpresw)
             SIZECORR=EXP(-2.0/lambda_local/SX) 
             DSTAR=((4.0*DSH*DSH/3.0)+(2.0*DSH*(SX-DSH)))/(SX-DSH+1.0)
             DSTAR=DSTAR+(gdoin*alfa*wpar_ice)
             tao=NHET_*SIZECORR*SX*sqrt(DSTAR)/(SX+1.0)/NSTAR
          else !monodisperse approximation
             DSH=SX-sh_mono
             DSTAR=((4.0*DSH*DSH/3.0)+(2.0*DSH*(SX-DSH)))/(SX-DSH+1.0)
             DLIM=-GAM+sqrt((GAM*GAM)+(2.0*DSTAR/g1/alfa/wpar_ice))
             tao=alfa*wpar_ice*(SX+1.0)/z/beta/SX
             tao=tao*((g1*DLIM)+g2)/DLIM/DLIM/nin_mono
          end if  
          !write(*,*) "findsmax:",Sx,tao
          FINDSMAX=1.0-tao
       end function FINDSMAX  
     end subroutine nice_param

     !------------------------------------------------------------------------------------------
     ! Function VPRESWATER. Calculates the saturated vapor pressure
     ! of water (Pa) according to Murphy & Koop (2005)
     ! T in K (173.15-373.15)
     !------------------------------------------------------------------------------------------

     real(r8) function VPRESWATER_ice(T)

       implicit none
       real(r8), intent(in) :: T
       real(r8) :: A(0:9) 
       data A /54.842763_r8, -6763.22_r8,   -4.21_r8, 0.000367_r8,   0.0415_r8, &
                   218.8_r8,   53.878_r8,-1331.22_r8, -9.44523_r8, 0.014025_r8/
       VPRESWATER_ice = A(0)+(A(1)/T)+(A(2)*dlog(T))+(A(3)*T) & 
                        +(tanh(A(4)*(T-A(5)))*((A(6)+(A(7)/T)) &
                        +(A(8)*dlog(T))+ (A(9)*T))) 
       VPRESWATER_ice= exp(VPRESWATER_ice)

     end function VPRESWATER_ice

     !------------------------------------------------------------------------------------------
     ! Function VPRESICE. Calculates the saturated vapor pressure
     ! of ice (pa) according to Murphy & Koop (2005)
     ! T in K (>110)
     !------------------------------------------------------------------------------------------

     real(r8) function VPRESICE(T)
       implicit none
       real(r8), intent(in) :: T
       real(r8) :: A(0:3)
       data A/9.550426_r8, -5723.265_r8, 3.53068_r8, -0.00728332_r8/
       VPRESICE = A(0)+(A(1)/T)+(A(2)*LOG(T))+(A(3)*T)
       VPRESICE = exp(VPRESICE)
     end function VPRESICE
     
     !------------------------------------------------------------------------------------------
     ! Function DHSUB. Calculates the latent heat of sublimation
     ! of ice (J/Kg) according to Murphy & Koop (2005)
     ! T in K (>30)
     !------------------------------------------------------------------------------------------

     real(r8) function DHSUB_ice(T)
       implicit none
       real(r8), intent(in) :: T
       real(r8)  :: A(0:4)
       DATA A/46782.5_r8, 35.8925_r8, -0.07414_r8, 541.5_r8, 123.75_r8/
       DHSUB_ice = A(0) + (A(1) * T) + (A(2)*T*T) + (A(3) * exp(-((T/ A(4))**2)))
       DHSUB_ice= 1000.*DHSUB_ice/18.
     end function DHSUB_ice
    
     !------------------------------------------------------------------------------------------
     ! Function ICEDENSITY. Calculates the DENSITY OF ICE
     ! of ice (Kg/m3) according to PK97 
     ! T in K (>30)
     !------------------------------------------------------------------------------------------

     real(r8) function DENSITYICE(T)
       implicit none
       real(r8), intent(in) :: T
       real(r8) :: A(0:2),  TTEMP
       DATA A/0.9167_r8, -1.75e-4_r8, -5.0e-7_r8/
       TTEMP=T-273
       DENSITYICE= 1000.*(A(0)+(A(1)*TTEMP)+(A(2)*TTEMP*TTEMP))
     end function DENSITYICE

     !------------------------------------------------------------------------------------------
     ! Function WATDENSITY. Calculates the DENSITY OF ICE
     ! of liquid water (Kg/m3) according to PK97 
     ! T in K (>240)
     !------------------------------------------------------------------------------------------

     real(r8) function WATDENSITY_ice(T)

       implicit none
       real(r8), intent(in) :: T
       real(r8) :: A(0:6),  TTEMP, WATDENSITY
       integer :: I
       DATA A / 0.99986_r8, 6.690e-5_r8, -8.486e-6_r8, 1.518e-7_r8, &  
               -6.9984e-9_r8, -3.6449e-10_r8, -7.497e-12_r8 /
       TTEMP=T-273.0
       if (TTEMP .le. -40.0) then
          TTEMP=-40.0
       end if
       WATDENSITY=A(6)*TTEMP 
       if (T .GE. 240.0) then 
          do I=5,1, -1
             WATDENSITY= (WATDENSITY+A(I))*(TTEMP)
          enddo
          WATDENSITY=WATDENSITY + A(0)
       else
          WATDENSITY=0.979
       end if 
       WATDENSITY=WATDENSITY*1000.0
       WATDENSITY_ice=WATDENSITY
     end function WATDENSITY_ice

     !------------------------------------------------------------------------------------------
     ! Subroutine PROPERTIES. Set physical an thermodynamic 
     ! properties at T and P 
     !------------------------------------------------------------------------------------------

     subroutine prop_ice(T_ice, P_ice, &
                         rv,dhs,vpresw,vpresi, &
                         denice,denwat,denair,diff, &
                         aircond,alfa,beta,shom, &
                         koft,dliq,g1,g2, &
                         doin,gdoin,z,lambda, &
                         kdust,si0dust,del0dust, &
                         Kbc,  si0bc,  del0bc  )
        implicit none
        real(r8), intent(in) :: T_ice, P_ice
        real(r8), intent(out):: rv,dhs,vpresw,vpresi
        real(r8), intent(out):: denice,denwat,denair,diff
        real(r8), intent(out):: aircond,alfa,beta,shom
        real(r8), intent(out):: koft,dliq,g1,g2
        real(r8), intent(out):: doin,gdoin,z,lambda
        real(r8), intent(out):: kdust,si0dust,del0dust
        real(r8), intent(out):: Kbc,  si0bc,  del0bc
        real(r8) ::  AUX 
        !real(r8) :: DHSUB_ice 
        !real(r8) :: VPRESWATER_ice, VPRESICE, DENSITYICE , WATDENSITY_ice
        !reak(r8) :: cubicint_ice
        real(r8) ::  AUX1, AUX2, SW, fice, mice, Tc, hdust, hbc
        real(r8) ::  b0, b1, b2, b3, x, T0bc, T0dust, gam, gamma
        integer i
        rv=rgas_ice/wmw_ice
        dhs=DHSUB_ice(T_ice)
       !DHSUB. Calculates the latent heat of sublimation
        vpresw=VPRESWATER_ice(T_ice)
       !VPRESWATER. Calculates the saturated vapor pressure of water (Pa)
        vpresi=VPRESICE(T_ice)
       !VPRESICE. Calculates the saturated vapor pressure of ice (pa)
        denice=DENSITYICE(T_ice)
       !ICEDENSITY. Calculates the DENSITY OF ICE  of ice (Kg/m3) 
        denwat=WATDENSITY_ice(T_ice)
       !WATDENSITY. Calculates the DENSITY OF ICE of liquid water (Kg/m3) 
        denair=P_ice*amw_ice/rgas_ice/T_ice
       !Kinetic properties of the bulk vapor (SI UNITS, Seinfel and Pandis, 1997)
        diff=(0.211*101325/P_ice)*((T_ice/273)**1.94)*1.e-4 !m^2/s
       ! diff_ice ::  water vapor diffusion coefficent
        AUX1=1.0e-3*(4.39+0.071*T_ice) !W/m
        ! thermal conductivity ::ka=10e-3*(4.39+0.071*T_in)
        !correcting Kair for size assuming D=1e-6 m  
        AUX2=(2.0*AUX1/(thaccom_ice*1.e-6*denair*cpa_ice))*((58.e-3*pi_ice/(rgas_ice*T_ice))**0.5)
        aircond=AUX1/(1.+AUX2)
        ! aricond :: thermal conductivity
        ! Physical constants
        AUX1=grav_ice*dhs/rv/T_ice/T_ice/cpa_ice
        AUX2=grav_ice*amw_ice/rgas_ice/T_ice
        alfa=AUX1-AUX2
        beta=amw_ice*P_ice/wmw_ice/vpresi
        gamma=1.5*dhs*dhs/rv/T_ice/T_ice/cpa_ice  !Correction for T>250 K
        beta=beta+gamma  !only for high T
        ! alfa beta (BN2008 eq.14)
        shom=2.349-(T_ice/259.) !threeshold Si according to Ren & McKenzie, 2005
        SW=shom*vpresi/vpresw
        shom=shom-1.
        koft=(0.0240*T_ice*T_ice)-(8.035*T_ice)+934.0 ! constant related to Jmax (Barahona and Nenes, 2008)
        !Calculate Dliq using an approximation derived from equilbrium calculations and the
        !expression proposed by Lewis (2008), 13, D03205, JGR 
        AUX1=(1.0/(1.0-SW))-1.1764
!       write(0,*) 'czhao check dliq',sw,aux1,ddry_ice,T_ice,shom,vpresi,vpresw
!       print*, 'czhao check dliq',sw,aux1,ddry_ice,T_ice,shom,vpresi,vpresw
        !czhao avoid negative AUX1
        if (aux1.lt.0) then 
        dliq=0.0
        else 
        dliq=ddry_ice*0.9344*(AUX1**(1./3.)) !
        endif
        ! ddry_ice    = 4.0e-8   !Droplet mean geometric diameter, m
        ! calculate average G    
        AUX1=denice*rv*T_ice/vpresi/diff
        AUX2=dhs*denice/aircond/T_ice
        AUX2=AUX2*((dhs/rv/T_ice)-1.0)
        g1=(AUX1+AUX2)/4.0
        g2=denice*rv*T_ice/2.0/vpresi/depcoef_ice
        g2=g2*((2.0*pi_ice/rv/T_ice)**0.5)   
        ! g1 g2 (ice growth, BN2008 eq.6)
        doin=1.0e-6 !assumed IN diameter
        gdoin=(g1*0.5*doin*doin)+(g2*doin)
        z=denice*pi_ice/2.0/denair
        gam=g2/g1
        lambda=1.0/sqrt(alfa*g1*gam*gam)  !divided by sqrt(wparcel) in niceparam 
        !!============Parameters needed for IN spectra=========
        !for CNT spectra !!!!!!!!!!!!!!!!!!!!!
        mice = 0.96 !compatibility parameter dust
        fice=0.25*((mice*mice*mice)-(3*mice)+2)
        kdust=koft*fice
        mice = 0.76 !compatibility parameter black carbon
        fice=0.25*((mice*mice*mice)-(3*mice)+2)
        kbc=koft*fice
        ! for Phillips, et. al. 2008 spectra PDA08!!!!!!!!!!!!!!!!!!!!!
        Tc=T_ice-273
        hdust=0.15
        T0dust=-40.
        b0=-1.0261
        b1=3.1656e-3
        b2=5.3938e-4
        b3=8.2584e-6  
        x=b0+(b1*Tc)+(b2*Tc*Tc)+(b3*Tc*Tc*Tc)
        si0dust=1.+(10.**x)
        del0dust=cubicint_ice(Tc, T0dust, T0dust+5, 1.0_r8, hdust)
        hbc=0.
        T0bc=-50.
        b0=0.5652 
        b1=1.085e-2
        b2=-3.118e-5
        si0bc=b0+(b1*Tc)+(b2*Tc*Tc)
        del0bc=cubicint_ice(Tc, T0bc, T0bc+5, 1.0_r8, hbc)
    end subroutine prop_ice 
  
!*************************************************************
!   Subroutine gauspdf (normalized over the width of the distribution).  
!************************************************************       
     subroutine gausspdf(x, dp, &
                         miuv,sigmav,normv )
       implicit none 
       real(r8), intent(in) :: x
       real(r8), intent(in) :: miuv,sigmav,normv 
       real(r8), intent(out) :: dp
       dp=EXP(-0.5*(x-miuv)*(x-miuv)/sigmav/sigmav) &
          /sigmav/sqrt(2.*pi_ice)/normv 
     endsubroutine gausspdf
     
!*************************************************************
!   Function cubicint_ice (cubic interpolation between y1 and y2 within a and b).  
!************************************************************  
    real(r8) function cubicint_ice(y, y1, y2, a, b)
       implicit none
       real(r8), intent(in) :: y, y1, y2, a, b   
       real(r8) :: A_, B_, a0, a1, a2, a3, d, AUX
       if (y .le. y1) then
          d=a
          goto 5065
       end if 
       if (y .ge. y2) then
          d=b
          goto 5065
       end if 
       AUX=y2-y1      
       A_=6*(a-b)/(AUX*AUX*AUX)
       B_=a+(A_*(y1*y1*y1)/6.)-(A_*(y1*y1)*y2*0.5)
       a0=B_
       a1=A_*y1*y2
       a2=-A_*(y1+y2)*0.5
       a3=A_/3.
       d=a0+(a1*y)+(a2*y*y)+(a3*y*y*y)    
5065   cubicint_ice=d
    end function cubicint_ice

!*************************************************************
!   Function dcubicint_ice (used in the PDA08 spectrum).  
!************************************************************  
     real(r8) function dcubicint_ice(y, y1, y2, a, b)
       implicit none
       real(r8), intent(in) :: y, y1, y2, a, b   
       real(r8) :: A_, a0, a1, a2, a3, d, AUX
       if (y .le. y1) then
          d=0
          goto 5066
       end if 
       if (y .ge. y2) then
          d=0
          goto 5066
       end if 
       AUX=y2-y1      
       A_=6*(a-b)/(AUX*AUX*AUX)
       a1=A_*y1*y2
       a2=-A_*(y1+y2)*0.5
       a3=A_/3.
       d=(a1)+(2*a2*y)+(3*a3*y*y)     
5066   dcubicint_ice=d
     end function dcubicint_ice 
     
!*************************************************************
! Function PDG07 (simplified ice nucleation 
!                     spectra according to Phillips et. al. 2007).  
! si is supersaturation wrt ice and T is in K 
!************************************************************ 
     real(r8) function PDG07_ice(si, T)     
       real(r8), intent(in) :: si, T
       real(r8) :: N 
       if (T .le. 243)then
         N=1000*exp(-0.388)*(exp(3.88*si)-1)
       else
         N=60*exp(-0.639)*(exp(12.96*si)-1)     
       end if
       PDG07_ice=N 
     end function PDG07_ice

!*************************************************************
! Subroutine INSPEC_ice
!  Provides the Ice Nuclei concentration (m-3) 
! and the characteristic freezing threeshold, DSh (Barahona & Nenes 2009), at given 
! si and T. The variable typeofspec_ice (integer) has the values
! 1 Meyers et. al. 1992
! 2  Phillips et. al. 2007
! 3  CNT Barahona & Nenes 2009
! 4  Phillips et. al. 2008 (simplifed) 
! si is supersaturation wrt ice and T is in K 
!************************************************************   
     subroutine INSPEC_ice(si, T, N, Dsh,&
                           shdust,ndust,effdust,kdust,si0dust,del0dust, &
                           shbc,nbc,effbc,kbc,si0bc,del0bc,&
                           typeofspec,vpresi,vpresw)
       implicit none
       real(r8), intent(in) :: si, T
       real(r8), intent(in) :: shdust,ndust,effdust,kdust,si0dust,del0dust
       real(r8), intent(in) :: shbc,nbc,effbc,kbc,si0bc,del0bc
       integer , intent(in) :: typeofspec
       real(r8), intent(in) :: vpresi,vpresw
       real(r8), intent(out) :: N, Dsh
       !real(r8) :: PDG07_ice
       real(r8) :: Nd, Nbc_, aux, Si_, SW, del1, ddel1 
       real(r8) :: fc, delw0, ddelw0, SW0, Hdust, Hbc 
       real(r8) :: Nbase, dNd, dNbc, dNbase, dH
       real(r8) ::  dfc
       !real(r8) :: dcubicint_ice,cubicint_ice
       select case  (typeofspec)
       case(1) 
          N=1000*exp(-0.639)*(exp(12.96*si)-1)
          Dsh=1/12.96
       case(2)
          N=PDG07_ice(si, T)
          if (T .le. 243)then
             Dsh=1/3.88
          else
             Dsh=1/12.96
          end if
       case(3)
          !dust contribution
          if (si .le.shdust)    then 
             Nd=(si/shdust)*ndust*effdust*exp(-kdust*(shdust-si))
             dNd=Nd*((1/si)+kdust)
          else
             Nd=ndust*effdust
             dNd=0.
          end if
          !soot/bc contribution
          if (si .le.shbc) then  
             Nbc_=(si/shbc)*nbc*effbc*exp(-kbc*(shbc-si))
             dNbc=Nbc_*((1/si)+kbc)
          else
             Nbc_=nbc*effbc
             dNbc=0
          end if
          N=Nd+Nbc_
          if ((dNd+dNbc) .gt. 0.) then
             Dsh=N/(dNd+dNbc)
          else
             Dsh=si
          end if 
       case(4)
          Si_=si+1
          SW=Si_*vpresi/vpresw
          SW0=0.97
          delw0=cubicint_ice(SW, SW0, 1.0_r8, 0.0_r8, 1.0_r8)
          ddelw0=dcubicint_ice(SW, SW0, 1.0_r8, 0.0_r8, 1.0_r8)
          Nbase=PDG07_ice(si, T)
          if (T .le. 243)then
             dNbase=3.88*Nbase
          else
             dNbase=12.96*Nbase
          end if
          !dust contribution
          del1=cubicint_ice(Si_, si0dust, si0dust+0.1, 0.0_r8, 1.0_r8)
          ddel1=dcubicint_ice(Si_, si0dust, si0dust+0.1, 0.0_r8, 1.0_r8)
          fc=0.5*del0dust*del1
          dfc=0.5*del0dust*ddel1
          Hdust=fc+((1.-fc)*delw0) 
          dH=(dfc*(1.-delw0))+(ddelw0*(1.-fc))
          if (Hdust .gt. 1.) then 
             Hdust=1.
             dH=0.
          end if
          !aux=exp(-2.*Hdust*Nbase/3./7.92e4) !++sxj 7.92e4 needs to be changed to 3.3e5 
          !Phillips et al 2008  BASE_DUST_OMEGA = 2.0e-6 (in table 2 ) changed from 5.0e-7
          aux=exp(-2.*Hdust*Nbase/3./3.3e5)
          Nd=ndust*(1.-aux)
          !dNd=(ndust-Nd)*((dH*Nbase)+(Hdust*dNbase))*2./3./7.92e4  !++sxj 7.92e4 needs to be changed to 3.3e5  
          !Phillips et al 2008  BASE_DUST_OMEGA = 2.0e-6 (in table 2 ) changed from 5.0e-7
          dNd=(ndust-Nd)*((dH*Nbase)+(Hdust*dNbase))*2./3./3.3e5
          !soot/bc contribution
          del1=cubicint_ice(Si_, si0bc, si0bc+0.1, 0.0_r8, 1.0_r8)
          ddel1=dcubicint_ice(Si_, si0bc, si0bc+0.1, 0.0_r8, 1.0_r8)
          fc=0.5*del0bc*del1
          Hbc=fc+((1.-fc)*delw0)
          !++sxj
          !dfc=0.5*del0dust*ddel1  !!!!!!! bug  "del0dust" should be revised to "del0bc"
          dfc=0.5*del0bc*ddel1 
          !write(*,*) "del0dust-del0bc",del0dust,del0bc
          !--sxj
          dH=(dfc*(1.-delw0))+(ddelw0*(1.-fc))
          if (Hbc .gt. 1.) then 
             Hbc=1.
             dH=0.
          end if
          aux=exp(-Hbc*Nbase/3./1.04e6) 
          Nbc_=nbc*(1.-aux)
          dNbc=(nbc-Nbc_)*((dH*Nbase)+(Hbc*dNbase))*1./3./1.04e6
          N=Nd+Nbc_
          if ((dNd+dNbc) .gt. 0.) then
             Dsh=N/(dNd+dNbc)
          else
             Dsh=si
          end if 
       case default 
          N=0.
          Dsh=si
       end select
       if (Dsh .ge. si) then
          Dsh=si
       end if
    end subroutine INSPEC_ice  

   endmodule module_icenucleation_bn

    


 
   

  
