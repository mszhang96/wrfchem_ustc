MODULE module_ctrans_kfeta

!************************************************************************
! This computer software was prepared by Battelle Memorial Institute,
! hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with
! the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
! CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
! LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to Compute Convective Transport following KFeta scheme 
! * Author: Chun Zhao 
! Last update: July 2011
!
! Contact:
! Chun Zhao, PhD
! Scientist
! Pacific Northwest National Laboratory
! P.O. Box 999, MSIN K9-30
! Richland, WA, 99352
! Phone: (509) 371-6372
! Email: Chun.Zhao@pnnl.gov
!
! Please report any bugs or problems to Chun Zhao
!
! Terms of Use:
!  1) Users are requested to consult the primary author prior to
!     modifying this module or incorporating it or its submodules in
!     another code. This is meant to ensure that the any linkages and/or
!     assumptions will not adversely affect the operation of this module.
!  2) The source code in this module is intended for research and
!     educational purposes. Users are requested to contact the primary
!     author regarding the use of the code for any commercial application.
!  3) Users preparing publications resulting from the usage of this code
!     are requested to cite one or more of the references below
!     (depending on the application) for proper acknowledgement.
!
! References: 
! * Zhao, C., Y. Wang, Y. Choi, T. Zeng: 
!   Summertime impact of convective transport and lightning NOx production
!   over North America: modeling dependence on meteorological simulations,
!   ACP, 9, 4315-4327, 2009. 
!
! Additional information:
! *  www.pnl.gov/atmospheric/research/wrf-chem
!
!************************************************************************
!===========================================================================
!CUMULUS CONVECTION PARAMETERIZATION USING NEW-KFETA SCHEME 
!
! New Kain-Fritsch Scheme: The scheme is based on Kain and Fritsh(1990) and 
! Kain and Fritsch(1993), but has been modified based on testing within the Eta
! model. As with the original KF scheme, it utilizes a simple cloud model with 
! moist updrafts and downdrafts, including the effects of detrainment, 
! entrainment, and relatively simple microphysics. It differs from the original
! KF scheme in the following ways:
!   1.A minimum entrainment rate is imposed to suppress widespread convection
!   in marginally unstable, relatively dry environments.
!   2.Shallow (non-precipitating) convection is allowed for any updraft that 
!   does not reach minimum cloud depth for precipitating clouds; this minimum 
!   depth varies as a function of cloud-base temperature.
!   3.The entrainment rate is allowed to vary as a function of low-level
!   convergence.
!   4.Downdraft changes:
!     ---Source layer is the entire 150-200 mb deep layer just above cloud base.
!     ---Mass flux is specified as a fraction of updraft mass flux at cloud
!     base. Fraction is a function of source layer RH rather than wind shear 
!     or other parameters, i.e., old precipitation efficiency relationship not
!     used.
!     ---Detrainment is specified to occur in updraft source layer and below.
!                         
!                                                   ----czhao  03/2006 
!IN THIS SUBROUTINE, IT INCLUDES FOLLOWING PARTS:
!  + CUMULUS SCAVENGING FOR SOLUBLE GASES  
!  + NEW_KFETA SCHEME
!                                         ---czhao    03/2006
!=============================================================================

USE module_dep_simple, only: is_aerosol
USE module_mosaic_wetscav, only:E_ICE,COMPUTE_L2G
USE module_state_description, only:param_first_scalar,   &
                p_so2, p_h2o2, p_sulf, p_msa,p_msa_a02,p_msa_cw02,   &
                p_hno3, p_hcl, p_nh3,p_hcho,  &
                RADM2SORG_AQ,RACMSORG_AQ,CBMZSORG_AQ,RADM2SORG,  &
                CBMZ_MOSAIC_4BIN_AQ,CBMZ_MOSAIC_8BIN_AQ, &
                CBMZ_MOSAIC_DMS_4BIN_AQ,CBMZ_MOSAIC_DMS_8BIN_AQ, &
                CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE,CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE, &
                CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN 

CONTAINS

!-------------------------------------------------------------
   SUBROUTINE kfeta_cup_transport(DT,itimestep,DX,           &
              rho_phy,alt,chem,                                  &
              U,V,t_phy,moist,dz8w,                          &
              p_phy,z,                          &
              UMF3D,UER3D,UDR3D,DMF3D,DER3D,DDR3D,WU3D,      & 
              LDB2D,LET2D,KSL2D,LFS2D,KLC2D,                 & 
              CUTOP,CUBOT,NCA2D,ISHALL2D,CHECK2D,  &
              num_moist,numgas,num_chem,chemopt,scavopt,     &
              ids,ide, jds,jde, kds,kde,                            &
              ims,ime, jms,jme, kms,kme,                            &
              its,ite, jts,jte, kts,kte                             )

!-------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------
! UMF3D      kg  s-1   CLOUD BASE MASS FLUX (UPDRAFT)
! DMF3D      kg  s-1   DOWNDRAFT FLUX(negative value)
! UER3D      kg  s-1   ENTRAINMENT RATE (UPDRAFT)
! UDR3D      kg  s-1   DETRAINMENT RATE (UPDRAFT)
! DER3D      kg  s-1   ENTRAINMENT RATE (DOWNDRAFT,negative)
! DDR3D      kg  s-1   DETRAINMENT RATE (DOWNDRAFT)
! WU3D       m   s-1   UPDRAFT VELOCITY  
! LDB2D      unitless     DOWNDRAFT END LEVEL
! LFS2D      unitless     DOWNDRAFT ORIGINATING LEVEL
! KSL2D      unitless     UPDRAFT SOURCE TOP LAYER
! KLC2D      unitless     UPDRAFT SOURCE BOT LAYER
! LET2D      unitless     TOPEST CLOUD DETRAINMENT LAYER
! CUTOP      unitless     CLOUD TOP LEVEL
! CUBOT      unitless     CLOUD BOT LEVEL
! ISHALL2D   unitless     SWITCH FOR DEEP OR SHALLOW CONVECTION
! CHECK2D    unitless     check if deep convection is really occur
! NCA2D      unitless     convection relaxation counter 
!-------------------------------------------------------------

   INTEGER,      INTENT(IN   ) ::                               &
                                  numgas,chemopt,scavopt,       &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte,    &
                                  num_chem,num_moist   
   INTEGER,      INTENT(IN   ) :: ITIMESTEP
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme,num_moist )         ,    &
          INTENT(IN   ) ::                              moist
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                         U,    &   
                         V,    &
                     t_phy,    &
                     z,        &       !heights of model levels in m
                     p_phy,    &
                      dz8w,    &       !layer thickness in m
                      alt,     &       !inversed dry air density in m3/kg
                   rho_phy             !moisture air density in kg/m3
   REAL, INTENT(IN   ) :: DT, DX
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme, num_chem ),    &
         INTENT(INOUT) ::               chem

   REAL, DIMENSION (ims:ime,kms:kme,jms:jme), INTENT(IN) ::    &
     UMF3D,UER3D,UDR3D,DMF3D,DER3D,DDR3D,WU3D
   INTEGER, DIMENSION(ims:ime,jms:jme),INTENT(IN) ::    &
     LDB2D,LET2D,KSL2D,LFS2D,KLC2D,CUTOP,CUBOT,ISHALL2D,CHECK2D

   REAL, DIMENSION(ims:ime,jms:jme),INTENT(IN) :: NCA2D

!  local variable

   REAL, DIMENSION( ims:ime , kms:kme , jms:jme) ::  & 
         con_trans_ten    !tendency difference due to cumulus convection

   REAL DTIN
   REAL,DIMENSION(kms:kme)  ::  AIRDENS_GAS   ! in mole/m3
   REAL,DIMENSION(kms:kme)  ::  AIRDENS_AER   ! in kg/m3
   REAL,DIMENSION(kms:kme)  ::  AIRDENS   ! in kg/m3
   REAL,DIMENSION(kms:kme)  ::  air_ratio   ! dry_air to moisture_air ratio 
!  REAL,DIMENSION(kms:kme)  ::     &
!    FLUXUP0,FLUXDOWN0,UENTR0,UDETR0,DENTR0,DDETR0
   REAL,DIMENSION(kms:kme)  ::     &
     FLUXUP,FLUXDOWN,UENTR,UDETR,DENTR,DDETR,DIFMASS,SUBMASS
   REAL,DIMENSION(kms:kme)  ::    TRACUP,TRACDW,TRACSUB1,TRACSUB2
   REAL   TRACSUM
   REAL UAIRMASS,UENMASS,UDEMASS,DENMASS,DDEMASS,TRACERMASS1,TRACERMASS2, &
    TRACERMASS,VOLUME,DENSITY,AIRMASS,AIRMASS1,AIRMASS2,FLUXTRA,FLUXRAT
   REAL,DIMENSION(kms:kme)  :: XS0  !ORIGINAL MIXING RATIO
   REAL,DIMENSION(kms:kme)  :: XS1  !NEW MIXING RATIO AFTER TRANSPORT
   REAL,DIMENSION(kms:kme)  :: TXS  !TEMPORARY MIXING RATIO
   REAL, PARAMETER :: AVOGADRO = 6.0225E23  ! molecules/mole
   REAL, PARAMETER :: AIRMW = 28.98  ! g/mole
   INTEGER I, J, K,  nv 
   INTEGER ZTOP              ! CLOUD TOP LAYER
   INTEGER ZB                ! CLOUD BOT LAYER
   INTEGER ZLC               ! UPDRAFT ORIGINATING LAYER
   INTEGER ZLDB              ! DOWNDRAFT END LAYER   
   INTEGER ZLFS              ! DOWNDRAFT ORIGINATING LAYER
   INTEGER ZLET              ! TOPEST CLOUD DETRAINMENT LAYER
   INTEGER ZKSL              ! THE UPDRAFT SOURCE TOP LAYER 
   INTEGER ISHALL            ! SWITCH FOR DEEP OR SHALLOW CONVECTION
   REAL   SUMASS1,SUMASS2,TOTRACER  ! CHECK MASS BALANCE 
   REAL,DIMENSION(kms:kme,1:2)  :: SCAVF  !updraft and downdraft scavenging coefficient
   REAL DX2
   REAL,DIMENSION(kms:kme)  :: DZQ,TMP,PMID
   REAL,DIMENSION(kms:kme)  :: VUD,VDD 
   REAL,DIMENSION(kms:kme,num_chem)  :: sol_f 
   REAL,DIMENSION(num_chem) ::  sol_f0 
   REAL    :: epsilc
   epsilc=1.e-30

      DTIN=DT
      dx2=dx*dx

!       con_trans_ten(its:ite,kts:kte,jts:jte)=0. ! tendency difference due to cumulus convection

      DO 1000 J=jts,jte
      DO 1000 I=its,ite

          !Obtain the soluble information for scavenging 
          do k=kts,kte
             sol_f(k,1:num_chem)=1.0
             sol_f0(1:num_chem)=1.0
           soluble_select: SELECT CASE(chemopt)
             CASE (RADM2SORG_AQ,RACMSORG_AQ,CBMZSORG_AQ,RADM2SORG)
                call sorgam_soluble(sol_f0(1:num_chem),num_chem,chem(i,k,j,1:num_chem))
                sol_f(k,1:num_chem)=sol_f0(1:num_chem)
             CASE (CBMZ_MOSAIC_4BIN_AQ,CBMZ_MOSAIC_8BIN_AQ, &
                   CBMZ_MOSAIC_DMS_4BIN_AQ,CBMZ_MOSAIC_DMS_8BIN_AQ, &
                   CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE,CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE, &
                   CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN)
                call mosaic_soluble(sol_f0(1:num_chem),num_chem,chem(i,k,j,1:num_chem))
                sol_f(k,1:num_chem)=sol_f0(1:num_chem)
             CASE DEFAULT
           END SELECT soluble_select
          enddo !k

      DO 900 nv=2,num_chem

         do k=kts,kte
         if(.not.is_aerosol(nv))then !gases in ppmv
         xs0(k)=1.e-6*max(epsilc,chem(i,k,j,nv))  ! ppmv to mole/mole
         xs1(k)=xs0(k)
         txs(k)=xs1(k)
         else !aerosol in ug/kg dry-air 
         xs0(k)=max(epsilc,chem(i,k,j,nv))
         xs1(k)=xs0(k)
         txs(k)=xs1(k)
         endif
         airdens_gas(k)=rho_phy(i,k,j)*1000.0/airmw !kg/m3 to mole/m3
         airdens_aer(k)=1.0/alt(i,k,j) !dry air density kg/m3
         air_ratio(k)=(1.0/alt(i,k,j))/rho_phy(i,k,j) !dry air to moist air ratio
         dzq(k)=dz8w(i,k,j)
         enddo

         !----------------------------------------------------------------   
         !BEFORE CONVECTION,CHECK TRACER MASS BALANCE
             SUMASS1=0.
             DO K=kts,kte
              if(.not.is_aerosol(nv))then !gases in mole/mole 
               SUMASS1=SUMASS1+dzq(k)*airdens_gas(k)*DX2*TXS(K)
               IF (SUMASS1 .GT. 1.e20) then  ! larger than 1.e10 mole
               WRITE(0,*) 'TRACER MASS > 1.e20 mole',nv,I,J,K
               PRINT *, 'TRACER MASS > 1.e20 mole',nv,I,J,K
               WRITE(0,*) SUMASS1
               PRINT *, SUMASS1
               ENDIF 
              else   ! aerosol in ug/kg-dryair 
               SUMASS1=SUMASS1+dzq(k)*airdens_aer(k)*DX2*TXS(K)
               IF (SUMASS1 .GT. 1.e30) then  ! larger than 1.e30 ug 
               WRITE(0,*) 'TRACER MASS > 1.e30 ug',nv,I,J,K
               PRINT *, 'TRACER MASS > 1.e30 ug',nv,I,J,K
               WRITE(0,*) SUMASS1
               PRINT *, SUMASS1
               ENDIF
              endif
             ENDDO

   !----------------------------------------------------------------   
   !0 Pre-initialize 
   !----------------------------------------------------------------   
             do k=kts,kte
             if(.not.is_aerosol(nv))then !for gas, in mole/m3
             AIRDENS(K)=AIRDENS_GAS(K)
             else  !for aerosol, in kg/m3 dry-air 
             AIRDENS(K)=AIRDENS_AER(K)
             endif
             TMP(K)=T_PHY(I,K,J)
             PMID(K)=P_PHY(I,K,J)
             enddo

             do k=kts,kte
             SCAVF(k,1:2)=0.0
             TRACUP(k)=0.
             TRACDW(k)=0.
             TRACSUB1(k)=0.
             TRACSUB2(k)=0.
             SUBMASS(k)=0.
             enddo

             ZTOP=CUTOP(I,J)
             ZB=CUBOT(I,J)
             ZLC=KLC2D(I,J)
             ZLFS =LFS2D(I,J)
             ZLDB =LDB2D(I,J)
             ZLET =LET2D(I,J)
             ZKSL= KSL2D(I,J)
             ISHALL=ISHALL2D(I,J)

             do k=kts,kte
              if(.not.is_aerosol(nv))then !gases in mole/mole 
              FLUXUP(k)=UMF3D(I,k,j)*1.e3/AIRMW !kg/s to mole/s
              FLUXDOWN(k)=-1.0*DMF3D(I,k,j)*1.e3/AIRMW !kg/s to mole/s
              !IF (FLUXDOWN(K).LT.0) FLUXDOWN(K)=0.0
              UENTR(k)=UER3D(I,k,j)*1.e3/AIRMW !kg/s to mole/s
              UDETR(k)=UDR3D(I,k,j)*1.e3/AIRMW !kg/s to mole/s
              DENTR(k)=-1.0*DER3D(I,k,j)*1.e3/AIRMW !kg/s to mole/s
              DDETR(k)=DDR3D(I,k,j)*1.e3/AIRMW !kg/s to mole/s
              else 
              FLUXUP(k)=UMF3D(I,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
              FLUXDOWN(k)=-1.0*DMF3D(I,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
              !IF (FLUXDOWN(K).LT.0) FLUXDOWN(K)=0.0
              UENTR(k)=UER3D(I,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
              UDETR(k)=UDR3D(I,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
              DENTR(k)=-1.0*DER3D(I,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
              DDETR(k)=DDR3D(I,k,j)*air_ratio(k) !kg/s air to kg/s dry-air
              endif
              IF (FLUXUP(K).GT.0) THEN 
              VUD(K)=MIN(50.0,MAX(WU3D(I,K,J),0.0))
              VDD(K)=MIN(50.0,MAX(VUD(K)*FLUXDOWN(k)/FLUXUP(k),0.0))
              ELSE 
              VUD(K)=0.0
              VDD(K)=0.0
              ENDIF
             enddo

             !make convective output consistent
             IF (MAXVAL(ABS(FLUXUP(kts:kte)))==0 ) THEN
               ZTOP=0 
               ZB=ZTOP
             ENDIF
             IF (MAXVAL(ABS(FLUXDOWN(kts:kte)))==0) THEN
               ZLFS=0
               ZLDB=ZLFS
             ENDIF

!------------------------------------------------------------------
!...DO NOT ALLOW ANY CLOUD FROM THIS LAYER IF:
!
!...            1.) if there is no CAPE, or 
!...            2.) cloud top is at model level just above LCL, or
!...            3.) cloud top is within updraft source layer, or
!...            4.) cloud-top detrainment layer begins within 
!...                updraft source layer.
!------------------------------------------------------------------

   IF((ZB.LT.ZTOP).AND.(CHECK2D(I,J).NE.0))THEN !ALLOWED CONVECTION BLOCK

   !================================================================   
   !1 Start deep convective transport 
   !================================================================   
   !note:
   !ishall: SWITCH FOR DEEP CONVECTION OR SHALLOW CONVECTION
   !FOR DEEP AND SHALLOW CONVECTION, THE DIFFERENCE IS JUST FOR DOWNDRAFT
   !FLUX. THE CALCULATION OF UPDRAFT IS SAME. SHALLOW CONVECTION HAS NO DOWNDRAFT.
   !----------------------------------------------------------------   

     IF (ISHALL==0) THEN  

     !----------------------------------------------------------------   
     !1.0 Calculate the scavenging coefficient for deep convective transport 
     !----------------------------------------------------------------   
         IF (scavopt==1) THEN 

         CALL CUSCAV_KF(I,J,nv,dzq,dx2,airdens,vud,vdd,tmp,PMID,         &
                        num_chem,chem(i,kms:kme,j,1:num_chem),SCAVF,  &
                        ids,ide, jds,jde, kds,kde,                       &
                        ims,ime, jms,jme, kms,kme,                       &
                        its,ite, jts,jte, kts,kte                        ) 
          !apply the internal-mixing soluble fraction for scavenging aerosol 
          if(is_aerosol(nv))then
          do k=kts,kte
           SCAVF(k,1:2)=SCAVF(k,1:2)*sol_f(k,nv)
          enddo
          endif
          !set downdraft scavenging zero for now
          do k=kts,kte
           SCAVF(k,2)=0.0
          enddo
         ENDIF

     !----------------------------------------------------------------   
     !1.1 COMPUTE UPDRAFT EFFECT IN DEEP CONVECTION                  
     !----------------------------------------------------------------   
         do k=kts,kte
         DIFMASS(k) = 0.0
         SUBMASS(k)=0.
         enddo

         DO K=ZLC+1,ZTOP-1
          IF (FLUXUP(K) .EQ.0 ) THEN
            WRITE(0,*) 'WRONG FLUXUP',I,J,K,ZTOP,ZLC
            PRINT *, 'WRONG FLUXUP',I,J,K,ZTOP,ZLC
            STOP
          ENDIF
          !fluxup ajust for mass consistence
          UENMASS=UENTR(K)*DTIN ! mole of entrain air
          UDEMASS=UDETR(K)*DTIN ! mole of detrain air
  
          IF(UENMASS.NE.UDEMASS) THEN 
          IF((FLUXUP(K)-FLUXUP(K-1))/(UENMASS-UDEMASS).GT.1.e-5.AND.  &
             (FLUXUP(K)-FLUXUP(K-1))/(UENMASS-UDEMASS).LT.0 ) THEN
                WRITE(0,*) 'UPDRAFT AJUSTMENT AT',I,J,K
                PRINT *, 'UPDRAFT AJUSTMENT AT',I,J,K
                WRITE(0,*) FLUXUP(K),FLUXUP(K+1),UENMASS,UDEMASS
                PRINT *, FLUXUP(K),FLUXUP(K+1),UENMASS,UDEMASS
                FLUXUP(K)=FLUXUP(K-1)+UENMASS-UDEMASS
          ENDIF
          ENDIF
         ENDDO !K

     !----------------------------------------------------------------   
     !1.1.1 calculate tracer tendence in LC layer
     !----------------------------------------------------------------   
         !gases in mole/mole or aerosol in ug/kg-dry air
         FLUXTRA=FLUXUP(ZLC)*TXS(ZLC)*DTIN  ! mole or ug of tracer
         UAIRMASS=FLUXUP(ZLC)*DTIN  ! mole or kg of air mass
         VOLUME=dzq(zlc)*DX2     !m3
         DENSITY=airdens(k)  !mole or kg /m3
         TRACUP(ZLC)=-FLUXTRA
         AIRMASS=VOLUME*DENSITY-UAIRMASS !mole or kg
         DIFMASS(ZLC)=-UAIRMASS  !save for subsidence
         FLUXRAT=TXS(ZLC)*(1.0-scavf(zlc,1)) ! mole/mole or ug/kg dry-air

     !----------------------------------------------------------------   
     !1.1.2 calculate other under cloud top layers
     !----------------------------------------------------------------   
         DO K=ZLC+1,ZTOP-1
           UAIRMASS=FLUXUP(K-1)*DTIN ! mole or kg of air mass from lower level
           UENMASS=UENTR(K)*DTIN ! mole or kg of entrain air
           UDEMASS=UDETR(K)*DTIN ! mole or kg of detrain air
           VOLUME=DZQ(K)*DX2 ! m3
           DENSITY=airdens(K) ! mole or kg /m3
           AIRMASS=VOLUME*DENSITY-UENMASS+UDEMASS  !mole or kg
           DIFMASS(K)=-UENMASS+UDEMASS  !save for subsidence
           FLUXTRA=UAIRMASS*FLUXRAT+UENMASS*TXS(K)*(1.-scavf(k,1))-UDEMASS*FLUXRAT !mole gas or ug aerosol
           TRACUP(K)=UDEMASS*FLUXRAT-UENMASS*TXS(K) !mole gas or ug aerosol
           FLUXRAT=FLUXTRA*(1.-scavf(k,1))/(FLUXUP(K)*DTIN) !mole/mole or ug/kg dry-air
         ENDDO

     !----------------------------------------------------------------   
     !1.1.3 calculate the top layer of cloud, no entrainment
     !----------------------------------------------------------------   
           UAIRMASS=FLUXUP(ZTOP-1)*DTIN  !mole or kg airmass from lower level
           VOLUME=DZQ(ZTOP)*DX2 
           DENSITY=airdens(ZTOP)
           AIRMASS=VOLUME*DENSITY+UAIRMASS ! mole or kg
           DIFMASS(ZTOP)=UAIRMASS !mole or kg dry-air
           SUBMASS(ZTOP)=DIFMASS(ZTOP) !mole or kg
           TRACUP(ZTOP)=UAIRMASS*FLUXRAT !mole gas or ug aerosol

     !----------------------------------------------------------------   
     !1.1.4 Check mass balance for updraft
     !----------------------------------------------------------------   
           TRACSUM=0.0
           DO K=ZLC,ZTOP
            TRACSUM=TRACSUM+TRACUP(K)
           ENDDO

           IF( (TRACSUM/SUMASS1) .GT. 1.e-5) THEN
            WRITE(0,*) 'TRACER MASS IS NOT BALANCE DURING UPDRAFT!!!'
            PRINT *, 'TRACER MASS IS NOT BALANCE DURING UPDRAFT!!!'
            WRITE(0,*) TRACSUM,I,J,nv
            PRINT *, TRACSUM,I,J,nv
            WRITE(0,*) TRACSUM/SUMASS1
            PRINT *, TRACSUM/SUMASS1
            STOP
           ENDIF

!          DO K=KTS,KTE
!           VOLUME=DZQ(K)*DX2
!           DENSITY=airdens(K)
!           AIRMASS=VOLUME*DENSITY+DIFMASS(K) ! mole or kg air
!           TRACERMASS=VOLUME*DENSITY*TXS(K)+TRACUP(K) ! mole gas or ug aerosol
!          ENDDO

     !------------------------------------------------------------
     !1.1.5 Calculate tracer concentration due to updraft in deep convection
     !------------------------------------------------------------
          DO K=ZLC,ZTOP
            !gases in mole/mole and aeroso in ug/kg dry-air 
            VOLUME=DZQ(K)*DX2
            DENSITY=airdens(K)
            AIRMASS=VOLUME*DENSITY+DIFMASS(K) ! mole or kg air
            TRACERMASS=VOLUME*DENSITY*TXS(K)+TRACUP(K) ! mole gas or ug aerosol
            XS1(K)=TRACERMASS/AIRMASS ! mole/mole or ug/kg
            TXS(K)=XS1(K)
          ENDDO

     !----------------------------------------------------------------   
     !1.1.6 Check mass balance for column 
     !----------------------------------------------------------------   
          SUMASS2=0.
          DO K=kts,kte
            !gas in mole/mole and aerosol in ug/kg
           if (.not.is_aerosol(nv)) then 
            SUMASS2=SUMASS2+XS1(K)*(dzq(k)*airdens(k)*DX2+DIFMASS(K))
            IF (SUMASS2 .GT. 1.e20) then  ! larger than 1.e10 mole
              WRITE(0,*) 'TRACER MASS2 > 1.e20 mole',nv,I,J,K
              PRINT *, 'TRACER MASS2 > 1.e20 mole',nv,I,J,K
              WRITE(0,*) SUMASS2
              PRINT *, SUMASS2
            ENDIF
           else   ! aerosol in ug/kg-dryair 
            SUMASS2=SUMASS2+XS1(K)*(dzq(k)*airdens(k)*DX2+DIFMASS(K))
            IF (SUMASS2 .GT. 1.e30) then  ! larger than 1.e30 ug 
              WRITE(0,*) 'TRACER MASS2 > 1.e30 ug',nv,I,J,K
              PRINT *, 'TRACER MASS2 > 1.e30 ug',nv,I,J,K
              WRITE(0,*) SUMASS2
              PRINT *, SUMASS2
            ENDIF
           endif
          ENDDO

          IF ( ((SUMASS2-SUMASS1)/SUMASS1) .GT. 1.e-5) THEN
            WRITE(0,*) 'TRACER MASS IS NOT BALANCE AFTER UPDRAFT!!!!!'
            PRINT *, 'TRACER MASS IS NOT BALANCE AFTER UPDRAFT!!!!!'
            WRITE(0,*) SUMASS1,SUMASS2
            PRINT *, SUMASS1,SUMASS2
            WRITE(0,*) (SUMASS2-SUMASS1)/SUMASS1
            PRINT *, (SUMASS2-SUMASS1)/SUMASS1
     
            TOTRACER=0.
            DO k=kts,kte
            TOTRACER=TOTRACER+TRACUP(K)
            write(0,*)nv,is_aerosol(nv),K,TRACUP(K)
            PRINT*,nv,is_aerosol(nv),K,TRACUP(K)
            ENDDO
            WRITE(0,*) 'TOTAL TRACER TENDENCE',I,J,ISHALL,TOTRACER
            PRINT *, 'TOTAL TRACER TENDENCE',I,J,ISHALL,TOTRACER
            STOP
          ENDIF

     !----------------------------------------------------------------   
     !1.2 CALCULATE EFFECT OF COMPENSATING DUE TO UPDRAFT IN DEEP CONVECTION 
     !----------------------------------------------------------------   
 
     !-----------------------------------------------   
     !1.2.1 calculate the top layer of cloud
     !-----------------------------------------------   
            TRACSUB1(ZTOP)=SUBMASS(ZTOP)*TXS(ZTOP)
            XS1(ZTOP)=XS1(ZTOP)
            TXS(ZTOP)= XS1(ZTOP)

     !-----------------------------------------------   
     !1.2.2 calclulate iner layers between cloud top and original level
     !-----------------------------------------------   
          DO K=ZTOP-1,ZLC+1,-1
            !a. mixing with upper layer downward flux
            VOLUME=DZQ(K)*DX2  ! m3
            DENSITY=airdens(K) ! mole or kg/m3
            AIRMASS=VOLUME*DENSITY+DIFMASS(K)+SUBMASS(K+1) ! mole or kg
            !TRACERMASS=(VOLUME*DENSITY+DIFMASS(K))*TXS(K)+TRACSUB1(K+1)*(1.-scavf(k,1))
            TRACERMASS=(VOLUME*DENSITY+DIFMASS(K))*TXS(K)+TRACSUB1(K+1) !mole or ug
            XS1(K)=TRACERMASS/AIRMASS !mole/mole or ug/kg
            TXS(K)=XS1(K)
            !b. consider the downward flux in the layer
            AIRMASS=VOLUME*DENSITY
            SUBMASS(K)=DIFMASS(K)+SUBMASS(K+1)
            !check submass 
            IF (SUBMASS(K).LT.0) THEN
              WRITE(0,*) 'ERROR SUBSIDE MASS FLUX'
              PRINT*,'ERROR SUBSIDE MASS FLUX'
              WRITE(0,*) I,J,K,SUBMASS(K),DIFMASS(K),SUBMASS(K+1)
              PRINT*,I,J,K,SUBMASS(K),DIFMASS(K),SUBMASS(K+1)
              STOP
            ENDIF
            TRACSUB1(K) =SUBMASS(K)*TXS(K)
            TRACERMASS=TRACERMASS-TRACSUB1(K)
            XS1(K)=TRACERMASS/AIRMASS
            TXS(K)=XS1(K)
          ENDDO

     !-----------------------------------
     !1.2.3 calculate original updraft level
     !-----------------------------------
            VOLUME=DZQ(ZLC)*DX2 !m3
            DENSITY=airdens(ZLC) !mole/m3 or kg/m3
            AIRMASS=VOLUME*DENSITY !mole or kg
            !TRACERMASS=(VOLUME*DENSITY+DIFMASS(ZLC))*TXS(ZLC)+TRACSUB1(ZLC+1)*(1.-scavf(k,1))
            TRACERMASS=(VOLUME*DENSITY+DIFMASS(ZLC))*TXS(ZLC)+TRACSUB1(ZLC+1) !mole or ug
            XS1(ZLC)=TRACERMASS/AIRMASS
            TXS(ZLC)=XS1(ZLC)
 
     !--------------------------------------
     !1.2.4 Check mass balance for column 
     !--------------------------------------
            SUMASS2=0.
            DO K=kts,kte
            SUMASS2=SUMASS2+XS1(K)*(DZQ(K)*airdens(K)*DX2)
            ENDDO

            IF ( ((SUMASS2-SUMASS1)/SUMASS1) .GT. 1.e-5) THEN
             WRITE(0,*) 'TRACER MASS IS NOT BALANCE AFTER SUBSIDENCE!!!!!'
             PRINT *, 'TRACER MASS IS NOT BALANCE AFTER SUBSIDENCE!!!!!'
             WRITE(0,*) SUMASS1, SUMASS2,(SUMASS2-SUMASS1)/SUMASS1
             PRINT *, SUMASS1, SUMASS2,(SUMASS2-SUMASS1)/SUMASS1

             TOTRACER=0.
             DO K=kts,kte
             TOTRACER=TOTRACER+TRACSUB1(K)
             write(0,*) k,tracsub1(k)
             print*, k,tracsub1(k)
             ENDDO
             WRITE(0,*) 'TOTAL TRACER TENDENCE',I,J,ISHALL,nv,is_aerosol(nv),TOTRACER
             PRINT *, 'TOTAL TRACER TENDENCE',I,J,ISHALL,nv,is_aerosol(nv),TOTRACER
             STOP
            ENDIF

     !----------------------------------------------------------------   
     !1.3 COMPUTE DOWNDRAFT EFFECT IN DEEP CONVECTION 
     !----------------------------------------------------------------   

         IF(ZLFS.GT.ZLDB ) THEN  ! ALLOWED DOWNDRAFT FOR DEEP CONVEC.   

     !----------------------------------------------------------------   
     !1.3.0 initializing 
     !----------------------------------------------------------------   
          do k=kts,kte
          DIFMASS(k) = 0.0
          SUBMASS(k) = 0.0
          enddo

          DO K=ZLFS-1,ZLDB+1,-1
           IF (FLUXDOWN(K).EQ.0 ) THEN
            WRITE(0,*) 'WRONG FLUXDOWN',I,J,K,ZLFS,ZLDB
            PRINT *, 'WRONG FLUXDOWN',I,J,K,ZLFS,ZLDB
            STOP
           ENDIF
           !fluxdown ajust for mass consistence
           DENMASS=DENTR(K)*DTIN ! mole of entrain air
           DDEMASS=DDETR(K)*DTIN ! mole of detrain air
           IF(DENMASS.NE.DDEMASS) THEN
           IF((FLUXDOWN(K-1)-FLUXDOWN(K))/(DENMASS-DDEMASS).GT.1.e-5.AND.  &
             (FLUXDOWN(K-1)-FLUXDOWN(K))/(DENMASS-DDEMASS).LT.0 ) THEN
                WRITE(0,*) 'DOWNDRAFT AJUSTMENT AT',I,J,K
                PRINT *, 'DOWNDRAFT AJUSTMENT AT',I,J,K
                WRITE(0,*) FLUXDOWN(K),FLUXDOWN(K+1),DENMASS,DDEMASS
                PRINT *, FLUXDOWN(K),FLUXDOWN(K+1),DENMASS,DDEMASS
                FLUXDOWN(K-1)=FLUXDOWN(K)+DENMASS-DDEMASS
           ENDIF
           ENDIF
          ENDDO !K


     !----------------------------------------------------------------   
     !1.3.1 calculate the LFS layer tracer tendence 
     !----------------------------------------------------------------   
            UAIRMASS=FLUXDOWN(ZLFS)*DTIN  ! mole or kg of air mass
            VOLUME=DZQ(ZLFS)*DX2 !m3
            DENSITY=airdens(ZLFS) !mole or kg /m3
            AIRMASS=VOLUME*DENSITY-UAIRMASS !mole or kg
            DIFMASS(ZLFS)=-UAIRMASS !mole or kg, save for uplifting
            FLUXRAT=TXS(ZLFS)*(1.-scavf(zlfs,2))
            TRACDW(ZLFS)=-UAIRMASS*TXS(ZLFS) !mole or ug of tracer

     !----------------------------------------------------------------   
     !1.3.2 calculate layer between LFS AND LDB 
     !----------------------------------------------------------------   
          DO K=ZLFS-1,ZLDB+1,-1
            UAIRMASS=FLUXDOWN(K+1)*DTIN ! mole or kg airmass from upper level
            DENMASS=DENTR(K)*DTIN ! mole or kg entrain air
            DDEMASS=DDETR(K)*DTIN ! mole or kg detrain air
            VOLUME=DZQ(K)*DX2  ! m3
            DENSITY=airdens(K)  !mole or kg /m3
            AIRMASS=VOLUME*DENSITY-DENMASS+DDEMASS !mole or kg
            DIFMASS(K)=-DENMASS+DDEMASS !save for uplifting
            FLUXTRA=UAIRMASS*FLUXRAT+DENMASS*TXS(K)*(1.-scavf(k,2))-DDEMASS*FLUXRAT !mole gas or ug aerosol
            TRACDW(K)=DDEMASS*FLUXRAT-DENMASS*TXS(K) !mole gas or ug aerosol
            FLUXRAT=(FLUXTRA/(FLUXDOWN(K)*DTIN))*(1.-scavf(k,2))
          ENDDO

     !----------------------------------------------------------------   
     !1.3.3 calculate LDB layer
     !----------------------------------------------------------------   
            UAIRMASS=FLUXDOWN(ZLDB+1)*DTIN
            VOLUME=DZQ(ZLDB)*DX2
            DENSITY=airdens(ZLDB)
            AIRMASS=VOLUME*DENSITY+UAIRMASS
            DIFMASS(ZLDB)=UAIRMASS
            SUBMASS(ZLDB)=DIFMASS(ZLDB)
            TRACDW(ZLDB)=UAIRMASS*FLUXRAT 

     !----------------------------------------------------------------   
     !1.3.4 Check mass balance for downdraft 
     !----------------------------------------------------------------   
            TRACSUM=0.0
            DO K=ZLDB,ZLFS
             TRACSUM=TRACSUM+TRACDW(K)
            ENDDO

            IF( (TRACSUM/SUMASS1).GT.1.e-5) THEN
             WRITE(0,*) 'TRACER MASS IS NOT BALANCE DURING DOWNDRAFT!!!'
             PRINT *, 'TRACER MASS IS NOT BALANCE DURING DOWNDRAFT!!!'
             WRITE(0,*) TRACSUM,I,J,nv
             PRINT *, TRACSUM,I,J,nv
             WRITE(0,*) TRACSUM/SUMASS1
             PRINT *, TRACSUM/SUMASS1
             STOP
            ENDIF

     !----------------------------------------------------------------   
     !1.3.5 Calculate tracer concentration due to downdraft in deep convection 
     !----------------------------------------------------------------   
           DO K=ZLFS,ZLDB,-1
             VOLUME=DZQ(K)*DX2 !m3
             DENSITY=airdens(K) !mole or kg /m3
             AIRMASS=VOLUME*DENSITY+DIFMASS(K) ! mole or kg
             TRACERMASS=VOLUME*DENSITY*TXS(K)+TRACDW(K) ! mole or ug
             XS1(K)=TRACERMASS/AIRMASS ! mole/mole or ug/kg
             TXS(K)=XS1(K)
           ENDDO

     !----------------------------------------------------------------   
     !1.3.6 Check mass balance for column 
     !----------------------------------------------------------------   
           SUMASS2=0.
           DO K=kts,kte
             SUMASS2=SUMASS2+XS1(K)*(DZQ(K)*airdens(K)*DX2+DIFMASS(K))
           ENDDO

           IF ( ((SUMASS2-SUMASS1)/SUMASS1) .GT. 1.e-5) THEN
             WRITE(0,*) 'TRACER MASS IS NOT BALANCE AFTER DOWNDRAFT!!!!'
             PRINT *, 'TRACER MASS IS NOT BALANCE AFTER DOWNDRAFT!!!!'
             WRITE(0,*) SUMASS1, SUMASS2,(SUMASS2-SUMASS1)/SUMASS1
             PRINT *, SUMASS1, SUMASS2,(SUMASS2-SUMASS1)/SUMASS1

             TOTRACER=0.
             DO K=kts,kte
             TOTRACER=TOTRACER+TRACDW(K)
             WRITE(0,*) K,TRACDW(K)
             PRINT*,K,TRACDW(K)
             ENDDO
             WRITE(0,*) 'TOTAL TRACER TENDENCE',I,J,nv,is_aerosol(nv),ISHALL, TOTRACER
             PRINT *, 'TOTAL TRACER TENDENCE',I,J,nv,is_aerosol(nv),ISHALL, TOTRACER
             STOP
           ENDIF

     !----------------------------------------------------------------   
     !1.4 CALCULATE THE EFFECT OF COMPENSATING DUE TO DOWNDRAFT IN DEEP CONVECTION 
     !----------------------------------------------------------------   

     !----------------------------------------------------------------   
     !1.4.1 calculate the LDB layer 
     !----------------------------------------------------------------   
            TRACSUB2(ZLDB)=SUBMASS(ZLDB)*TXS(ZLDB)
            XS1(ZLDB)=XS1(ZLDB)
            TXS(ZLDB)=TXS(ZLDB)

     !----------------------------------------------------------------   
     !1.4.2 calculate the iner layer between cloud top and LFS layer 
     !----------------------------------------------------------------   
          DO K=ZLDB+1,ZLFS-1
            !a. mixing with lower layer upward flux
            VOLUME=DZQ(K)*DX2
            DENSITY=airdens(K)
            AIRMASS=VOLUME*DENSITY+DIFMASS(K)+SUBMASS(K-1) ! mole or kg
            !TRACERMASS=(VOLUME*DENSITY+DIFMASS(K))*TXS(K)+TRACSUB2(K-1)*(1.-scavf(k,2))
            TRACERMASS=(VOLUME*DENSITY+DIFMASS(K))*TXS(K)+TRACSUB2(K-1) !mole or ug
            XS1(K)=TRACERMASS/AIRMASS !mole/mole or ug/kg
            TXS(K)=XS1(K)
            !b. consider the upward flux in the layer
            AIRMASS=VOLUME*DENSITY
            SUBMASS(K)=DIFMASS(K)+SUBMASS(K-1)
            ! check submass 
            IF (SUBMASS(K).LT.0) THEN
              WRITE(0,*) 'ERROR UPLIFT MASS FLUX'
              PRINT*,'ERROR UPLIFT MASS FLUX'
              WRITE(0,*) I,J,K,NV,SUBMASS(K),DIFMASS(K),SUBMASS(K-1)
              PRINT*,I,J,K,NV,SUBMASS(K),DIFMASS(K),SUBMASS(K-1)
              STOP
            ENDIF
            TRACSUB2(K)=SUBMASS(K)*TXS(K)
            TRACERMASS=TRACERMASS-TRACSUB2(K)
            XS1(K)=TRACERMASS/AIRMASS
            TXS(K)=XS1(K)
          ENDDO

     !----------------------------------------------------------------   
     !1.4.3 calculate the LFS layer
     !----------------------------------------------------------------   
            VOLUME=DZQ(ZLFS)*DX2
            DENSITY=airdens(ZLFS)
            AIRMASS=VOLUME*DENSITY ! mole or kg 
            !TRACERMASS=(VOLUME*DENSITY+DIFMASS(ZLFS))*TXS(ZLFS)+TRACSUB2(ZLFS-1)*(1.-scavf(k,2))
            TRACERMASS=(VOLUME*DENSITY+DIFMASS(ZLFS))*TXS(ZLFS)+TRACSUB2(ZLFS-1) !mole or ug
            XS1(ZLFS)=TRACERMASS/AIRMASS
            TXS(ZLFS)=XS1(ZLFS)

     !----------------------------------------------------------------   
     !1.4.4 check mass balance for column 
     !----------------------------------------------------------------   
            SUMASS2=0.
            DO K=kts,kte
            SUMASS2=SUMASS2+XS1(K)*(DZQ(K)*airdens(K)*DX2)
            ENDDO

            IF ( ((SUMASS2-SUMASS1)/SUMASS1) .GT. 1.e-5) THEN
             WRITE(0,*) 'TRACER MASS IS NOT BALANCE AFTER UPLIFTING!!!!'
             PRINT *, 'TRACER MASS IS NOT BALANCE AFTER UPLIFTING!!!!'
             WRITE(0,*) SUMASS1, SUMASS2, (SUMASS2-SUMASS1)/SUMASS1
             PRINT *, SUMASS1, SUMASS2, (SUMASS2-SUMASS1)/SUMASS1

             TOTRACER=0.
             DO K=kts,kte
             TOTRACER=TOTRACER+TRACSUB2(K)
             WRITE(0,*) K,TRACSUB2(K)
             PRINT*,K,TRACSUB2(K)
             ENDDO
             WRITE(0,*) 'TOTAL TRACER TENDENCE',I,J,ISHALL,nv,is_aerosol(nv),TOTRACER
             PRINT *, 'TOTAL TRACER TENDENCE',I,J,ISHALL,nv,is_aerosol(nv),TOTRACER
             STOP
            ENDIF

         ENDIF !LFS, DOWNDRAFT

     ENDIF ! ISHALL, for deep convection

   !================================================================   
   !2 Start shallow convective transport 
   !================================================================   
     IF ( ISHALL==1) THEN  ! FOR SHALLOW CONVECTION

     !----------------------------------------------------------------   
     !2.1 calculate tracer tendence between LC layer and cloud top 
     !----------------------------------------------------------------   
     !Note:the mixing between each two layers every time 
     !----------------------------------------------------------------   

           DO K=ZLC,ZTOP-1
             UAIRMASS=FLUXUP(K)*DTIN  ! mole or kg
             AIRMASS1=DZQ(K)*DX2*airdens(K) !mole or kg, current layer
             AIRMASS2=DZQ(K+1)*DX2*airdens(K+1) !mole or kg, upper layer
             TRACERMASS1=AIRMASS1*TXS(K)       & ! current layer
                        -UAIRMASS*TXS(K)       & ! upward flux
                        +UAIRMASS*TXS(K+1)       ! downward flux
             TRACERMASS2=AIRMASS2*TXS(K+1)     & ! upper layer
                        -UAIRMASS*TXS(K+1)     & ! downward flux
                        +UAIRMASS*TXS(K)         ! upward flux
             XS1(K)=TRACERMASS1/AIRMASS1   ! mole/mole or ug/kg
             XS1(K+1)=TRACERMASS2/AIRMASS2
             TXS(K:K+1)=XS1(K:K+1)
           ENDDO

     !----------------------------------------------------------------   
     !2.2 check tracer mass balance 
     !----------------------------------------------------------------   
             SUMASS2=0.
             DO K=kts,kte
              SUMASS2=SUMASS2+XS1(K)*(DZQ(K)*airdens(K)*DX2) ! mole or ug
             ENDDO

             IF ( (SUMASS2-SUMASS1).GT.SUMASS1*1.e-5) THEN
              WRITE(0,*) 'TRACER MASS IS NOT BALANCE AFTER SHALLOW SUBSIDE!'
              PRINT *,'TRACER MASS IS NOT BALANCE AFTER SHALLOW SUBSIDE!'
              WRITE(0,*) SUMASS1, SUMASS2
              PRINT *, SUMASS1, SUMASS2
              STOP
             ENDIF

     ENDIF   !ISHALL, for shallow convection

   !================================================================   
   !3 Check tracer mass balance due to convective transport 
   !================================================================   
          SUMASS2=0.
          DO K=kts,kte
            SUMASS2=SUMASS2+XS1(K)*(dzq(k)*airdens(k)*DX2)
           !gas in mole/mole and aerosol in ug/kg
           if (.not.is_aerosol(nv)) then
            IF (SUMASS2 .GT. 1.e20) then  ! larger than 1.e10 mole
              WRITE(0,*) 'TRACER MASS2 > 1.e20 mole',nv,I,J,K,SUMASS2
              PRINT *, 'TRACER MASS2 > 1.e20 mole',nv,I,J,K,SUMASS2
            ENDIF
           else   ! aerosol in ug/kg-dryair 
            IF (SUMASS2 .GT. 1.e30) then  ! larger than 1.e30 ug 
              WRITE(0,*) 'TRACER MASS2 > 1.e30 ug',nv,I,J,K,SUMASS2
              PRINT *, 'TRACER MASS2 > 1.e30 ug',nv,I,J,K,SUMASS2
            ENDIF
           endif
          ENDDO


          IF ( ((SUMASS2-SUMASS1)/SUMASS1) .GT. 1.e-5) THEN
            WRITE(0,*) 'TRACER MASS IS NOT BALANCE AFTER CONVECTION!!!!!'
            PRINT *, 'TRACER MASS IS NOT BALANCE AFTER CONVECTION!!!!!'
            WRITE(0,*) SUMASS1, SUMASS2, (SUMASS2-SUMASS1)/SUMASS1
            PRINT *, SUMASS1, SUMASS2, (SUMASS2-SUMASS1)/SUMASS1

            TOTRACER=0.
            DO K=ZLC,ZTOP
             TOTRACER=TOTRACER+TRACUP(K)+TRACSUB1(K)+TRACDW(K)+TRACSUB2(K)
             WRITE(0,*) K,TRACUP(K),TRACSUB1(K),TRACDW(K),TRACSUB2(K)
             PRINT*,K,TRACUP(K),TRACSUB1(K),TRACDW(K),TRACSUB2(K)
            ENDDO
            WRITE(0,*) 'TOTAL TRACER TENDENCE',I,J,ISHALL,nv,TOTRACER
            PRINT *, 'TOTAL TRACER TENDENCE',I,J,ISHALL,nv,TOTRACER
            STOP
          ENDIF


   ENDIF !CHECK2D, for convective transport allowed

        DO K=KTS,KTE
         if(.not.is_aerosol(nv))then !gases in ppmv
         chem(i,k,j,nv)=xs1(k)*1.e6  ! mole/mole to ppmv
         else !aerosol in ug/kg dry-air 
         chem(i,k,j,nv)=xs1(k) 
         endif
        ENDDO


 900   CONTINUE  !num_chem


1000   CONTINUE  !I,J

   END SUBROUTINE kfeta_cup_transport 

!===========================================================================
 subroutine sorgam_soluble(sol_f0,num_chem,chem0)

 !calculate the soluble fraction of the internal-mixing aerosol based on hygroscopicity
 !this is not accurate, just for temporary application !!!! 

 USE module_data_sorgam

 IMPLICIT NONE

   real, DIMENSION(num_chem ),intent(inout)  :: sol_f0 
   integer, intent(in)  :: num_chem
   REAL, DIMENSION(num_chem ),INTENT(IN) :: chem0
   integer n,m,l,lmass
   real sumhygro,sumvol,hygro

          do n=1,ntype_aer
          do m=1,nsize_aer(n)
           sumhygro=0.
           sumvol=0. 
           do l=1,ncomp_aer(n)
             lmass=massptr_aer(l,m,n,ai_phase)
             sumhygro = sumhygro+hygro_aer(l,n)*chem0(lmass)/dens_aer(l,n)
             sumvol = sumvol+chem0(lmass)/dens_aer(l,n)
           enddo ! comp
           hygro=sumhygro/sumvol
           do l=1,ncomp_aer(n)
             lmass=massptr_aer(l,m,n,ai_phase)
             sol_f0(lmass)=min(max(hygro/1.0,0.0),1.0)
           enddo
          enddo
          enddo

 end subroutine sorgam_soluble

!===========================================================================
 subroutine mosaic_soluble(sol_f0,num_chem,chem0)
 
 !calculate the soluble fraction of the internal-mixing aerosol based on hygroscopicity
 !this is not accurate, just for temporary application !!!! 
 
 USE module_data_mosaic_asect
 
 IMPLICIT NONE

   real, DIMENSION(num_chem ),intent(inout)  :: sol_f0
   integer, intent(in)  :: num_chem
   REAL, DIMENSION(num_chem ),INTENT(IN) :: chem0
   integer n,m,l,lmass
   real sumhygro,sumvol,hygro

          do n=1,ntype_aer
          do m=1,nsize_aer(n)
           sumhygro=0.
           sumvol=0.
           do l=1,ncomp_aer(n)
             lmass=massptr_aer(l,m,n,ai_phase)
             sumhygro = sumhygro+hygro_aer(l,n)*chem0(lmass)/dens_aer(l,n)
             sumvol = sumvol+chem0(lmass)/dens_aer(l,n)
           enddo ! comp
           hygro=sumhygro/sumvol
           do l=1,ncomp_aer(n)
             lmass=massptr_aer(l,m,n,ai_phase)
             sol_f0(lmass)=min(max(hygro/1.0,0.0),1.0)
           enddo
          enddo
          enddo

 end subroutine mosaic_soluble

!===========================================================================
 subroutine CUSCAV_KF (II,JJ,N,DZQ,dx2,airdens,vud,vdd,tmp,pmid,        &
                       num_chem,chem0,F,   &
                       ids,ide, jds,jde, kds,kde,                       &
                       ims,ime, jms,jme, kms,kme,                       &
                       its,ite, jts,jte, kts,kte                        )

!******************************************************************************
!  Subroutine CUSCAV_KF computes F, the fraction of soluble tracer lost by 
!  scavenging in convective cloud updrafts. 
!
!  Arguments as Input:
!  (1 ) N    (INTEGER) : Tracer number
!
!  Arguments as Output:
!  (2 ) F    (REAL)  : Fraction of tracer scavenged in cloud updraft [0-1]
!
   IMPLICIT NONE


   ! Arguments
   INTEGER,INTENT(IN)  ::    &
                          II,JJ,N, &
                          ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte
   REAL,INTENT(IN),DIMENSION(kms:kme)  :: DZQ,TMP,PMID,AIRDENS
   REAL,INTENT(IN)  :: DX2
   REAL,INTENT(IN),DIMENSION(kms:kme)  :: VUD,VDD
   integer, intent(in)  :: num_chem

   REAL,INTENT(INOUT),DIMENSION(kms:kme,num_chem) :: chem0
   REAL,INTENT(INOUT),DIMENSION(kms:kme,2) :: F

   ! Local variables 
   INTEGER            :: I, J, K,NN
   REAL               :: L2G, I2G, C_TOT, F_L, F_I, KR, SO2LOSS

   ! Kc is the conversion rate from cloud condensate to precip [s^-1]
   REAL, PARAMETER    :: KC        = 5.e-3

   ! CONV_H2O2 = 0.6 * SQRT( 1.9 ), used for the ice to gas ratio for H2O2
   ! 0.6 is ( sticking  coeff H2O2  / sticking  coeff  water )
   ! 1.9 is ( molecular weight H2O2 / molecular weight water )
   REAL*8, PARAMETER    :: CONV_H2O2 = 8.27042925126d-1

   ! CONV_NH3 = 0.6 * SQRT( 0.9 ), used for the ice to gas ratio for NH3
   ! 0.6 is ( sticking  coeff  NH3 / sticking  coeff  water )
   ! 0.9 is ( molecular weight NH3 / molecular weight water )
   REAL*8, PARAMETER    :: CONV_NH3  = 5.69209978831d-1

   REAL csr_conv(3) ! for three modes, aerosol activation rate
   REAL CSR

   REAL, DIMENSION (kms:kme) :: CLDLIQ,CLDICE 
   REAL C_H2O

   csr_conv      = (/0.60, 0.99, 0.99/) ! see module_mosaic_wetscav.F

   !calculate the updraft velocity (m/s)
   !This is grid-average updraft velocity, 
   !Now use real updraft velocity from convective scheme
!  DO K = KTS,KTE
!    !fluxup: mole or kg/s  VUD,VDD: m/s
!    VUD(k) = FLUXUP(K)/dx2/airdens(k)
!    VDD(k) = FLUXDOWN(K)/dx2/airdens(k)
!    VUD(k) = 5.0 ! over ocean, or 10.0 over land 
!    VDD(k) = 5.0 
!  ENDDO

   !-------------------------------------------------------
   !calculate the convective cloud property
   ! CLDLIQ, the cloud liquid water content [cm3 H2O/cm3 air],
   !         is a function of the local Kelvin temperature:
   ! CLDLIQ = 2e-6                    [     T >= 268 K    ]
   ! CLDLIQ = 2e-6 * ((T - 248) / 20) [ 248 K < T < 268 K ]
   ! CLDLIQ = 0                       [     T <= 248 K    ]
   !-------------------------------------------------------
   DO K = KTS,KTE
    IF ( TMP(K) >= 268.0 ) THEN
      CLDLIQ(K) = 2e-6
    ELSE IF ( TMP(K) > 248.0 .and. TMP(K) < 268.0 ) THEN
      CLDLIQ(K) = 2e-6 * ( ( TMP(K) - 248.0 ) / 20.0 )
    ELSE
      CLDLIQ(K) = 0.0
    ENDIF
   !-------------------------------------------------------
   ! CLDICE, the cloud ice content [cm3 ice/cm3 air] is given by:
   ! CLDICE = 2e-6 - CLDLIQ
   !-------------------------------------------------------
      CLDICE(K) = 2e-6 - CLDLIQ(K)
   ENDDO

   !-------------------------------
   ! COMPUTE_F begins here!
   !-------------------------------
      IF (is_aerosol(n)) THEN 
       CSR=CSR_CONV(2)  ! now assume aerosol has the same activation rate for different size
                        ! in future, it should pass massptr_aer here to decide the size 
       CALL F_AEROSOL( KC,DZQ,CSR,F,VUD,VDD,      &
                       ids,ide, jds,jde, kds,kde,       &
                       ims,ime, jms,jme, kms,kme,       &
                       its,ite, jts,jte, kts,kte        )
      !-------------------------------
      ! HNO3 (gases), treated as aerosol
      !-------------------------------
      ELSE IF ( N == p_hno3 ) THEN
       CSR=CSR_CONV(2)  ! assume HNO3 as accumulation mode
       CALL F_AEROSOL( KC,DZQ,CSR,F,VUD,VDD,      &
                       ids,ide, jds,jde, kds,kde,       &
                       ims,ime, jms,jme, kms,kme,       &
                       its,ite, jts,jte, kts,kte        )
      !-------------------------------
      ! H2O2 (liquid & ice phases)
      !-------------------------------
      ELSE IF ( N == p_h2o2 ) THEN
         ! No scavenging at the surface
         F(kts,1:2) = 0.0
         ! Start scavenging at level 2
         DO K = kts+1, kte 
            ! Compute ice to gas ratio for H2O2 by co-condensation
            ! (Eq. 9, Jacob et al, 2000)
            c_h2o= E_ICE(TMP(K)) / (pmid(k)*0.01)  ! hPa 
            IF ( C_H2O > 0.0 ) THEN
               I2G = ( CLDICE(K) / C_H2O ) * CONV_H2O2
            ELSE
               I2G = 0.0
            ENDIF
            ! Compute liquid to gas ratio for H2O2, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            CALL COMPUTE_L2G( 8.3d4,-7.4d3, TMP(K), CLDLIQ(K), L2G )
            ! Fraction of H2O2 in liquid & ice phases
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1.0 + L2G + I2G
            F_L   = L2G / C_TOT
            F_I   = I2G / C_TOT
            ! Compute the rate constant Kr.  The retention factor for 
            ! liquid H2O2 is 0.05 for 248 K < T < 268 K and 1.0 for 
            ! T >= 268 K. (Eq. 1, Jacob et al, 2000)
            IF ( TMP(K) >= 268d0 ) THEN
               KR = KC * ( F_L + F_I )
            ELSE IF ( TMP(K) > 248d0  .and. TMP(K) < 268d0 ) THEN
               KR = KC * ( ( 5d-2 * F_L ) + F_I )
            ELSE
               KR = KC * F_I
            ENDIF
            ! Compute F, the fraction of scavenged H2O2.
            ! (Eq. 2, Jacob et al, 2000)
            IF (VUD(K).GT.0) THEN 
             F(k,1) = 1.0 - EXP( -KR * DZQ(K) / Vud(K) )
            ELSE
             F(k,1) = 0.0
            ENDIF
            IF (VDD(K).GT.0) THEN 
             F(k,2) = 1.0 - EXP( -KR * DZQ(K) / Vdd(K) )
            ELSE  
             F(k,2) = 0.0
            ENDIF
         ENDDO

         DO K = KTS,KTE
          IF (II==20.AND.JJ==10.and.F(K,1).GT.0.1) THEN
!          WRITE(0,*) 'CZHAO CHECK F1_H2O2',K,F(K,1),VUD(K)
!          PRINT*, 'CZHAO CHECK F1_H2O2',K,F(K,1),VUD(K)
          ENDIF   
         ENDDO

      !-------------------------------
      ! CH2O (liquid phase only)
      !-------------------------------
      ELSE IF ( N == p_hcho ) THEN
         ! No scavenging at the surface
         F(kts,1:2) = 0.0
         ! Apply scavenging in levels 2 and higher
         DO K = KTS+1, KTE 
            ! Compute liquid to gas ratio for CH2O, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            CALL COMPUTE_L2G( 3.0d3, -7.2d3, TMP(K), CLDLIQ(K), L2G )
            ! Fraction of CH2O in liquid phase 
            ! NOTE: CH2O does not exist in the ice phase!
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1.0 + L2G
            F_L   = L2G / C_TOT
            ! Compute the rate constant K.  The retention factor 
            ! for liquid CH2O is 0.0 for T <= 248K and 0.02 for 
            ! 248 K < T < 268 K. (Eq. 1, Jacob et al, 2000)
            IF ( TMP(K) >= 268d0 ) THEN
               KR = KC * F_L
            ELSE IF ( TMP(K) > 248d0 .and. TMP(K) < 268d0 ) THEN
               KR = KC * ( 2d-2 * F_L )
            ELSE
               KR = 0.0
            ENDIF
            ! F is the fraction of CH2O scavenged out of the updraft
            ! (Eq. 2, Jacob et al, 2000)
            IF (VUD(K).GT.0) THEN
            F(K,1) = 1.0 - EXP( -KR * DZQ(K) / Vud(K) )
            ELSE 
             F(k,1) = 0.0
            ENDIF
            IF (VDD(K).GT.0) THEN
             F(k,2) = 1.0 - EXP( -KR * DZQ(K) / Vdd(K) )
            ELSE
             F(k,2) = 0.0
            ENDIF
         ENDDO

      !------------------------------
      ! SO2 (treat as aerosol)
      !------------------------------
      ELSE IF ( N == p_so2 ) THEN 
         ! Compute fraction of SO2 scavenged
         CSR=CSR_CONV(2)  ! now assume aerosol has the same activation rate for different size
                          ! in future, it should pass massptr_aer here to decide the size 
         CALL F_AEROSOL(KC,DZQ,CSR,F,VUD,VDD,      &
                       ids,ide, jds,jde, kds,kde,       &
                       ims,ime, jms,jme, kms,kme,       &
                       its,ite, jts,jte, kts,kte        )
         !==============================================================
         ! Coupled full chemistry/aerosol simulation:
         ! Use the wet scavenging formula of Chin et al [1996], 
         ! such that a soluble fraction of SO2 is limited by the
         ! availability of H2O2 in the precipitating grid box. 
         ! Scavenge the soluble SO2 at the same rate as the sulfate.
         ! Update H2O2_sav and SO2_sav for use in RAINOUT, WASHOUT
         !==============================================================
         DO K = KTS+1, KTE 
            ! Make sure to deplete H2O2s the same as SO2s. 
            ! (dkh, rjp, bmy, 11/17/05)
            IF ( chem0(k,p_so2) > 1.e-16 ) THEN
               ! Limit F
               SO2LOSS      = MIN( chem0(k,p_so2), chem0(k,p_h2o2) )
               F(K,1)     = F(K,1) * SO2LOSS / chem0(k,p_so2) 
               F(K,1)     = MAX(F(K,1), 0.0)
               F(K,2)     = F(K,2) * SO2LOSS / chem0(k,p_so2) 
               F(K,2)     = MAX(F(K,2), 0.0)
               ! Update saved H2O2 concentration
               chem0(k,p_h2o2) = chem0(k,p_h2o2) - (chem0(k,p_h2o2) * F(k,1))
               chem0(k,p_h2o2) = MAX(chem0(k,p_h2o2), 1.e-16)
            ELSE
               ! Set F = 0 if SO2s < EPSILON
               F(K,1:2)     = 0.0
            ENDIF
         ENDDO

      !-------------------------------
      ! SO4 (gaseous aerosol), assume HCL the same as SO4
      !-------------------------------
      ELSE IF ( N == p_sulf.or.N == p_hcl) THEN
         CSR=CSR_CONV(2)  ! now assume aerosol has the same activation rate for different size
                          ! in future, it should pass massptr_aer here to decide the size 
         CALL F_AEROSOL(KC,DZQ,CSR,F,VUD,VDD,      &
                       ids,ide, jds,jde, kds,kde,       &
                       ims,ime, jms,jme, kms,kme,       &
                       its,ite, jts,jte, kts,kte        )

      !-------------------------------
      ! MSA (treat as aerosol)
      !------------------------------- 
      ELSE IF ( N == p_msa ) THEN 
         CSR=CSR_CONV(2)  ! now assume aerosol has the same activation rate for different size
         CALL F_AEROSOL(KC,DZQ,CSR,F,VUD,VDD,      &
                       ids,ide, jds,jde, kds,kde,       &
                       ims,ime, jms,jme, kms,kme,       &
                       its,ite, jts,jte, kts,kte        )

      !-------------------------------
      ! NH3 (liquid & ice phases)
      !-------------------------------
      ELSE IF ( N == p_nh3 ) THEN 
         ! No scavenging at surface
         F(kts,1:2) = 0.0
         ! Start scavenging at level 2
         DO K = KTS+1, KTE 
            ! Compute ice to gas ratio for NH3 by co-condensation
            ! (Eq. 9, Jacob et al, 2000)
            c_h2o= E_ICE(TMP(K)) / (pmid(k)*0.01)  ! hPa 
            IF ( C_H2O > 0.0 ) THEN
               I2G = ( CLDICE(K) / C_H2O ) * CONV_NH3
            ELSE
               I2G = 0.0
            ENDIF
            ! Compute liquid to gas ratio for NH3, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            CALL COMPUTE_L2G( 3.3d6, -4.1d3, TMP(K),CLDLIQ(K),L2G )
            ! Fraction of NH3 in liquid & ice phases
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1.0 + L2G + I2G
            F_L   = L2G / C_TOT
            F_I   = I2G / C_TOT
            ! Compute the rate constant K.  The retention factor  
            ! for liquid NH3 is 0.0 for T <= 248 K and 0.05 for 
            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
            IF ( TMP(K) >= 268d0 ) THEN
               KR = KC * ( F_L + F_I )
            ELSE IF ( TMP(K) > 248d0 .and. TMP(K) < 268d0 ) THEN
               KR = KC * ( ( 5d-2 * F_L ) + F_I )
            ELSE
               KR = KC * F_I
            ENDIF
            ! F is the fraction of NH3 scavenged out of the updraft
            ! (Eq. 2, Jacob et al, 2000)
            IF (VUD(K).GT.0) THEN
             F(K,1) = 1.0 - EXP( -KR * DZQ(K) / Vud(K) )
            ELSE
             F(k,1) = 0.0
            ENDIF
            IF (VDD(K).GT.0) THEN
             F(k,2) = 1.0 - EXP( -KR * DZQ(K) / Vdd(K) )
            ELSE
             F(k,2) = 0.0
            ENDIF
         ENDDO

      ENDIF

 end subroutine CUSCAV_KF
!===========================================================================

 SUBROUTINE F_AEROSOL( KC,DZQ,CSR,F,VUD,VDD,         &
                       ids,ide, jds,jde, kds,kde,       &
                       ims,ime, jms,jme, kms,kme,       &
                       its,ite, jts,jte, kts,kte        )
!
!******************************************************************************
!  Subroutine F_AEROSOL returns the fraction of aerosol scavenged in updrafts
!
!  Arguments as Input:
!  (1 ) KC (REAL) : Conversion rate from cloud condensate to precip [s^-1]
!
!  Arguments as Output:
!  (2 ) F  (REAL) : Fraction of aerosol scavenged in updrafts [unitless]
!******************************************************************************
!
   IMPLICIT NONE


   ! Arguments
   INTEGER,INTENT(IN)  ::    &
                          ids,ide, jds,jde, kds,kde,    &
                          ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte
   REAL, INTENT(IN)  :: KC
   REAL, INTENT(IN)  :: CSR !aerosol activation rate 
   REAL, INTENT(IN)  :: DZQ(kms:kme)
   REAL, INTENT(IN)  :: VUD(kms:kme)
   REAL, INTENT(IN)  :: VDD(kms:kme)
   REAL, INTENT(INOUT) :: F(kms:kme,2)

   ! Local variables
   INTEGER           :: I, J, K
   REAL              :: KR

   !=================================================================
   ! F_AEROSOL begins here!
   ! we set K = Kc*CSR, and compute F accordingly (cf Jacob et al 2000 )    
   !=================================================================
   KR=KC*CSR

   ! Turn off scavenging in the first level by setting F = 0
   F(kts,1:2) = 0.0

   ! Apply scavenging in levels 2 and higher
   DO k = kts+1, kte 
     ! (Eq. 2, Jacob et al, 2000, with K = Kc)
     IF (VUD(K).GT.0) THEN 
      F(k,1) = 1.0 - EXP( -KR * DZQ(K) / Vud(K) )
     ELSE 
      F(k,1) = 0.0
     ENDIF
     IF (VDD(K).GT.0) THEN 
      F(k,2) = 1.0 - EXP( -KR * DZQ(K) / Vdd(K) )
     ELSE  
      F(k,2) = 0.0
     ENDIF
   ENDDO

 END SUBROUTINE F_AEROSOL

END MODULE module_ctrans_kfeta

