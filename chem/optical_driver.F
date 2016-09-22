!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! MOSAIC module: see module_mosaic_driver.F for references and terms of use
!**********************************************************************************  
!
! WRF-chem V3.0 : Original version of optical_driver written by Jerome Fast (PNNL)
!                 and James Barnard (PNNL)
!
!WRF:MODEL_LAYER:CHEMISTRY
!
      SUBROUTINE optical_driver(id,curr_secs,dtstep,config_flags,haveaer,&
               chem,dz8w,alt,relhum,                                     &
               h2oai,h2oaj,h2oai_cz,h2oaj_cz,h2oac_cz,                &
               tauaersw,absaersw,extaersw,gaersw,waersw,            &
               bscoef1,bscoef2,bscoef3,bscoef4,                          &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                      &
               totoa_a01,totoa_a02,totoa_a03,totoa_a04,                  &
               totoa_a05,totoa_a06,totoa_a07,totoa_a08,                  &
               tauaerlw,  & 
               tauaersw_dust,absaersw_dust,gaersw_dust,waersw_dust,                    &
               tauaersw_bc,absaersw_bc,gaersw_bc,waersw_bc,                    &
               tauaersw_sulf,absaersw_sulf,gaersw_sulf,waersw_sulf,            &
               tauaersw_oc,absaersw_oc,gaersw_oc,waersw_oc,                    &
               tauaersw_other,absaersw_other,gaersw_other,waersw_other,                    &
               tauaerlw_dust,tauaerlw_bc,tauaerlw_sulf,tauaerlw_oc,tauaerlw_other,  & 
               ids,ide, jds,jde, kds,kde,                                &
               ims,ime, jms,jme, kms,kme,                                &
               its,ite, jts,jte, kts,kte                                 )

!------------------------------------------------------------------------
   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_optical_averaging
   USE module_data_mosaic_therm, only: nbin_a
   USE module_data_rrtmgaeropt, only: nswbands,nlwbands 
   USE module_peg_util, only:  peg_error_fatal, peg_message
   use infnan,                 only: inf
   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) :: id,                                  &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL,         INTENT(IN   ) :: dtstep
!
! array that holds all advected chemical species
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) ::  chem
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(IN ) ::  relhum, dz8w, alt

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(INOUT) :: h2oai, h2oaj,h2oai_cz,h2oaj_cz,h2oac_cz
                         
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(IN) :: totoa_a01, totoa_a02, totoa_a03, totoa_a04,  &
                       totoa_a05,totoa_a06,totoa_a07,totoa_a08


!
! arrays that hold the aerosol optical properties
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         INTENT(INOUT ) ::                                             &
           bscoef1, bscoef2, bscoef3, bscoef4                              
   !for rrtmg shortwave and longwave --czhao
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,num_tauaersw),INTENT(INOUT ) ::  &
         tauaersw,absaersw,extaersw,gaersw,waersw 
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,num_tauaerlw),INTENT(INOUT ) ::  &
         tauaerlw 

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,num_tauaersw),INTENT(INOUT ) ::  &
         tauaersw_dust,absaersw_dust,gaersw_dust,waersw_dust, & 
         tauaersw_sulf,absaersw_sulf,gaersw_sulf,waersw_sulf, & 
         tauaersw_other,absaersw_other,gaersw_other,waersw_other, & 
         tauaersw_oc,absaersw_oc,gaersw_oc,waersw_oc, & 
         tauaersw_bc,absaersw_bc,gaersw_bc,waersw_bc
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,num_tauaerlw),INTENT(INOUT ) ::  &
         tauaerlw_dust,tauaerlw_bc,tauaerlw_sulf,tauaerlw_oc,tauaerlw_other

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:4) ::  & 
         tauaersw0,absaersw0,extaersw0,gaersw0,waersw0,bscoefsw0 
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ),                  &  
         INTENT(INOUT ) ::                                             &
           l2aer, l3aer, l4aer, l5aer, l6aer, l7aer

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:nlwbands) ::  &
         extaerlw0,tauaerlw0 

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:4) ::  &
             tauaersw_bc0,absaersw_bc0,gaersw_bc0,waersw_bc0, &
             tauaersw_sulf0,absaersw_sulf0,gaersw_sulf0,waersw_sulf0, &
             tauaersw_other0,absaersw_other0,gaersw_other0,waersw_other0, &
             tauaersw_oc0,absaersw_oc0,gaersw_oc0,waersw_oc0, &
             tauaersw_dust0,absaersw_dust0,gaersw_dust0,waersw_dust0
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:nlwbands) ::  &
             tauaerlw_bc0,tauaerlw_dust0,tauaerlw_sulf0,tauaerlw_oc0,tauaerlw_other0

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   character*100 msg
   integer lunerr 

   LOGICAL, INTENT(IN) :: haveaer
!         
! local variables
!
      logical processingAerosols
      integer nbin_o
      integer option_method, option_mie
      integer i,j,k

!-----------------------------------------------------------------
! compute only if simulating aerosols and aer_ra_feedback=1

!czhao change here, so even aer_ra_feedback==0, AOD is still calculated
!  IF (config_flags%aer_ra_feedback .eq. 0) THEN
!        call wrf_debug(15,'no feedback, return from optical driver')
!   return
!  ENDIF
   select case (config_flags%chem_opt)
   case ( RADM2SORG,           RADM2SORG_KPP,      RADM2SORG_AQ, RADM2SORG_AQCHEM, &
          GOCART_SIMPLE,       RACMSORG_KPP,       RACMSORG_AQ,  RACMSORG_AQCHEM_KPP, &
          RACM_ESRLSORG_AQCHEM_KPP, RACM_SOA_VBS_KPP,                                      &
          GOCARTRACM_KPP,      GOCARTRADM2_KPP,    GOCARTRADM2,  &
          RACM_ESRLSORG_KPP,   MOZCART_KPP,                      &
          CBMZ_MOSAIC_4BIN,    CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_KPP,   &
          CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, CBMZSORG, CBMZSORG_AQ, &
          CBMZ_MOSAIC_DMS_4BIN,    CBMZ_MOSAIC_DMS_8BIN,   &
          CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &    
          CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE, CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE, &    
          CBMZ_MOSAIC_4BIN_VBS2_KPP,SAPRC99_MOSAIC_4BIN_VBS2_KPP,MOZART_MOSAIC_4BIN_VBS0_KPP , &
          CBMZ_CAM_MAM3_NOAQ,CBMZ_CAM_MAM7_NOAQ,  CBMZ_CAM_MAM3_AQ,  &
          CBMZ_CAM_MAM7_AQ, CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP,SAPRC99_MOSAIC_4BIN_VBS9_KPP, &
          SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_KPP )!BSINGH(12/05/2013): Added for SAPRC 8 bin vbs and added non-aq on 04/03/2014
      processingAerosols = .true.
      call wrf_debug(15,'optical driver: process aerosols true')
   case default
      processingAerosols = .false.
      call wrf_debug(15,'optical driver: process aerosols false')
   end select

  if( processingAerosols ) then
!
! select aerosol optical property option
! VOLUME: volume averaging of refractive indicies
! * for MADE/SORGAM, assume same 8 size bins as MOSAIC by default
! SHELL: shell-core approach, placeholder
!
   select case (config_flags%chem_opt)
   case ( RADM2SORG, RACM_ESRLSORG_KPP, RADM2SORG_KPP, RADM2SORG_AQ, RADM2SORG_AQCHEM, &
          GOCARTRACM_KPP,      GOCARTRADM2_KPP,    GOCARTRADM2,      &
          GOCART_SIMPLE,       RACMSORG_KPP,       RACMSORG_AQ,      RACMSORG_AQCHEM_KPP, &
          RACM_ESRLSORG_AQCHEM_KPP, RACM_SOA_VBS_KPP,                                          &
          CBMZSORG,            CBMZSORG_AQ,        MOZCART_KPP,       &
          CBMZ_CAM_MAM3_NOAQ,  CBMZ_CAM_MAM7_NOAQ,  CBMZ_CAM_MAM3_AQ,  CBMZ_CAM_MAM7_AQ)
     nbin_o = 8
   case (CBMZ_MOSAIC_4BIN,    CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_KPP,  &
         CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
         CBMZ_MOSAIC_DMS_4BIN,    CBMZ_MOSAIC_DMS_8BIN,   &
         CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &    
         CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE, CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE, &    
         CBMZ_MOSAIC_4BIN_VBS2_KPP,SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
         MOZART_MOSAIC_4BIN_VBS0_KPP, CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP, &
         SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_KPP )!BSINGH(12/05/2013): Added for SAPRC 8 bin vbs and added non-aq on 04/03/2014
     nbin_o = nbin_a
   end select
!
     call wrf_debug(15,'optical averaging')
     aer_op_opt_select: SELECT CASE(config_flags%aer_op_opt)
     CASE (VOLUME_APPROX)
       option_method=1
       option_mie=1
     CASE (MAXWELL_APPROX)
       option_method=2
       option_mie=1
     CASE (VOLUME_EXACT)
       option_method=1
       option_mie=2
     CASE (MAXWELL_EXACT)
       option_method=2
       option_mie=2
     CASE (SHELL_EXACT)
       option_method=3
       option_mie=2
     CASE DEFAULT
        if( config_flags%aer_op_opt > 0 ) then
           call wrf_message('WARNING: Invalid aer_op_opt. Defaulting to VOLUME_APPROX.')
           option_method=1
           option_mie=1
        end if
     END SELECT aer_op_opt_select

     if( config_flags%aer_op_opt > 0 ) then
        call wrf_debug(15,'optical driver: call optical averaging')
!        lunerr=-1
!        write( msg, '(a, 6i4)' )	&
!                 'jdf ', ids, ide, jds, jde, kds, kde
!                 call peg_message( lunerr, msg )
!        write( msg, '(a, 6i4)' )	&
!                 'jdf ', ims, ime, jms, jme, kms, kme
!                 call peg_message( lunerr, msg )
!        write( msg, '(a, 6i4)' )	&
!                 'jdf ', its, ite, jts, jte, kts, kte
!                 call peg_message( lunerr, msg )

        !BSINGH(PNNL)- The followingvariables had an undefined behavior at the boundaries.
        tauaersw0(its:ite,kts:kte,jts:jte,1:4) = 0.0 
        extaersw0(its:ite,kts:kte,jts:jte,1:4) = 0.0
        gaersw0(its:ite,kts:kte,jts:jte,1:4)   = 0.0 
        waersw0(its:ite,kts:kte,jts:jte,1:4)   = 0.0 
        bscoefsw0(its:ite,kts:kte,jts:jte,1:4) = 0.0
        !long wave
        extaerlw0(its:ite,kts:kte,jts:jte,1:nlwbands) = 0.0 
        tauaerlw0(its:ite,kts:kte,jts:jte,1:nlwbands) = 0.0

        call optical_averaging(id,curr_secs,dtstep,config_flags,     &
             nbin_o,haveaer,option_method,option_mie,chem,dz8w,alt,  &
             relhum,h2oai,h2oaj,h2oai_cz,h2oaj_cz,h2oac_cz,         &
!czhao       tauaer1,tauaer2,tauaer3,tauaer4,                        &
!            gaer1,gaer2,gaer3,gaer4,                                &
!            waer1,waer2,waer3,waer4,                                &
!            bscoef1,bscoef2,bscoef3,bscoef4,                        &
             tauaersw0,absaersw0,extaersw0,gaersw0,waersw0,bscoefsw0,               &
             l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                    &
             tauaerlw0,extaerlw0,                                      &
             tauaersw_bc0,absaersw_bc0,gaersw_bc0,waersw_bc0, &
             tauaersw_dust0,absaersw_dust0,gaersw_dust0,waersw_dust0,   &
             tauaersw_sulf0,absaersw_sulf0,gaersw_sulf0,waersw_sulf0,   &
             tauaersw_oc0,absaersw_oc0,gaersw_oc0,waersw_oc0,   &
             tauaersw_other0,absaersw_other0,gaersw_other0,waersw_other0,   &
             tauaerlw_bc0,tauaerlw_dust0,tauaerlw_sulf0,tauaerlw_oc0,tauaerlw_other0,  &
             totoa_a01,totoa_a02,totoa_a03,totoa_a04,                &
             totoa_a05,totoa_a06,totoa_a07,totoa_a08,                &
             ids,ide, jds,jde, kds,kde,                              &
             ims,ime, jms,jme, kms,kme,                              &
             its,ite, jts,jte, kts,kte                               )
             !short wave
             tauaersw(its:ite,kts:kte,jts:jte,2:num_tauaersw)=tauaersw0(its:ite,kts:kte,jts:jte,1:4)
             absaersw(its:ite,kts:kte,jts:jte,2:num_tauaersw)=absaersw0(its:ite,kts:kte,jts:jte,1:4)
             extaersw(its:ite,kts:kte,jts:jte,2:num_extaersw)=extaersw0(its:ite,kts:kte,jts:jte,1:4)
             gaersw(its:ite,kts:kte,jts:jte,2:num_gaersw)=gaersw0(its:ite,kts:kte,jts:jte,1:4)
             waersw(its:ite,kts:kte,jts:jte,2:num_waersw)=waersw0(its:ite,kts:kte,jts:jte,1:4)
             bscoef1(its:ite,kts:kte,jts:jte)=bscoefsw0(its:ite,kts:kte,jts:jte,1)
             bscoef2(its:ite,kts:kte,jts:jte)=bscoefsw0(its:ite,kts:kte,jts:jte,2)
             bscoef3(its:ite,kts:kte,jts:jte)=bscoefsw0(its:ite,kts:kte,jts:jte,3)
             bscoef4(its:ite,kts:kte,jts:jte)=bscoefsw0(its:ite,kts:kte,jts:jte,4)

             tauaersw_dust(its:ite,kts:kte,jts:jte,2:num_tauaersw)=tauaersw_dust0(its:ite,kts:kte,jts:jte,1:4)
             absaersw_dust(its:ite,kts:kte,jts:jte,2:num_tauaersw)=absaersw_dust0(its:ite,kts:kte,jts:jte,1:4)
             gaersw_dust(its:ite,kts:kte,jts:jte,2:num_gaersw)=gaersw_dust0(its:ite,kts:kte,jts:jte,1:4)
             waersw_dust(its:ite,kts:kte,jts:jte,2:num_waersw)=waersw_dust0(its:ite,kts:kte,jts:jte,1:4)
             tauaersw_bc(its:ite,kts:kte,jts:jte,2:num_tauaersw)=tauaersw_bc0(its:ite,kts:kte,jts:jte,1:4)
             absaersw_bc(its:ite,kts:kte,jts:jte,2:num_tauaersw)=absaersw_bc0(its:ite,kts:kte,jts:jte,1:4)
             gaersw_bc(its:ite,kts:kte,jts:jte,2:num_gaersw)=gaersw_bc0(its:ite,kts:kte,jts:jte,1:4)
             waersw_bc(its:ite,kts:kte,jts:jte,2:num_waersw)=waersw_bc0(its:ite,kts:kte,jts:jte,1:4)
             tauaersw_sulf(its:ite,kts:kte,jts:jte,2:num_tauaersw)=tauaersw_sulf0(its:ite,kts:kte,jts:jte,1:4)
             absaersw_sulf(its:ite,kts:kte,jts:jte,2:num_tauaersw)=absaersw_sulf0(its:ite,kts:kte,jts:jte,1:4)
             gaersw_sulf(its:ite,kts:kte,jts:jte,2:num_gaersw)=gaersw_sulf0(its:ite,kts:kte,jts:jte,1:4)
             waersw_sulf(its:ite,kts:kte,jts:jte,2:num_waersw)=waersw_sulf0(its:ite,kts:kte,jts:jte,1:4)
             tauaersw_other(its:ite,kts:kte,jts:jte,2:num_tauaersw)=tauaersw_other0(its:ite,kts:kte,jts:jte,1:4)
             absaersw_other(its:ite,kts:kte,jts:jte,2:num_tauaersw)=absaersw_other0(its:ite,kts:kte,jts:jte,1:4)
             gaersw_other(its:ite,kts:kte,jts:jte,2:num_gaersw)=gaersw_other0(its:ite,kts:kte,jts:jte,1:4)
             waersw_other(its:ite,kts:kte,jts:jte,2:num_waersw)=waersw_other0(its:ite,kts:kte,jts:jte,1:4)
             tauaersw_oc(its:ite,kts:kte,jts:jte,2:num_tauaersw)=tauaersw_oc0(its:ite,kts:kte,jts:jte,1:4)
             absaersw_oc(its:ite,kts:kte,jts:jte,2:num_tauaersw)=absaersw_oc0(its:ite,kts:kte,jts:jte,1:4)
             gaersw_oc(its:ite,kts:kte,jts:jte,2:num_gaersw)=gaersw_oc0(its:ite,kts:kte,jts:jte,1:4)
             waersw_oc(its:ite,kts:kte,jts:jte,2:num_waersw)=waersw_oc0(its:ite,kts:kte,jts:jte,1:4)

             !long wave
             tauaerlw(its:ite,kts:kte,jts:jte,2:num_tauaerlw)=tauaerlw0(its:ite,kts:kte,jts:jte,1:nlwbands)
             tauaerlw_dust(its:ite,kts:kte,jts:jte,2:num_tauaerlw)=tauaerlw_dust0(its:ite,kts:kte,jts:jte,1:nlwbands)
             tauaerlw_bc(its:ite,kts:kte,jts:jte,2:num_tauaerlw)=tauaerlw_bc0(its:ite,kts:kte,jts:jte,1:nlwbands)
             tauaerlw_sulf(its:ite,kts:kte,jts:jte,2:num_tauaerlw)=tauaerlw_sulf0(its:ite,kts:kte,jts:jte,1:nlwbands)
             tauaerlw_other(its:ite,kts:kte,jts:jte,2:num_tauaerlw)=tauaerlw_other0(its:ite,kts:kte,jts:jte,1:nlwbands)
             tauaerlw_oc(its:ite,kts:kte,jts:jte,2:num_tauaerlw)=tauaerlw_oc0(its:ite,kts:kte,jts:jte,1:nlwbands)
     else
        !If aer_op_opt==0 then the optical arrays are already set to
        !zero in chemics_init so there will not be a problem if the
        !user has selected aer_ra_feedback=1.
     end if
!
   endif

   return

END SUBROUTINE optical_driver
