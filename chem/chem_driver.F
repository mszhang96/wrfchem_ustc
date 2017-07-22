
!
!   WRF-chem V3.1 : Original version of chem_driver written by Georg Grell (ESRL/GSD)
!                   Further developments, bugfixes and improvements  by
!                   William Gustafson (PNNL), Marc Salzmann (GFDL), and Georg Grell
! 10/12/2011 - Ravan Ahmadov (NOAA) updated to include the RACM_SOA_VBS option
!
#if ( NMM_CORE == 1 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NCEP_MESO:MEDIATION_LAYER:SOLVER
!
!-----------------------------------------------------------------------
#include "../dyn_nmm/nmm_loop_basemacros.h"
#include "../dyn_nmm/nmm_loop_macros.h"
!-----------------------------------------------------------------------
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine chem_driver ( grid , config_flags   &
 
#include "dummy_new_args.inc"
 
                 )
!----------------------------------------------------------------------
  USE module_domain , only : domain
  USE module_configure
#if ( EM_CORE == 1 )
  USE module_driver_constants
  USE module_machine
  USE module_tiles
#endif
  USE module_dm
  USE module_model_constants
  USE module_state_description
#if ( NMM_CORE == 1 )
  USE MODULE_PHYSICS_CALLS 
#endif
  USE module_data_radm2
  USE module_data_sorgam
  USE module_radm
  USE module_dep_simple
  USE module_bioemi_simple
  USE module_phot_mad
  USE module_ftuv_driver, only : ftuv_timestep_init
  USE module_aerosols_sorgam
  USE module_chem_utilities
  USE module_gocart_so2so4
  USE module_aer_opt_out,only: aer_opt_out
  USE module_ctrans_grell
  USE module_ctrans_kfeta
! USE module_aerosols_soa_vbs, only:
  USE module_data_soa_vbs, only: ldrog_vbs
!
  USE module_chem_cup, only: chem_cup_driver !BSINGH - For WRFCuP scheme
  USE module_dry_dep_driver
  USE module_emissions_driver
  USE module_input_tracer, only: set_tracer
  USE module_wetscav_driver, only: wetscav_driver
  USE module_wetdep_ls, only:wetdep_ls
  USE module_input_chem_data, only: last_chem_time, &
#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
                                     chem_dbg,      &
#endif
                                     get_last_gas,mozcart_lbc_set
  USE module_upper_bc_driver, only: upper_bc_driver
  USE module_tropopause,      only: tropopause_driver
  USE modal_aero_data, only: ntot_amode_cam_mam => ntot_amode
  USE module_cam_support, only: gas_pcnst => gas_pcnst_modal_aero,gas_pcnst_pos => gas_pcnst_modal_aero_pos, &
       pcnst =>pcnst_runtime, numgas_mam, cam_mam_aerosols
  USE module_cu_camzm_driver, only: zm_conv_tend_2
  USE module_cam_mam_gas_wetdep_driver, only: cam_mam_gas_wetdep_driver

  IMPLICIT NONE


  !BSINGH(PNNL)- Lahey compiler forces to declare the following interface
  interface
     SUBROUTINE sum_pm_driver ( config_flags,                            &
          alt,         chem,        h2oaj,        h2oai, h2oac,                &
            pm2_5_dry, pm2_5_water,                                    &
            pm2_5_dry_dust,pm2_5_dry_oc,pm2_5_dry_ec,                  &
            pm2_5_dry_seas,pm2_5_dry_so4f,pm2_5_dry_no3f,pm2_5_dry_nh4f,    &
            pm10,                                                      &
            !czhao
            dgnuc,dgacc,dgcor,                                      & 
          hoa_a01,     hoa_a02,     hoa_a03,      hoa_a04,               &
          hoa_a05,     hoa_a06,     hoa_a07,      hoa_a08,               & !BSINGH Added 4 more bins for each species for SAPRC 8 bin version
          bboa_a01,    bboa_a02,    bboa_a03,     bboa_a04,              &
          bboa_a05,    bboa_a06,    bboa_a07,     bboa_a08,              &
          soa_a01,     soa_a02,     soa_a03,      soa_a04,               &
          soa_a05,     soa_a06,     soa_a07,      soa_a08,               &
          bbsoa_a01,   bbsoa_a02,   bbsoa_a03,    bbsoa_a04,             &
          bbsoa_a05,   bbsoa_a06,   bbsoa_a07,    bbsoa_a08,             &
          hsoa_a01,    hsoa_a02,    hsoa_a03,     hsoa_a04,              &
          hsoa_a05,    hsoa_a06,    hsoa_a07,     hsoa_a08,              &
          biog_a01,    biog_a02,    biog_a03,     biog_a04,              &
          biog_a05,    biog_a06,    biog_a07,     biog_a08,              &
          asmpsoa_a01, asmpsoa_a02, asmpsoa_a03, asmpsoa_a04,            & !BSINGH - Not adding 5-8 bins for asmpsoa, as it is not req. for 8 bin SAPRC
          arosoa_a01,  arosoa_a02,  arosoa_a03,  arosoa_a04,             &
          arosoa_a05,  arosoa_a06,  arosoa_a07,  arosoa_a08,             &
          totoa_a01,   totoa_a02,   totoa_a03,   totoa_a04,              &
          totoa_a05,   totoa_a06,   totoa_a07,   totoa_a08,              &
          hsoa_c,      hsoa_o,      bbsoa_c,     bbsoa_o,                &
          biog_v1,     biog_v2,     biog_v3,     biog_v4,                &
          ant_v1,      ant_v2,      ant_v3,      ant_v4,                 &
          smpa_v1,     smpbb_v1,                                         &
          !BSINGH - Added cw aerosols(for VBS)
          hoa_cw01,    hoa_cw02,    hoa_cw03,    hoa_cw04,               &
          hoa_cw05,    hoa_cw06,    hoa_cw07,    hoa_cw08,               &
          bboa_cw01,   bboa_cw02,   bboa_cw03,   bboa_cw04,              &
          bboa_cw05,   bboa_cw06,   bboa_cw07,   bboa_cw08,              &
          soa_cw01,    soa_cw02,    soa_cw03,    soa_cw04,               &
          soa_cw05,    soa_cw06,    soa_cw07,    soa_cw08,               &
          bbsoa_cw01,  bbsoa_cw02,  bbsoa_cw03,  bbsoa_cw04,             &
          bbsoa_cw05,  bbsoa_cw06,  bbsoa_cw07,  bbsoa_cw08,             &
          biog_cw01,   biog_cw02,   biog_cw03,   biog_cw04,              &
          biog_cw05,   biog_cw06,   biog_cw07,   biog_cw08,              &
          hsoa_cw01,   hsoa_cw02,   hsoa_cw03,   hsoa_cw04,              &
          hsoa_cw05,   hsoa_cw06,   hsoa_cw07,   hsoa_cw08,              &
          arosoa_cw01, arosoa_cw02, arosoa_cw03, arosoa_cw04,            &
          arosoa_cw05, arosoa_cw06, arosoa_cw07, arosoa_cw08,            &
          totoa_cw01,  totoa_cw02,  totoa_cw03,  totoa_cw04,             &
          totoa_cw05,  totoa_cw06,  totoa_cw07,  totoa_cw08,             &        
          hsoa_cw_c,   hsoa_cw_o,   bbsoa_cw_c,  bbsoa_cw_o,             &
          biog_cw_v1,                                                    &
          ant_cw_v1,                                                     &
          !BSINGH -ENDS
          ids,ide, jds,jde, kds,kde,                                     &
          ims,ime, jms,jme, kms,kme,                                     &
          its,ite, jts,jte, kts,kte                                      )

       
       
       USE module_configure
       USE module_aerosols_sorgam, only: sum_pm_sorgam
       USE module_mosaic_driver, only: sum_pm_mosaic,sum_pm_mosaic_vbs2,sum_pm_mosaic_vbs0,sum_vbs9,sum_vbs2,sum_vbs0
       USE module_gocart_aerosols, only: sum_pm_gocart
       USE module_aerosols_soa_vbs, only: sum_pm_soa_vbs
       
       IMPLICIT NONE
       
       INTEGER,      INTENT(IN   )    ::                                   &
            ids,ide, jds,jde, kds,kde,       &
            ims,ime, jms,jme, kms,kme,       &
            its,ite, jts,jte, kts,kte
       
       REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
            INTENT(IN ) :: chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(IN ) ::                                             &
           dgnuc, dgacc,dgcor
       
       REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
            INTENT(IN ) :: alt
       REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
            OPTIONAL,                                                     &
            INTENT(IN ) :: h2oaj,h2oai,h2oac
       
       REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         OPTIONAL,                                                     &
         INTENT(INOUT) :: &
                        pm2_5_dry,pm2_5_water,pm2_5_dry_dust,pm2_5_dry_oc,&
                        pm2_5_dry_ec,pm2_5_dry_seas,pm2_5_dry_so4f, &
                        pm2_5_dry_no3f,pm2_5_dry_nh4f, pm10,           &             
             hoa_a01,     hoa_a02,     hoa_a03,      hoa_a04,          &
             hoa_a05,     hoa_a06,     hoa_a07,      hoa_a08,          &!BSINGH Added 4 more bins for each species for SAPRC 8 bin version             
             bboa_a01,    bboa_a02,    bboa_a03,     bboa_a04,         &
             bboa_a05,    bboa_a06,    bboa_a07,     bboa_a08,         &             
             soa_a01,     soa_a02,     soa_a03,      soa_a04,          &
             soa_a05,     soa_a06,     soa_a07,      soa_a08,          &            
             bbsoa_a01,   bbsoa_a02,   bbsoa_a03,    bbsoa_a04,        &
             bbsoa_a05,   bbsoa_a06,   bbsoa_a07,    bbsoa_a08,        &             
             hsoa_a01,    hsoa_a02,    hsoa_a03,     hsoa_a04,         &
             hsoa_a05,    hsoa_a06,    hsoa_a07,     hsoa_a08,         &             
             biog_a01,    biog_a02,    biog_a03,     biog_a04,         &
             biog_a05,    biog_a06,    biog_a07,     biog_a08,         &             
             arosoa_a01,  arosoa_a02,  arosoa_a03,   arosoa_a04,       &
             arosoa_a05,  arosoa_a06,  arosoa_a07,   arosoa_a08,       &             
             totoa_a01,   totoa_a02,   totoa_a03,    totoa_a04,        &
             totoa_a05,   totoa_a06,   totoa_a07,    totoa_a08,        &             
             hsoa_c,      hsoa_o,      bbsoa_c,      bbsoa_o,          &
             biog_v1,     biog_v2,     biog_v3,      biog_v4,          &
             ant_v1,      ant_v2,      ant_v3,       ant_v4,           &
             smpa_v1,                                                  &
             smpbb_v1,                                                 &
             asmpsoa_a01, asmpsoa_a02, asmpsoa_a03, asmpsoa_a04,       &!BSINGH - Not adding 5-8 bins for asmpsoa, as it is not req. for 8 bin SAPRC
             !BSINGH - Added cw aerosols(for VBS)
             hoa_cw01,    hoa_cw02,    hoa_cw03,    hoa_cw04,          &
             hoa_cw05,    hoa_cw06,    hoa_cw07,    hoa_cw08,          &
             bboa_cw01,   bboa_cw02,   bboa_cw03,   bboa_cw04,         &
             bboa_cw05,   bboa_cw06,   bboa_cw07,   bboa_cw08,         &
             soa_cw01,    soa_cw02,    soa_cw03,    soa_cw04,          &
             soa_cw05,    soa_cw06,    soa_cw07,    soa_cw08,          &
             bbsoa_cw01,  bbsoa_cw02,  bbsoa_cw03,  bbsoa_cw04,        &
             bbsoa_cw05,  bbsoa_cw06,  bbsoa_cw07,  bbsoa_cw08,        &
             biog_cw01,   biog_cw02,   biog_cw03,   biog_cw04,         &
             biog_cw05,   biog_cw06,   biog_cw07,   biog_cw08,         &
             hsoa_cw01,   hsoa_cw02,   hsoa_cw03,   hsoa_cw04,         &
             hsoa_cw05,   hsoa_cw06,   hsoa_cw07,   hsoa_cw08,         &
             arosoa_cw01, arosoa_cw02, arosoa_cw03, arosoa_cw04,       &
             arosoa_cw05, arosoa_cw06, arosoa_cw07, arosoa_cw08,       &
             totoa_cw01,  totoa_cw02,  totoa_cw03,  totoa_cw04,        &
             totoa_cw05,  totoa_cw06,  totoa_cw07,  totoa_cw08,        &        
             hsoa_cw_c,   hsoa_cw_o,   bbsoa_cw_c,  bbsoa_cw_o,        &
             biog_cw_v1,                                               &
             ant_cw_v1                                                 
             !BSINGH -ENDS

       TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
end SUBROUTINE sum_pm_driver


SUBROUTINE photolysis_driver (id,curr_secs,ktau,dtstep,            &
               config_flags,haveaer,                                     &
               gmt,julday,t_phy,moist,aerwrf,p8w,t8w,p_phy,              &
               chem,rho_phy,dz8w,xlat,xlong,z_at_w,gd_cloud,gd_cloud2,   &
               ph_macr,ph_o31d,ph_o33p,ph_no2,ph_no3o2,                  &
               ph_no3o,ph_hno2,ph_hno3,ph_hno4,ph_h2o2,                  &
               ph_ch2or,ph_ch2om,ph_ch3cho,ph_ch3coch3,                  &
               ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,                       &
               ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,             &
               ph_hcochob, ph_n2o5,ph_o2,ph_n2o,                         &
               ph_pan,ph_mpan,ph_acetol,ph_gly,                          &
               ph_bigald,ph_mek,ph_c2h5ooh,ph_c3h7ooh,ph_pooh,           &
               ph_rooh,ph_xooh,ph_isopooh,ph_alkooh,                     &
               ph_mekooh,ph_tolooh,ph_terpooh,ph_mvk,                    &
               ph_glyald,ph_hyac,                                        &
               nref0, nw0, tuv_jmax0,                                    &
               ph_radfld, ph_adjcoe, ph_prate,                           &
               wc, zref,                                                 &
               tauaer1,tauaer2,tauaer3,tauaer4,                          &
               gaer1,gaer2,gaer3,gaer4,                                  &
               waer1,waer2,waer3,waer4,                                  &
               bscoef1,bscoef2,bscoef3,bscoef4,                          &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                      &
               pm2_5_dry,pm2_5_water,uvrad,ivgtyp,                       &
               ids,ide, jds,jde, kds,kde,                                &
               ims,ime, jms,jme, kms,kme,                                &
               its,ite, jts,jte, kts,kte                                 )

!------------------------------------------------------------------------
   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_phot_mad
   USE module_phot_fastj
   USE module_ftuv_driver

   INTEGER,      INTENT(IN   ) :: id,julday,                           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) :: ktau
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL,         INTENT(IN   ) :: dtstep,gmt
!
! advected moisture variables
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),            &
         INTENT(IN ) ::                                   moist
!
!  aerosol interaction
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::                                             &
               pm2_5_dry,pm2_5_water, aerwrf
!
! arrays that hold the photolysis rates
!

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::                                             &
           ph_macr,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,    &
           ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,        &
           ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,             &
           ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,    &
           ph_n2o5,ph_o2,ph_n2o,ph_pan,ph_mpan,ph_acetol,ph_gly,       &
           ph_bigald,ph_mek,ph_c2h5ooh,ph_c3h7ooh,ph_pooh,ph_rooh,     &
           ph_xooh,ph_isopooh,ph_alkooh,ph_mekooh,ph_tolooh,           &
           ph_terpooh,ph_mvk,ph_glyald,ph_hyac

   INTEGER,      INTENT(IN   ) :: nref0, nw0, tuv_jmax0 
   real, dimension( ims:ime, nref0, jms:jme, nw0 ),                    &
            intent(out  ) :: ph_radfld
   real, dimension( ims:ime, nref0, jms:jme, tuv_jmax0 ),              &
            intent(out  ) :: ph_adjcoe
   real, dimension( ims:ime, nref0, jms:jme, tuv_jmax0 ),              &
            intent(out  ) :: ph_prate
   real, dimension(nw0),                                               &
            intent(out  ) :: wc
   real, dimension(nref0),                                             &
            intent(out  ) :: zref
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &  
         OPTIONAL,                                                     &
         INTENT(INOUT ) ::                                             &
           gd_cloud,gd_cloud2
!
! arrays that hold the aerosol optical properties
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(IN ) ::                                                &
           tauaer1,tauaer2,tauaer3,tauaer4,                            &
           gaer1,gaer2,gaer3,gaer4,                                    &
           waer1,waer2,waer3,waer4,                                    &
           bscoef1,bscoef2,bscoef3,bscoef4
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ),                  &
         INTENT(IN ) ::                                                &
           l2aer,l3aer,l4aer,l5aer,l6aer,l7aer
!
! array that holds all advected chemical species
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(INOUT ) ::                                   chem
!
! on input from meteorological part of model
!
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,           &
          INTENT(IN   ) ::                                             &
                                                      t_phy,           &
                                                      p_phy,           &
                                                      dz8w,            &
                                              t8w,p8w,z_at_w ,         &
                                                    rho_phy
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,           &
          INTENT(INOUT   ) ::  uvrad
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,           &
          INTENT(IN   ) ::                                             &
                                                     xlat,             &
                                                     xlong
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   LOGICAL, INTENT(IN)  :: haveaer
   integer, INTENT(IN ) :: ivgtyp(ims:ime,jms:jme)
   character(len=132)   :: dbg_msg

END SUBROUTINE photolysis_driver


subroutine mechanism_driver(id,curr_secs,ktau,dtstep,ktauc,dtstepc,&
               config_flags,gmt,julday,t_phy,moist,p8w,t8w,gd_cldfr,   &
               p_phy,chem,rho_phy,dz8w,dx,g,z,z_at_w,xlat,xlong,vdrog3,&
               vcsulf_old,vcso2_old,vch2o2_old,ttday,tcosz,            &
               ph_macr,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,&
               ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,    &
               ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,         &
               ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,&
               ph_n2o5,ph_o2,backg_oh,backg_h2o2,backg_no3,            &
               addt,addx,addc,etep,oltp,olip,cslp,limp,hc5p,hc8p,tolp, &
               xylp,apip,isop,hc3p,ethp,o3p,tco3,mo2,o1d,olnn,rpho,xo2,&
               ketp,olnd,                                              &
               ids,ide, jds,jde, kds,kde,                              &
               ims,ime, jms,jme, kms,kme,                              &
               its,ite, jts,jte, kts,kte                               )
!----------------------------------------------------------------------
  USE module_configure
  USE module_state_description
  USE module_data_radm2
  USE module_data_sorgam
  USE module_radm
  USE module_gocart_chem
  USE module_aerosols_sorgam
  USE module_cbmz, only:        cbmz_driver
   IMPLICIT NONE


   INTEGER,      INTENT(IN   ) :: id,julday,                           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) :: ktau,ktauc
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL,         INTENT(IN   ) :: dtstep,dtstepc,gmt,dx,g
!
! advected moisture variables
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),            &
         INTENT(IN ) ::                                   moist
!
! advected chemical species
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, jms:jme ),                       &
         INTENT(IN ) ::                                             &
          xlat,xlong,ttday,tcosz

!
! arrays that hold the photolysis rates
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::                                             &
           ph_macr,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,    &
           ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho,        &
           ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho,             &
           ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,    &
           ph_n2o5,ph_o2
!
! RACM radicals
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(INOUT ) ::                                             &
               addt,addx,addc,etep,oltp,olip,cslp,limp,hc5p,hc8p,tolp, &
               xylp,apip,isop,hc3p,ethp,o3p,tco3,mo2,o1d,olnn,rpho,xo2,&
               ketp,olnd


!
! on input from meteorological model 
!
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,           &
          INTENT(IN   ) ::                                             &
                                                      t_phy,           &
                                                      p_phy,           &
                                                    dz8w,              &
                        z    ,backg_oh,backg_h2o2,backg_no3,           &
                                              t8w,p8w,z_at_w ,         &
                                                    rho_phy
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,           &
          OPTIONAL,                                                    &
          INTENT(IN   ) ::  gd_cldfr
! ..
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
!
! for interaction of aerosls/chemistry (MADE/SORGAM only)
!
      real ,      INTENT(INOUT   ) ::                                  &
                      vdrog3(ims:ime,kms:kme-0,jms:jme,ldrog)
      real ,      INTENT(INOUT   ) ::                                  &
                      vcsulf_old(ims:ime,kms:kme,jms:jme),             &
                      vcso2_old(ims:ime,kms:kme,jms:jme),              &
                      vch2o2_old(ims:ime,kms:kme,jms:jme)

END subroutine mechanism_driver


  end interface
  
  
   !  Input data.

   TYPE(domain) , TARGET          :: grid
   !
   !  Definitions of dummy arguments to solve
# include <dummy_new_decl.inc>
#if ( EM_CORE == 1 )
# define NO_I1_OLD
#endif
#if ( NMM_CORE == 1 )
# ifdef DM_PARALLEL
      INCLUDE "mpif.h"
# endif
#endif

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER                     :: ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  ips,ipe, jps,jpe, kps,kpe,    &
                                  its,ite, jts,jte, kts,kte
! ..
! .. Local Scalars ..
   INTEGER :: stepave,i,j,k,l,numgas,nv,n, nr,ktau,k_start,k_end,idf,jdf,kdf
   INTEGER :: ijulian
!  REAL    :: convtrans_avglen_m

! ................................................................
! ..
!
! necessary for aerosols (module dependent)
!
#if ( NMM_CORE == 1 )
      real, dimension(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32) ::vcsulf_old,vcso2_old
      real, dimension(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,ldrog) ::vdrog3

      real, dimension(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,ldrog_vbs) ::vdrog3_vbs 
      
      real, dimension(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32) ::n2o5_het
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32) ::              &
                                                              p_phy,u_phy,v_phy                   &
                                                             ,t_phy,dz8w,t8w,p8w                  &
                                                             ,rho,rri,z_at_w,vvel,zmid
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32) :: pbl_h
      REAL,DIMENSION(grid%sm33:grid%em33-1) :: QL,TL
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32) :: REXNSFC,FACTRS                &
                                        ,TOT,TSFC
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: moist_trans
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: scalar_trans
#endif
#if ( EM_CORE == 1 )
      real, dimension(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33) ::vcsulf_old,vcso2_old,vch2o2_old
      real, dimension(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,ldrog) ::vdrog3

      real, dimension(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,ldrog_vbs) ::vdrog3_vbs

      real, dimension(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33) ::n2o5_het 
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33) ::              &
                                                              p_phy,u_phy,v_phy                   &
                                                             ,t_phy,dz8w,t8w,p8w                  &
                                                             ,rho,rri,z_at_w,vvel,zmid,rh
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33) :: pbl_h
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33, 5) :: seasin,dustin
      REAL,DIMENSION(grid%sm32:grid%em32-1) :: QL,TL
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33) :: REXNSFC,FACTRS                &
                                        ,TOT,TSFC

      ! temporary arrays for old chemistry values
      REAL,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_chem_ct) :: chem_old
      INTEGER,DIMENSION(num_chem_ct) :: chem_ct_indices

#endif
! Variables for adaptive time steps...
      TYPE(WRFU_TimeInterval) :: tmpTimeInterval
      REAL(KIND=8) :: curr_secs
      REAL(KIND=8) :: real_time_r8                                 !ext. function in adapt_timestep_em.F
      LOGICAL      :: adapt_step_flag, do_chemstep, do_photstep

      REAL :: DAYI,DPL,FICE,FRAIN,HOUR,PLYR          &
     &       ,QI,QR,QW,RADT,TIMES,WC,TDUM,WMSK,RWMSK
 

      INTEGER                         :: ij 
      INTEGER                         :: im , num_3d_m , ic , num_3d_c, num_3d_s
      INTEGER                         :: ijds, ijde
      INTEGER                         :: ksubt

      REAL :: chem_minval, dtstepc

      REAL, DIMENSION(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, numgas_mam) :: gas_aqfrac ! fraction of gas that is in cloud water
      
      REAL, DIMENSION( grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, ntot_amode_cam_mam ) :: &
         wetdens_ap  ! modal diameters and wet densities for cam_mam packages

      REAL, DIMENSION( grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33 ) :: &
         del_h2so4_gasprod  ! change to h2so4 during gas-phase chemistry, cam_mam packages

      !Balwinder.Singh@pnnl.gov:dvmrdt_sv13d,dvmrcwdt_sv13d are the tendencies which are passsed on from the CAM-MAM cloud chemistry
      !to gasaerexch subroutine in cam_mam_aerchem_driver
      REAL, DIMENSION( grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33,gas_pcnst_pos) :: dvmrdt_sv13d,dvmrcwdt_sv13d 

      ! temporary for dry deposition budget -czhao
      REAL, DIMENSION( grid%sm31:grid%em31, grid%sm33:grid%em33, 8,3) :: &
      temporary,temporary_cw
      REAL, DIMENSION( grid%sm31:grid%em31, grid%sm33:grid%em33, 8) :: &
      deposition_dryrate

      !for snow albedo
      REAL, DIMENSION( grid%sm31:grid%em31, grid%sm33:grid%em33, 4,5) :: &
      tend_rate_dep,snowchem_old

      LOGICAL :: haveaer
      CHARACTER (LEN=1000) :: msg !BSINGH:01/31/2013: For message string
      CHARACTER (LEN=256) :: current_date_char !shc
      integer :: current_month                 !shc

      !Chun Zhao check the process time 
      real :: startt, finisht

! ..
! .. Intrinsic Functions ..
      INTRINSIC max, min
!<DESCRIPTION>
!<pre>
! chem_driver is the main driver for handling chemistry related tasks
! for a particular timestep. chem_driver is a mediation-layer routine ->
! DM and SM calls are made where needed for parallel processing.
!
! The main sections of chem_driver are:
!
! (1) Initialization of meteorology variables as needed for chemistry
!
! (2) Calls to drivers for the various chemistry tasks:
!        emissions_driver
!        photolysis_driver
!        dry_dep_driver
!        grelldrvct (convective tracer transport)
!        mechanism_driver (gases)
!        cloud_chem_driver
!        aerosols_driver
!        wetscav_driver
!        sum_pm_driver
!
! Handling of tile indices in chem_driver is as close as possible to
! what is done in solve_em. For subroutines called from chem_driver,
! the its:ite, jts:jte, and kts:kte variables represent the extent of
! the domain that each processor should loop over. For example, a do
! loop in the vertical for the chem array should go from kts to kte.
! For the EM core, kte=kde-1. For the NMM core, kte=kde-2.
!
! Note that the tile indices for the chemistry initialization differ
! from the integration loop indices in that the initializataion routines
! use kte=kde. Go figure, this is how the met. folks set things up.
!
!</pre>
!</DESCRIPTION>

! ..

! Number of levels to exclude from the chem calculations counting from
! the model top.
!  ksubt=0

!
! Setup the adaptive timestep for the chem routines. Most of this follows
! what is in solve_em, except for the call to adjust time_step.
!
#if ( NMM_CORE == 1)
  !NMM defaults to the old step counting methodology in physics so we
  !will do the same here in chemistry. In theory, adapt_step_flag can
  !probably be set to TRUE for NMM too using the curr_secs calculated
  !with ktau, but I do not have input files to test NMM. (wig, 12-May-2008)
  adapt_step_flag = .FALSE.
  KTAU=GRID%NTSD
  curr_secs = (ktau-1)*real(grid%dt,8) !I think this breaks around 68 yrs w/ i4
  ijulian=ifix(grid%julian)
#endif
#if ( EM_CORE == 1 )
  !The necesssary variables exist for the EM core and using the adaptive
  !techniques will work even with a constant time step. In fact, they
  !prevent issues with restarts and changed time steps. So, we will
  !always use them with the EM core.
  adapt_step_flag = .TRUE.
  ktau = grid%itimestep
  tmpTimeInterval = domain_get_time_since_sim_start(grid)
  curr_secs = real_time_r8(tmpTimeInterval)
  ijulian=ifix(grid%julian)
#endif

  do_photstep = .false.
  IF ( ktau==1 ) then
     do_photstep = .true.
  ELSE IF ( adapt_step_flag ) THEN
     IF ( (grid%photdt<=0) .or. &
          ( curr_secs+real(grid%dt,8)+0.01 >= &
          ( INT( curr_secs/real(grid%photdt*60.,8)+1,8 )*real(grid%photdt*60.,8) ) ) &
          ) then
          !NOTE: The 0.01 added to the LHS of the conditional is to compensate
          !      for occasional round off errors that prevented the >= from
          !      ever testing as true. This adjustment has been added to the
          !      other checks within the chem directory as well. wig, 30-Sep-2008
        do_photstep = .true.
     ENDIF
  ELSE IF ( (MOD(ktau,grid%stepphot)==0) .or. (grid%stepphot==1) ) THEN
     do_photstep = .true.
  ENDIF

#if (NMM_CORE == 1)
  if( ktau==1 ) then
     dtstepc=grid%dt
  else
     dtstepc=grid%dt*float(grid%stepchem)
  end if
#endif
#if (EM_CORE == 1)
  if( ktau==1 ) then
     dtstepc = grid%dt
  else
     tmpTimeInterval = domain_get_current_time(grid) - last_chem_time(grid%id)
     dtstepc = real(real_time_r8(tmpTimeInterval),4)
  end if
#endif
     
  ! initializing diagnostics and macros

#if (EM_CORE == 1)
  if( ktau==1 ) then
     grid%conv_ct(:,:,:,:)   = 0.
     grid%chem_ct(:,:,:,:)   = 0.
     grid%vmix_ct(:,:,:,:)   = 0.  
  endif
  if(config_flags%chemdiag == USECHEMDIAG)then
  ! modify tendency list here
      chem_ct_indices(p_chem_co)   = p_co
      chem_ct_indices(p_chem_o3)   = p_o3
      chem_ct_indices(p_chem_no)   = p_no
      chem_ct_indices(p_chem_no2)  = p_no2
      chem_ct_indices(p_chem_hno3) = p_hno3
      chem_ct_indices(p_chem_iso)  = p_iso
      chem_ct_indices(p_chem_ho)   = p_ho
      chem_ct_indices(p_chem_ho2)  = p_ho2
  endif
#endif
  
  do_chemstep = .false.
  IF ( ktau==1 ) then
     do_chemstep = .true.
     grid%ktauc = 1
  ELSE IF ( adapt_step_flag ) THEN
     IF ( (grid%chemdt<=0) .or. &
          ( curr_secs+real(grid%dt,8)+0.01 >= &
          ( INT( curr_secs/real(grid%chemdt*60.,8)+1,8 )*real(grid%chemdt*60.,8) ) ) &
          ) then
        do_chemstep = .true.
        grid%ktauc = grid%ktauc+1
        last_chem_time(grid%id) = domain_get_current_time( grid )
        call WRFU_TimeGet( last_chem_time(grid%id),         &
                           YY = grid%last_chem_time_year,   &
                           MM = grid%last_chem_time_month,  &
                           DD = grid%last_chem_time_day,    &
                           H  = grid%last_chem_time_hour,   &
                           M  = grid%last_chem_time_minute, &
                           S  = grid%last_chem_time_second  )
     ENDIF
  ELSE IF ( (MOD(ktau,grid%stepchem)==0) .or. (grid%stepchem==1) ) THEN
     do_chemstep = .true.
     grid%ktauc=max(ktau/grid%stepchem,1)
  ENDIF

! ..
!
! Some stuff for convective transport...
!
!  convtrans_avglen_m = 30.  !Averaging time for convective transport in min.
! stepave=convtrans_avglen_m*60./grid%dt

  CALL get_ijk_from_grid (  grid ,                   &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            ips, ipe, jps, jpe, kps, kpe    )

! following two lines needed for MEGAN
  CALL domain_clock_get( grid, current_timestr=current_date_char )
  read(current_date_char(6:7),FMT='(I2)') current_month
!
!initialize
!
  seasin(:,:,:)=0.
  dustin(:,:,:)=0.

  if(config_flags%cu_diag == 0 ) grid%raincv_b(:,:) = grid%raincv(:,:)

#if ( NMM_CORE == 1 )
!***  IN NMM SET CONTROLS FOR TILES TO PATCHES
!
!-----------------------------------------------------------------------
      IDF=IDE-1
      JDF=JDE-1
      KDF=KDE-2 ! DO NOT do chem at the top level to mimic what used to be done (also prevents a solver failure at kde-1 for MADE/SORGAM)
      ITS=IPS
      ITE=MIN(IPE,IDF)
      JTS=JPS
      JTE=MIN(JPE,JDF)
      KTS=KPS
      KTE=MIN(KPE,KDF)

#endif


  num_3d_m        = num_moist
  num_3d_c        = num_chem
  num_3d_s        = num_scalar
  numgas          = get_last_gas(config_flags%chem_opt)



#if ( EM_CORE == 1 )

   !  Compute these starting and stopping locations for each tile and number of tiles.
  CALL set_tiles ( grid , ids , ide , jds , jde , ips , ipe , jps , jpe )
  k_start         = kps
  k_end           = kpe

#endif
  ijds = min(ids, jds)
  ijde = max(ide, jde)
#if ( NMM_CORE ==1) 
  allocate(moist_trans(ims:ime,kms:kme,jms:jme,num_3d_m))
  allocate(scalar_trans(ims:ime,kms:kme,jms:jme,num_3d_s))
      DO l=1,num_3d_m
       DO j=jts,jte
         DO k=kts,kte
             DO i=its,ite
                moist_trans(i,k,j,l)=moist(i,j,k,l)
             ENDDO
         ENDDO
       ENDDO
       ENDDO
       DO l=1,num_3d_s
       DO j=jts,jte
         DO k=kts,kte
             DO i=its,ite
               scalar_trans(i,k,j,l)=scalar(i,j,k,l)
             ENDDO
         ENDDO
       ENDDO
       ENDDO
#endif
   chem_minval = epsilc !chem_minval can be case dependant and set below...
   chem_select: SELECT CASE(config_flags%chem_opt)
     CASE (RADM2)
       CALL wrf_debug(15,'calling radm2 from chem_driver')
       haveaer = .false.
     CASE (RADM2_KPP)
       CALL wrf_debug(15,'calling radm2_kpp from chem_driver')
       haveaer = .false.
     CASE (CRIMECH_KPP)
       CALL wrf_debug(15,'calling crimech_kpp from chem_driver')
       haveaer = .false.
     CASE (RADM2SORG)
       CALL wrf_debug(15,'calling radm2sorg aerosols driver from chem_driver')
       haveaer = .true.
     CASE (RADM2SORG_KPP)
       CALL wrf_debug(15,'calling radm2sorg aerosols driver from chem_driver')
       haveaer = .false.
     CASE (RADM2SORG_AQ)
       CALL wrf_debug(15,'calling radm2sorg_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (RACMSORG_AQ)
       CALL wrf_debug(15,'calling racmsorg_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (RADM2SORG_AQCHEM)
       CALL wrf_debug(15,'calling radm2sorg_aqchem aerosols driver from chem_driver')
       haveaer = .true.
     CASE (RACMSORG_AQCHEM_KPP)
       CALL wrf_debug(15,'calling racmsorg_aqchem_kpp aerosols driver from chem_driver')
       haveaer = .true.
     CASE (RACM_ESRLSORG_AQCHEM_KPP)
       CALL wrf_debug(15,'calling racm_esrlsorg_aqchem_kpp aerosols driver from chem_driver')
       haveaer = .true.
     CASE (RACM_KPP)
       CALL wrf_debug(15,'calling racm_kpp from chem_driver')
     CASE (RACMPM_KPP)
       CALL wrf_debug(15,'calling racmpm_kpp from chem_driver')
       haveaer = .false.
     CASE (RACM_MIM_KPP)
       CALL wrf_debug(15,'calling racm_mim_kpp from chem_driver')
       haveaer = .false.
     CASE (RACM_ESRLSORG_KPP)
       CALL wrf_debug(15,'calling racmsorgesrl_kpp aerosols driver from chem_driver')
       haveaer = .false.
     CASE (RACMSORG_KPP)
       CALL wrf_debug(15,'calling racmsorg_kpp aerosols driver from chem_driver')
       haveaer = .false.
     CASE (RACM_SOA_VBS_KPP)
       CALL wrf_debug(15,'calling racm_soa_vbs_kpp aerosols driver from chem_driver')
       haveaer = .false.
     CASE (GOCART_SIMPLE)
       CALL wrf_debug(15,'calling only gocart aerosols driver from chem_driver')
       haveaer = .false.
     CASE (GOCARTRACM_KPP)
       CALL wrf_debug(15,'calling gocart and racm driver from chem_driver')
       haveaer = .false.
     CASE (GOCARTRADM2_KPP)
       CALL wrf_debug(15,'calling gocart and radmkpp driver from chem_driver')
       haveaer = .false.
     CASE (GOCARTRADM2)
       CALL wrf_debug(15,'calling gocart and radm driver from chem_driver')
       haveaer = .false.
     CASE (SAPRC99_KPP)
       CALL wrf_debug(15,'calling saprc99_kpp from chem_driver')
       haveaer = .false.
     CASE (CBMZ_MOSAIC_4BIN_VBS2_KPP)
       CALL wrf_debug(15,'calling cbmz_mosaic_4bin_vbs2_kpp from chem_driver')
       haveaer = .false.
     CASE (SAPRC99_MOSAIC_4BIN_VBS2_KPP)
       CALL wrf_debug(15,'calling saprc99_mosaic_4bin_vbs2_kpp from chem_driver')
       haveaer = .false.
     CASE (MOZART_MOSAIC_4BIN_VBS0_KPP) 
       CALL wrf_debug(15,'calling mozart_mosaic_4bin_vbs0_kpp from chem_driver')
       haveaer = .true.
       !BSINGH -  Added case statement for SAPRC 8 bin
     CASE (SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP)
       CALL wrf_debug(15,'calling saprc99_mosaic_8bin_vbs2_aq_kpp from chem_driver')
       haveaer = .false.
    CASE (SAPRC99_MOSAIC_8BIN_VBS2_KPP) !BSINGH(04/04/2014)
       CALL wrf_debug(15,'calling saprc99_mosaic_8bin_vbs2_kpp from chem_driver')
       haveaer = .false.
       !BSINGH - ENDS

     CASE (CBMZSORG)
       CALL wrf_debug(15,'calling cbmzsorg aerosols from chem_driver')
       haveaer = .true.
     CASE (CBMZSORG_AQ)
       CALL wrf_debug(15,'calling cbmzsorg_aq aerosols from chem_driver')
       haveaer = .true.
     CASE (CBMZ)
       CALL wrf_debug(15,'calling cbmz from chem_driver')
       haveaer = .false.
     CASE (CBMZ_BB)
       CALL wrf_debug(15,'calling cbmz_bb from chem_driver')
       haveaer = .false.
     CASE (CBMZ_BB_KPP)
       CALL wrf_debug(15,'calling cbmz_bb_kpp from chem_driver')
     CASE (CBMZ_SORGAM_KPP)
       CALL wrf_debug(15,'calling cbmz_sorgam_kpp from chem_driver')
       haveaer = .false.
     CASE (CBMZ_MOSAIC_KPP)
       CALL wrf_debug(15,'calling cbmz_mosaic_kpp from chem_driver')
       haveaer = .false.
     CASE (CBMZ_MOSAIC_8BIN_VBS9_KPP)
       CALL wrf_debug(15,'calling cbmz_8bin_vbs9_kpp from chem_driver')
       haveaer = .false.
     CASE (CBMZ_MOSAIC_4BIN_VBS9_KPP)
       CALL wrf_debug(15,'calling cbmz_4bin_vbs9_kpp from chem_driver')
       haveaer = .false.
     CASE (CBMZ_MOSAIC_4BIN)
       CALL wrf_debug(15,'calling cbmz_mosaic_4bin aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_8BIN)
       CALL wrf_debug(15,'calling cbmz_mosaic_8bin aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_4BIN_AQ)
       CALL wrf_debug(15,'calling cbmz_mosaic_4bin_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_8BIN_AQ)
       CALL wrf_debug(15,'calling cbmz_mosaic_8bin_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_DMS_4BIN)
       CALL wrf_debug(15,'calling cbmz_mosaic_dms_4bin aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_DMS_8BIN)
       CALL wrf_debug(15,'calling cbmz_mosaic_dms_8bin aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_DMS_4BIN_AQ)
       CALL wrf_debug(15,'calling cbmz_mosaic_dms_4bin_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_DMS_8BIN_AQ)
       CALL wrf_debug(15,'calling cbmz_mosaic_dms_8bin_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE)
       CALL wrf_debug(15,'calling cbmz_mosaic_dms_4bin_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE)
       CALL wrf_debug(15,'calling cbmz_mosaic_dms_8bin_aq aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CRI_MOSAIC_8BIN_AQ_KPP)
       CALL wrf_debug(15,'calling cri_mosaic_8bin_aq_kpp aerosols driver from chem_driver')
       haveaer = .true.
     CASE (CRI_MOSAIC_4BIN_AQ_KPP)
       CALL wrf_debug(15,'calling cri_mosaic_4bin_aq_kpp aerosols driver from chem_driver')
       haveaer = .true.
     CASE (MOZART_KPP)
       CALL wrf_debug(15,'calling mozart driver from chem_driver')
     CASE (MOZCART_KPP)
       CALL wrf_debug(15,'calling mozcart driver from chem_driver')
    CASE (CHEM_TRACER,CHEM_TRACE2)
       CALL wrf_debug(15,'tracer mode: only doing emissions and dry dep in chem_driver')
    CASE (CHEM_VOLC)
       CALL wrf_debug(15,'Full Volcanic Ash mode: doing emissions (SO2 + ASH), settling, and subgrid transport in chem_driver')
    CASE (CHEM_VOLC_4BIN)
       CALL wrf_debug(15,'4bin Volcanic Ash mode: doing emissions (ASH), settling, and subgrid transport in chem_driver')
    CASE (CHEM_VASH)
       CALL wrf_debug(15,'Volcanic Ash mode: only doing emissions, settling, and subgrid transport in chem_driver')
    CASE (DUST)
       CALL wrf_debug(15,'Dust only mode: only doing emissions, settling, and subgrid transport chem_driver')
    CASE (CO2_TRACER,GHG_TRACER)
      CALL wrf_debug(15,'Greenhouse gas mode: fluxes and transport of GHG')
    CASE DEFAULT
       if(config_flags%tracer_opt > 0 )then
       CALL wrf_debug(15,'only doing tracer transport in chem_driver')
       else
       CALL wrf_debug(15,'calling chem_opt=? from chem_driver')
       endif
   END SELECT chem_select                              
   tracer_select: SELECT CASE(config_flags%tracer_opt)
    CASE (TRACER_SMOKE)
       CALL wrf_debug(15,'tracer mode: 1 tracer for fires')
    CASE (TRACER_TEST1)
       CALL wrf_debug(15,'tracer mode: 8 tracers')
    CASE (TRACER_TEST2)
       CALL wrf_debug(15,'tracer mode: 8 tracers')
    CASE (TRACER_TEST3)
       CALL wrf_debug(15,'tracer mode: 10 tracers')
     CASE DEFAULT
       CALL wrf_debug(15,'calling chem_opt=? from chem_driver')
   END SELECT tracer_select

! initialize cam_mam local arrays
   if ((config_flags%chem_opt == CBMZ_CAM_MAM3_NOAQ) .or. &
       (config_flags%chem_opt == CBMZ_CAM_MAM3_AQ  ) .or. &
       (config_flags%chem_opt == CBMZ_CAM_MAM7_NOAQ) .or. &
       (config_flags%chem_opt == CBMZ_CAM_MAM7_AQ  )) then
      grid%dgnum4d(:,:,:,:) = 0.0 !Balwinder.Singh@pnnl: dgnum is now defined in Registry as dgnum4d
      grid%dgnumwet4d(:,:,:,:) = 0.0 !Balwinder.Singh@pnnl: dgnumwet is now defined in Registry as dgnumwet4d
      wetdens_ap(:,:,:,:) = 0.0
      !BSINGH: 01/31/2013: Sanity check for numgas_mam
      if(numgas_mam < numgas) then
         write(msg,*)'CHEM_DRIVER - NUMGAS_MAM is should be equal to numgas (check chemics_init.F), numgas_mam=',numgas_mam,' and numgas=',numgas
         call wrf_error_fatal( msg )
      endif
      !BSINGH:02/01/2013 - Sanity check for cam_mam_aerols, it should be true (updated in module_cam_mam_init.F)
      if(.NOT.cam_mam_aerosols) then
         write(msg,*)'CHEM_DRIVER - cam_mam_aerosol should be TRUE (check module_physics_init.F), module_cam_mam_aerosol=',cam_mam_aerosols
         call wrf_error_fatal(msg)
      endif
   end if

!
!
!
#if ( NMM_CORE == 1 )
  k_start         = kts
  k_end           = min(kpe,kde-1)
! this should be in seperate routine!!!!!!
      GRID%SIGMA=1
      grid%HYDRO=.FALSE.
      its=max(its,MYIS1)
      jts=max(jts,MYJS2)
      ite=min(ite,MYIE1)
      jte=min(jte,MYJE2)
      DO J=jts,jte
      DO I=its,ite
        pbl_h(i,j)=grid%pblh(i,j)
!
!       PDSL=PD(I,J)*RES(I,J)
!-----------------------------------------------------------------------
!*** LONG AND SHORTWAVE FLUX AT GROUND SURFACE
!-----------------------------------------------------------------------
        IF(grid%CZMEAN(I,J)>0.) THEN
          FACTRS(I,J)=grid%CZEN(I,J)/grid%CZMEAN(I,J)
        ELSE
          FACTRS(I,J)=0.
        ENDIF
        grid%GSW(I,J)=(grid%RSWIN(I,J)-grid%RSWOUT(I,J))*grid%HBM2(I,J)*FACTRS(I,J)
        P8W(I,KTE+1,J)=grid%PT
        grid%XLAT(I,J)=grid%GLAT(I,J)/DEGRAD
        grid%XLONG(I,J)=grid%GLON(I,J)/DEGRAD
        grid%XLAND(I,J)=grid%SM(I,J)+1.
        grid%PSFC(i,j)=grid%PD(I,J)+grid%PDTOP+grid%PT
        grid%UST(I,J)=grid%USTAR(I,J)
        REXNSFC(I,J)=(grid%PSFC(i,j)*1.E-5)**CAPA
        TSFC(I,J)=grid%THS(I,J)*REXNSFC(I,J)
        grid%TSK(I,J)=TSFC(I,J)

        T8W(I,1,J)=TSFC(I,J) 
        P8W(I,KTS,J)=grid%ETA1(KTS)*grid%PDTOP+grid%ETA2(KTS)*grid%PDSL(i,j)+grid%PT
!       
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        z_at_w(i,kts,j)=grid%fis(i,j)/g
        DO K=KTS,KTE+1
           vvel(i,k,j)=grid%w(i,j,k)
           DPL=grid%DETA1(K)*grid%PDTOP+grid%DETA2(K)*grid%PDSL(i,j)
          QL(K)=AMAX1(grid%Q(I,J,K),EPSQ)
          PLYR=grid%AETA1(K)*grid%PDTOP+grid%AETA2(K)*grid%PDSL(i,j)+grid%PT
          TL(K)=grid%T(I,J,K)
!
! here rri is inverse density!
!
          RHO(I,K,J)=PLYR/(R_D*TL(K)*(1.+P608*QL(K)))
          RRI(I,K,J)=1./RHO(i,k,j)
          T_PHY(I,K,J)=TL(K)
          moist_trans(I,K,J,P_QV)=QL(K)/(1.-QL(K))
          P8W(I,K+1,J)=grid%ETA1(K+1)*grid%PDTOP+grid%ETA2(K+1)*grid%PDSL(i,j)+grid%PT
          P_PHY(I,K,J)=PLYR
          DZ8W(I,K,J)=TL(K)*(P608*QL(K)+1.)*R_D                           &
     &                 *(P8W(I,K,J)-P8W(I,K+1,J))                       &
     &                 /(P_PHY(I,K,J)*G)
          if(K.gt.kts)then
            Z_AT_W(i,k,j)=Z_AT_W(I,k-1,j)+DZ8W(I,K-1,J)
            ZMID(I,K-1,J)=.5*(Z_AT_W(I,K-1,J)+Z_AT_W(I,K,J))
          endif

         ENDDO
!         
        DO K=KTS+1,KTE+1
          T8W(I,K,J)=0.5*(TL(K-1)+TL(K))
        ENDDO 
        T8W(I,KTE+2,J)=-1.E20
        ZMID(I,KTE+1,J)=Z_AT_W(I,KTE+1,J)
!         
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  COMPUTE VELOCITY COMPONENTS AT MASS POINTS
!
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(i,j,k,rwmsk,wmsk)
      DO J=MYJS1_P1,MYJE1_P1
!
        DO K=KTS,KTE
          DO I=MYIS_P1,MYIE_P1
            WMSK=grid%VTM(I+grid%IHE(J),J,K)+grid%VTM(I+grid%IHW(J),J,K)                    &
     &          +grid%VTM(I,J+1,K)+grid%VTM(I,J-1,K)
            IF(WMSK>0.)THEN
              RWMSK=1./WMSK
              U_PHY(I,K,J)=(grid%U(I+grid%IHE(J),J,K)*grid%VTM(I+grid%IHE(J),J,K)           &
     &                     +grid%U(I+grid%IHW(J),J,K)*grid%VTM(I+grid%IHW(J),J,K)       &
     &                     +grid%U(I,J+1,K)*grid%VTM(I,J+1,K)                 &
     &                     +grid%U(I,J-1,K)*grid%VTM(I,J-1,K))*RWMSK
              V_PHY(I,K,J)=(grid%V(I+grid%IHE(J),J,K)*grid%VTM(I+grid%IHE(J),J,K)           &
     &                     +grid%V(I+grid%IHW(J),J,K)*grid%VTM(I+grid%IHW(J),J,K)       &
     &                     +grid%V(I,J+1,K)*grid%VTM(I,J+1,K)                 &
     &                     +grid%V(I,J-1,K)*grid%VTM(I,J-1,K))*RWMSK
            ELSE
              U_PHY(I,K,J)=0.
              V_PHY(I,K,J)=0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO     
#endif

      do nv=1,num_chem
         do j=jps,jpe
            do k=kps,kpe
               do i=ips,ipe
                  chem(i,k,j,nv)=max(chem(i,k,j,nv),chem_minval)
               enddo
            enddo
         enddo
      enddo
      select case (config_flags%chem_opt)
      case (RADM2SORG, RADM2SORG_AQ, RADM2SORG_AQCHEM, RADM2SORG_KPP, &
            RACM_ESRLSORG_KPP,RACMSORG_AQ,RACMSORG_KPP, RACMSORG_AQCHEM_KPP, RACM_ESRLSORG_AQCHEM_KPP, RACM_SOA_VBS_KPP)
         do j=jps,jpe
            do k=kps,kpe
               do i=ips,ipe
                  if(chem(i,k,j,p_nu0).lt.1.e07) then
                     chem(i,k,j,p_nu0)=1.e7
                  endif
               enddo
            enddo
         enddo

! Special treatment of CH4 in SAPRC99
      case (CBMZ_MOSAIC_4BIN_VBS2_KPP, SAPRC99_MOSAIC_4BIN_VBS2_KPP,SAPRC99_KPP,SAPRC99_MOSAIC_4BIN_VBS9_KPP, &
           SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_KPP) !BSINGH -  Added for SAPRC 8 bin and non-aq on (04/04/2014)
         CALL wrf_debug ( 15 , ' fixing ch4 conc using co conc' )
         do j=jps,jpe
            do k=kps,kpe
               do i=ips,ipe
                  chem(i,k,j,p_ch4)=1.74
               enddo
            enddo
         enddo
      end select

      vdrog3=0.
#if ( EM_CORE == 1 )
      do j=jps,min(jde-1,jpe)
         do k=kps,kpe
            do i=ips,min(ide-1,ipe)
              vvel(i,k,j)=grid%w_2(i,k,j)
              zmid(i,k,j)=grid%z(i,k,j)
            enddo
         enddo
      enddo
      do j=jps,min(jde-1,jpe)
         do k=kps,min(kde-1,kpe)
            do i=ips,min(ide-1,ipe)
              rri(i,k,j)=grid%alt(i,k,j)
            enddo
         enddo
      enddo
      do j=jps,min(jde-1,jpe)
         do i=ips,min(ide-1,ipe)
            pbl_h(i,j)=grid%pblh(i,j)
         enddo
      enddo

!------------------------------------------------------------------------
! setup ftuv column density timing
!------------------------------------------------------------------------
      if( config_flags%chem_opt == MOZART_KPP .or. &                               
          config_flags%chem_opt == MOZCART_KPP .or. &
          config_flags%chem_opt == MOZART_MOSAIC_4BIN_VBS0_KPP) then
        if( config_flags%phot_opt == FTUV ) then
          CALL ftuv_timestep_init( grid%id, grid%julday )
        endif
      endif

!------------------------------------------------------------------------
! Main chemistry tile loop
!------------------------------------------------------------------------
     
!$OMP PARALLEL DO   &
!$OMP PRIVATE ( ij, its, ite, jts, jte )
   chem_tile_loop_1: DO ij = 1 , grid%num_tiles
       its = grid%i_start(ij) 
       ite = min(grid%i_end(ij),ide-1)
       jts = grid%j_start(ij)
       jte = min(grid%j_end(ij),jde-1) 

       kts=k_start
       kte=min(k_end,kde-1)

#endif
! 
! no time average available in first half hour
!   
!      if( config_flags%chem_conv_tr>0)then
!          call convtrans_prep(grid%gd_cloud,grid%gd_cloud2,grid%gd_cloud_a,&
!              grid%gd_cloud_b,grid%raincv,grid%raincv_a,grid%raincv_b,     &
!              grid%gd_cloud2_a,grid%gd_cloud2_b,convtrans_avglen_m,stepave,&
!              adapt_step_flag,curr_secs,grid%stepave_count,                &
!              ktau,grid%dt,                                                &
!              config_flags%cu_rad_feedback,  config_flags%cu_physics,      &
!              ids,ide, jds,jde, kds,kde,                                   &
!              ims,ime, jms,jme, kms,kme,                                   &
!              its,ite, jts,jte,kts,kte                                     )
!      endif ! chem_conv_tr

!
!
#if ( EM_CORE == 1 )
         CALL wrf_debug ( 15 , ' call chem_prep' )
         CALL chem_prep ( config_flags,                                               &
                         grid%u_2, grid%v_2, grid%p, grid%pb,             &
                         grid%alt,grid%ph_2, grid%phb, grid%t_2,          &
                         moist, num_3d_m, rho,                                        &
                         p_phy,  u_phy, v_phy,                                        &
                         p8w, t_phy, t8w, grid%z, z_at_w,                          &
                         dz8w, rh, grid%fnm, grid%fnp,                              &
                         ids, ide, jds, jde, kds, kde,                                &
                         ims, ime, jms, jme, kms, kme,                                &
                         its,ite,jts,jte,                                             &
                         k_start, k_end                                               )
#endif

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
       if( (its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and.                          &
           (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and.                          &
           (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)  ) then
          call wrf_debug(15,"calling chem_dbg at top of chem_driver")
          call chem_dbg(CHEM_DBG_I,CHEM_DBG_J,CHEM_DBG_K,grid%dt,ktau,                &
               dz8w,t_phy,p_phy,rho,chem,emis_ant,                       &
               ids,ide, jds,jde, kds,kde,                                             &
               ims,ime, jms,jme, kms,kme,                                             &
               its,ite, jts,jte,kts,kte,                                              &
               config_flags%kemit                                                     )
       end if
#endif

!-- set upper boundary condition
!   if( config_flags%have_bcs_upper )then
!      call wrf_debug(15,'set upper boundary condition')
!      call tropopause_driver( grid%id, grid%dt, current_date_char,         &
!                              t_phy, p_phy, p8w, zmid, z_at_w,             &
!                              grid%tropo_lev, grid%tropo_p,  grid%tropo_z, &
!                              ids, ide, jds, jde, kds, kde,                &
!                              ims, ime, jms, jme, kms, kme,                &
!                              its, ite, jts, jte, kts, kte                 )
!      call upper_bc_driver  ( grid%id, grid%dt, current_date_char, &
!                              chem, p_phy, p8w, grid%tropo_lev,    &
!                              ids,ide, jds,jde, kds,kde,           &
!                              ims,ime, jms,jme, kms,kme,           &
!                              its,ite, jts,jte, kts,kte            )
!   endif


    !czhao-initialize for budget analysis
    if(ktau==1)then
     grid%tmp_budget(its:ite, kts:kte, jts:jte,:)=0.
     grid%tend_tran(its:ite, kts:kte, jts:jte,:)=0.
     grid%tend_dryd(its:ite, kts:kte, jts:jte,:)=0.
     grid%tend_acti(its:ite, kts:kte, jts:jte,:)=0.
     grid%tend_wetdcv(its:ite, kts:kte, jts:jte,:)=0.
     grid%tend_wetdls(its:ite, kts:kte, jts:jte,:)=0.
    end if  !ktau

  snowchem_old(its:ite,jts:jte,:,:)=0.0
  call budget_calc(1,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
       grid%tmp_budget,grid%tend_tran,  &
       tend_rate_dep,snowchem_old,  &
       ims,ime, jms,jme, kms,kme,                              &
       its,ite, jts,jte, kts,kte                                )

!--- emissions

! call budget_calc(0,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
!      grid%tmp_budget,grid%tend_emis, &
!      tend_rate_dep,snowchem_old,  &
!      ims,ime, jms,jme, kms,kme,                              &
!      its,ite, jts,jte, kts,kte                                )

    if(config_flags%emiss_inpt_opt > 0 .or. config_flags%dust_opt > 0)then
      call wrf_debug(15,'calling emissions driver')
 
      !czhao
      call cpu_time(startt)

      call emissions_driver(grid%id,ktau,grid%dt,grid%DX,                                  &
              adapt_step_flag, curr_secs,                                                  &
              grid%plumerisefire_frq,grid%stepfirepl,                                      &
              grid%bioemdt,grid%stepbioe,                                                  &
              config_flags,                                                                &
#if (NMM_CORE == 1)
              grid%gmt,ijulian,rri,t_phy,moist_trans,p8w,t8w,u_phy,v_phy,vvel,         &
#endif
#if (EM_CORE == 1 )
              grid%gmt,ijulian,rri,t_phy,moist,p8w,t8w,u_phy,v_phy,vvel,               &
#endif
              grid%e_bio,p_phy,chem,rho,dz8w,grid%ne_area,emis_ant,emis_vol,grid%tsk,      &
              grid%erod,g,emis_seas,emis_dust,tracer,grid%sst,  &
              ebu , ebu_in,grid%mean_fct_agtf,grid%mean_fct_agef,grid%mean_fct_agsv,       &
              grid%mean_fct_aggr,grid%firesize_agtf, &
              grid%firesize_agef,grid%firesize_agsv,grid%firesize_aggr,                    &
              grid%u10,grid%v10,grid%ivgtyp,grid%isltyp,grid%gsw,grid%vegfra,grid%rmol,    &
              grid%ust,grid%znt,grid%dms_0,grid%erup_beg,grid%erup_end,                    &
              grid%xland,grid%xlat,grid%xlong,                                             &
              z_at_w,zmid,grid%smois,dustin,seasin,                                        &
              !additional for DEAD dust emission scheme -czhao
              grid%landusef,grid%soilctop,grid%num_land_cat,grid%num_soil_cat, &
              grid%sh2o,grid%snowc,grid%lai,grid%watsat2d,    & 
              !additional for DEAD dust emission scheme -czhao
              !additional for subgrid winds -czhao
              grid%hfx,grid%tke_pbl,grid%pblh,grid%psfc,grid%rainncv,grid%raincv,            &
              grid%sigma_dry,grid%sigma_wet,grid%sigma_tke,grid%k_shape,grid%c_scale,           &
              !additional for subgrid winds -czhao
              grid%sebio_iso,grid%sebio_oli,grid%sebio_api,grid%sebio_lim,                 &
              grid%sebio_xyl,grid%sebio_hc3,grid%sebio_ete,grid%sebio_olt,                 &
              grid%sebio_ket,grid%sebio_ald,grid%sebio_hcho,grid%sebio_eth,                &
              grid%sebio_ora2,grid%sebio_co,grid%sebio_nr,                                 &
              grid%sebio_sesq,grid%sebio_mbo,                                              &
              grid%noag_grow,grid%noag_nongrow,grid%nononag,grid%slai,                     &
              grid%ebio_iso,grid%ebio_oli,grid%ebio_api,grid%ebio_lim,grid%ebio_xyl,       &
              grid%ebio_hc3,grid%ebio_ete,grid%ebio_olt,grid%ebio_ket,grid%ebio_ald,       &
              grid%ebio_hcho,grid%ebio_eth,grid%ebio_ora2,grid%ebio_co,grid%ebio_nr,       &
              grid%ebio_no,grid%ebio_sesq,grid%ebio_mbo,                                   &
              grid%ebio_c10h16,grid%ebio_tol,grid%ebio_bigalk,                             &
              grid%ebio_ch3oh,grid%ebio_acet,grid%ebio_nh3,grid%ebio_no2,                  &
              grid%ebio_c2h5oh,grid%ebio_ch3cooh,grid%ebio_mek,grid%ebio_bigene,           &
              grid%ebio_c2h6,grid%ebio_c2h4,grid%ebio_c3h6,grid%ebio_c3h8,grid%ebio_so2,   &
              grid%ebio_dms,                                                               &
              grid%ebio_alk3, grid%ebio_alk4, grid%ebio_alk5, grid%ebio_ole1, grid%ebio_ole2,            &    
              grid%ebio_aro1, grid%ebio_aro2, grid%ebio_ccho, grid%ebio_meoh,             &    
              grid%ebio_ethene, grid%ebio_hcooh, grid%ebio_terp, grid%ebio_bald,          &    
              grid%ebio_cco_oh, grid%ebio_rco_oh,                                         &    
              grid%clayfrac,grid%sandfrac,grid%dust_alpha,grid%dust_gamma,grid%dust_smtune,&
              grid%snowh,grid%zs,                                                          &
              grid%ebio_c5h8,grid%ebio_apinene,grid%ebio_bpinene,grid%ebio_toluene,        &
              grid%ebio_ch3cho,grid%ebio_ch3co2h,grid%ebio_tbut2ene,                      &
              grid%ebio_c2h5cho,grid%ebio_nc4h10,					                      &
#if (NMM_CORE == 1)
              grid%T2,grid%RSWIN,                                                         &
#endif
#if (EM_CORE == 1 )
              grid%T2,grid%swdown,                                                         &
#endif
              grid%nmegan,grid%EFmegan,                                                    &
              grid%msebio_isop,                                                            &
              grid%mlai,                                                                   &
              grid%pftp_bt, grid%pftp_nt, grid%pftp_sb, grid%pftp_hb,                      &
              grid%mtsa,                                                                   &
              grid%mswdown,                                                                &
              grid%mebio_isop,grid%mebio_apin,grid%mebio_bpin, grid%mebio_bcar,            &
              grid%mebio_acet,grid%mebio_mbo,grid%mebio_no,                                &
              current_month,                                                               &
         ! stuff for LNOx emissions
             grid%ht, grid%refl_10cm, grid%ic_flashrate, grid%cg_flashrate,                 &
         ! stuff for aircraft emissions
              emis_aircraft,grid%e_dust,grid%e_seasalt,grid%e_dms,                         &
         ! stuff for the GHG fluxes
              vprm_in,grid%rad_vprm,grid%lambda_vprm,                                      &
              grid%alpha_vprm,grid%resp_vprm,grid%xtime,                                   &
              grid%TSLB, wet_in,grid%RAINC,grid%RAINNC,                                    &
              grid%potevp,grid%SFCEVP,grid%LU_INDEX,                                       &
              grid%biomt_par,grid%emit_par,grid%ebio_co2oce,                               &
              eghg_bio,                                                                    &

              ids,ide, jds,jde, kds,kde,                                                   &
              ims,ime, jms,jme, kms,kme,                                                   &
              its,ite,jts,jte,kts,kte)

            call cpu_time(finisht)
            print '("Chun Zhao check: emissions_driver, Time = ",f6.3," seconds.")',finisht-startt


              if( config_flags%chem_opt == MOZART_KPP .or. &
                  config_flags%chem_opt == MOZCART_KPP .or. &
                  config_flags%chem_opt == MOZART_MOSAIC_4BIN_VBS0_KPP ) then
                 call mozcart_lbc_set( chem, num_chem, grid%id, &
                                       ims, ime, jms, jme, kms, kme,    &                  
                                       its, ite, jts, jte, kts )                           
              end if

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
       if( (its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and.                          &
           (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and.                          &
           (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)  ) then
          call wrf_debug(15,'calling chem_dbg after emissions_driver')
          call chem_dbg(CHEM_DBG_I,CHEM_DBG_J,CHEM_DBG_K,grid%dt,ktau,                &
               dz8w,t_phy,p_phy,rho,chem,emis_ant,                                    &
               ids,ide, jds,jde, kds,kde,                                             &
               ims,ime, jms,jme, kms,kme,                                             &
               its,ite, jts,jte, kts, kte,                                            &
               config_flags%kemit                                                     )
       end if
#endif 
     endif

! call budget_calc(1,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
!      grid%tmp_budget,grid%tend_emis,  &
!      tend_rate_dep,snowchem_old,  &
!      ims,ime, jms,jme, kms,kme,                              &
!      its,ite, jts,jte, kts,kte                                )

!
! calculate aerosol optical properities
!
      !czhao initialize totoa
      if (ktau==1.or.config_flags%restart) then 
      grid%totoa_a01(its:ite,kts:kte,jts:jte)=0.0
      grid%totoa_a02(its:ite,kts:kte,jts:jte)=0.0
      grid%totoa_a03(its:ite,kts:kte,jts:jte)=0.0
      grid%totoa_a04(its:ite,kts:kte,jts:jte)=0.0
      endif

      if( do_photstep .and. &
           config_flags%chem_opt /= CHEM_TRACER .and. &
           config_flags%chem_opt /= CHEM_VASH .and. &
           config_flags%chem_opt /= CHEM_VOLC .and. &
           config_flags%chem_opt /= CHEM_VOLC_4BIN .and. &
           config_flags%chem_opt /= DUST .and. &
           config_flags%chem_opt /= CHEM_TRACE2 .and. &
           config_flags%chem_opt /= CO2_TRACER  .and. &
           config_flags%chem_opt /= GHG_TRACER ) then
         call wrf_debug(15,'calling optical driver')
         call optical_driver (grid%id,curr_secs,grid%dt,config_flags,haveaer,         &
              chem,dz8w,rri,rh,                                                       &
              grid%h2oai,grid%h2oaj,grid%h2oai_adjust,grid%h2oaj_adjust,grid%h2oac_adjust, &
              !czhao
              tauaersw,absaersw,extaersw,gaersw,waersw,     &
              grid%bscoef1,grid%bscoef2,grid%bscoef3,grid%bscoef4,                    &
              grid%l2aer,grid%l3aer,grid%l4aer,grid%l5aer,grid%l6aer,grid%l7aer,      &
              grid%totoa_a01,grid%totoa_a02,grid%totoa_a03,grid%totoa_a04,            &
              grid%totoa_a05,grid%totoa_a06,grid%totoa_a07,grid%totoa_a08,            &
              tauaerlw,    & 
              tauaersw_dust,absaersw_dust,gaersw_dust,waersw_dust,                    &
              tauaersw_bc,absaersw_bc,gaersw_bc,waersw_bc,                    &
              tauaersw_sulf,absaersw_sulf,gaersw_sulf,waersw_sulf,                    &
              tauaersw_oc,absaersw_oc,gaersw_oc,waersw_oc,                    &
              tauaersw_other,absaersw_other,gaersw_other,waersw_other,                    &
              tauaerlw_dust,tauaerlw_bc,tauaerlw_sulf,tauaerlw_oc,tauaerlw_other,  &
              !czhao
              ids,ide, jds,jde, kds,kde,                                              &
              ims,ime, jms,jme, kms,kme,                                              &
              its,ite, jts,jte, kts,kte)
      endif

!
! calculate photolysis rates
!
      if( do_photstep .and. &
           config_flags%chem_opt /= CHEM_TRACER .and. &
           config_flags%chem_opt /= CHEM_VASH .and. &
           config_flags%chem_opt /= CHEM_VOLC .and. &
           config_flags%chem_opt /= CHEM_VOLC_4BIN .and. &
           config_flags%chem_opt /= DUST .and. &
           config_flags%chem_opt /= CHEM_TRACE2 .and. &
           config_flags%chem_opt /= CO2_TRACER  .and. &
           config_flags%chem_opt /= GHG_TRACER  ) then
         call wrf_debug(15,'calling photolysis driver')
         call photolysis_driver (grid%id,curr_secs,ktau,grid%dt,                      &
              config_flags,haveaer,                                                   &
#if (NMM_CORE == 1)
              grid%gmt,ijulian,t_phy,moist_trans,grid%aerwrf,p8w,t8w,p_phy,       &
#endif
#if (EM_CORE == 1)
              grid%gmt,ijulian,t_phy,moist,grid%aerwrf,p8w,t8w,p_phy,             &
#endif
              chem,rho,dz8w,grid%xlat,grid%xlong,                                     &
              z_at_w,                                                                 &
              grid%gd_cloud_b,grid%gd_cloud2_b,                                       &
              grid%ph_macr,grid%ph_o31d,grid%ph_o33p,grid%ph_no2,grid%ph_no3o2,       &
              grid%ph_no3o,grid%ph_hno2,grid%ph_hno3,grid%ph_hno4,grid%ph_h2o2,       &
              grid%ph_ch2or,grid%ph_ch2om,grid%ph_ch3cho,grid%ph_ch3coch3,            &
              grid%ph_ch3coc2h5,grid%ph_hcocho,grid%ph_ch3cocho,                      &
              grid%ph_hcochest,grid%ph_ch3o2h,grid%ph_ch3coo2h,grid%ph_ch3ono2,       &
              grid%ph_hcochob,grid%ph_n2o5,grid%ph_o2,grid%ph_n2o,                    &
              grid%ph_pan,grid%ph_mpan,grid%ph_acetol,grid%ph_gly,                    &
              grid%ph_open,grid%ph_mek,grid%ph_etooh,grid%ph_prooh,grid%ph_pooh,      &
              grid%ph_acetp,grid%ph_xooh,grid%ph_isooh,grid%ph_alkooh,                &
              grid%ph_mekooh,grid%ph_tolooh,grid%ph_terpooh,grid%ph_mvk,              &
              grid%ph_glyald,grid%ph_hyac,                                            &
              config_flags%track_tuv_lev,                                             &
              config_flags%track_rad_num,                                             &
              config_flags%track_tuv_num,                                             &
              grid%radfld,grid%adjcoe,grid%phrate,                                    &
              grid%track_wc,grid%track_zref,                                          &
              tauaersw(:,:,:,2),tauaersw(:,:,:,3),tauaersw(:,:,:,4),tauaersw(:,:,:,5),&
              gaersw(:,:,:,2),gaersw(:,:,:,3),gaersw(:,:,:,4),gaersw(:,:,:,5),        &
              waersw(:,:,:,2),waersw(:,:,:,3),waersw(:,:,:,4),waersw(:,:,:,5),        &
              grid%bscoef1,grid%bscoef2,grid%bscoef3,grid%bscoef4,                    &
              grid%l2aer,grid%l3aer,grid%l4aer,grid%l5aer,grid%l6aer,grid%l7aer,      &
              grid%pm2_5_dry,grid%pm2_5_water,grid%uvrad,grid%ivgtyp,                 &
              ids,ide, jds,jde, kds,kde,                                              &
              ims,ime, jms,jme, kms,kme,                                              &
              its,ite,jts,jte,kts,kte)

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
         if( (its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
             (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
             (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)  ) then
          call wrf_debug(15,'calling chem_dbg after photolysis_driver')
          call chem_dbg(CHEM_DBG_I,CHEM_DBG_J,CHEM_DBG_K,grid%dt,ktau,                &
               dz8w,t_phy,p_phy,rho,chem,emis_ant,                                    &
               ids,ide, jds,jde, kds,kde,                                             &
               ims,ime, jms,jme, kms,kme,                                             &
               its,ite, jts,jte, kts,kte,                                             &
               config_flags%kemit,                                                    &
               grid%ph_macr,grid%ph_o31d,grid%ph_o33p,grid%ph_no2,grid%ph_no3o2,      &
               grid%ph_no3o,grid%ph_hno2,grid%ph_hno3,grid%ph_hno4,grid%ph_h2o2,      &
               grid%ph_ch2or,grid%ph_ch2om,grid%ph_ch3cho,grid%ph_ch3coch3,           &
               grid%ph_ch3coc2h5,grid%ph_hcocho,grid%ph_ch3cocho,                     &
               grid%ph_hcochest,grid%ph_ch3o2h,grid%ph_ch3coo2h,                      &
               grid%ph_ch3ono2,grid%ph_hcochob,grid%ph_n2o5,grid%ph_o2                &
                                                                                      )
       end if
#endif
      endif


!
! do vertical mixing with dry deposition
!
! save old concentrations for vertical mixing tendencies
#if (EM_CORE == 1)
   DO nv=PARAM_FIRST_SCALAR,num_chem_ct
     chem_old(its:ite,kts:kte,jts:jte,nv) = chem(its:ite,kts:kte,jts:jte,chem_ct_indices(nv))
   ENDDO
#endif

!budget for dry deposition
call budget_calc2(0,config_flags%chem_opt,chem,rho,  &
     grid%tmp_budget,grid%tend_dryd,grid%tend_acti,temporary,temporary_cw, &
     deposition_dryrate,dz8w,grid%dt,  &
     grid%drydep_flux, &
     ims,ime, jms,jme, kms,kme,                              &
     its,ite, jts,jte, kts,kte                                )

      if (config_flags%vertmix_onoff>0) then
         if (ktau.gt.2) then
            call wrf_debug(15,'calling dry_deposition_driver')
            call dry_dep_driver(grid%id,curr_secs,ktau,grid%dt,config_flags,          &
#if (NMM_CORE == 1)
                 grid%gmt,ijulian,t_phy,moist_trans,scalar_trans,p8w,t8w,vvel,    &
#endif
#if (EM_CORE == 1)
                 grid%gmt,ijulian,t_phy,moist,scalar,p8w,t8w,vvel,                &
#endif
!czhao modify to use exch_m for diffusion
!cz                rri,p_phy,chem,tracer,rho,dz8w,rh,grid%exch_h,grid%hfx,grid%dx,      & 
                 rri,p_phy,chem,tracer,rho,dz8w,rh,grid%exch_m,grid%hfx,grid%dx,      & 
                 grid%cldfra, grid%cldfra_old,grid%raincv_b,seasin,dustin,            &
	         grid%ccn1, grid%ccn2, grid%ccn3, grid%ccn4, grid%ccn5, grid%ccn6,    &
                 grid%qndropsource,grid%ivgtyp,grid%tsk,grid%gsw,grid%vegfra,pbl_h,   &
                 grid%rmol,grid%ust,grid%znt,grid%xlat,grid%xlong,                    &
                 zmid,z_at_w,grid%xland,grid%ash_fall,                                &
                 grid%h2oaj,grid%h2oai,grid%nu3,grid%ac3,grid%cor3,grid%asulf,        &
                 grid%ahno3,grid%anh3,grid%cvaro1,grid%cvaro2,grid%cvalk1,grid%cvole1,&
                 grid%cvapi1,grid%cvapi2,grid%cvlim1,grid%cvlim2,grid%dep_vel_o3,     &
                 !czhao 
                 grid%dep_vel_nnuc,grid%dep_vel_mnuc,grid%dep_vel_nacc, &
                 grid%dep_vel_macc,grid%dep_vel_ncor,grid%dep_vel_mcor,  &
                 grid%dep_vel_a01,grid%dep_vel_a02,grid%dep_vel_a03,grid%dep_vel_a04, & 
                 grid%dep_vel_a05,grid%dep_vel_a06,grid%dep_vel_a07,grid%dep_vel_a08, &
                 emis_ant,ebu_in,                     &
                 config_flags%sf_urban_physics,numgas,current_month,dvel,grid%snowh,grid%is_CAMMGMP_used,grid%f_ice_phy, grid%f_rain_phy,  & !BSINGH:01/31/2013: Added is_CAMMGMP_used for MAM drydep
                 !czhao
                 grid%nihf3d,grid%niimm3d,grid%nidep3d,grid%nimey3d,        &
                 grid%rndst1,grid%rndst2,grid%rndst3,grid%rndst4,                     & 
                 grid%nacon1,grid%nacon2,grid%nacon3,grid%nacon4,                     & 
                 grid%dgnuc,grid%dgacc,grid%dgcor,                                    &
                 grid%tke_pbl,grid%relhum3d,grid%wsub3d,    &
                 grid%dep_vel,grid%num_vert_mix,                                                   &
                 ids,ide, jds,jde, kds,kde,                                           &
                 ims,ime, jms,jme, kms,kme,                                           &
                 its,ite,jts,jte,kts,kte)

            if( config_flags%chem_opt == MOZART_KPP .or. &
                config_flags%chem_opt == MOZCART_KPP .or. &
                config_flags%chem_opt == MOZART_MOSAIC_4BIN_VBS0_KPP ) then
               call mozcart_lbc_set( chem, num_chem, grid%id, &
                                     ims, ime, jms, jme, kms, kme,    & 
                                     its, ite, jts, jte, kts )
            end if

         end if


! accumulate vertical mixing tendencies
#if (EM_CORE == 1)
   DO nv=PARAM_FIRST_SCALAR,num_chem_ct
      grid%vmix_ct(its:ite,kts:kte,jts:jte,nv) = grid%vmix_ct(its:ite,kts:kte,jts:jte,nv) + &
                                                        (chem(its:ite,kts:kte,jts:jte,chem_ct_indices(nv)) - &
                                                     chem_old(its:ite,kts:kte,jts:jte,nv))
   ENDDO
#endif

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
       if( (its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
           (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
           (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)  ) then
          call wrf_debug(15,'calling chem_dbg after dry_deposition_driver')
          call chem_dbg(CHEM_DBG_I,CHEM_DBG_J,CHEM_DBG_K,grid%dt,ktau,                &
               dz8w,t_phy,p_phy,rho,chem,emis_ant,                                    &
               ids,ide, jds,jde, kds,kde,                                             &
               ims,ime, jms,jme, kms,kme,                                             &
               its,ite, jts,jte, kts, kte,                                            &
               config_flags%kemit,                                                    &
               grid%ph_macr,grid%ph_o31d,grid%ph_o33p,grid%ph_no2,grid%ph_no3o2,      &
               grid%ph_no3o,grid%ph_hno2,grid%ph_hno3,grid%ph_hno4,grid%ph_h2o2,      &
               grid%ph_ch2or,grid%ph_ch2om,grid%ph_ch3cho,grid%ph_ch3coch3,           &
               grid%ph_ch3coc2h5,grid%ph_hcocho,grid%ph_ch3cocho,                     &
               grid%ph_hcochest,grid%ph_ch3o2h,grid%ph_ch3coo2h,                      &
               grid%ph_ch3ono2,grid%ph_hcochob,grid%ph_n2o5,grid%ph_o2                &
                                                                                      )
       end if
#endif
      end if


!budget for dry deposition
!assign dry deposition rate for budget 
aer_budget_drydep1: SELECT CASE(config_flags%chem_opt)
CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_KPP, &
      RACMSORG_AQ,RACMSORG_KPP,RACM_ESRLSORG_KPP,CBMZSORG,CBMZSORG_AQ)
deposition_dryrate(its:ite,jts:jte,1)=grid%dep_vel_mnuc(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,2)=grid%dep_vel_macc(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,3)=grid%dep_vel_mcor(its:ite,jts:jte)
do nv=1,3
do j=jts,jte
do i=its,ite
if (deposition_dryrate(i,j,nv).lt.0) deposition_dryrate(i,j,nv)=0.0
enddo
enddo
enddo
CASE (CBMZ_MOSAIC_8BIN,  &
      CBMZ_MOSAIC_8BIN_AQ,CBMZ_MOSAIC_DMS_8BIN, &
      CBMZ_MOSAIC_DMS_8BIN_AQ,CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE)
deposition_dryrate(its:ite,jts:jte,1)=grid%dep_vel_a01(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,2)=grid%dep_vel_a02(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,3)=grid%dep_vel_a03(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,4)=grid%dep_vel_a04(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,5)=grid%dep_vel_a05(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,6)=grid%dep_vel_a06(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,7)=grid%dep_vel_a07(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,8)=grid%dep_vel_a08(its:ite,jts:jte)
do nv=1,8
do j=jts,jte
do i=its,ite
if (deposition_dryrate(i,j,nv).lt.0) deposition_dryrate(i,j,nv)=0.0
enddo
enddo
enddo
CASE (CBMZ_MOSAIC_4BIN,CBMZ_MOSAIC_4BIN_AQ,  &
      CBMZ_MOSAIC_DMS_4BIN,CBMZ_MOSAIC_DMS_4BIN_AQ,CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE)
deposition_dryrate(its:ite,jts:jte,1)=grid%dep_vel_a01(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,2)=grid%dep_vel_a02(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,3)=grid%dep_vel_a03(its:ite,jts:jte)
deposition_dryrate(its:ite,jts:jte,4)=grid%dep_vel_a04(its:ite,jts:jte)
do nv=1,4
do j=jts,jte
do i=its,ite
if (deposition_dryrate(i,j,nv).lt.0) deposition_dryrate(i,j,nv)=0.0
enddo
enddo
enddo
CASE DEFAULT
  call wrf_debug(15,'aerosol budget before drydep NOT CALLING')
END SELECT aer_budget_drydep1

grid%drydep_flux(its:ite,jts:jte,:,:)=0.0
call budget_calc2(1,config_flags%chem_opt,chem,rho,  &
     grid%tmp_budget,grid%tend_dryd,grid%tend_acti,temporary,temporary_cw, &
     deposition_dryrate,dz8w,grid%dt,    &
     grid%drydep_flux, &
     ims,ime, jms,jme, kms,kme,                              &
     its,ite, jts,jte, kts,kte                                )

!
!   convective transport/wet deposition
!

call budget_calc(0,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
     grid%tmp_budget,grid%tend_wetdcv,  &
     tend_rate_dep,snowchem_old,  &
     ims,ime, jms,jme, kms,kme,                              &
     its,ite, jts,jte, kts,kte                                )

     if( config_flags%cu_physics>0 .and. config_flags%chem_conv_tr>0 &
          .and. config_flags%cu_physics/=kfcupscheme ) then            !BSINGH - For WRFCuP scheme
        call wrf_debug(15,'calling conv transport for chemical species')
        if(config_flags%chem_opt >0 )then

           if (config_flags%chem_conv_tr==1) then ! original Grell convection, no wet scavenging
             if (config_flags%cu_physics.ne.3.and.config_flags%cu_physics.ne.5) then 
             call wrf_error_fatal("Convective transport is not consistent with Convective scheme used in physics!!!")
             endif

        ! save old concentrations for convective tendencies
#if (EM_CORE == 1)
        DO nv=PARAM_FIRST_SCALAR,num_chem_ct
           chem_old(its:ite,kts:kte,jts:jte,nv) = chem(its:ite,kts:kte,jts:jte,chem_ct_indices(nv))
        ENDDO
#endif

           call grelldrvct(grid%DT,ktau,grid%DX,rho,grid%RAINCV_B,chem,  &
#if (NMM_CORE == 1)
           U_phy,V_phy,t_phy,moist_trans,dz8w,         &
#endif
#if (EM_CORE == 1)
           U_phy,V_phy,t_phy,moist,dz8w,               &
#endif
             p_phy,XLV,CP,G,r_v,                                                      &
             z_at_w,grid%cu_co_ten,                                                   &
             grid%wd_no3_cu,grid%wd_so4_cu,                                 &
             grid%k22_shallow,grid%kbcon_shallow,grid%ktop_shallow,grid%xmb_shallow,  &
             config_flags%ishallow,num_moist,numgas,num_chem,config_flags%chem_opt,0, &
             config_flags%conv_tr_wetscav,config_flags%conv_tr_aqchem,                &
             ids,ide, jds,jde, kds,kde,                                               &
             ims,ime, jms,jme, kms,kme,                                               &
             its,ite,jts,jte,kts,k_end)
           endif 
           if (config_flags%chem_conv_tr==2) then ! Grell convection, include wet scavenging
             if (config_flags%cu_physics.ne.3) then 
             call wrf_error_fatal("Convective transport is not consistent with Convective scheme used in physics!!!")
             endif
           !czhao
!          call grell_cup_transport
           endif
           if (config_flags%chem_conv_tr==3) then ! GD3 convection, include wet scavenging
             if (config_flags%cu_physics.ne.5) then 
             call wrf_error_fatal("Convective transport is not consistent with Convective scheme used in physics!!!")
             endif
           !czhao
!          call grell3_cup_transport
           endif
           if (config_flags%chem_conv_tr==4) then ! KF-eta convection, include wet scavenging
             if (config_flags%cu_physics.ne.1) then 
             call wrf_error_fatal("Convective transport is not consistent with Convective scheme used in physics!!!")
             endif
           !czhao
           call kfeta_cup_transport(grid%DT,ktau,grid%DX,rho,grid%alt,chem,  &
#if (NMM_CORE == 1)
           U_phy,V_phy,t_phy,moist_trans,dz8w,         &
#endif
#if (EM_CORE == 1)
           U_phy,V_phy,t_phy,moist,dz8w,               &
#endif
           p_phy,z_at_w,  &
           grid%UMF_KFETA,grid%UER_KFETA,     &
           grid%UDR_KFETA,grid%DMF_KFETA,     &
           grid%DER_KFETA,grid%DDR_KFETA,     & 
           grid%WU_KFETA,                     & 
           grid%LDB_KFETA,grid%LET_KFETA,     & 
           grid%KSL_KFETA,grid%LFS_KFETA,     &
           grid%KLC_KFETA,grid%CUTOP_KFETA,   &
           grid%CUBOT_KFETA,grid%NCA,         &
           grid%ISHALL_KFETA,grid%CHECK_KFETA,  &
           num_moist,numgas,num_chem,         &
           config_flags%chem_opt,config_flags%chem_conv_scav, &
           ids,ide, jds,jde, kds,kde,                           &
           ims,ime, jms,jme, kms,kme,                           &
           its,ite,jts,jte,kts,kte)
           endif
        
             if( config_flags%chem_opt == MOZART_KPP .or. &
                 config_flags%chem_opt == MOZCART_KPP .or. &
                 config_flags%chem_opt == MOZART_MOSAIC_4BIN_VBS0_KPP ) then
                     call mozcart_lbc_set( chem, num_chem, grid%id, &
                                 ims, ime, jms, jme, kms, kme,    &
                                 its, ite, jts, jte, kts )
             end if
        ! accumulate convective tendencies
#if (EM_CORE == 1)
        DO nv=PARAM_FIRST_SCALAR,num_chem_ct
          grid%conv_ct(its:ite,kts:kte,jts:jte,nv) = grid%conv_ct(its:ite,kts:kte,jts:jte,nv) + &
                                                            (chem(its:ite,kts:kte,jts:jte,chem_ct_indices(nv)) - &
                                                         chem_old(its:ite,kts:kte,jts:jte,nv))
        ENDDO
#endif
        endif !chem_opt
        if (config_flags%tracer_opt > 0)then
        call wrf_debug(15,'calling conv transport for tracers')
        call grelldrvct(grid%DT,ktau,grid%DX,                    &
             rho,grid%RAINCV_B,tracer,                               &
#if (NMM_CORE == 1)
             U_phy,V_phy,t_phy,moist_trans,dz8w,                                      &
#endif         
#if (EM_CORE == 1)
             U_phy,V_phy,t_phy,moist,dz8w,                                            &
#endif
             p_phy,XLV,CP,G,r_v,                                                      &
             z_at_w, grid%cu_co_ten,                                                  &
             grid%wd_no3_cu,grid%wd_so4_cu,                                 &
             grid%k22_shallow,grid%kbcon_shallow,grid%ktop_shallow,grid%xmb_shallow,  &
             config_flags%ishallow,num_moist,0,num_tracer,0,config_flags%tracer_opt,  &
             config_flags%conv_tr_wetscav,config_flags%conv_tr_aqchem,                &
             ids,ide, jds,jde, kds,kde,                                               &
             ims,ime, jms,jme, kms,kme,                                               &
             its,ite,jts,jte,kts,k_end)


          end if  ! tracer_opt
     end if !chem_conv_tr

!
!
!
!
     n2o5_het(its:ite,kts:kte,jts:jte)=0.
! Calculate rate of n2o5 hydrolysis 
     call wrf_debug(15,'calling calc_het_n2o5')

grid%wetdep_flux(its:ite,jts:jte,1:4,:)=0.0
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcphi)=0.0 !bcphi
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocphi)=0.0 !ocphi
  call budget_calc(1,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
       grid%tmp_budget,grid%tend_wetdcv, &
       tend_rate_dep,snowchem_old,  &
       ims,ime, jms,jme, kms,kme,                              &
       its,ite, jts,jte, kts,kte                                )
!add convective scavenging
do nv=1,5 
do n=1,4 
do j=jts,jte
do i=its,ite
if (tend_rate_dep(i,j,n,nv).gt.0) tend_rate_dep(i,j,n,nv)=0.0
enddo
enddo
enddo
enddo
!dust
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_dust)= &
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_dust)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,1)
!bcphi-wet from bccw-wet
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_bcphi)= &
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_bcphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,3)
!bcphi-dry from bc-wet
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcphi)= &
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,2)
!ocphi-wet from occw-wet
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_ocphi)= &
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_ocphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,5)
!ocphi-dry from oc-wet
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocphi)= &
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,4)

!
! For the chemistry tracer mode, only emissions and vertical mixing are done.
! So, finish any remaining tiles and then skip to the end of chem_driver.
!

!Chemistry time step check!
if( do_chemstep .and.                           &
          config_flags%chem_opt /= CHEM_TRACER .and. &
          config_flags%chem_opt /= CHEM_VASH .and. &
          config_flags%chem_opt /= CHEM_VOLC .and. &
          config_flags%chem_opt /= CHEM_VOLC_4BIN .and. &
          config_flags%chem_opt /= DUST .and. &
          config_flags%chem_opt /= CHEM_TRACE2 .and. &
          config_flags%chem_opt /= CO2_TRACER  .and. &
          config_flags%chem_opt /= GHG_TRACER ) then

!
! chemical mechanisms
!

! save old concentrations for chemistry tendencies
#if (EM_CORE == 1)
   DO nv=PARAM_FIRST_SCALAR,num_chem_ct
      chem_old(its:ite,kts:kte,jts:jte,nv) = chem(its:ite,kts:kte,jts:jte,chem_ct_indices(nv))
   ENDDO
#endif
        if ( cam_mam_aerosols ) &
           del_h2so4_gasprod(:,:,:) = chem(:,:,:,p_sulf)

        if(config_flags%gaschem_onoff>0)then

          call mechanism_driver(grid%id,curr_secs,ktau,grid%dt,grid%ktauc,dtstepc,config_flags, &
#if (NMM_CORE == 1)
              grid%gmt,ijulian,t_phy,moist_trans,p8w,t8w,                         &
#endif
#if (EM_CORE == 1)
              grid%gmt,ijulian,t_phy,moist,p8w,t8w,grid%gd_cldfr,                 &
#endif
              p_phy,chem,rho,dz8w,grid%dx,g,                                          &
              zmid,z_at_w,grid%xlat,grid%xlong,                                       &
              vdrog3,vcsulf_old,vcso2_old,vch2o2_old,grid%ttday,grid%tcosz,           &
              grid%ph_macr,grid%ph_o31d,grid%ph_o33p,grid%ph_no2,grid%ph_no3o2,       &
              grid%ph_no3o,grid%ph_hno2,grid%ph_hno3,grid%ph_hno4,grid%ph_h2o2,       &
              grid%ph_ch2or,grid%ph_ch2om,grid%ph_ch3cho,grid%ph_ch3coch3,            &
              grid%ph_ch3coc2h5,grid%ph_hcocho,grid%ph_ch3cocho,grid%ph_hcochest,     &
              grid%ph_ch3o2h,grid%ph_ch3coo2h,grid%ph_ch3ono2,grid%ph_hcochob,        &
              grid%ph_n2o5,grid%ph_o2,grid%backg_oh,grid%backg_h2o2,grid%backg_no3,   &
              grid%addt,grid%addx,grid%addc,grid%etep,                                &
              grid%oltp,grid%olip,grid%cslp,grid%limp,grid%hc5p,grid%hc8p,grid%tolp,  &
              grid%xylp,grid%apip,grid%isop,grid%hc3p,grid%ethp,grid%o3p,grid%tco3,   &
              grid%mo2,grid%o1d,grid%olnn,grid%rpho,grid%xo2,                         &
              grid%ketp,grid%olnd,                                                    &
              ids,ide, jds,jde, kds,kde,                                              &
              ims,ime, jms,jme, kms,kme,                                              &
              its,ite,jts,jte,kts,kte        )
!

#ifdef WRF_KPP

         if( config_flags%chem_opt == CBMZ_MOSAIC_4BIN_VBS2_KPP .or. &
             config_flags%chem_opt == SAPRC99_MOSAIC_4BIN_VBS2_KPP .or. &
            config_flags%chem_opt ==SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP ) then
          do k=kts,kte
           do i=its,ite
            do j=jts,jte
              chem(i,k,j,p_psd1)=0.0
              chem(i,k,j,p_psd2)=0.0
            enddo
           enddo
          enddo
         endif

   CALL wrf_debug(15,'calling kpp_mechanism_driver')

CALL kpp_mechanism_driver (chem,                                                      &
   grid%id,dtstepc,config_flags,                                                      &
   p_phy,t_phy,rho,                                                                   &
#if (NMM_CORE == 1)
   moist_trans,                                                                       &
#endif
#if (EM_CORE == 1)
   moist,                                                                             &
#endif
   vdrog3, ldrog, vdrog3_vbs, ldrog_vbs,                                              &
!
#include <call_to_kpp_mech_drive.inc>
!
   ids,ide, jds,jde, kds,kde,                                                         &
   ims,ime, jms,jme, kms,kme,                                                         &
   its,ite,jts,jte,kts,kte)
          if( config_flags%chem_opt == MOZART_KPP .or. &
              config_flags%chem_opt == MOZCART_KPP .or. &
              config_flags%chem_opt == MOZART_MOSAIC_4BIN_VBS0_KPP ) then
             call mozcart_lbc_set( chem, num_chem, grid%id, &
                                   ims, ime, jms, jme, kms, kme,    &
                                    its, ite, jts, jte, kts )
          end if
    if(config_flags%chem_opt == 301 ) then
       chem(its:ite,kts:kte,jts:jte,p_sulf)=vcsulf_old(its:ite,kts:kte,jts:jte)
       chem(its:ite,kts:kte,jts:jte,p_so2)=vcso2_old(its:ite,kts:kte,jts:jte)
!      chem(its:ite,kts:kte,jts:jte,p_h2o2)=vch2o2_old(its:ite,kts:kte,jts:jte)
   endif

   IF(config_flags%conv_tr_aqchem == 0 ) THEN
      so2so4_selecta: SELECT CASE(config_flags%chem_opt)
      CASE (RADM2SORG,RADM2SORG_KPP,RACMSORG_KPP,RACM_SOA_VBS_KPP)
         CALL wrf_debug(15,'gocart so2-so4 conversion')
         CALL  so2so4(0,chem,p_so2,p_sulf,p_h2o2,p_QC,T_PHY,MOIST,           &
              grid%gd_cloud_b, grid%gd_cldfr,                                    &
              NUM_CHEM,NUM_MOIST,                                                &
              ids,ide, jds,jde, kds,kde,                                         &
              ims,ime, jms,jme, kms,kme,                                         &
              its,ite, jts,jte, kts,kte                                          )
      CASE DEFAULT
         CALL wrf_debug(15,'no gocart so2-so4 conversion')
      END SELECT so2so4_selecta
   else IF(config_flags%conv_tr_aqchem == 1 ) THEN
      so2so4_selectb: SELECT CASE(config_flags%chem_opt)
      CASE (RADM2SORG,RADM2SORG_KPP,RACMSORG_KPP,RACM_SOA_VBS_KPP)
         CALL wrf_debug(15,'gocart so2-so4 conversion')
         CALL  so2so4(1,chem,p_so2,p_sulf,p_h2o2,p_QC,T_PHY,MOIST,           &
              grid%gd_cloud_b, grid%gd_cldfr,                                    &
              NUM_CHEM,NUM_MOIST,                                                &
              ids,ide, jds,jde, kds,kde,                                         &
              ims,ime, jms,jme, kms,kme,                                         &
              its,ite, jts,jte, kts,kte                                          )
      CASE DEFAULT
         CALL wrf_debug(15,'no gocart so2-so4 conversion')
      END SELECT so2so4_selectb

   ENDIF

!
#endif
#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
       if( (its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and.                          &
           (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and.                          &
           (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)  ) then
          call wrf_debug(15,'calling chem_dbg after mechanism_driver')
          call chem_dbg(CHEM_DBG_I,CHEM_DBG_J,CHEM_DBG_K,grid%dt,ktau,                &
               dz8w,t_phy,p_phy,rho,chem,emis_ant,                                    &
               ids,ide, jds,jde, kds,kde,                                             &
               ims,ime, jms,jme, kms,kme,                                             &
               its,ite, jts,jte, kts,kte,                                             &
               config_flags%kemit,                                                    &
               grid%ph_macr,grid%ph_o31d,grid%ph_o33p,grid%ph_no2,grid%ph_no3o2,      &
               grid%ph_no3o,grid%ph_hno2,grid%ph_hno3,grid%ph_hno4,grid%ph_h2o2,      &
               grid%ph_ch2or,grid%ph_ch2om,grid%ph_ch3cho,grid%ph_ch3coch3,           &
               grid%ph_ch3coc2h5,grid%ph_hcocho,grid%ph_ch3cocho,                     &
               grid%ph_hcochest,grid%ph_ch3o2h,grid%ph_ch3coo2h,                      &
               grid%ph_ch3ono2,grid%ph_hcochob,grid%ph_n2o5,grid%ph_o2                &
                                                                                      )
       end if
#endif
        endif !gaschem_onoff


! accumulate gas phase chemistry tendencies
#if (EM_CORE == 1)
   DO nv=PARAM_FIRST_SCALAR,num_chem_ct
      grid%chem_ct(its:ite,kts:kte,jts:jte,nv) = grid%chem_ct(its:ite,kts:kte,jts:jte,nv) + &
                                                        (chem(its:ite,kts:kte,jts:jte,chem_ct_indices(nv)) - &
                                                     chem_old(its:ite,kts:kte,jts:jte,nv))
   ENDDO
#endif
if ( cam_mam_aerosols ) &
           del_h2so4_gasprod(:,:,:) = chem(:,:,:,p_sulf) - del_h2so4_gasprod(:,:,:)

! initialize gas_aqfrac if either cldchem or wetscav is on
        if ( (config_flags%cldchem_onoff > 0) .or.                                    &
             (config_flags%wetscav_onoff > 0) ) then
            gas_aqfrac(its:ite,kts:kte,jts:jte,:) = 0.0
        end if
        
        !BSINGH -  The following call is for wet deposition of gases for MAM
        !aerosol model ONLY. In other WRF aerosol models wet deposition of gasses
        !is handled after dry deposition. This call doesn't belong here as there should 
        !be a formal driver for gas wet deposition and this call should be called
        !from that driver (something for future mods)

        !**IMPORTANT-WARNING**: We are going to operate on only 6 gases (following MOSAIC -
        !module_mosaic_wetscav.F). MOSAIC actually operates upon 7 gases but MAM
        !doesn't have msa, therefore MAM will operate on only 6 gases        
        !so2,h2o2,h2so4,hno3,hcl,nh3 (msa not included)
        
        !For adding more gasses, look into phys/module_cam_support.F
        if (((config_flags%chem_opt == CBMZ_CAM_MAM3_NOAQ) .or. &
             (config_flags%chem_opt == CBMZ_CAM_MAM3_AQ  ) .or. &
             (config_flags%chem_opt == CBMZ_CAM_MAM7_NOAQ) .or. &
             (config_flags%chem_opt == CBMZ_CAM_MAM7_AQ  )).and. &
             (config_flags%cu_physics == CAMZMSCHEME)) then
           !BSINGH - Wetdep for gases (MAM only)
           
           call cam_mam_gas_wetdep_driver(                &
                !Intent in-outs
                chem,                                     &
                !Intent ins
                dtstepc, config_flags, grid%ht,grid%XLAT, &
                grid%nevapr3d, grid%rprdsh,               &
                grid%rprddp3d, grid%prain3d, grid%z,      & 
                p_phy, t_phy, grid%alt, moist, scalar,    &
                ids,ide, jds,jde, kds,kde,                &
                ims,ime, jms,jme, kms,kme,                &
                its,ite, jts,jte, kts,kte                 )
           
        endif
       
        !BSINGH - For WRFCuP scheme
        
        !
        !   now do convective cloud processing for cup cumulus scheme
        !
        !   issues for future work
        !      convective cloud processing should be done at same place (in the code)
        !         and same time frequency for cup and grell schemes,
        !         but doing it here is ok for initial testing
        !      there should be a "convchem_driver" or "ctrans_driver" routine
        !         to interface to the various options
        !
        if( config_flags%chem_conv_tr>0 .and. &
             config_flags%cu_physics==kfcupscheme ) then
           
           call chem_cup_driver(                                                        &
                grid%id, ktau, grid%ktauc, grid%dt, dtstepc, config_flags,              &
                t_phy, p_phy, rho, rri, dz8w, zmid, z_at_w,                             &
#if (NMM_CORE == 1)
                moist_trans, grid%cldfra, grid%ph_no2,                                  &
#endif
#if (EM_CORE == 1)
                moist, grid%cldfra, grid%ph_no2,                                        &
#endif
                chem, grid%chem_cupflag,                                                &
                grid%cupflag, grid%shall, grid%tcloud_cup, grid%nca, grid%wact_cup,     &
                grid%cldfra_cup, grid%updfra_cup, grid%qc_ic_cup, grid%qc_iu_cup,       &
                grid%mfup_cup, grid%mfup_ent_cup, grid%mfdn_cup, grid%mfdn_ent_cup,     &
                grid%fcvt_qc_to_pr_cup, grid%fcvt_qc_to_qi_cup, grid%fcvt_qi_to_pr_cup, &
                !BSINGH(12/03/2013): Commenting out *_ams_* variables and replacing them with 
                !*_ic_cup* variables
                !grid%so4_a_ams_ic_cup, grid%so4_cw_ams_ic_cup,                         &
                !grid%nh4_a_ams_ic_cup, grid%nh4_cw_ams_ic_cup,                         &
                !grid%no3_a_ams_ic_cup, grid%no3_cw_ams_ic_cup,                         &
                !grid%oa_a_ams_ic_cup,  grid%oa_cw_ams_ic_cup,                          &
                grid%co_a_ic_cup,       grid%hno3_a_ic_cup,                             &
                grid%so4_a_1to4_ic_cup, grid%so4_cw_1to4_ic_cup,                        &
                grid%nh4_a_1to4_ic_cup, grid%nh4_cw_1to4_ic_cup,                        &
                grid%no3_a_1to4_ic_cup, grid%no3_cw_1to4_ic_cup,                        &
                grid%oa_a_1to4_ic_cup,  grid%oa_cw_1to4_ic_cup,                         &
                grid%oin_a_1to4_ic_cup, grid%oin_cw_1to4_ic_cup,                        &
                grid%bc_a_1to4_ic_cup,  grid%bc_cw_1to4_ic_cup,                         &
                grid%na_a_1to4_ic_cup,  grid%na_cw_1to4_ic_cup,                         &
                grid%cl_a_1to4_ic_cup,  grid%cl_cw_1to4_ic_cup,                         &
                grid%water_1to4_ic_cup,                                                 &
                grid%so4_a_5to6_ic_cup, grid%so4_cw_5to6_ic_cup,                        &
                grid%nh4_a_5to6_ic_cup, grid%nh4_cw_5to6_ic_cup,                        &
                grid%no3_a_5to6_ic_cup, grid%no3_cw_5to6_ic_cup,                        &
                grid%oa_a_5to6_ic_cup,  grid%oa_cw_5to6_ic_cup,                         &
                grid%oin_a_5to6_ic_cup, grid%oin_cw_5to6_ic_cup,                        &
                grid%bc_a_5to6_ic_cup,  grid%bc_cw_5to6_ic_cup,                         &
                grid%na_a_5to6_ic_cup,  grid%na_cw_5to6_ic_cup,                         &
                grid%cl_a_5to6_ic_cup,  grid%cl_cw_5to6_ic_cup,                         &
                grid%water_5to6_ic_cup,                                                 & 
                ids,ide, jds,jde, kds,kde,                                              &
                ims,ime, jms,jme, kms,kme,                                              &
                its,ite, jts,jte, kts,kte                                               )
           
        else
           grid%chem_cupflag = 0
        endif
        
        
        !BSINGH - Ends
        


!
!   now do cloud chemistry
!
! call budget_calc(0,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
!      grid%tmp_budget,grid%tend_cldchem,  &
!      tend_rate_dep,snowchem_old,  &
!      ims,ime, jms,jme, kms,kme,                              &
!      its,ite, jts,jte, kts,kte                                )

        if (config_flags%cldchem_onoff > 0) then

        call cloudchem_driver(                                                        &
               grid%id, ktau, grid%ktauc, grid%dt, dtstepc, config_flags,             &
               t_phy, p_phy, rho, rri, dz8w,                                          &
               p8w,grid%prain3d,scalar,grid%dvmrdt_sv13d,grid%dvmrcwdt_sv13d, grid%f_ice_phy,   & !Balwinder.Singh@pnnl.gov: Variables required for CAM-MAM cloud chemistry
               grid%f_rain_phy,grid%cldfrai, grid%cldfral,                            &
#if (NMM_CORE == 1)
               moist_trans, grid%cldfra, grid%ph_no2,                                 &
#endif
#if (EM_CORE == 1)
               moist, grid%cldfra, grid%cldfra_mp_all, grid%ph_no2,                   &
#endif
               chem, gas_aqfrac, numgas_mam,grid%is_CAMMGMP_used,                     &!BSINGH:01/31/2013: Added is_CAMMGMP_used for CAM_MAM_cloudchem
               ids,ide, jds,jde, kds,kde,                                             &
               ims,ime, jms,jme, kms,kme,                                             &
               its,ite, jts,jte, kts,kte                                              )

       endif

! call budget_calc(1,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
!      grid%tmp_budget,grid%tend_cldchem,  &
!      tend_rate_dep,snowchem_old,  &
!      ims,ime, jms,jme, kms,kme,                              &
!      its,ite, jts,jte, kts,kte                                )


!------------
!   now do aerosols
!
! call budget_calc(0,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
!      grid%tmp_budget,grid%tend_aerchem,  &
!      tend_rate_dep,snowchem_old,  &
!      ims,ime, jms,jme, kms,kme,                              &
!      its,ite, jts,jte, kts,kte                                )


	if(config_flags%aerchem_onoff>0)then

        call aerosols_driver (grid%id,curr_secs,ktau,grid%dt,grid%ktauc,              &
             config_flags,dtstepc,grid%dx,                                            &
#if (NMM_CORE==1)
              rri,t_phy,moist_trans,grid%aerwrf,p8w,t8w,                              &
#endif
#if (EM_CORE == 1)
              rri,t_phy,moist,grid%aerwrf,p8w,t8w,                                    &
#endif
              p_phy,chem,rho,dz8w, rh,                                                & 
              zmid,z_at_w,pbl_h,grid%cldfra,grid%cldfra_mp_all,grid%vbs_nbin,                                              &
              grid%h2oaj,grid%h2oai,grid%h2oaj_adjust,grid%h2oai_adjust,grid%h2oac_adjust,      &
              grid%nu3,grid%ac3,grid%cor3,grid%asulf,           &
              grid%ahno3,grid%anh3,grid%cvaro1,grid%cvaro2,grid%cvalk1,grid%cvole1,   &
              grid%cvapi1,grid%cvapi2,grid%cvlim1,grid%cvlim2,vcsulf_old,             &
              vdrog3,vdrog3_vbs,grid%br_rto,grid%dgnum4d,grid%dgnumwet4d,wetdens_ap,  &
              del_h2so4_gasprod,grid%dvmrdt_sv13d,grid%dvmrcwdt_sv13d,                &
              grid%is_CAMMGMP_used,                                                   &!BSINGH:01/31/2013: Added is_CAMMGMP_used for cam_mam_aerchem_driver
              ids,ide, jds,jde, kds,kde,                                              &
              ims,ime, jms,jme, kms,kme,                                              &
              its,ite,jts,jte,kts,kte,                         & 
              !czhao 
              grid%dgnuc,grid%dgacc,grid%dgcor  )

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
       if( (its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and.                          &
           (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and.                          &
           (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)  ) then
          call wrf_debug(15,'calling chem_dbg after aerosols_driver')
          call chem_dbg(CHEM_DBG_I,CHEM_DBG_J,CHEM_DBG_K,grid%dt,ktau,                &
               dz8w,t_phy,p_phy,rho,chem,emis_ant,                                    &
               ids,ide, jds,jde, kds,kde,                                             &
               ims,ime, jms,jme, kms,kme,                                             &
               its,ite, jts,jte, kts, kte,                                            &
               config_flags%kemit,                                                    &
               grid%ph_macr,grid%ph_o31d,grid%ph_o33p,grid%ph_no2,grid%ph_no3o2,      &
               grid%ph_no3o,grid%ph_hno2,grid%ph_hno3,grid%ph_hno4,grid%ph_h2o2,      &
               grid%ph_ch2or,grid%ph_ch2om,grid%ph_ch3cho,grid%ph_ch3coch3,           &
               grid%ph_ch3coc2h5,grid%ph_hcocho,grid%ph_ch3cocho,                     &
               grid%ph_hcochest,grid%ph_ch3o2h,grid%ph_ch3coo2h,                      &
               grid%ph_ch3ono2,grid%ph_hcochob,grid%ph_n2o5,grid%ph_o2                &
                                                                                      )
       end if
#endif
       endif


! call budget_calc(1,config_flags%chem_opt,chem,rho,dz8w,grid%dt,  &
!      grid%tmp_budget,grid%tend_aerchem,  &
!      tend_rate_dep,snowchem_old,  &
!      ims,ime, jms,jme, kms,kme,                              &
!      its,ite, jts,jte, kts,kte                                )


!
!   now do wet removal
call budget_calc(0,config_flags%chem_opt,chem,rho,dz8w,dtstepc, &
     grid%tmp_budget,grid%tend_wetdls, &
     tend_rate_dep,snowchem_old,  &
     ims,ime, jms,jme, kms,kme,                              &
     its,ite, jts,jte, kts,kte                                )

	if (config_flags%wetscav_onoff > 0) then
        call wetscav_driver (grid%id,ktau,grid%dt,grid%ktauc,config_flags,dtstepc,    &
#if (NMM_CORE == 1)
              rri,t_phy,moist_trans,p8w,t8w,                                          &
#endif
#if (EM_CORE == 1)
              rri,t_phy,moist,p8w,t8w,                                                &
#endif
              grid%dx, grid%dy,                                                       &
              p_phy,chem,rho,grid%cldfra,grid%cldfra2,                                &
              grid%rainprod,grid%evapprod,grid%hno3_col_mdel,                         &
!czhao              grid%qlsink,grid%precr,grid%preci,grid%precs,grid%precg,                &
              grid%qlsink,grid%ilsink,grid%precr,grid%preci,grid%precs,grid%precg,                &
              grid%r2evap3d,grid%i2evap3d,                &
              gas_aqfrac, numgas_mam,dz8w,                                         &
              grid%h2oaj,grid%h2oai,grid%nu3,grid%ac3,grid%cor3,                      &
              grid%asulf,grid%ahno3,grid%anh3,grid%cvaro1,grid%cvaro2,                &
              grid%cvalk1,grid%cvole1,grid%cvapi1,grid%cvapi2,                        &
              grid%cvlim1,grid%cvlim2,                                                &
              grid%wd_no3_sc,grid%wd_so4_sc,                                          &
              grid%qv_b4mp, grid%qc_b4mp, grid%qi_b4mp, grid%qs_b4mp,                 &
!======================================================================================
!Variables required for CAM_MAM_WETSCAV
              grid%p_hyd,scalar,grid%dgnum4d,grid%dgnumwet4d,grid%dlf,grid%dlf2,      &
              grid%qme3d,grid%prain3d,grid%nevapr3d,grid%rate1ord_cw2pr_st3d,         &
              grid%shfrc3d,grid%cmfmc,grid%cmfmc2,grid%evapcsh,grid%icwmrsh,          &
              grid%rprdsh,grid%evapcdp3d,grid%icwmrdp3d,grid%rprddp3d,grid%fracis3d,  &               
              grid%f_ice_phy,grid%f_rain_phy,grid%cldfrai,grid%cldfral,               &
              grid%cldfra_mp_all,grid%is_CAMMGMP_used,                                &
!======================================================================================
              ids,ide, jds,jde, kds,kde,                                              &
              ims,ime, jms,jme, kms,kme,                                              &
              its,ite, jts,jte, kts,kte)                                              
       
       endif

       !BSINGH - zm_conv_tend_2 call HAS to be called after wetscavenging call as it uses 'fracis3d'
       !variable from wetscav driver.
       if (((config_flags%chem_opt == CBMZ_CAM_MAM3_NOAQ) .or. &
            (config_flags%chem_opt == CBMZ_CAM_MAM3_AQ  ) .or. &
            (config_flags%chem_opt == CBMZ_CAM_MAM7_NOAQ) .or. &
            (config_flags%chem_opt == CBMZ_CAM_MAM7_AQ  )).and. &
            (config_flags%cu_physics == CAMZMSCHEME)) then
          !BSINGH - CALL ZM cumulus scheme for convective transport of the aerosols
          call zm_conv_tend_2(grid%itimestep, grid%dt, p8w, grid%fracis3d, grid%dp3d,    &
               grid%du3d, grid%ed3d, grid%eu3d, grid%md3d, grid%mu3d, grid%dsubcld2d,    &
               grid%ideep2d, grid%jt2d, grid%maxg2d, grid%lengath2d, moist, scalar, chem,&
               ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte)     
       endif

call budget_calc(1,config_flags%chem_opt,chem,rho,dz8w,dtstepc, &
     grid%tmp_budget,grid%tend_wetdls, &
     tend_rate_dep,snowchem_old,  &
     ims,ime, jms,jme, kms,kme,                              &
     its,ite, jts,jte, kts,kte                                )
!add large-scale wet scavenging 
do nv=1,5
do n=1,4
do j=jts,jte
do i=its,ite
if (tend_rate_dep(i,j,n,nv).gt.0) tend_rate_dep(i,j,n,nv)=0.0
enddo
enddo
enddo
enddo
!dust
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_dust)= &
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_dust)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,1)
!bcphi-wet from bccw wet
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_bcphi)= &
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_bcphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,3)
!bcphi-dry from bc wet
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcphi)= &
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,2)
!ocphi-wet from occw wet
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_ocphi)= &
grid%wetdep_flux(its:ite,jts:jte,1:4,p_wetdep_ocphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,5)
!ocphi-dry from oc wet
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocphi)= &
grid%drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocphi)+ &
(-1.)*tend_rate_dep(its:ite,jts:jte,1:4,4)

!-----
end if !Chemistry time step check!
!
! do so2 to sulf conversion for full volc case
! since we do not have h2o2 as a variable, pass in p_h2o2 as zero
! will have to use backgrund value
!
        if(config_flags%chem_opt == CHEM_VOLC)then
          CALL wrf_debug(15,'gocart so2-so4 conversion')
          CALL  so2so4(0,chem,p_so2,p_sulf,p_h2o2,p_QC,T_PHY,MOIST,           &
         grid%gd_cloud_b, grid%gd_cldfr,                                    &
         NUM_CHEM,NUM_MOIST,                                                &
         ids,ide, jds,jde, kds,kde,                                         &
         ims,ime, jms,jme, kms,kme,                                         &
         its,ite, jts,jte, kts,kte                                          )
        endif
!
!   now do wet removal; first LS if there is no explicit aqeous phase
!
        if(config_flags%wetscav_onoff<0)then
           call wrf_debug(15,'calculate LS wet deposition')
           call wetdep_ls(grid%dt,chem,grid%rainncv,moist,rho,num_moist, &
                num_chem,numgas,dz8w,vvel,grid%chem_opt,             &
                ids,ide, jds,jde, kds,kde,                                        &
                ims,ime, jms,jme, kms,kme,                                        &
                its,ite, jts,jte, kts,kte                                         )
         endif
!
! Sum up the aerosol mass for radiation and diagnostic purposes. Unlike
! aerosol_driver, which is called every dtchem, this must be done every
! time step because of emissions and deposition.
!
call budget_calc(0,config_flags%chem_opt,chem,rho,dz8w,grid%dt, &
     grid%tmp_budget,grid%tend_tran, &
     tend_rate_dep,snowchem_old,  &
     ims,ime, jms,jme, kms,kme,                              &
     its,ite, jts,jte, kts,kte                                )

      call sum_pm_driver( config_flags,                                               &
        rri,              chem,             grid%h2oaj,        grid%h2oai,  grid%h2oac_adjust,          &
           grid%pm2_5_dry, grid%pm2_5_water,                                          &
           grid%pm2_5_dry_dust,grid%pm2_5_dry_oc,grid%pm2_5_dry_ec,                   & 
           grid%pm2_5_dry_seas,grid%pm2_5_dry_so4f,grid%pm2_5_dry_no3f,grid%pm2_5_dry_nh4f,   & 
           grid%pm10,                                                                 &
           !czhao 
           grid%dgnuc,grid%dgacc,grid%dgcor,   & 
        grid%hoa_a01,     grid%hoa_a02,     grid%hoa_a03,      grid%hoa_a04,          &
        grid%hoa_a05,     grid%hoa_a06,     grid%hoa_a07,      grid%hoa_a08,          & !BSINGH(12/03/2013) Added 4 more bins for each species for SAPRC 8 bin version
        grid%bboa_a01,    grid%bboa_a02,    grid%bboa_a03,     grid%bboa_a04,         &
        grid%bboa_a05,    grid%bboa_a06,    grid%bboa_a07,     grid%bboa_a08,         &
        grid%soa_a01,     grid%soa_a02,     grid%soa_a03,      grid%soa_a04,          &
        grid%soa_a05,     grid%soa_a06,     grid%soa_a07,      grid%soa_a08,          &
        grid%bbsoa_a01,   grid%bbsoa_a02,   grid%bbsoa_a03,    grid%bbsoa_a04,        &
        grid%bbsoa_a05,   grid%bbsoa_a06,   grid%bbsoa_a07,    grid%bbsoa_a08,        &
        grid%hsoa_a01,    grid%hsoa_a02,    grid%hsoa_a03,     grid%hsoa_a04,         &
        grid%hsoa_a05,    grid%hsoa_a06,    grid%hsoa_a07,     grid%hsoa_a08,         &
        grid%biog_a01,    grid%biog_a02,    grid%biog_a03,     grid%biog_a04,         &
        grid%biog_a05,    grid%biog_a06,    grid%biog_a07,     grid%biog_a08,         &
        grid%asmpsoa_a01, grid%asmpsoa_a02, grid%asmpsoa_a03,  grid%asmpsoa_a04,      & !BSINGH(12/03/2013) - Not adding 5-8 bins for asmpsoa, as it is not req. for 8 bin SAPRC
        grid%arosoa_a01,  grid%arosoa_a02,  grid%arosoa_a03,   grid%arosoa_a04,       &
        grid%arosoa_a05,  grid%arosoa_a06,  grid%arosoa_a07,   grid%arosoa_a08,       &
        grid%totoa_a01,   grid%totoa_a02,   grid%totoa_a03,    grid%totoa_a04,        &
        grid%totoa_a05,   grid%totoa_a06,   grid%totoa_a07,    grid%totoa_a08,        &
        grid%hsoa_c,      grid%hsoa_o,      grid%bbsoa_c,      grid%bbsoa_o,          &
        grid%biog_v1,     grid%biog_v2,     grid%biog_v3,      grid%biog_v4,          &
        grid%ant_v1,      grid%ant_v2,      grid%ant_v3,       grid%ant_v4,           &
        grid%smpa_v1,     grid%smpbb_v1,                                              &
        !BSINGH(12/03/2013) - Added cw aerosols(for VBS)
        grid%hoa_cw01,    grid%hoa_cw02,    grid%hoa_cw03,    grid%hoa_cw04,          &
        grid%hoa_cw05,    grid%hoa_cw06,    grid%hoa_cw07,    grid%hoa_cw08,          &
        grid%bboa_cw01,   grid%bboa_cw02,   grid%bboa_cw03,   grid%bboa_cw04,         &
        grid%bboa_cw05,   grid%bboa_cw06,   grid%bboa_cw07,   grid%bboa_cw08,         &
        grid%soa_cw01,    grid%soa_cw02,    grid%soa_cw03,    grid%soa_cw04,          &
        grid%soa_cw05,    grid%soa_cw06,    grid%soa_cw07,    grid%soa_cw08,          &
        grid%bbsoa_cw01,  grid%bbsoa_cw02,  grid%bbsoa_cw03,  grid%bbsoa_cw04,        &
        grid%bbsoa_cw05,  grid%bbsoa_cw06,  grid%bbsoa_cw07,  grid%bbsoa_cw08,        &
        grid%biog_cw01,   grid%biog_cw02,   grid%biog_cw03,   grid%biog_cw04,         &
        grid%biog_cw05,   grid%biog_cw06,   grid%biog_cw07,   grid%biog_cw08,         &
        grid%hsoa_cw01,   grid%hsoa_cw02,   grid%hsoa_cw03,   grid%hsoa_cw04,         &
        grid%hsoa_cw05,   grid%hsoa_cw06,   grid%hsoa_cw07,   grid%hsoa_cw08,         &
        grid%arosoa_cw01, grid%arosoa_cw02, grid%arosoa_cw03, grid%arosoa_cw04,       &
        grid%arosoa_cw05, grid%arosoa_cw06, grid%arosoa_cw07, grid%arosoa_cw08,       &
        grid%totoa_cw01,  grid%totoa_cw02,  grid%totoa_cw03,  grid%totoa_cw04,        &
        grid%totoa_cw05,  grid%totoa_cw06,  grid%totoa_cw07,  grid%totoa_cw08,        &        
        grid%hsoa_cw_c,   grid%hsoa_cw_o,   grid%bbsoa_cw_c,  grid%bbsoa_cw_o,        &
        grid%biog_cw_v1,                                                              &
        grid%ant_cw_v1,                                                               &
        !BSINGH(12/03/2013) -ENDS
        ids,ide, jds,jde, kds,kde,                                                    &
        ims,ime, jms,jme, kms,kme,                                                    &
        its,ite, jts,jte, kts,kte                                                     )

! Fill top level to prevent spurious interpolation results (no extrapolation)
      do nv=1,num_chem
         do j=jts,jte
            do i=its,ite
                  chem(i,k_end,j,nv)=chem(i,kte,j,nv)
            enddo
         enddo
      enddo
      call wrf_debug(15,'done tileloop in chem_driver')
   if( grid%OPT_PARS_OUT == 1) then
      call wrf_debug(15,'calculate optical output stuff')
    call aer_opt_out(TAUAER300=tauaersw(:,:,:,2), TAUAER400=tauaersw(:,:,:,3) & 
     &    ,TAUAER600=tauaersw(:,:,:,4), TAUAER999=tauaersw(:,:,:,5) &
     &    ,GAER300=gaersw(:,:,:,2), GAER400=gaersw(:,:,:,3), GAER600=gaersw(:,:,:,4), GAER999=gaersw(:,:,:,5) & 
     &    ,WAER300=waersw(:,:,:,2), WAER400=waersw(:,:,:,3), WAER600=waersw(:,:,:,4), WAER999=waersw(:,:,:,5) & 
          ,ext_coeff=grid%ext_coef,bscat_coeff=grid%bscat_coef,asym_par=grid%asym_par   &
          ,num_ext_coef=num_ext_coef,num_bscat_coef=num_bscat_coef,num_asym_par=num_asym_par    &
     &    ,dz8w=dz8w &
     &    ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &    ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &    ,its=its,ite=ite,jts=jts,jte=jte,kts=kts, kte=kte)


   endif
   tracer2: SELECT CASE(config_flags%tracer_opt)
    CASE (TRACER_TEST1, TRACER_TEST2, TRACER_TEST3)
       CALL wrf_debug(15,'tracer mode: reset some tracers')
       call set_tracer(grid%dt,ktau,pbl_h,tracer,t_phy,         &
                        config_flags%tracer_opt,num_tracer,     &
                       zmid,grid%ht,ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,       & ! memory dims
                               its,ite, jts,jte, kts,kte )
   END SELECT tracer2

!-- set upper boundary condition
    if( config_flags%have_bcs_upper )then
       call wrf_debug(15,'set upper boundary condition')
       call tropopause_driver( grid%id, grid%dt, current_date_char,         &
                               t_phy, p_phy, p8w, zmid, z_at_w,             &
                               grid%tropo_lev, grid%tropo_p,  grid%tropo_z, &
                               ids, ide, jds, jde, kds, kde,                &
                               ims, ime, jms, jme, kms, kme,                &
                               its, ite, jts, jte, kts, kte                 )
       call upper_bc_driver  ( grid%id, grid%dt, current_date_char, &
                               chem, p_phy, p8w, grid%tropo_lev,    &
                               ids,ide, jds,jde, kds,kde,           &
                               ims,ime, jms,jme, kms,kme,           &
                               its,ite, jts,jte, kts,kte            )
    endif


# if ( EM_CORE == 1 ) 
   END DO chem_tile_loop_1
#endif
#if (NMM_CORE==1) 
       DO l=1,num_3d_m
       DO k=kts,kte   
         DO j=jts,jte
             DO i=its,ite
                moist(i,j,k,l)=moist_trans(i,k,j,l)
             ENDDO
         ENDDO
       ENDDO
       ENDDO
       DO l=1,num_3d_s
       DO k=kts,kte   
         DO j=jts,jte
             DO i=its,ite
               scalar(i,j,k,l)=scalar_trans(i,k,j,l)
             ENDDO
         ENDDO
       ENDDO
       ENDDO
    deallocate(moist_trans)
    deallocate(scalar_trans)
#endif

    END subroutine chem_driver


!=================================================================================
!budget calculation
!=================================================================================
subroutine budget_calc(iswitch,chem_opt,chem,rho_phy,dz8w,dt, &
                       tmp_chem,tend_chem, &
                       tend_rate_dep,snowchem_old, &
                       ims,ime,jms,jme,kms,kme,                  &
                       its,ite,jts,jte,kts,kte                                )

USE module_state_description
USE module_data_mosaic_asect, only: lptr_dust_aer,lptr_tr1r1_aer,lptr_tr1r2_aer,lptr_tr1r3_aer,lptr_tr1r4_aer,nsize_aer

IMPLICIT NONE

integer,intent(in) :: iswitch,chem_opt

INTEGER,INTENT(IN) :: ims,ime, jms,jme, kms,kme,    &
                      its,ite, jts,jte, kts,kte

REAL,DIMENSION( ims:ime , kms:kme , jms:jme ),INTENT(IN   ) :: dz8w,rho_phy
REAL,      INTENT(IN   ) :: dt

REAL, DIMENSION( ims:ime,kms:kme,jms:jme,num_chem),INTENT(IN) :: chem
REAL, DIMENSION( ims:ime,kms:kme,jms:jme,num_tmp_budget),INTENT(INOUT) :: tmp_chem
REAL, DIMENSION( ims:ime,kms:kme,jms:jme,num_tend_tran),INTENT(INOUT) :: tend_chem
!REAL, DIMENSION( ims:ime,jms:jme,4,3),INTENT(INOUT) :: tend_rate_dep,snowchem_old
REAL, DIMENSION( ims:ime,jms:jme,4,5),INTENT(INOUT) :: tend_rate_dep,snowchem_old

integer :: nv,kk,iphase,isize,itype,l,lcw,p1st,ibin,bsize,esize
real :: massfrac

itype=1
p1st=param_first_scalar

tend_rate_dep(its:ite,jts:jte,:,:)=0.0
aer_budget_select: SELECT CASE(chem_opt)
!tendency (as rain) needs to be accumulated since output is not as frequent as timestep (e.g., 1 h v.s. 150 s)

CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_KPP, &
      RACMSORG_AQ,RACMSORG_KPP,RACM_ESRLSORG_KPP,CBMZSORG,CBMZSORG_AQ)

if (iswitch.eq.0) then !before the process

tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
    chem(its:ite, kts:kte, jts:jte,p_dustj)+chem(its:ite, kts:kte, jts:jte,p_soila)
if (chem_opt.ne.RADM2SORG)  & 
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
    chem(its:ite, kts:kte, jts:jte,p_dustcwj)+chem(its:ite, kts:kte, jts:jte,p_soilcw)

snowchem_old(its:ite,jts:jte,:,:)=0.0
do kk=kts,kte 
if (chem_opt.ne.RADM2SORG) then
!dust
snowchem_old(its:ite,jts:jte,1,1)=snowchem_old(its:ite,jts:jte,1,1)+&
 (chem(its:ite, kk, jts:jte,p_dustj)+chem(its:ite, kk, jts:jte,p_dustcwj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2,1)=0.0
snowchem_old(its:ite,jts:jte,3,1)=snowchem_old(its:ite,jts:jte,3,1)+&
 (chem(its:ite, kk, jts:jte,p_soila)+chem(its:ite, kk, jts:jte,p_soilcw))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,4,1)=0.0
!bc
snowchem_old(its:ite,jts:jte,1,2)=snowchem_old(its:ite,jts:jte,1,2)+&
 (chem(its:ite, kk, jts:jte,p_eci)+chem(its:ite, kk, jts:jte,p_ecj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,2)=0.0
!bccw
snowchem_old(its:ite,jts:jte,1,3)=snowchem_old(its:ite,jts:jte,1,3)+&
 (chem(its:ite, kk, jts:jte,p_eccwi)+chem(its:ite, kk, jts:jte,p_eccwj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,3)=0.0
!oc
snowchem_old(its:ite,jts:jte,1,4)=snowchem_old(its:ite,jts:jte,1,4)+&
 (chem(its:ite, kk, jts:jte,p_orgpai)+chem(its:ite, kk, jts:jte,p_orgpaj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,4)=0.0
!occw
snowchem_old(its:ite,jts:jte,1,5)=snowchem_old(its:ite,jts:jte,1,5)+&
 (chem(its:ite, kk, jts:jte,p_orgpacwi)+chem(its:ite, kk, jts:jte,p_orgpacwj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,5)=0.0
else 
!dust
snowchem_old(its:ite,jts:jte,1,1)=snowchem_old(its:ite,jts:jte,1,1)+&
  chem(its:ite, kk, jts:jte,p_dustj)*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2,1)=0.0
snowchem_old(its:ite,jts:jte,3,1)=snowchem_old(its:ite,jts:jte,3,1)+&
  chem(its:ite, kk, jts:jte,p_soila)*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,4,1)=0.0
!bc
snowchem_old(its:ite,jts:jte,1,2)=snowchem_old(its:ite,jts:jte,1,2)+&
  (chem(its:ite, kk, jts:jte,p_eci)+chem(its:ite, kk, jts:jte,p_ecj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,2)=0.0
!bccw
snowchem_old(its:ite,jts:jte,1:4,3)=0.0
!oc
snowchem_old(its:ite,jts:jte,1,4)=snowchem_old(its:ite,jts:jte,1,4)+&
  (chem(its:ite, kk, jts:jte,p_orgpai)+chem(its:ite, kk, jts:jte,p_orgpaj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,4)=0.0
!occw
snowchem_old(its:ite,jts:jte,1:4,5)=0.0
endif

enddo

endif

if (iswitch.eq.1) then !after the process

tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)=tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
    ((chem(its:ite, kts:kte, jts:jte,p_dustj)+chem(its:ite, kts:kte, jts:jte,p_soila))-&
      tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dust))

if (chem_opt.ne.RADM2SORG) then 
tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)=tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
    ((chem(its:ite, kts:kte, jts:jte,p_dustcwj)+chem(its:ite, kts:kte, jts:jte,p_soilcw))-&
      tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw))
endif

tend_rate_dep(its:ite,jts:jte,:,:)=0.0
do kk=kts,kte

if (chem_opt.ne.RADM2SORG) then
!dust
tend_rate_dep(its:ite,jts:jte,1,1)=tend_rate_dep(its:ite,jts:jte,1,1)+&
 (chem(its:ite, kk, jts:jte,p_dustj)+chem(its:ite, kk, jts:jte,p_dustcwj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2,1)=0.0
tend_rate_dep(its:ite,jts:jte,3,1)=tend_rate_dep(its:ite,jts:jte,3,1)+&
 (chem(its:ite, kk, jts:jte,p_soila)+chem(its:ite, kk, jts:jte,p_soilcw))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,4,1)=0.0
!bc
tend_rate_dep(its:ite,jts:jte,1,2)=tend_rate_dep(its:ite,jts:jte,1,2)+&
 (chem(its:ite, kk, jts:jte,p_eci)+chem(its:ite, kk, jts:jte,p_ecj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,2)=0.0
!bccw
tend_rate_dep(its:ite,jts:jte,1,3)=tend_rate_dep(its:ite,jts:jte,1,3)+&
 (chem(its:ite, kk, jts:jte,p_eccwi)+chem(its:ite, kk, jts:jte,p_eccwj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,3)=0.0
!oc
tend_rate_dep(its:ite,jts:jte,1,4)=tend_rate_dep(its:ite,jts:jte,1,4)+&
 (chem(its:ite, kk, jts:jte,p_orgpai)+chem(its:ite, kk, jts:jte,p_orgpacwi)+chem(its:ite, kk, jts:jte,p_orgpaj)+chem(its:ite, kk, jts:jte,p_orgpacwj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,4)=0.0
!occw
tend_rate_dep(its:ite,jts:jte,1,5)=tend_rate_dep(its:ite,jts:jte,1,5)+&
 (chem(its:ite, kk, jts:jte,p_orgpacwi)+chem(its:ite, kk, jts:jte,p_orgpacwj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,5)=0.0
else 
!dust
tend_rate_dep(its:ite,jts:jte,1,1)=tend_rate_dep(its:ite,jts:jte,1,1)+&
  chem(its:ite, kk, jts:jte,p_dustj)*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2,1)=0.0
tend_rate_dep(its:ite,jts:jte,3,1)=tend_rate_dep(its:ite,jts:jte,3,1)+&
  chem(its:ite, kk, jts:jte,p_soila)*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,4,1)=0.0
!bc
tend_rate_dep(its:ite,jts:jte,1,2)=tend_rate_dep(its:ite,jts:jte,1,2)+&
  (chem(its:ite, kk, jts:jte,p_eci)+chem(its:ite, kk, jts:jte,p_ecj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,2)=0.0
!bccw
tend_rate_dep(its:ite,jts:jte,1:4,3)=0.0
!oc
tend_rate_dep(its:ite,jts:jte,1,4)=tend_rate_dep(its:ite,jts:jte,1,4)+&
  (chem(its:ite, kk, jts:jte,p_orgpai)+chem(its:ite, kk, jts:jte,p_orgpaj))*rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,4)=0.0
!occw
tend_rate_dep(its:ite,jts:jte,1:4,5)=0.0
endif

enddo
tend_rate_dep(its:ite,jts:jte,:,:)=(tend_rate_dep(its:ite,jts:jte,:,:)-snowchem_old(its:ite,jts:jte,:,:))/dt
endif

CASE (CBMZ_MOSAIC_8BIN, &
      CBMZ_MOSAIC_8BIN_AQ,CBMZ_MOSAIC_DMS_8BIN, &
      CBMZ_MOSAIC_DMS_8BIN_AQ,CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE)

if (iswitch.eq.0) then !before the process

tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)=0.0
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)=0.0
do isize = 1, nsize_aer(itype) 
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then 
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then 
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= & 
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ & 
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then 
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
enddo !isize 

snowchem_old(its:ite,jts:jte,:,:)=0.0
do kk=kts,kte
!dust
 do ibin=1,4
  if (ibin.eq.1) then 
   bsize=1
   esize=4
  endif
  if (ibin.eq.2) then 
   bsize=5
   esize=6
  endif
  if (ibin.eq.3) then 
   bsize=7
   esize=7
  endif
  if (ibin.eq.4) then 
   bsize=8
   esize=8
  endif
  do isize = bsize, esize 
   l=lptr_dust_aer(isize,itype,1) !dry phase
   lcw=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif 
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif 
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif 
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif 
  enddo !isize
 enddo !ibin
!bc
snowchem_old(its:ite,jts:jte,1,2)=snowchem_old(its:ite,jts:jte,1,2)+&
 (chem(its:ite, kk, jts:jte,p_bc_a01)+chem(its:ite, kk, jts:jte,p_bc_a02)+&
  chem(its:ite, kk, jts:jte,p_bc_a03)+chem(its:ite, kk, jts:jte,p_bc_a04)+&
  chem(its:ite, kk, jts:jte,p_bc_a05)+chem(its:ite, kk, jts:jte,p_bc_a06)+&
  chem(its:ite, kk, jts:jte,p_bc_a07)+chem(its:ite, kk, jts:jte,p_bc_a08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,2)=0.0
!bccw
snowchem_old(its:ite,jts:jte,1,3)=snowchem_old(its:ite,jts:jte,1,3)+&
 (chem(its:ite, kk, jts:jte,p_bc_cw01)+chem(its:ite, kk, jts:jte,p_bc_cw02)+&
  chem(its:ite, kk, jts:jte,p_bc_cw03)+chem(its:ite, kk, jts:jte,p_bc_cw04)+&
  chem(its:ite, kk, jts:jte,p_bc_cw05)+chem(its:ite, kk, jts:jte,p_bc_cw06)+&
  chem(its:ite, kk, jts:jte,p_bc_cw07)+chem(its:ite, kk, jts:jte,p_bc_cw08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,3)=0.0
!oc
snowchem_old(its:ite,jts:jte,1,4)=snowchem_old(its:ite,jts:jte,1,4)+&
 (chem(its:ite, kk, jts:jte,p_oc_a01)+chem(its:ite, kk, jts:jte,p_oc_a02)+&
  chem(its:ite, kk, jts:jte,p_oc_a03)+chem(its:ite, kk, jts:jte,p_oc_a04)+&
  chem(its:ite, kk, jts:jte,p_oc_a05)+chem(its:ite, kk, jts:jte,p_oc_a06)+&
  chem(its:ite, kk, jts:jte,p_oc_a07)+chem(its:ite, kk, jts:jte,p_oc_a08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,4)=0.0
!occw
snowchem_old(its:ite,jts:jte,1,5)=snowchem_old(its:ite,jts:jte,1,5)+&
 (chem(its:ite, kk, jts:jte,p_oc_cw01)+chem(its:ite, kk, jts:jte,p_oc_cw02)+&
  chem(its:ite, kk, jts:jte,p_oc_cw03)+chem(its:ite, kk, jts:jte,p_oc_cw04)+&
  chem(its:ite, kk, jts:jte,p_oc_cw05)+chem(its:ite, kk, jts:jte,p_oc_cw06)+&
  chem(its:ite, kk, jts:jte,p_oc_cw07)+chem(its:ite, kk, jts:jte,p_oc_cw08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,5)=0.0
enddo !kk

endif !iswitch

if (iswitch.eq.1) then !after the process

do isize = 1, nsize_aer(itype) 
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then 
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then 
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
enddo !isize
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)- &
   tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dust)

   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)- &
   tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dustcw)

tend_rate_dep(its:ite,jts:jte,:,:)=0.0
do kk=kts,kte
!dust
 do ibin=1,4
  if (ibin.eq.1) then
   bsize=1
   esize=4
  endif
  if (ibin.eq.2) then
   bsize=5
   esize=6
  endif
  if (ibin.eq.3) then
   bsize=7
   esize=7
  endif
  if (ibin.eq.4) then
   bsize=8
   esize=8
  endif
  do isize = bsize, esize
   l=lptr_dust_aer(isize,itype,1) !dry phase
   lcw=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
  enddo !isize
 enddo !ibin
!bc
tend_rate_dep(its:ite,jts:jte,1,2)=tend_rate_dep(its:ite,jts:jte,1,2)+&
 (chem(its:ite, kk, jts:jte,p_bc_a01)+chem(its:ite, kk, jts:jte,p_bc_a02)+&
  chem(its:ite, kk, jts:jte,p_bc_a03)+chem(its:ite, kk, jts:jte,p_bc_a04)+&
  chem(its:ite, kk, jts:jte,p_bc_a05)+chem(its:ite, kk, jts:jte,p_bc_a06)+&
  chem(its:ite, kk, jts:jte,p_bc_a07)+chem(its:ite, kk, jts:jte,p_bc_a08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,2)=0.0
!bccw
tend_rate_dep(its:ite,jts:jte,1,3)=tend_rate_dep(its:ite,jts:jte,1,3)+&
 (chem(its:ite, kk, jts:jte,p_bc_cw01)+chem(its:ite, kk, jts:jte,p_bc_cw02)+&
  chem(its:ite, kk, jts:jte,p_bc_cw03)+chem(its:ite, kk, jts:jte,p_bc_cw04)+&
  chem(its:ite, kk, jts:jte,p_bc_cw05)+chem(its:ite, kk, jts:jte,p_bc_cw06)+&
  chem(its:ite, kk, jts:jte,p_bc_cw07)+chem(its:ite, kk, jts:jte,p_bc_cw08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,3)=0.0
!oc
tend_rate_dep(its:ite,jts:jte,1,4)=tend_rate_dep(its:ite,jts:jte,1,4)+&
 (chem(its:ite, kk, jts:jte,p_oc_a01)+chem(its:ite, kk, jts:jte,p_oc_a02)+&
  chem(its:ite, kk, jts:jte,p_oc_a03)+chem(its:ite, kk, jts:jte,p_oc_a04)+&
  chem(its:ite, kk, jts:jte,p_oc_a05)+chem(its:ite, kk, jts:jte,p_oc_a06)+&
  chem(its:ite, kk, jts:jte,p_oc_a07)+chem(its:ite, kk, jts:jte,p_oc_a08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,4)=0.0
!occw
tend_rate_dep(its:ite,jts:jte,1,5)=tend_rate_dep(its:ite,jts:jte,1,5)+&
 (chem(its:ite, kk, jts:jte,p_oc_cw01)+chem(its:ite, kk, jts:jte,p_oc_cw02)+&
  chem(its:ite, kk, jts:jte,p_oc_cw03)+chem(its:ite, kk, jts:jte,p_oc_cw04)+&
  chem(its:ite, kk, jts:jte,p_oc_cw05)+chem(its:ite, kk, jts:jte,p_oc_cw06)+&
  chem(its:ite, kk, jts:jte,p_oc_cw07)+chem(its:ite, kk, jts:jte,p_oc_cw08))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,5)=0.0
enddo
tend_rate_dep(its:ite,jts:jte,:,:)=(tend_rate_dep(its:ite,jts:jte,:,:)-snowchem_old(its:ite,jts:jte,:,:))/dt

endif

CASE (CBMZ_MOSAIC_4BIN,CBMZ_MOSAIC_4BIN_AQ,  &
      CBMZ_MOSAIC_DMS_4BIN,CBMZ_MOSAIC_DMS_4BIN_AQ,CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE)

if (iswitch.eq.0) then !before the process

tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)=0.0
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)=0.0
do isize = 1, nsize_aer(itype) 
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
enddo !isize 

snowchem_old(its:ite,jts:jte,:,:)=0.0
do kk=kts,kte
!dust
 do ibin=1,4
  if (ibin.eq.1) then
   bsize=1
   esize=2
   massfrac=1.0
  endif
  if (ibin.eq.2) then
   bsize=3
   esize=3
   massfrac=1.0
  endif
  if (ibin.eq.3) then
   bsize=4
   esize=4
   massfrac=0.5
  endif
  if (ibin.eq.4) then
   bsize=4
   esize=4
   massfrac=0.5
  endif
  do isize = bsize, esize
   l=lptr_dust_aer(isize,itype,1) !dry phase
   lcw=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   snowchem_old(its:ite,jts:jte,ibin,1)=snowchem_old(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
  enddo !isize
 enddo !ibin
!bc
snowchem_old(its:ite,jts:jte,1,2)=snowchem_old(its:ite,jts:jte,1,2)+&
 (chem(its:ite, kk, jts:jte,p_bc_a01)+chem(its:ite, kk, jts:jte,p_bc_a02)+&
  chem(its:ite, kk, jts:jte,p_bc_a03)+chem(its:ite, kk, jts:jte,p_bc_a04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,2)=0.0
!bccw
snowchem_old(its:ite,jts:jte,1,3)=snowchem_old(its:ite,jts:jte,1,3)+&
 (chem(its:ite, kk, jts:jte,p_bc_cw01)+chem(its:ite, kk, jts:jte,p_bc_cw02)+&
  chem(its:ite, kk, jts:jte,p_bc_cw03)+chem(its:ite, kk, jts:jte,p_bc_cw04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,3)=0.0
!oc
snowchem_old(its:ite,jts:jte,1,4)=snowchem_old(its:ite,jts:jte,1,4)+&
 (chem(its:ite, kk, jts:jte,p_oc_a01)+chem(its:ite, kk, jts:jte,p_oc_a02)+&
  chem(its:ite, kk, jts:jte,p_oc_a03)+chem(its:ite, kk, jts:jte,p_oc_a04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,4)=0.0
!occw
snowchem_old(its:ite,jts:jte,1,5)=snowchem_old(its:ite,jts:jte,1,5)+&
 (chem(its:ite, kk, jts:jte,p_oc_cw01)+chem(its:ite, kk, jts:jte,p_oc_cw02)+&
  chem(its:ite, kk, jts:jte,p_oc_cw03)+chem(its:ite, kk, jts:jte,p_oc_cw04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
snowchem_old(its:ite,jts:jte,2:4,5)=0.0

enddo !kk

endif !iswitch

if (iswitch.eq.1) then !after the process
do isize = 1, nsize_aer(itype)
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   endif
enddo !isize
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)- &
   tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dust)

   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tend_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)- &
   tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dustcw)

tend_rate_dep(its:ite,jts:jte,:,:)=0.0
do kk=kts,kte
!dust
 do ibin=1,4
  if (ibin.eq.1) then
   bsize=1
   esize=2
   massfrac=1.
  endif
  if (ibin.eq.2) then
   bsize=3
   esize=3
   massfrac=1.
  endif
  if (ibin.eq.3) then
   bsize=4
   esize=4
   massfrac=0.5
  endif
  if (ibin.eq.4) then
   bsize=4
   esize=4
   massfrac=0.5
  endif
  do isize = bsize, esize
   l=lptr_dust_aer(isize,itype,1) !dry phase
   lcw=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tend_rate_dep(its:ite,jts:jte,ibin,1)=tend_rate_dep(its:ite,jts:jte,ibin,1)+ &
   (chem(its:ite, kk, jts:jte,l)+chem(its:ite, kk, jts:jte,lcw))*massfrac* &
   rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
   endif
  enddo !isize
 enddo !ibin
!bc
tend_rate_dep(its:ite,jts:jte,1,2)=tend_rate_dep(its:ite,jts:jte,1,2)+&
 (chem(its:ite, kk, jts:jte,p_bc_a01)+chem(its:ite, kk, jts:jte,p_bc_a02)+&
  chem(its:ite, kk, jts:jte,p_bc_a03)+chem(its:ite, kk, jts:jte,p_bc_a04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,2)=0.0
!bccw
tend_rate_dep(its:ite,jts:jte,1,3)=tend_rate_dep(its:ite,jts:jte,1,3)+&
 (chem(its:ite, kk, jts:jte,p_bc_cw01)+chem(its:ite, kk, jts:jte,p_bc_cw02)+&
  chem(its:ite, kk, jts:jte,p_bc_cw03)+chem(its:ite, kk, jts:jte,p_bc_cw04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,3)=0.0
!oc
tend_rate_dep(its:ite,jts:jte,1,4)=tend_rate_dep(its:ite,jts:jte,1,4)+&
 (chem(its:ite, kk, jts:jte,p_oc_a01)+chem(its:ite, kk, jts:jte,p_oc_a02)+&
  chem(its:ite, kk, jts:jte,p_oc_a03)+chem(its:ite, kk, jts:jte,p_oc_a04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,4)=0.0
tend_rate_dep(its:ite,jts:jte,1,5)=tend_rate_dep(its:ite,jts:jte,1,5)+&
 (chem(its:ite, kk, jts:jte,p_oc_cw01)+chem(its:ite, kk, jts:jte,p_oc_cw02)+&
  chem(its:ite, kk, jts:jte,p_oc_cw03)+chem(its:ite, kk, jts:jte,p_oc_cw04))*&
  rho_phy(its:ite,kk,jts:jte)*1.e-9*dz8w(its:ite,kk,jts:jte)
tend_rate_dep(its:ite,jts:jte,2:4,5)=0.0
enddo
tend_rate_dep(its:ite,jts:jte,:,:)=(tend_rate_dep(its:ite,jts:jte,:,:)-snowchem_old(its:ite,jts:jte,:,:))/dt

endif


CASE DEFAULT
 !call wrf_debug(15,'aerosol budget NOT CALLING')
END SELECT aer_budget_select

END subroutine budget_calc

!--------------------------------------------------------
subroutine budget_calc2(iswitch,chem_opt,chem,rho_phy,   &
                        tmp_chem,tend_dryd,tend_acti,temporary,temporary_cw, &
                        deposition_dryrate,dz8w,dt,  &
                        drydep_flux, &
                       ims,ime,jms,jme,kms,kme,                  &
                       its,ite,jts,jte,kts,kte                                )

USE module_state_description
USE module_data_mosaic_asect, only: lptr_dust_aer,lptr_tr1r1_aer,lptr_tr1r2_aer,lptr_tr1r3_aer,lptr_tr1r4_aer,nsize_aer
use module_peg_util, only:  peg_message, peg_error_fatal

IMPLICIT NONE

integer,intent(in) :: iswitch,chem_opt

INTEGER,INTENT(IN) :: ims,ime, jms,jme, kms,kme,    &
                      its,ite, jts,jte, kts,kte

REAL,DIMENSION( ims:ime , kms:kme , jms:jme ),INTENT(IN   ) :: dz8w,rho_phy
REAL,      INTENT(IN   ) :: dt
REAL, DIMENSION( ims:ime,jms:jme,8),INTENT(IN) :: deposition_dryrate 

REAL, DIMENSION( ims:ime,kms:kme,jms:jme,num_chem),INTENT(IN) :: chem
REAL, DIMENSION( ims:ime,kms:kme,jms:jme,num_tmp_budget),INTENT(INOUT) :: tmp_chem
REAL, DIMENSION( ims:ime,kms:kme,jms:jme,num_tend_dryd),INTENT(INOUT) :: tend_dryd,tend_acti
REAL, DIMENSION( ims:ime,jms:jme,4,num_drydep_flux),INTENT(INOUT) :: drydep_flux 
REAL, DIMENSION( ims:ime,jms:jme,8,3),INTENT(INOUT) :: temporary,temporary_cw 

integer :: i,j,nv,n,iphase,isize,itype,l,lcw,p1st,ibin,bsize,esize
real :: massfrac,l2,l3,l4
integer lunerr

itype=1
p1st=param_first_scalar

aer_budget_select: SELECT CASE(chem_opt)
!tendency (as rain) needs to be accumulated since output is not as frequent as timestep (e.g., 1 h v.s. 150 s)

CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_KPP, &
      RACMSORG_AQ,RACMSORG_KPP,RACM_ESRLSORG_KPP,CBMZSORG,CBMZSORG_AQ)

if (iswitch.eq.0) then !before the process

tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
    chem(its:ite, kts:kte, jts:jte,p_dustj)+chem(its:ite, kts:kte, jts:jte,p_soila)
if (chem_opt.ne.RADM2SORG) then
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
    chem(its:ite, kts:kte, jts:jte,p_dustcwj)+chem(its:ite, kts:kte, jts:jte,p_soilcw)
endif

!dust
temporary(its:ite,  jts:jte,1,1)=0.0
temporary(its:ite,  jts:jte,2,1)=chem(its:ite, kts, jts:jte,p_dustj)
temporary(its:ite,  jts:jte,3,1)=chem(its:ite, kts, jts:jte,p_soila)
if (chem_opt.ne.RADM2SORG) then
temporary_cw(its:ite,  jts:jte,1,1)=0.0
temporary_cw(its:ite,  jts:jte,2,1)=chem(its:ite, kts, jts:jte,p_dustj)
temporary_cw(its:ite,  jts:jte,3,1)=chem(its:ite, kts, jts:jte,p_soilcw)
endif
!bc
temporary(its:ite,  jts:jte,1,2)=chem(its:ite, kts, jts:jte,p_eci)
temporary(its:ite,  jts:jte,2,2)=chem(its:ite, kts, jts:jte,p_ecj)
temporary(its:ite,  jts:jte,3,2)=0.0
if (chem_opt.ne.RADM2SORG) then
temporary_cw(its:ite,  jts:jte,1,2)=chem(its:ite, kts, jts:jte,p_eccwi)
temporary_cw(its:ite,  jts:jte,2,2)=chem(its:ite, kts, jts:jte,p_eccwj)
temporary_cw(its:ite,  jts:jte,3,2)=0.0
endif
!oc
temporary(its:ite,  jts:jte,1,3)=chem(its:ite, kts, jts:jte,p_orgpai)
temporary(its:ite,  jts:jte,2,3)=chem(its:ite, kts, jts:jte,p_orgpaj)
temporary(its:ite,  jts:jte,3,3)=0.0
if (chem_opt.ne.RADM2SORG) then
temporary_cw(its:ite,  jts:jte,1,3)=chem(its:ite, kts, jts:jte,p_orgpacwi)
temporary_cw(its:ite,  jts:jte,2,3)=chem(its:ite, kts, jts:jte,p_orgpacwj)
temporary_cw(its:ite,  jts:jte,3,3)=0.0
endif
    
endif

if (iswitch.eq.1) then !after the process

tend_dryd(its:ite, kts, jts:jte, p_tend_dryd_dust)= &
     tend_dryd(its:ite, kts, jts:jte, p_tend_dryd_dust)+ &
(-1.)*(deposition_dryrate(its:ite,jts:jte,2)/dz8w(its:ite,kts,jts:jte)*dt*temporary(its:ite,jts:jte,2,1)+ &
       deposition_dryrate(its:ite,jts:jte,3)/dz8w(its:ite,kts,jts:jte)*dt*temporary(its:ite,jts:jte,3,1))

if (chem_opt.ne.RADM2SORG) then
tend_dryd(its:ite, kts, jts:jte, p_tend_dryd_dustcw)= &
     tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)+ &
(-1.)*(0.1/dz8w(its:ite,kts,jts:jte)*dt*temporary_cw(its:ite,jts:jte,2,1)+ &
       0.1/dz8w(its:ite,kts,jts:jte)*dt*temporary_cw(its:ite,jts:jte,3,1))
endif

tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
     tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
    ((chem(its:ite, kts:kte, jts:jte,p_dustj)+chem(its:ite, kts:kte, jts:jte,p_soila))-&
      tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dust))

if (chem_opt.ne.RADM2SORG) then
tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
     tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
    ((chem(its:ite, kts:kte, jts:jte,p_dustj)+chem(its:ite, kts:kte, jts:jte,p_soilcw))-&
      tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw))
endif

!tendency for activation should be added by the loss due to dry depostion at the surface layer
tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)=&
     tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)+&
    (deposition_dryrate(its:ite,jts:jte,2)/dz8w(its:ite,kts,jts:jte)*dt*temporary(its:ite,jts:jte,2,1)+ &
     deposition_dryrate(its:ite,jts:jte,3)/dz8w(its:ite,kts,jts:jte)*dt*temporary(its:ite,jts:jte,3,1))

if (chem_opt.ne.RADM2SORG) then
tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)=&
     tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)+&
    (0.1/dz8w(its:ite,kts,jts:jte)*dt*temporary_cw(its:ite,jts:jte,2,1)+ &
     0.1/dz8w(its:ite,kts,jts:jte)*dt*temporary_cw(its:ite,jts:jte,3,1))
endif

!calculate the instaneous dry deposition fluxes(positive) (kg/m2/s)
if (chem_opt.ne.RADM2SORG) then
!dust
drydep_flux(its:ite,jts:jte,1,p_drydep_dust)= deposition_dryrate(its:ite,jts:jte,2)*temporary(its:ite,jts:jte,2,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9+0.1*temporary_cw(its:ite,jts:jte,2,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9
drydep_flux(its:ite,jts:jte,2,p_drydep_dust)= 0.0 
drydep_flux(its:ite,jts:jte,3,p_drydep_dust)= deposition_dryrate(its:ite,jts:jte,3)*temporary(its:ite,jts:jte,3,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9+0.1*temporary_cw(its:ite,jts:jte,3,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9 
drydep_flux(its:ite,jts:jte,4,p_drydep_dust)= 0.0 
!bc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcpho)=0.0
do n=1,2 
drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)+&
 deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,2)*rho_phy(its:ite,kts,jts:jte)*1.e-9+0.1*temporary_cw(its:ite,jts:jte,n,2)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo
!oc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocpho)=0.0
do n=1,2 
drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)+&
 deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,3)*rho_phy(its:ite,kts,jts:jte)*1.e-9+0.1*temporary_cw(its:ite,jts:jte,n,3)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo
else 
!dust
drydep_flux(its:ite,jts:jte,1,p_drydep_dust)= deposition_dryrate(its:ite,jts:jte,2)*temporary(its:ite,jts:jte,2,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9
drydep_flux(its:ite,jts:jte,2,p_drydep_dust)= 0.0
drydep_flux(its:ite,jts:jte,3,p_drydep_dust)= deposition_dryrate(its:ite,jts:jte,3)*temporary(its:ite,jts:jte,3,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9
drydep_flux(its:ite,jts:jte,4,p_drydep_dust)= 0.0
!bc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcpho)=0.0
do n=1,2
drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)+&
 deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,2)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo
!oc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocpho)=0.0
do n=1,2
drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)+&
 deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,3)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo

endif


endif

CASE (CBMZ_MOSAIC_8BIN, &
      CBMZ_MOSAIC_8BIN_AQ,CBMZ_MOSAIC_DMS_8BIN, &
      CBMZ_MOSAIC_DMS_8BIN_AQ,CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE)

if (iswitch.eq.0) then !before the process

tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)=0.0
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)=0.0
temporary(its:ite,jts:jte,:,:)=0.0
temporary_cw(its:ite,jts:jte,:,:)=0.0
do isize = 1, nsize_aer(itype)
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
enddo !isize 

!dust
!bc
temporary(its:ite,  jts:jte,1,2)=chem(its:ite,kts,jts:jte,p_bc_a01)
temporary(its:ite,  jts:jte,2,2)=chem(its:ite,kts,jts:jte,p_bc_a02)
temporary(its:ite,  jts:jte,3,2)=chem(its:ite,kts,jts:jte,p_bc_a03)
temporary(its:ite,  jts:jte,4,2)=chem(its:ite,kts,jts:jte,p_bc_a04)
temporary(its:ite,  jts:jte,5,2)=chem(its:ite,kts,jts:jte,p_bc_a05)
temporary(its:ite,  jts:jte,6,2)=chem(its:ite,kts,jts:jte,p_bc_a06)
temporary(its:ite,  jts:jte,7,2)=chem(its:ite,kts,jts:jte,p_bc_a07)
temporary(its:ite,  jts:jte,8,2)=chem(its:ite,kts,jts:jte,p_bc_a08)
temporary_cw(its:ite,  jts:jte,1,2)=chem(its:ite,kts,jts:jte,p_bc_cw01)
temporary_cw(its:ite,  jts:jte,2,2)=chem(its:ite,kts,jts:jte,p_bc_cw02)
temporary_cw(its:ite,  jts:jte,3,2)=chem(its:ite,kts,jts:jte,p_bc_cw03)
temporary_cw(its:ite,  jts:jte,4,2)=chem(its:ite,kts,jts:jte,p_bc_cw04)
temporary_cw(its:ite,  jts:jte,5,2)=chem(its:ite,kts,jts:jte,p_bc_cw05)
temporary_cw(its:ite,  jts:jte,6,2)=chem(its:ite,kts,jts:jte,p_bc_cw06)
temporary_cw(its:ite,  jts:jte,7,2)=chem(its:ite,kts,jts:jte,p_bc_cw07)
temporary_cw(its:ite,  jts:jte,8,2)=chem(its:ite,kts,jts:jte,p_bc_cw08)
!oc
temporary(its:ite,  jts:jte,1,3)=chem(its:ite,kts,jts:jte,p_oc_a01)
temporary(its:ite,  jts:jte,2,3)=chem(its:ite,kts,jts:jte,p_oc_a02)
temporary(its:ite,  jts:jte,3,3)=chem(its:ite,kts,jts:jte,p_oc_a03)
temporary(its:ite,  jts:jte,4,3)=chem(its:ite,kts,jts:jte,p_oc_a04)
temporary(its:ite,  jts:jte,5,3)=chem(its:ite,kts,jts:jte,p_oc_a05)
temporary(its:ite,  jts:jte,6,3)=chem(its:ite,kts,jts:jte,p_oc_a06)
temporary(its:ite,  jts:jte,7,3)=chem(its:ite,kts,jts:jte,p_oc_a07)
temporary(its:ite,  jts:jte,8,3)=chem(its:ite,kts,jts:jte,p_oc_a08)
temporary_cw(its:ite,  jts:jte,1,3)=chem(its:ite,kts,jts:jte,p_oc_cw01)
temporary_cw(its:ite,  jts:jte,2,3)=chem(its:ite,kts,jts:jte,p_oc_cw02)
temporary_cw(its:ite,  jts:jte,3,3)=chem(its:ite,kts,jts:jte,p_oc_cw03)
temporary_cw(its:ite,  jts:jte,4,3)=chem(its:ite,kts,jts:jte,p_oc_cw04)
temporary_cw(its:ite,  jts:jte,5,3)=chem(its:ite,kts,jts:jte,p_oc_cw05)
temporary_cw(its:ite,  jts:jte,6,3)=chem(its:ite,kts,jts:jte,p_oc_cw06)
temporary_cw(its:ite,  jts:jte,7,3)=chem(its:ite,kts,jts:jte,p_oc_cw07)
temporary_cw(its:ite,  jts:jte,8,3)=chem(its:ite,kts,jts:jte,p_oc_cw08)

endif

if (iswitch.eq.1) then !after the process

do isize = 1, nsize_aer(itype)
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)+ &
(-1.)*deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)

      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   l2=lptr_tr1r2_aer(isize,itype,1) !dry phase
   l3=lptr_tr1r3_aer(isize,itype,1) !dry phase
   l4=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
     if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)+ &
(-1.)*deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)

      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)+ &
(-1.)*0.1*temporary_cw(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   l2=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   l3=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   l4=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)+ &
(-1.)*0.1*temporary_cw(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
enddo !isize

!activation tendency
tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)- &
tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dust)

tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)- &
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)

!tendency for activation should be added by the loss due to dry depostion at the surface layer
do isize = 1, nsize_aer(itype)
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)+&
    deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
    dt/dz8w(its:ite,kts,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   l2=lptr_tr1r2_aer(isize,itype,1) !dry phase
   l3=lptr_tr1r3_aer(isize,itype,1) !dry phase
   l4=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
    if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)+&
    deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
    dt/dz8w(its:ite,kts,jts:jte)
   endif
 
   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)+&
    0.1*temporary_cw(its:ite,jts:jte,isize,1)*&
    dt/dz8w(its:ite,kts,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   l2=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   l3=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   l4=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
    if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)+&
    0.1*temporary_cw(its:ite,jts:jte,isize,1)* &
    dt/dz8w(its:ite,kts,jts:jte)
   endif
enddo !isize

!calculate the instaneous dry deposition fluxes(positive) (kg/m2/s)
!dust
drydep_flux(its:ite,jts:jte,:,p_drydep_dust)=0.0
do ibin=1,4
  if (ibin.eq.1) then
   bsize=1
   esize=4
   massfrac=1.
  endif
  if (ibin.eq.2) then
   bsize=5
   esize=6
   massfrac=1.
  endif
  if (ibin.eq.3) then
   bsize=7
   esize=7
   massfrac=1.0
  endif
  if (ibin.eq.4) then
   bsize=8
   esize=8
   massfrac=1.0
  endif
  do isize = bsize, esize
   l=lptr_dust_aer(isize,itype,1) !dry phase
   lcw=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)= &
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)+ &
   deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9 + &
   0.1*temporary_cw(its:ite,jts:jte,isize,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   l2=lptr_tr1r2_aer(isize,itype,1) !dry phase
   l3=lptr_tr1r3_aer(isize,itype,1) !dry phase
   l4=lptr_tr1r4_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
    if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)= &
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)+ &
   deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9 + &
   0.1*temporary_cw(its:ite,jts:jte,isize,1)*rho_phy(its:ite,kts,jts:jte)*1.e-9
   endif
  enddo !isize
enddo !ibin
!bc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcpho)=0.0
do n=1,8 
drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)+ &
deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,2)*rho_phy(its:ite,kts,jts:jte)*1.e-9+&
0.1*temporary_cw(its:ite,jts:jte,n,2)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo
!oc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocpho)=0.0
do n=1,8
drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)+ &
deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,3)*rho_phy(its:ite,kts,jts:jte)*1.e-9+&
0.1*temporary_cw(its:ite,jts:jte,n,3)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo

endif

CASE (CBMZ_MOSAIC_4BIN,CBMZ_MOSAIC_4BIN_AQ,  &
      CBMZ_MOSAIC_DMS_4BIN,CBMZ_MOSAIC_DMS_4BIN_AQ,CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE)

if (iswitch.eq.0) then !before the process

tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)=0.0
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)=0.0
temporary(its:ite,jts:jte,:,:)=0.0 
temporary_cw(its:ite,jts:jte,:,:)=0.0
do isize = 1, nsize_aer(itype)
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dust)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary(its:ite,jts:jte,isize,1)= &
   temporary(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)= &
   tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)+ &
   chem(its:ite,kts:kte,jts:jte,l)
   temporary_cw(its:ite,jts:jte,isize,1)= &
   temporary_cw(its:ite,jts:jte,isize,1)+ &
   chem(its:ite,kts,jts:jte,l)
   endif
enddo !isize 

!dust
!bc
temporary(its:ite,  jts:jte,1,2)=chem(its:ite,kts,jts:jte,p_bc_a01)
temporary(its:ite,  jts:jte,2,2)=chem(its:ite,kts,jts:jte,p_bc_a02)
temporary(its:ite,  jts:jte,3,2)=chem(its:ite,kts,jts:jte,p_bc_a03)
temporary(its:ite,  jts:jte,4,2)=chem(its:ite,kts,jts:jte,p_bc_a04)
temporary_cw(its:ite,  jts:jte,1,2)=chem(its:ite,kts,jts:jte,p_bc_cw01)
temporary_cw(its:ite,  jts:jte,2,2)=chem(its:ite,kts,jts:jte,p_bc_cw02)
temporary_cw(its:ite,  jts:jte,3,2)=chem(its:ite,kts,jts:jte,p_bc_cw03)
temporary_cw(its:ite,  jts:jte,4,2)=chem(its:ite,kts,jts:jte,p_bc_cw04)
!oc
temporary(its:ite,  jts:jte,1,3)=chem(its:ite,kts,jts:jte,p_oc_a01)
temporary(its:ite,  jts:jte,2,3)=chem(its:ite,kts,jts:jte,p_oc_a02)
temporary(its:ite,  jts:jte,3,3)=chem(its:ite,kts,jts:jte,p_oc_a03)
temporary(its:ite,  jts:jte,4,3)=chem(its:ite,kts,jts:jte,p_oc_a04)
temporary_cw(its:ite,  jts:jte,1,3)=chem(its:ite,kts,jts:jte,p_oc_cw01)
temporary_cw(its:ite,  jts:jte,2,3)=chem(its:ite,kts,jts:jte,p_oc_cw02)
temporary_cw(its:ite,  jts:jte,3,3)=chem(its:ite,kts,jts:jte,p_oc_cw03)
temporary_cw(its:ite,  jts:jte,4,3)=chem(its:ite,kts,jts:jte,p_oc_cw04)

endif

if (iswitch.eq.1) then !after the process

do isize = 1, nsize_aer(itype)
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)+ &
(-1.)*deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)

      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   l2=lptr_tr1r2_aer(isize,itype,1) !dry phase
   l3=lptr_tr1r3_aer(isize,itype,1) !dry phase
   l4=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
     if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dust)+ &
(-1.)*deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)

      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)+ &
(-1.)*0.1*temporary_cw(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   l2=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   l3=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   l4=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)= &
      tend_dryd(its:ite, kts, jts:jte,p_tend_dryd_dustcw)+ &
(-1.)*0.1*temporary_cw(its:ite,jts:jte,isize,1)* &
      dt/dz8w(its:ite,kts,jts:jte)
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
   l=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
      tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)+ &
      chem(its:ite,kts:kte,jts:jte,l)
   endif
enddo !isize

!activation tendency
tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)= &
tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dust)- &
tmp_chem(its:ite, kts:kte,jts:jte,p_tmp_dust)

tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)= &
tend_acti(its:ite, kts:kte, jts:jte,p_tend_acti_dustcw)- &
tmp_chem(its:ite, kts:kte, jts:jte,p_tmp_dustcw)

!tendency for activation should be added by the loss due to dry depostion at the surface layer
do isize = 1, nsize_aer(itype)
   l=lptr_dust_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)+&
    deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
    dt/dz8w(its:ite,kts,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   l2=lptr_tr1r2_aer(isize,itype,1) !dry phase
   l3=lptr_tr1r3_aer(isize,itype,1) !dry phase
   l4=lptr_tr1r4_aer(isize,itype,1) !dry phase
   if (l.ge.p1st) then
    if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dust)+&
    deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)* &
    dt/dz8w(its:ite,kts,jts:jte)
   endif

   l=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)+&
    0.1*temporary_cw(its:ite,jts:jte,isize,1)*&
    dt/dz8w(its:ite,kts,jts:jte)
   endif
   l=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   l2=lptr_tr1r2_aer(isize,itype,2) !cloud phase
   l3=lptr_tr1r3_aer(isize,itype,2) !cloud phase
   l4=lptr_tr1r4_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
    if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)=&
    tend_acti(its:ite, kts, jts:jte,p_tend_acti_dustcw)+&
    0.1*temporary_cw(its:ite,jts:jte,isize,1)* &
    dt/dz8w(its:ite,kts,jts:jte)
   endif
enddo !isize

!calculate the instaneous dry deposition fluxes(positive) (kg/m2/s)
!dust
drydep_flux(its:ite,jts:jte,:,p_drydep_dust)=0.0
do ibin=1,4
  if (ibin.eq.1) then
   bsize=1
   esize=2
   massfrac=1.
  endif
  if (ibin.eq.2) then
   bsize=3
   esize=3
   massfrac=1.
  endif
  if (ibin.eq.3) then
   bsize=4
   esize=4
   massfrac=0.5
  endif
  if (ibin.eq.4) then
   bsize=4
   esize=4
   massfrac=0.5
  endif
  do isize = bsize, esize
   l=lptr_dust_aer(isize,itype,1) !dry phase
   lcw=lptr_dust_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)= &
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)+ &
   deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)*massfrac*rho_phy(its:ite,kts,jts:jte)*1.e-9 + &
   0.1*temporary_cw(its:ite,jts:jte,isize,1)*massfrac*rho_phy(its:ite,kts,jts:jte)*1.e-9
   endif
   l=lptr_tr1r1_aer(isize,itype,1) !dry phase
   l2=lptr_tr1r2_aer(isize,itype,1) !dry phase
   l3=lptr_tr1r3_aer(isize,itype,1) !dry phase
   l4=lptr_tr1r4_aer(isize,itype,1) !dry phase
   lcw=lptr_tr1r1_aer(isize,itype,2) !cloud phase
   if (l.ge.p1st) then
    if (l2.le.p1st) call peg_error_fatal(lunerr, 'only one region??')
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)= &
   drydep_flux(its:ite,jts:jte,ibin,p_drydep_dust)+ &
   deposition_dryrate(its:ite,jts:jte,isize)*temporary(its:ite,jts:jte,isize,1)*massfrac*rho_phy(its:ite,kts,jts:jte)*1.e-9 + &
   0.1*temporary_cw(its:ite,jts:jte,isize,1)*massfrac*rho_phy(its:ite,kts,jts:jte)*1.e-9
   endif
  enddo !isize
enddo !ibin
!bc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_bcpho)=0.0
do n=1,4
drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_bcpho)+ &
deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,2)*rho_phy(its:ite,kts,jts:jte)*1.e-9+&
0.1*temporary_cw(its:ite,jts:jte,n,2)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo
!oc
drydep_flux(its:ite,jts:jte,1:4,p_drydep_ocpho)=0.0
do n=1,4
drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)=drydep_flux(its:ite,jts:jte,1,p_drydep_ocpho)+ &
deposition_dryrate(its:ite,jts:jte,n)*temporary(its:ite,jts:jte,n,3)*rho_phy(its:ite,kts,jts:jte)*1.e-9+&
0.1*temporary_cw(its:ite,jts:jte,n,3)*rho_phy(its:ite,kts,jts:jte)*1.e-9
enddo

endif


CASE DEFAULT
 !call wrf_debug(15,'aerosol budget NOT CALLING')
END SELECT aer_budget_select

END subroutine budget_calc2

