MODULE module_cam_mam_addemiss
!WRF:MODEL_LAYER:CHEMICS

!----------------------------------------------------------------------
! module_cam_mam_addemiss.F
!
! adapted from module_mosaic_addemiss.F in nov, 2010 by r.c.easter
!
!----------------------------------------------------------------------

#include "MODAL_AERO_CPP_DEFINES.h"

! rce 2005-feb-18 - one fix for indices of volumcen_sect, [now (isize,itype)]
! rce 2005-jan-14 - added subr mosaic_seasalt_emiss (and a call to it)
! rce 2004-dec-03 - many changes associated with the new aerosol "pointer"
!     variables in module_data_mosaic_asect


   private
   public :: cam_mam_addemiss


CONTAINS



!----------------------------------------------------------------------
   subroutine cam_mam_addemiss( id, dtstep, u10, v10, alt, dz8w, xland,    &
               config_flags, chem, slai, ust, smois, ivgtyp, isltyp,       &
               emis_ant,ebio_iso,ebio_olt,ebio_oli,rho_phy,                &
               dust_emiss_active, seasalt_emiss_active,sstemp,             &
               hfx,pblh,psfc,rainncv,raincv,tke_pbl,t2,                                       &
               sigma_dry,sigma_wet,sigma_tke,k_shape,c_scale,           &
            !czhao 
            u_phy,v_phy,g,emis_dust_factor,  &
            num_land,landf,num_soil,soilf,sh2o,snowc,xlai,   &
            watsat2d,sf_surface_physics,e_dust,e_seasalt,      &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
!
! adds emissions for cam-mam aerosol species
! (i.e., emissions tendencies over time dtstep are applied 
! to the aerosol mass and number mixing ratios)
!
! currently this routine
! - applies emissions held in the emis_ant array
! - does online sea-salt emissions using a modified version of the
!   mosaic sea-salt emissions routine
! - does online dust emissions using a modified version of the
!   mosaic dust emissions routine
!
! still to do
! - for emis_ant emissions, implement emitted sizes similar to those used in cam
! - implement the cam sea-salt and dust emissions routines
!

   USE module_configure, only:  grid_config_rec_type
   USE module_state_description
!  USE module_state_description, only:  num_chem, param_first_scalar,   &
!     emiss_inpt_default, emiss_inpt_pnnl_rs, emiss_inpt_pnnl_cm,num_emis_ant

   USE module_data_cam_mam_asect, only:  ai_phase,   &
         dens_so4_aer, dens_nh4_aer, dens_pom_aer,   &
         dens_bc_aer, dens_dust_aer, dens_seas_aer,   &
         lptr_so4_aer, lptr_nh4_aer, lptr_pom_aer,   &
         lptr_bc_aer, lptr_dust_aer, lptr_seas_aer,   &
         maxd_asize, maxd_atype, nsize_aer, ntype_aer, numptr_aer

   USE modal_aero_data, only:  &
         modeptr_accum, modeptr_aitken, modeptr_coarse,  &
         modeptr_coardust, modeptr_finedust, modeptr_pcarbon, &
         modeptr_fineseas,  modeptr_coarseas

   USE module_cam_mam_init, only:  &
         pom_emit_1p4_factor, so4_emit_1p2_factor

   USE module_data_sorgam, only:  &
         dgvem_i, dgvem_j, dgvem_c, sgem_i, sgem_j, sgem_c

  USE module_data_gocart_dust, only: porosity


   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   INTEGER, INTENT(IN) ::    dust_emiss_active, seasalt_emiss_active

  INTEGER,      INTENT(IN   ) :: num_land,num_soil

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme )         ,                         &
          INTENT(IN   ) ::   u10, v10, xland, slai, ust
   INTEGER,  DIMENSION( ims:ime , jms:jme )      ,                         &
          INTENT(IN   ) ::   ivgtyp, isltyp

   REAL, DIMENSION( ims:ime, jms:jme),               &
         INTENT(IN) ::                                   sstemp

  real, dimension (ims:ime, jms:jme ) ,                               &
        intent(in) ::                                                  &
        hfx,pblh,psfc,rainncv,raincv,t2

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
         INTENT(IN ) :: tke_pbl

   REAL,DIMENSION(ims:ime,jms:jme),INTENT(INOUT) :: sigma_dry,sigma_wet,sigma_tke
   REAL,DIMENSION(ims:ime,jms:jme),INTENT(INOUT) :: k_shape,c_scale

   REAL, DIMENSION( ims:ime, jms:jme), INTENT(INOUT) ::  e_seasalt,e_dust

! trace species mixing ratios (aerosol mass = ug/kg-air; number = #/kg-air)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem
!
! aerosol emissions arrays ((ug/m3)*m/s)
!
!   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                          &
   REAL, DIMENSION( ims:ime,   1:config_flags%kemit, jms:jme,num_emis_ant ),            &
         INTENT(IN ) ::                                                    &
                         emis_ant
!jdf
   REAL,  DIMENSION( ims:ime , jms:jme )      ,                            &
          INTENT(IN   ) ::   ebio_iso, ebio_olt, ebio_oli
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::   rho_phy
!jdf

! 1/(dry air density) and layer thickness (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

   !czhao add for gocart dust
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::                                                 &
                           u_phy,v_phy
   REAL, INTENT(IN   ) ::  g
   REAL, INTENT(IN) :: emis_dust_factor

  REAL, DIMENSION(ims:ime,num_land,jms:jme),               &
         INTENT(IN) ::  landf
  REAL, DIMENSION(ims:ime,num_soil,jms:jme),               &
         INTENT(IN) ::  soilf 
  REAL, DIMENSION(ims:ime,jms:jme),               &
         INTENT(IN) :: snowc,xlai
  REAL, DIMENSION(ims:ime,jms:jme),               &
         INTENT(IN) :: watsat2d
  INTEGER, INTENT(IN) :: sf_surface_physics

  REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,      &
          INTENT(IN   ) ::   smois,sh2o

! local variables

!    anth_emiss_size_opt = 1 --> sizes follow those used in cam5
!    anth_emiss_size_opt = 2 --> sizes follow those used in sorgam
   integer, parameter :: anth_emiss_size_opt = 1   ! user can change this value
!    dust_emiss_size_opt = 1 --> sizes follow those used in cam5
!    dust_emiss_size_opt = 2 --> sizes follow those used in sorgam
   integer, parameter :: dust_emiss_size_opt = 1   ! user can change this value
!    seas_emiss_size_opt = 1 --> sizes follow those used in cam5
!    seas_emiss_size_opt = 2 --> sizes follow those used in sorgam
   integer, parameter :: seas_emiss_size_opt = 1   ! user can change this value

   integer :: i, iphase, isize, itype
   integer :: itype_accum, itype_aitken, itype_pcarbon, &
              itype_finedust, itype_coardust
   integer :: j, k
   integer :: l, lemit, lemit_type, lnumb
   integer :: m, n
   integer :: p1st

   real, parameter :: pi = 3.1415926536
   real :: dp_anth_emiss_aitken, dp_anth_emiss_accum, dp_anth_emiss_coarse
   real :: dlo_seas_emiss( maxd_asize, maxd_atype ), &
           dhi_seas_emiss( maxd_asize, maxd_atype )
   real :: fact_mass, factbb_mass, factbb_numb
   real :: tmp_dens, tmp_dp, tmp_vol1p
   real :: total_soag, fact_soag

   character(len=100) :: msg


	iphase = ai_phase


!
! code applies emissions in the emis_ant array,
! then does "online" sea-salt and dust emissions
!
! currently the code only works with
! config_flags%emiss_inpt_opt = emiss_inpt_default OR emiss_inpt_pnnl_rs
!

	emiss_inpt_select_1: SELECT CASE( config_flags%emiss_inpt_opt )

!jdf
	    CASE( emiss_inpt_default, emiss_inpt_pnnl_rs, emiss_inpt_pnnl_mam )
		write(*,*) 'cam_mam_addemiss DOING emiss - emiss_inpt_opt =', &
			config_flags%emiss_inpt_opt
!jdf
	    CASE DEFAULT
		write(*,*) 'cam_mam_addemiss SKIPPING emiss - emiss_inpt_opt =', &
			config_flags%emiss_inpt_opt
		return

	END SELECT emiss_inpt_select_1


	p1st = param_first_scalar

	itype_accum    = modeptr_accum
	itype_aitken   = modeptr_aitken
#if ( defined MODAL_AERO_3MODE )
	itype_pcarbon  = modeptr_accum
	itype_finedust = modeptr_accum
	itype_coardust = modeptr_coarse
#elif ( defined MODAL_AERO_7MODE )
    if ( (config_flags%chem_opt /= CBMZ_CAM_MAM7_NOAQ) .and. &
         (config_flags%chem_opt /= CBMZ_CAM_MAM7_AQ  ) ) then
       call wrf_error_fatal( 'emission package for MAM7 is not implemented yet' )
    end if
#else
	itype_pcarbon  = modeptr_pcarbon
	itype_finedust = modeptr_finedust
	itype_coardust = modeptr_coardust
#endif

! set emitted sizes
	if (anth_emiss_size_opt == 1) then
	    ! volume-mean diameter for emitted particles (cm)
	    dp_anth_emiss_aitken = 0.0504e-4
	    dp_anth_emiss_accum  = 0.134e-4
	    dp_anth_emiss_coarse = 2.06e-4
	else if (anth_emiss_size_opt == 2) then
	    ! volume-mean diameter for emitted particles (cm)
	    dp_anth_emiss_aitken = 1.0e2 * dgvem_i / exp( 1.5 * (log(sgem_i)**2) )
	    dp_anth_emiss_accum  = 1.0e2 * dgvem_j / exp( 1.5 * (log(sgem_j)**2) )
	    dp_anth_emiss_coarse = 1.0e2 * dgvem_c / exp( 1.5 * (log(sgem_c)**2) )
	else
	    write(msg,'(2a,i7)') 'subr cam_mam_addemiss', &
		' - illegal anth_emiss_size_opt = ', anth_emiss_size_opt
	    call wrf_error_fatal( msg )
	end if


!jdf add case statements here to switch between 2 emission file options
	emiss_inpt_select_2: SELECT CASE( config_flags%emiss_inpt_opt )

	CASE( emiss_inpt_pnnl_rs )
!jdf
	do 1900 lemit_type = 1, 11

	iphase = 1
	isize = 1

! following if/then/else blocks determine
!    which chem species will receive a particular emis_ant(...,p_e_...)
!    additional mass factor for consistency with cam_mam assumptions
!    particle size and density for converting mass to number emissions
!
! oc (org), ec, and pm25 - both the "i" and "j" go to the same
!    cam_mam mode
! pm25 and pm_10 - go to the cam_mam dust species
! no3 - ignore cam_mam does not treat no3
!
	if (lemit_type == 1) then
	    lemit = p_e_so4i
	    itype = itype_aitken
	    l = lptr_so4_aer(isize,itype,iphase)
	    tmp_dens = dens_so4_aer
	    factbb_mass = so4_emit_1p2_factor
	    tmp_dp = dp_anth_emiss_aitken
	else if (lemit_type == 2) then
	    lemit = p_e_so4j
	    itype = itype_accum
	    l = lptr_so4_aer(isize,itype,iphase)
	    tmp_dens = dens_so4_aer
! so4_emit_1p2_factor is 115/96 if so4 (in 3-mode) is being treated 
! as ammonium bisulfate, or 1.0 if treated as pure sulfate
	    factbb_mass = so4_emit_1p2_factor
	    tmp_dp = dp_anth_emiss_accum

	else if (lemit_type == 3) then
	    lemit = p_e_no3i
	    cycle   ! cam_mam does not treat no3
	else if (lemit_type == 4) then
	    lemit = p_e_no3j
	    cycle   ! cam_mam does not treat no3

	else if (lemit_type == 5) then
	    lemit = p_e_orgi
	    itype = itype_pcarbon
	    l = lptr_pom_aer(isize,itype,iphase)
	    tmp_dens = dens_pom_aer
! pom_emit_1p4_factor is 1.4 if pom emissions are carbon-only emissions,
! or 1.0 if they are organic matter emissions
	    factbb_mass = pom_emit_1p4_factor
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 6) then
	    lemit = p_e_orgj
	    itype = itype_pcarbon
	    l = lptr_pom_aer(isize,itype,iphase)
	    tmp_dens = dens_pom_aer
	    factbb_mass = pom_emit_1p4_factor
	    tmp_dp = dp_anth_emiss_accum

	else if (lemit_type == 7) then
	    lemit = p_e_eci
	    itype = itype_pcarbon
	    l = lptr_bc_aer(isize,itype,iphase)
	    tmp_dens = dens_bc_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 8) then
	    lemit = p_e_ecj
	    itype = itype_pcarbon
	    l = lptr_bc_aer(isize,itype,iphase)
	    tmp_dens = dens_bc_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum

	else if (lemit_type == 9) then
	    lemit = p_e_pm25i
	    itype = itype_finedust
	    l = lptr_dust_aer(isize,itype,iphase)
	    tmp_dens = dens_dust_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 10) then
	    lemit = p_e_pm25j
	    itype = itype_finedust
	    l = lptr_dust_aer(isize,itype,iphase)
	    tmp_dens = dens_dust_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_accum
	else if (lemit_type == 11) then
	    lemit = p_e_pm_10
	    itype = itype_coardust
	    l = lptr_dust_aer(isize,itype,iphase)
	    tmp_dens = dens_dust_aer
	    factbb_mass = 1.0
	    tmp_dp = dp_anth_emiss_coarse
	else
	    cycle
	end if

	if ((l < p1st) .or. (l > num_chem)) cycle

	lnumb = numptr_aer(isize,itype,iphase)
	if ((lnumb < p1st) .or. (lnumb > num_chem)) lnumb = -999888777

	tmp_vol1p = (pi/6.0)*(tmp_dp**3)   ! mean volume of emitted particles (cm3)
! factbb_numb convert mass emissions (ug/m2/s) to number emissions (#/m2/s)
	factbb_numb = 1.0e-6/(tmp_dens*tmp_vol1p)

	do j = jts, jte
	do k =   1, min(config_flags%kemit,kte)
	do i = its, ite

! emissions fluxes are ug/m2/s
! mult by tmp_fact_mass to get mixing ratio change (ug/kg) over dtstep
	fact_mass = (dtstep/dz8w(i,k,j))*alt(i,k,j)*factbb_mass

	chem(i,k,j,l) = chem(i,k,j,l) + emis_ant(i,k,j,lemit)*fact_mass

	if (lnumb > 0) then
	chem(i,k,j,lnumb) = chem(i,k,j,lnumb) + &
		emis_ant(i,k,j,lemit)*fact_mass*factbb_numb
	end if

	end do ! i
	end do ! k
	end do ! j

1900	continue
!
! now do emissions for e_soag
! add biogenic isoprene for soag_isoprene
! no monoterpenes in CBMZ at present for soag_terpene
! lump anthropogenic ald, hc3, hc5, hc8, ket, oli, olt, and ora2 which is par in CBMZ for soag_bigalk
! lump biogenic olet and olei for soag_bigene
! lump anthropogenic tol and xyl for soag_toluene
!
        do j = jts, jte
        do k =   1, 1
        do i = its, ite
          fact_soag = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
          total_soag = (emis_ant(i,k,j,p_e_tol) * (92.14*0.15/12.0))              + &
                       (emis_ant(i,k,j,p_e_xyl) * (106.16*0.15/12.0))             + &
                       (emis_ant(i,k,j,p_e_hc5) * (77.6*0.05/12.0))               + &
                       (emis_ant(i,k,j,p_e_olt) * (72.34*0.05/12.0))              + &
                       (emis_ant(i,k,j,p_e_oli) * (75.78*0.05/12.0))              + &
                       (ebio_iso(i,j) * (68.11*0.04/12.0))
!                      0.4*emis_ant(i,k,j,p_e_ald) + 2.9*emis_ant(i,k,j,p_e_hc3)  + &
!                      4.8*emis_ant(i,k,j,p_e_hc5) + 7.9*emis_ant(i,k,j,p_e_hc8)  + &
!                      0.9*emis_ant(i,k,j,p_e_ket) + 2.8*emis_ant(i,k,j,p_e_oli)  + &
!                      1.8*emis_ant(i,k,j,p_e_olt) + 1.0*emis_ant(i,k,j,p_e_ora2) + &
!                      (ebio_iso(i,j) * (68.11*0.04/12.0))                        + &
!                      ebio_olt(i,j)                                              + &
!                      ebio_oli(i,j) 
          chem(i,k,j,p_soag) = chem(i,k,j,p_soag) + total_soag*fact_soag
        end do ! i
        end do ! k
        end do ! j
        do j = jts, jte
        do k =   2, min(config_flags%kemit,kte)
        do i = its, ite
          fact_soag = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
          total_soag = (emis_ant(i,k,j,p_e_tol) * (92.14*0.15/12.0))              + &
                       (emis_ant(i,k,j,p_e_xyl) * (106.16*0.15/12.0))             + &
                       (emis_ant(i,k,j,p_e_hc5) * (77.6*0.05/12.0))               + &
                       (emis_ant(i,k,j,p_e_olt) * (72.34*0.05/12.0))              + &
                       (emis_ant(i,k,j,p_e_oli) * (75.78*0.05/12.0))
!                      emis_ant(i,k,j,p_e_xyl)                                    + &
!                      0.4*emis_ant(i,k,j,p_e_ald) + 2.9*emis_ant(i,k,j,p_e_hc3)  + &
!                      4.8*emis_ant(i,k,j,p_e_hc5) + 7.9*emis_ant(i,k,j,p_e_hc8)  + &
!                      0.9*emis_ant(i,k,j,p_e_ket) + 2.8*emis_ant(i,k,j,p_e_oli)  + &
!                      1.8*emis_ant(i,k,j,p_e_olt) + 1.0*emis_ant(i,k,j,p_e_ora2)
          chem(i,k,j,p_soag) = chem(i,k,j,p_soag) + total_soag*fact_soag
        end do ! i
        end do ! k
        end do ! j

	CASE( emiss_inpt_pnnl_mam )

	do 1910 lemit_type = 1, 14

	iphase = 1
        isize = 1
	if (lemit_type == 1) then
	    lemit = p_e_so4i
	    itype = itype_aitken
	    l = lptr_so4_aer(isize,itype,iphase)
	    factbb_mass = so4_emit_1p2_factor
	else if (lemit_type == 2) then
	    lemit = p_e_so4j
	    itype = itype_accum
	    l = lptr_so4_aer(isize,itype,iphase)
! so4_emit_1p2_factor is 115/96 if so4 (in 3-mode) is being treated 
! as ammonium bisulfate, or 1.0 if treated as pure sulfate
	    factbb_mass = so4_emit_1p2_factor
	else if (lemit_type == 3) then
	    lemit = p_e_orgj
	    itype = itype_pcarbon
	    l = lptr_pom_aer(isize,itype,iphase)
	    factbb_mass = pom_emit_1p4_factor
	else if (lemit_type == 4) then
	    lemit = p_e_ecj
	    itype = itype_pcarbon
	    l = lptr_bc_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 5) then
	    lemit = p_e_dust_a1
	    itype = itype_finedust
	    l = lptr_dust_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 6) then
	    lemit = p_e_dust_a3
	    itype = itype_coardust
	    l = lptr_dust_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 7) then
	    lemit = p_e_ncl_a1
	    itype = itype_accum
	    l = lptr_seas_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 8) then
	    lemit = p_e_ncl_a2
	    itype = itype_aitken
	    l = lptr_seas_aer(isize,itype,iphase)
	    factbb_mass = 1.0
	else if (lemit_type == 9) then
	    lemit = p_e_ncl_a3
	    itype = itype_coardust
	    l = lptr_seas_aer(isize,itype,iphase)
	    factbb_mass = 1.0
!PMA
!Currently my pre-processing script put num_a1 and num_a2 emissions from 
!sea salt and dust into so4j_num and so4i_num 
!num_a3 emissions are solely from dust and sea salt
!

        else if (lemit_type == 10) then
            lemit = p_e_so4i_num
            itype = itype_aitken
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 11) then
            lemit = p_e_so4j_num
            itype = itype_accum
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 12) then
            lemit = p_e_orgj_num
            itype = itype_accum
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 13) then
            lemit = p_e_ecj_num
            itype = itype_accum
            l = numptr_aer(isize,itype,iphase)
        else if (lemit_type == 14) then
            lemit = p_e_num_a3
            itype = itype_coardust
            l = numptr_aer(isize,itype,iphase)
        else
	    cycle
	end if

	if ((l < p1st) .or. (l > num_chem)) cycle

	do j = jts, jte
	do k =   1, min(config_flags%kemit,kte)
	do i = its, ite

! emissions fluxes are ug/m2/s
! mult by tmp_fact_mass to get mixing ratio change (ug/kg) over dtstep
	fact_mass = (dtstep/dz8w(i,k,j))*alt(i,k,j)*factbb_mass
!PMA 
        factbb_numb = (dtstep/dz8w(i,k,j))*alt(i,k,j)
        if (lemit_type < 10) then 
	chem(i,k,j,l) = chem(i,k,j,l) + emis_ant(i,k,j,lemit)*fact_mass
        else
        chem(i,k,j,l) = chem(i,k,j,l) + emis_ant(i,k,j,lemit)*factbb_numb
        end if

	end do ! i
	end do ! k
	end do ! j

1910	continue
!
! units for e_soag* assumed to be ug/m2/s for now
!
	do j = jts, jte
	do k =   1, min(config_flags%kemit,kte)
	do i = its, ite
	  fact_soag = (dtstep/dz8w(i,k,j))*alt(i,k,j)/1000.0
          total_soag = emis_ant(i,k,j,p_e_soag_bigalk)   + &
                       emis_ant(i,k,j,p_e_soag_bigene)   + &
                       emis_ant(i,k,j,p_e_soag_isoprene) + &
                       emis_ant(i,k,j,p_e_soag_terpene)  + &
                       emis_ant(i,k,j,p_e_soag_toluene)
	  chem(i,k,j,p_soag) = chem(i,k,j,p_soag) + total_soag*fact_soag
	end do ! i
	end do ! k
	end do ! j


	END SELECT emiss_inpt_select_2
!jdf end of add case statements here to switch between 2 emission file options

!
! seasalt emissions
!

! set emitted sizes
	dlo_seas_emiss(:,:) = 0.0
	dhi_seas_emiss(:,:) = 0.0

	if (seas_emiss_size_opt == 1) then
#if ( defined MODAL_AERO_3MODE )
! cut sizes are 0.02, 0.08, 1.0, 10.0 um
	    dlo_seas_emiss(1,modeptr_aitken  ) = 0.02e-4
	    dhi_seas_emiss(1,modeptr_aitken  ) = 0.08e-4
	    dlo_seas_emiss(1,modeptr_accum   ) = 0.08e-4
	    dhi_seas_emiss(1,modeptr_accum   ) = 1.00e-4
	    dlo_seas_emiss(1,modeptr_coarse  ) = 1.00e-4
	    dhi_seas_emiss(1,modeptr_coarse  ) = 10.0e-4
#else
! cut sizes are 0.02, 0.08, 0.3, 1.0, 10.0 um
	    dlo_seas_emiss(1,modeptr_aitken  ) = 0.02e-4
	    dhi_seas_emiss(1,modeptr_aitken  ) = 0.08e-4
	    dlo_seas_emiss(1,modeptr_accum   ) = 0.08e-4
	    dhi_seas_emiss(1,modeptr_accum   ) = 0.30e-4
	    dlo_seas_emiss(1,modeptr_fineseas) = 0.30e-4
	    dhi_seas_emiss(1,modeptr_fineseas) = 1.00e-4
	    dlo_seas_emiss(1,modeptr_coarseas) = 1.00e-4
	    dhi_seas_emiss(1,modeptr_coarseas) = 10.0e-4
#endif

	else if (seas_emiss_size_opt == 2) then
#if ( defined MODAL_AERO_3MODE )
! cut sizes are 0.1, 1.25, 10.0 um and no aitken
	    dlo_seas_emiss(1,modeptr_accum   ) = 0.10e-4
	    dhi_seas_emiss(1,modeptr_accum   ) = 1.25e-4
	    dlo_seas_emiss(1,modeptr_coarse  ) = 1.25e-4
	    dhi_seas_emiss(1,modeptr_coarse  ) = 10.0e-4
#else
! cut sizes are 0.1, 0.3, 1.25, 10.0 um and no aitken
	    dlo_seas_emiss(1,modeptr_accum   ) = 0.10e-4
	    dhi_seas_emiss(1,modeptr_accum   ) = 0.30e-4
	    dlo_seas_emiss(1,modeptr_fineseas) = 0.30e-4
	    dhi_seas_emiss(1,modeptr_fineseas) = 1.25e-4
	    dlo_seas_emiss(1,modeptr_coarseas) = 1.25e-4
	    dhi_seas_emiss(1,modeptr_coarseas) = 10.0e-4
#endif

	else
	    write(msg,'(2a,i7)') 'subr cam_mam_addemiss', &
		' - illegal seas_emiss_size_opt = ', seas_emiss_size_opt
	    call wrf_error_fatal( msg )
	end if

! now do the sea-salt emissions
	if (config_flags%seas_opt == 2) then
	    call cam_mam_mosaic_seasalt_emiss(                             &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dlo_seas_emiss, dhi_seas_emiss,                             &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
	end if

        if (config_flags%seas_opt == 4) then
            call cam_mam_cam5_seasalt_emiss(                             &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dlo_seas_emiss, dhi_seas_emiss,e_seasalt,sstemp,            &
               hfx,pblh,psfc,rainncv,raincv,tke_pbl,t2,rho_phy,            &
               sigma_dry,sigma_wet,sigma_tke,k_shape,c_scale,         &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
        end if

!
! dust emissions
!
  !soil porosity calculation if necessary
  if (sf_surface_physics==CLMSCHEME) then
        do j=jts,jte
        do i=its,ite
          porosity(isltyp(i,j))=watsat2d(i,j)
        enddo
        enddo
  endif

! jdf: preliminary version that has not been made generic for situation
	if (config_flags%dust_opt == 2) then
            call cam_mam_mosaic_dust_emiss(                                &
               slai, ust, smois, ivgtyp, isltyp,                           &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dust_emiss_size_opt,                                        &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
	end if

        if (config_flags%dust_opt == 12 ) then
        !for now, the option is only tested with Noah land surface model
         call cam_mam_dust_DEADemis (                                    &
         dtstep,config_flags%num_soil_layers,alt,u_phy,v_phy,  &
         chem,rho_phy,dz8w,smois,ust,u10,v10,                  &
         ivgtyp,isltyp,xland,g,e_dust,       &
         num_land,landf,num_soil,soilf,sh2o,snowc,xlai,            &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
        endif

	return


   END subroutine cam_mam_addemiss



!----------------------------------------------------------------------
   subroutine cam_mam_mosaic_seasalt_emiss(                                &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dlo_seas_emiss, dhi_seas_emiss,                             &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
!
! adds seasalt emissions for mosaic aerosol species
! (i.e., seasalt emissions tendencies over time dtstep are applied 
! to the aerosol mixing ratios)
!

   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_cam_mam_asect, only:  ai_phase,   &
         lptr_seas_aer,   &
         maxd_asize, maxd_atype, nsize_aer, ntype_aer, numptr_aer
!  USE module_mosaic_addemiss, only:  seasalt_emitfactors_1bin

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland

! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

! emissions cut sizes (cm)
   REAL,  DIMENSION( maxd_asize, maxd_atype ), &
          INTENT(IN   ) ::  dlo_seas_emiss, dhi_seas_emiss

! local variables
	integer i, j, k, l, l_seas, n
	integer iphase, itype
	integer p1st

	real dum, dumdlo, dumdhi, dumoceanfrac, dumspd10
	real factaa, factbb

	real :: ssemfact_numb( maxd_asize, maxd_atype )
	real :: ssemfact_mass( maxd_asize, maxd_atype )

    write(*,*) 'in subr cam_mam_mosaic_seasalt_emiss'
    p1st = PARAM_FIRST_SCALAR

!   for now just do itype=1
	iphase = ai_phase

!   compute emissions factors for each size bin
!   (limit emissions to dp > 0.1 micrometer)
	ssemfact_mass(:,:) = 0.0
	ssemfact_numb(:,:) = 0.0
	do itype = 1, ntype_aer
	do n = 1, nsize_aer(itype)
	    dumdlo = max( dlo_seas_emiss(n,itype), 0.1e-4 )
	    dumdhi = max( dhi_seas_emiss(n,itype), 0.1e-4 )
	    if (dumdlo >= dumdhi) cycle
	    call seasalt_emitfactors_1bin( 1, dumdlo, dumdhi,   &
		ssemfact_numb(n,itype), dum, ssemfact_mass(n,itype) )

!   convert mass emissions factor from (g/m2/s) to (ug/m2/s)
	    ssemfact_mass(n,itype) = ssemfact_mass(n,itype)*1.0e6
	end do
	end do


!   loop over i,j and apply seasalt emissions
	k = kts
	do 1830 j = jts, jte
	do 1820 i = its, ite

    !Skip this point if over land. xland=1 for land and 2 for water.
    !Also, there is no way to differentiate fresh from salt water.
    !Currently, this assumes all water is salty.
	if( xland(i,j) < 1.5 ) cycle

    !wig: As far as I can tell, only real.exe knows the fractional breakdown
    !     of land use. So, in wrf.exe, dumoceanfrac will always be 1.
	dumoceanfrac = 1. !fraction of grid i,j that is salt water
	dumspd10 = dumoceanfrac* &
         ( (u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5*3.41) )

!   factaa is (s*m2/kg-air)
!   factaa*ssemfact_mass(n) is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
!   factaa*ssemfact_numb(n) is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
	factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)

	factbb = factaa * dumspd10

	do 1815 itype = 1, ntype_aer
	do 1810 n = 1, nsize_aer(itype)
	    if (ssemfact_mass(n,itype) <= 0.0) cycle

!   only apply emissions if bin has both na and cl species
	    l_seas = lptr_seas_aer(n,itype,iphase)
	    if (l_seas < p1st) cycle

	    chem(i,k,j,l_seas) = chem(i,k,j,l_seas) +   &
			factbb * ssemfact_mass(n,itype)

	    l = numptr_aer(n,itype,iphase)
	    if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
			factbb * ssemfact_numb(n,itype)

1810	continue
1815	continue

1820	continue
1830	continue

	return

   END subroutine cam_mam_mosaic_seasalt_emiss


   subroutine cam_mam_cam5_seasalt_emiss(                                &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dlo_seas_emiss, dhi_seas_emiss,e_seasalt,sstemp,            &
               hfx,pblh,psfc,rainncv,raincv,tke_pbl,t2,rho,                &
               sigma_dry2d,sigma_wet2d,sigma_tke2d,k_2d,c_2d,         &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )

!----------------------------------------------------------------------
! Routine to calculate seasalt emissions (original code is from CAM5 bulk_aero/progseasalts_intr.F90) 
! Chun.Zhao@pnnl.gov; 18-Dec-2012
!
! Derives from Martensson et al. (2003) (JGR-108, 4297,doi:10.1029/2002JD002263)
! valid from 20nm to ~2500nm dry diameter (based on lab experiment with artificial sea water)
!
! currently we recommend that it is combined with
! the parameterisation by Monahan et al. (1986) for
! dry diameters > 2-3 um even if it lacks
! temperature dependence (despite that Bowyer et
! al. (1990) found a similar dependency in the tank
! from Monahan et al. (1986))
!
!-----------------------------------------------------------------------

   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar
   USE module_data_cam_mam_asect, only:  ai_phase,   &
         lptr_seas_aer,   &
         maxd_asize, maxd_atype, nsize_aer, ntype_aer, numptr_aer

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland

  real, dimension (ims:ime, jms:jme ) ,                               &
        intent(in) ::                                                  &
        hfx,pblh,psfc,rainncv,raincv,t2

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
         INTENT(IN ) :: tke_pbl,rho

! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w


! emissions cut sizes (cm)
   REAL,  DIMENSION( maxd_asize, maxd_atype ), &
          INTENT(IN   ) ::  dlo_seas_emiss, dhi_seas_emiss

   REAL, DIMENSION( ims:ime, jms:jme),               &
         INTENT(IN) ::                                   sstemp

   REAL, DIMENSION( ims:ime, jms:jme), INTENT(INOUT) ::  e_seasalt 

! local variables

! use Ekman's ss
! here is used for Ekman's ss
   integer, parameter :: sections = 31 ! number of sections used to calculate fluxes 
! only use up to ~20um
   real, dimension(sections) :: Dg = (/  &
                                0.0020e-5, 0.0025e-5, 0.0032e-5,  &
                                0.0040e-5, 0.0051e-5, 0.0065e-5,  &
                                0.0082e-5, 0.0104e-5, 0.0132e-5,  &
                                0.0167e-5, 0.0211e-5, 0.0267e-5,  &
                                0.0338e-5, 0.0428e-5, 0.0541e-5,  &
                                0.0685e-5, 0.0867e-5, 0.1098e-5,  &
                                0.1389e-5, 0.1759e-5, 0.2226e-5,  &
                                0.2818e-5, 0.3571e-5, 0.4526e-5,  &
                                0.5735e-5, 0.7267e-5, 0.9208e-5,  &
                                1.1668e-5, 1.4786e-5, 1.8736e-5,  &
                                2.3742e-5 /)

   real :: fi(sections)
   real, dimension(sections) :: bm, rdry, rm
   real, dimension(4,sections) :: consta, constb  !constants for calculating emission polynomial
   integer, parameter:: nsst =3 !for MAM3
   real :: sst_sz_range_lo(nsst) = &
               (/0.08e-6, 0.02e-6, 1.0e-6/)  ! accu, aitken, coarse
   real :: sst_sz_range_hi(nsst) = &
               (/1.0e-6, 0.08e-6, 10.0e-6/)
   real :: seasflux_num0(nsst),seasflux_mass0(nsst)
   real :: seasflux_num(nsst),seasflux_mass(nsst)

   !parameter for sub-grid wind variability
   REAL,DIMENSION(ims:ime,jms:jme),INTENT(INOUT) :: sigma_dry2d,sigma_wet2d,sigma_tke2d 
   REAL,DIMENSION(ims:ime,jms:jme),INTENT(INOUT) :: k_2d,c_2d 
   real(kind=8) :: sigma,sigma_dry,sigma_wet,sigma_tke
   real(kind=8) :: k_shape, c_scale
   real(kind=8),parameter :: min_sigma=0.3
   real :: rain
   real, parameter :: gravity=9.81
   real, parameter :: cp=1.00464
   real, parameter :: min_wind=0.
   real, parameter :: max_wind=50.
   integer,parameter :: nwind=200
   real(kind=8) :: subwind(nwind),windpdf(nwind)

   real w10m,w10m_0,factaa
   integer ii,jj,i, j, k, l, l_seas, n,m,ibin
   integer iphase, itype
   integer p1st

  p1st = PARAM_FIRST_SCALAR
  iphase = ai_phase 
  k=kts
  do j=jts,jte
  do i=its,ite

   e_seasalt(i,j)=0.0


   !sea-salt only over water
   if(xland(i,j).gt.1.5)then
    
   !-----------------------------------------------------
   !1. Pre-setting inputs 
   !-----------------------------------------------------
      !1.0 calculate sub-grid variability of winds
      !dry sigma
      sigma_dry=0.0
      if (hfx(i,j).gt.0) sigma_dry=((hfx(i,j)*gravity*pblh(i,j))/(rho(i,k,j)*cp*t2(i,j)))**(1/3.0)
    
      !wet sigma 
      !rain should in unit of cm/day, from mm/time_step
      rain=(rainncv(i,j)+raincv(i,j))/10./dtstep*3600.*24.
      rain=rain/4.0
      sigma_wet=alog(1.0+6.69*rain-0.476*rain**2)
      !tke sigma
      sigma_tke = 0.4*sqrt(tke_pbl(i,k,j))
      if (sigma_wet.ge.0.3) sigma_tke=0.

      !sigma=sqrt(sigma_wet**2+sigma_dry**2+sigma_tke**2)
      !not use sigma_dry over ocean for now
      !total sigma
      sigma=REAL(sqrt(sigma_wet**2+sigma_tke**2),8)
      if (sigma.lt.min_sigma) sigma=min_sigma
      !write(0,*) 'czhao check sigma', i,j,sigma, sigma_wet,sigma_tke
      !print*, 'czhao check sigma', i,j,sigma, sigma_wet,sigma_tke

      w10m_0=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
      k_shape= REAL((w10m_0/sigma)**1.086,8)
      c_scale= REAL(REAL(w10m_0,8)/fgamma(1.+1./k_shape),8) 
      !write(0,*) 'czhao check k,c', k_shape,c_scale,w10m_0,fgamma(1.+1./k_shape) 
      !print*, 'czhao check k,c', k_shape,c_scale,w10m_0,fgamma(1.+1./k_shape)

     do ibin=1,nwind 
      subwind(ibin)=REAL(min_wind+ibin*(max_wind-min_wind)/nwind,8)
      !write(0,*) 'czhao check windpdf',ibin,subwind(ibin),k_shape/c_scale,subwind(ibin)/c_scale,k_shape-1,(subwind(ibin)/c_scale)**k_shape,(k_shape/c_scale)*(subwind(ibin)/c_scale)**(k_shape-1)*exp(-((subwind(ibin)/c_scale)**k_shape))
      !print*, 'czhao check windpdf',ibin,subwind(ibin),k_shape/c_scale,subwind(ibin)/c_scale,k_shape-1,(subwind(ibin)/c_scale)**k_shape,(k_shape/c_scale)*(subwind(ibin)/c_scale)**(k_shape-1)*exp(-((subwind(ibin)/c_scale)**k_shape))

      windpdf(ibin)=(k_shape/c_scale)*(subwind(ibin)/c_scale)**(k_shape-1)*exp(-((subwind(ibin)/c_scale)**k_shape))
     enddo
      windpdf=windpdf/sum(windpdf)

      seasflux_num=0.0
      seasflux_mass=0.0
     do ibin=1,nwind 

      w10m_0=subwind(ibin)

      !1.1 10-m wind speed
      ! we don't trust the u10,v10 values, for thin model surface layer
!     if(dz8w(i,k,j).le.10.) w10m=sqrt(u_phy(i,k,j)*u_phy(i,k,j)+v_phy(i,k,j)*v_phy(i,k,j))
      ! we need them to the 3.41 power, according to Gong et al., 1997:
      w10m_0=w10m_0**3.41

      ! Calculations of source strength and size distribution
      ! NB the 0.1 is the dlogDp we have to multiplie with to get the flux, but the value dependence
      ! of course on what dlogDp you have. You will also have to change the sections of Dg if you use
      ! a different number of size bins with different intervals.
      w10m=3.84e-6*w10m_0*0.1 ! whitecap area


      !1.2 calculate constants form emission polynomials
      ! use Ekman's ss
      rdry(:)=Dg(:)/2.   ! meter
      rm(:)=1.814*rdry(:)*1.e6   ! um
      bm(:)=(0.380-log10(rm(:)))/0.65  ! use in Manahan

      do m=1,sections
       if (m.le.9)then
         consta(1,m) = (-2.576)*10.**35*Dg(m)**4+5.932*10.**28  &
                   * Dg(m)**3+(-2.867)*10.**21*Dg(m)**2+(-3.003)  &
                   * 10.**13*Dg(m) + (-2.881)*10.**6
         constb(1,m) = 7.188*10.**37  &
                   * Dg(m)**4+(-1.616)*10.**31*Dg(m)**3+6.791*10.**23  &
                   * Dg(m)**2+1.829*10.**16*Dg(m)+7.609*10.**8
       elseif (m.ge.10.and.m.le.13)then
         consta(2,m) = (-2.452)*10.**33*Dg(m)**4+2.404*10.**27  &
                   * Dg(m)**3+(-8.148)*10.**20*Dg(m)**2+(1.183)*10.**14  &
                   * Dg(m)+(-6.743)*10.**6
         constb(2,m) = 7.368*10.**35  &
                   * Dg(m)**4+(-7.310)*10.**29*Dg(m)**3+ 2.528*10.**23  &
                   * Dg(m)**2+(-3.787)*10.**16*Dg(m)+ 2.279*10.**9
       elseif (m.ge.14.and.m.lt.22)then
         consta(3,m) = 1.085*10.**29*Dg(m)**4+(-9.841)*10.**23  &
                   * Dg(m)**3+(3.132)*10.**18*Dg(m)**2+(-4.165)*10.**12  &
                   * Dg(m)+(2.181)*10.**6
         constb(3,m) = (-2.859)*10.**31  &
                   * Dg(m)**4+(2.601)*10.**26*Dg(m)**3+(-8.297)*10.**20  &
                   * Dg(m)**2+(1.105)*10.**15*Dg(m)+(-5.800)*10.**8
       elseif (m.ge.22.and.m.le.40)then
         ! use monahan
         consta(4,m) = (1.373*rm(m)**(-3)*(1+0.057*rm(m)**1.05)  &
                   * 10**(1.19*exp(-bm(m)**2)))  &
                   * (rm(m)-rm(m-1))
       endif
      enddo


   !-----------------------------------------------------
   !2. Emissions 
   !-----------------------------------------------------
      !2.1 calculate number flux fi (#/m2/s)
      fi(:)=0.
      do m=1,sections
       if (m.le.9)then
         fi(m)=w10m*(sstemp(i,j)*consta(1,m)+constb(1,m))
       elseif (m.ge.10.and.m.le.13)then
         fi(m)=w10m*(sstemp(i,j)*consta(2,m)+constb(2,m))
       elseif (m.ge.14.and.m.lt.22)then
         fi(m)=w10m*(sstemp(i,j)*consta(3,m)+constb(3,m))
       elseif (m.ge.22.and.m.le.40)then
! use Monahan
         fi(m)=consta(4,m)*w10m_0
       endif
      enddo

      !2.2 calculate mass flux e_seasalt (kg/m2/s)
      seasflux_num0(:)=0.0
      seasflux_mass0(:)=0.0
      do m=1,nsst
         do n=1, sections
          if (Dg(n).ge.sst_sz_range_lo(m) .and. Dg(n).lt.sst_sz_range_hi(m)) then
           seasflux_num0(m)=seasflux_num0(m)+fi(n)*1.35  !++ ag: scale sea-salt
          endif
         enddo
         do n=1, sections
          if (Dg(n).ge.sst_sz_range_lo(m) .and. Dg(n).lt.sst_sz_range_hi(m)) then
           seasflux_mass0(m)=seasflux_mass0(m)+fi(n)*1.35  &   !++ ag: scale sea-salt
           *4./3.*3.1415926*rdry(n)**3*2.20*1.e3  ! should use dry size, convert from number to mass flux (kg/m2/s)
          endif
         enddo
      enddo

      seasflux_mass=seasflux_mass+windpdf(ibin)*seasflux_mass0
      seasflux_num=seasflux_num+windpdf(ibin)*seasflux_num0

    enddo !nwind


    sigma_dry2d(i,j)=sigma_dry
    sigma_wet2d(i,j)=sigma_wet
    sigma_tke2d(i,j)=sigma_tke
    k_2d(i,j)=k_shape
    c_2d(i,j)=c_scale

      e_seasalt(i,j)=sum(seasflux_mass(:)*1.e9)  !kg/m2/s to ug/m2/s

   !-----------------------------------------------------
   !3. Emissions 
   !-----------------------------------------------------
     !3.1 Factor converting emission to concentration
     !factaa is (s*m2/kg-air)
     !factaa*ssemfact_mass is (s*m2/kg-air)*(ug/m2/s) = ug/kg-air
     !factaa*ssemfact_numb is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
     factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)

     !3.2 Update 
     do itype = 1, ntype_aer
     do n = 1, nsize_aer(itype)  !nsize_aer should be 1 for MAM3
        if (n.gt.1) then 
          print*,'wrong size number for MAM3'
          stop
        endif
        if (seasflux_mass(itype) <= 0.0) cycle

        !only apply emissions if bin has both na and cl species
        l_seas = lptr_seas_aer(n,itype,iphase)
        if (l_seas < p1st) cycle
        chem(i,k,j,l_seas) = chem(i,k,j,l_seas) +   &
                      factaa * seasflux_mass(itype)*1.e9 ! 1.e9 is for kg/m2/s to ug/m2/s

        l = numptr_aer(n,itype,iphase)
        if (l >= p1st) chem(i,k,j,l) = chem(i,k,j,l) +   &
                      factaa * seasflux_num(itype)
     enddo !n
     enddo !itype

   endif !xland

  enddo !i
  enddo !j


        return

END subroutine cam_mam_cam5_seasalt_emiss


!c----------------------------------------------------------------------
!c   following is from gong06b.f in
!c	/net/cirrus/files1/home/rce/oldfiles1/box/seasaltg
!c----------------------------------------------------------------------
	subroutine seasalt_emitfactors_1bin( ireduce_smallr_emit,	&
      		dpdrylo_cm, dpdryhi_cm,	  &
                emitfact_numb, emitfact_surf, emitfact_mass )
!c
!c   computes seasalt emissions factors for a specifed 
!c   dry particle size range
!c	dpdrylo_cm  = lower dry diameter (cm)
!c	dpdryhi_cm  = upper dry diameter (cm)
!c
!c   number and mass emissions are then computed as
!c	number   emissions (#/m2/s)   == emitfact_numb * (spd10*3.41)
!c	dry-sfc  emissions (cm2/m2/s) == emitfact_surf * (spd10*3.41)
!c	dry-mass emissions (g/m2/s)   == emitfact_mass * (spd10*3.41)
!c
!c   where spd10 = 10 m windspeed in m/s
!c
!c   uses bubble emissions formula (eqn 5a) from 
!c	Gong et al. [JGR, 1997, p 3805-3818]
!c
!c   *** for rdry < rdry_star, this formula overpredicts emissions.
!c	A strictly ad hoc correction is applied to the formula,
!c	based on sea-salt size measurements of
!c	O'Dowd et al. [Atmos Environ, 1997, p 73-80]
!c
!c   *** the correction is only applied when ireduce_smallr_emit > 0
!c
	implicit none

!c   subr arguments
	integer ireduce_smallr_emit
	real dpdrylo_cm, dpdryhi_cm,				&
                emitfact_numb, emitfact_surf, emitfact_mass

!c   local variables
	integer isub_bin, nsub_bin

	real alnrdrylo
	real drydens, drydens_43pi_em12, x_4pi_em8
	real dum, dumadjust, dumb, dumexpb
	real dumsum_na, dumsum_ma, dumsum_sa
	real drwet, dlnrdry
	real df0drwet, df0dlnrdry, df0dlnrdry_star
	real relhum
	real rdry, rdrylo, rdryhi, rdryaa, rdrybb
	real rdrylowermost, rdryuppermost, rdry_star
	real rwet, rwetaa, rwetbb
	real rdry_cm, rwet_cm
	real sigmag_star
	real xmdry, xsdry

	real pi
	parameter (pi = 3.1415936536)

!c   c1-c4 are constants for seasalt hygroscopic growth parameterization
!c   in Eqn 3 and Table 2 of Gong et al. [1997]
	real c1, c2, c3, c4, onethird
	parameter (c1 = 0.7674)
	parameter (c2 = 3.079)
	parameter (c3 = 2.573e-11)
	parameter (c4 = -1.424)
	parameter (onethird = 1.0/3.0)


!c   dry particle density (g/cm3)
	drydens = 2.165
!c   factor for radius (micrometers) to mass (g)
	drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12
!c   factor for radius (micrometers) to surface (cm2)
	x_4pi_em8 = 4.0*pi*1.0e-8
!c   bubble emissions formula assume 80% RH
	relhum = 0.80

!c   rdry_star = dry radius (micrometers) below which the
!c   dF0/dr emission formula is adjusted downwards
	rdry_star = 0.1
	if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20
!c   sigmag_star = geometric standard deviation used for
!c   rdry < rdry_star
	sigmag_star = 1.9

!c   initialize sums
	dumsum_na = 0.0
	dumsum_sa = 0.0
	dumsum_ma = 0.0

!c   rdrylowermost, rdryuppermost = lower and upper 
!c   dry radii (micrometers) for overall integration
        rdrylowermost = dpdrylo_cm*0.5e4
        rdryuppermost = dpdryhi_cm*0.5e4

!c
!c   "section 1"
!c   integrate over rdry > rdry_star, where the dF0/dr emissions 
!c   formula is applicable
!c   (when ireduce_smallr_emit <= 0, rdry_star = -1.0e20,
!c   and the entire integration is done here)
!c
	if (rdryuppermost .le. rdry_star) goto 2000

!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = max( rdrylowermost, rdry_star )
        rdryhi = rdryuppermost

	nsub_bin = 1000

	alnrdrylo = log( rdrylo )
	dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

!c   compute rdry, rwet (micrometers) at lowest size
	rdrybb = exp( alnrdrylo )
	rdry_cm = rdrybb*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetbb = rwet_cm*1.0e4

	do 1900 isub_bin = 1, nsub_bin

!c   rdry, rwet at sub_bin lower boundary are those
!c   at upper boundary of previous sub_bin
	rdryaa = rdrybb
	rwetaa = rwetbb

!c   compute rdry, rwet (micrometers) at sub_bin upper boundary
	dum = alnrdrylo + isub_bin*dlnrdry
	rdrybb = exp( dum )

	rdry_cm = rdrybb*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetbb = rwet_cm*1.0e4

!c   geometric mean rdry, rwet (micrometers) for sub_bin
	rdry = sqrt(rdryaa * rdrybb)
	rwet = sqrt(rwetaa * rwetbb)
	drwet = rwetbb - rwetaa

!c   xmdry is dry mass in g
	xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
	xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumb is "B" in Gong's Eqn 5a
!c   df0drwet is "dF0/dr" in Gong's Eqn 5a
	dumb = ( 0.380 - log10(rwet) ) / 0.650
	dumexpb = exp( -dumb*dumb)
	df0drwet = 1.373 * (rwet**(-3.0)) * 			&
      		(1.0 + 0.057*(rwet**1.05)) * 			&
      		(10.0**(1.19*dumexpb))

	dumsum_na = dumsum_na + drwet*df0drwet
	dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
	dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

1900	continue


!c
!c   "section 2"
!c   integrate over rdry < rdry_star, where the dF0/dr emissions 
!c   formula is just an extrapolation and predicts too many emissions
!c
!c   1.  compute dF0/dln(rdry) = (dF0/drwet)*(drwet/dlnrdry) 
!c	at rdry_star
!c   2.  for rdry < rdry_star, assume dF0/dln(rdry) is lognormal,
!c	with the same lognormal parameters observed in 
!c	O'Dowd et al. [1997]
!c
2000	if (rdrylowermost .ge. rdry_star) goto 3000

!c   compute dF0/dln(rdry) at rdry_star
	rdryaa = 0.99*rdry_star
	rdry_cm = rdryaa*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetaa = rwet_cm*1.0e4

	rdrybb = 1.01*rdry_star
	rdry_cm = rdrybb*1.0e-4
	rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/		&
      		( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
	rwetbb = rwet_cm*1.0e4

	rwet = 0.5*(rwetaa + rwetbb)
	dumb = ( 0.380 - log10(rwet) ) / 0.650
	dumexpb = exp( -dumb*dumb)
	df0drwet = 1.373 * (rwet**(-3.0)) * 			&
      		(1.0 + 0.057*(rwet**1.05)) * 			&
      		(10.0**(1.19*dumexpb))

	drwet = rwetbb - rwetaa
	dlnrdry = log( rdrybb/rdryaa )
	df0dlnrdry_star = df0drwet * (drwet/dlnrdry)


!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = rdrylowermost
        rdryhi = min( rdryuppermost, rdry_star )

	nsub_bin = 1000

	alnrdrylo = log( rdrylo )
	dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

	do 2900 isub_bin = 1, nsub_bin

!c   geometric mean rdry (micrometers) for sub_bin
	dum = alnrdrylo + (isub_bin-0.5)*dlnrdry
	rdry = exp( dum )

!c   xmdry is dry mass in g
	xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
	xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumadjust is adjustment factor to reduce dF0/dr
	dum = log( rdry/rdry_star ) / log( sigmag_star )
	dumadjust = exp( -0.5*dum*dum )

	df0dlnrdry = df0dlnrdry_star * dumadjust

	dumsum_na = dumsum_na + dlnrdry*df0dlnrdry
	dumsum_ma = dumsum_ma + dlnrdry*df0dlnrdry*xmdry
	dumsum_sa = dumsum_sa + dlnrdry*df0dlnrdry*xsdry

2900	continue


!c
!c  all done
!c
3000	emitfact_numb = dumsum_na
	emitfact_mass = dumsum_ma
	emitfact_surf = dumsum_sa

	return
	end subroutine seasalt_emitfactors_1bin



!----------------------------------------------------------------------
   subroutine cam_mam_mosaic_dust_emiss(  slai,ust, smois, ivgtyp, isltyp, &
               id, dtstep, u10, v10, alt, dz8w, xland, config_flags, chem, &
               dust_emiss_size_opt,                                        &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
!
! adds dust emissions for mosaic aerosol species (i.e. emission tendencies
! over time dtstep are applied to the aerosol mixing ratios)
!
! This is a simple dust scheme based on Shaw et al. (2008) to appear in
! Atmospheric Environment, recoded by Jerome Fast
!
! NOTE: 
! 1) This version only works with the 8-bin version of MOSAIC.
! 2) Dust added to MOSAIC's other inorganic specie, OIN.  If Ca and CO3 are 
!    activated in the Registry, a small fraction also added to Ca and CO3.
! 3) The main departure from Shaw et al., is now alphamask is computed since
!    the land-use categories in that paper and in WRF differ.  WRF currently 
!    does not have that many land-use categories and adhoc assumptions had to
!    be made. This version was tested for Mexico in the dry season.  The main
!    land-use categories in WRF that are likely dust sources are grass, shrub,
!    and savannna (that WRF has in the desert regions of NW Mexico).  Having
!    dust emitted from these types for other locations and other times of the
!    year is not likely to be valid.
! 4) An upper bound on ustar was placed because the surface parameterizations
!    in WRF can produce unrealistically high values that lead to very high
!    dust emission rates.
! 5) Other departures' from Shaw et al. noted below, but are probably not as
!    important as 2) and 3).
!

   USE module_configure, only:  grid_config_rec_type
   USE module_state_description, only:  num_chem, param_first_scalar

   USE module_data_cam_mam_asect, only:  ai_phase,   &
         lptr_dust_aer, ntype_aer, numptr_aer

   USE modal_aero_data, only:  &
         modeptr_accum, modeptr_coarse, modeptr_coardust, modeptr_finedust

   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,                                      &
                                  dust_emiss_size_opt,                     &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   ) ::    dtstep

! 10-m wind speed components (m/s)
   REAL,  DIMENSION( ims:ime , jms:jme ),                                  &
          INTENT(IN   ) ::   u10, v10, xland, slai, ust
   INTEGER,  DIMENSION( ims:ime , jms:jme ),                               &
          INTENT(IN   ) ::   ivgtyp, isltyp

! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::   chem

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::   alt, dz8w

   REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,     &
          INTENT(IN   ) ::   smois

! local variables
        integer, parameter :: max_dust_mode = 2
        integer :: i, ii, j, k, l, l_dust, l_numb, n
        integer :: imode, iphase, isize, itype, izob
        integer :: p1st
        integer :: isize_dust_mode(max_dust_mode)
        integer :: itype_dust_mode(max_dust_mode)
        integer :: n_dust_mode

        real :: dum, dumdlo, dumdhi, dumlandfrac, dumspd10
        real :: factaa, factbb, factcc, fracoin, fracca, fracco3, fractot
        real :: ustart, ustar1, ustart0
        real :: alphamask, f8, f50, f51, f52, wetfactor, sumdelta, ftot
        real :: smois_grav, wp, pclay
        real :: beta(4,7)
        real :: gamma(4), delta(4)
        real :: sz(8)
        real :: dustflux, densdust
        real :: mass1part_mos8bin(8)
        real :: sz_wght_dust_mode(8,max_dust_mode)
        real :: dcen_mos8bin(8)

        character(len=100) :: msg

        write(*,*) 'in subr cam_mam_mosaic_dust_emiss'
        p1st = param_first_scalar


! the mass emission fractions [sz(1:8)] were derived for the 
!    mosaic 8-bin size bins
! the emissions for each mosaic bin are assigned to the cam-mam fine or coarse mode
!    (e.g, the fine mode may get bins 1-5 and the coarse mode bins 6-8)
! for calculating number emissions from mass emissions, use the mosaic bin
!    sizes (and corresponding 1-particle masses) in order to get number emissions
!    consistent with mosaic
!
! note:  the sz() values were recommend by Natalie Mahowold, so probably follow the
!    assumptions used in the "dead" dust module.  
!    applying those assumptions along with the cut sizes might be a better approach
!    than adapting the mosaic 8-bin sz() values

#if (defined MODAL_AERO)
! 3 modes -- dust is in accum and coarse modes
        n_dust_mode = 2
        itype_dust_mode(1) = modeptr_accum
        itype_dust_mode(2) = modeptr_coarse
        isize_dust_mode(:) = 1
        if (dust_emiss_size_opt == 1) then
! cam5 cut sizes for dust emiss are 0.1, 1.0, 10.0 um
! fine mode gets bins 1-4 and 60% of bin 5
           sz_wght_dust_mode(1:8,1) = (/ 1.,  1.,  1.,  1.,  0.6, 0.,  0.,  0.  /)
           sz_wght_dust_mode(1:8,2) = (/ 0.,  0.,  0.,  0.,  0.4, 1.,  1.,  1.  /)
        else if (dust_emiss_size_opt == 2) then
! sorgam cut sizes for dust emiss are 0.1, 2.5, 10.0 um
! fine mode gets bins 1-6
           sz_wght_dust_mode(1:8,1) = (/ 1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.  /)
           sz_wght_dust_mode(1:8,2) = (/ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.  /)
        else
            write(msg,'(2a,i7)') 'subr cam_mam_mosaic_dust_emiss', &
                ' - illegal dust_emiss_size_opt = ', dust_emiss_size_opt
            call wrf_error_fatal( msg )
        end if
#else
! 7 modes -- dust is in fine-dust and coarse-dust modes
        n_dust_mode = 2
        itype_dust_mode(1) = modeptr_finedust
        itype_dust_mode(2) = modeptr_coardust
        isize_dust_mode(:) = 1
        if (dust_emiss_size_opt == 1) then
! cam5 cut sizes for dust emiss are 0.1, 2.0, 10.0 um
! fine mode gets bins 1-5 and 60% of bin 6
! (for 7-mode cam, the fine-dust mode is larger than the accumulation mode,
!  so the cut size is larger than with 3-mode cam)
           sz_wght_dust_mode(1:8,1) = (/ 1.,  1.,  1.,  1.,  1.,  0.6, 0.,  0.  /)
           sz_wght_dust_mode(1:8,2) = (/ 0.,  0.,  0.,  0.,  0.,  0.4, 1.,  1.  /)
        else if (dust_emiss_size_opt == 2) then
! sorgam cut sizes for dust emiss are 0.1, 2.5, 10.0 um
! fine mode gets bins 1-6
           sz_wght_dust_mode(1:8,1) = (/ 1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.  /)
           sz_wght_dust_mode(1:8,2) = (/ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.  /)
        else
            write(msg,'(2a,i7)') 'subr cam_mam_mosaic_dust_emiss', &
                'illegal dust_emiss_size_opt = ', dust_emiss_size_opt
            call wrf_error_fatal( msg )
        end if
#endif

! bin center diameters for mosaic-8bin (cm)
        dcen_mos8bin(8) = 1.0e-4*(10.0/sqrt(2.0))
        do n = 7, 1, -1
           dcen_mos8bin(n) = dcen_mos8bin(n+1)*0.5
        end do
! mass1part_mos8bin is mass of a single particle in ug, 
! assuming density of dust ~2.5 g cm-3
        densdust=2.5
        do n = 1, 8
           mass1part_mos8bin(n)=0.523598*(dcen_mos8bin(n)**3)*densdust*1.0e06
        end do


!
! from Nickovic et al., JGR, 2001 and Shaw et al. 2007
! beta: fraction of clay, small silt, large silt, and sand correcsponding to Zobler class (7)
! beta (1,*) for 0.5-1 um
! beta (2,*) for 1-10 um
! beta (3,*) for 10-25 um
! beta (4,*) for 25-50 um
!
        beta(1,1)=0.12
        beta(2,1)=0.04
        beta(3,1)=0.04
        beta(4,1)=0.80
        beta(1,2)=0.34
        beta(2,2)=0.28
        beta(3,2)=0.28
        beta(4,2)=0.10
        beta(1,3)=0.45
        beta(2,3)=0.15
        beta(3,3)=0.15
        beta(4,3)=0.25
        beta(1,4)=0.12
        beta(2,4)=0.09
        beta(3,4)=0.09
        beta(4,4)=0.70
        beta(1,5)=0.40
        beta(2,5)=0.05
        beta(3,5)=0.05
        beta(4,5)=0.50
        beta(1,6)=0.34
        beta(2,6)=0.18
        beta(3,6)=0.18
        beta(4,6)=0.30
        beta(1,7)=0.22
        beta(2,7)=0.09
        beta(3,7)=0.09
        beta(4,7)=0.60
        gamma(1)=0.08
        gamma(2)=1.00
        gamma(3)=1.00
        gamma(4)=0.12
!
! * Mass fractions for each size bin. These values were recommended by 
!   Natalie Mahowold, with bins 7 and 8 the same as bins 3 and 4 from CAM.
! * Changed slightly since Natelie's estimates do not add up to 1.0
! * This would need to be made more generic for other bin sizes.
!       sz(1)=0
!       sz(2)=1.78751e-06
!       sz(3)=0.000273786
!       sz(4)=0.00847978
!       sz(5)=0.056055
!       sz(6)=0.0951896
!       sz(7)=0.17
!       sz(8)=0.67
!jdf    sz(1)=0.0
!jdf    sz(2)=0.0
!jdf    sz(3)=0.0005
!jdf    sz(4)=0.0095
!jdf    sz(5)=0.03
!jdf    sz(6)=0.10
!jdf    sz(7)=0.18
!jdf    sz(8)=0.68
        sz(1)=0.0
        sz(2)=0.0
        sz(3)=0.0005
        sz(4)=0.0095
        sz(5)=0.01
        sz(6)=0.06
        sz(7)=0.20
        sz(8)=0.72

        iphase = ai_phase

!   loop over i,j and apply dust emissions
        k = kts
        do 1830 j = jts, jte
        do 1820 i = its, ite

    if( xland(i,j) > 1.5 ) cycle

! compute wind speed anyway, even though ustar is used below

        dumlandfrac = 1.
        dumspd10=(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))**(0.5)
        if(dumspd10 >= 5.0) then
           dumspd10 = dumlandfrac* &
         ( dumspd10*dumspd10*(dumspd10-5.0))
         else
            dumspd10=0.
         endif

! part1 - compute vegetation mask
!
! * f8, f50, f51, f52 refer to vegetation classes from the Olsen categories
!   for desert, sand desert, grass semi-desert, and shrub semi-desert
! * in WRF, vegetation types 7, 8 and 10 are grassland, shrubland, and savanna
!   that are dominate types in Mexico and probably have some erodable surface
!   during the dry season
! * currently modified these values so that only a small fraction of cell
!   area is erodable
! * these values are highly tuneable!

         alphamask=0.001
         if (ivgtyp(i,j) .eq. 7) then
           f8=0.005
           f50=0.00
           f51=0.10
           f51=0.066
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 8) then
           f8=0.010
           f50=0.00
           f51=0.00
           f52=0.15
           f52=0.10
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif
         if (ivgtyp(i,j) .eq. 10) then
           f8=0.00
           f50=0.00
           f51=0.01
           f52=0.00
           alphamask=(f8+f50)*1.0+(f51+f52)*0.5
         endif

! part2 - zobler
! 
! * in Shaw's paper, dust is computed for 4 size ranges:
!   0.5-1 um 
!    1-10 um  
!   10-25 um  
!   25-50 um
! * Shaw's paper also accounts for sub-grid variability in soil
!   texture, but here we just assume the same soil texture for each
!   grid cell
! * since MOSAIC is currently has a maximum size range up to 10 um,
!   neglect upper 2 size ranges and lowest size range (assume small)
! * map WRF soil classes arbitrarily to Zolber soil textural classes
! * skip dust computations for WRF soil classes greater than 13, i.e. 
!   do not compute dust over water, bedrock, and other surfaces
! * should be skipping for water surface at this point anyway
!
         izob=0
         if(isltyp(i,j).eq.1) izob=1
         if(isltyp(i,j).eq.2) izob=1
         if(isltyp(i,j).eq.3) izob=4
         if(isltyp(i,j).eq.4) izob=2
         if(isltyp(i,j).eq.5) izob=2
         if(isltyp(i,j).eq.6) izob=2
         if(isltyp(i,j).eq.7) izob=7
         if(isltyp(i,j).eq.8) izob=2
         if(isltyp(i,j).eq.9) izob=6
         if(isltyp(i,j).eq.10) izob=5
         if(isltyp(i,j).eq.11) izob=2
         if(isltyp(i,j).eq.12) izob=3
         if(isltyp(i,j).ge.13) izob=0
         if(izob.eq.0) goto 1820
!
! part3 - dustprod
!
         do ii=1,4
           delta(ii)=0.0
         enddo
         sumdelta=0.0
         do ii=1,4
           delta(ii)=beta(ii,izob)*gamma(ii)
           if(ii.lt.4) then
             sumdelta=sumdelta+delta(ii)
           endif
         enddo
         do ii=1,4
           delta(ii)=delta(ii)/sumdelta
         enddo

! part4 - wetness
!
! * assume dry for now, have passed in soil moisture to this routine
!   but needs to be included here
! * wetfactor less than 1 would reduce dustflux
! * convert model soil moisture (m3/m3) to gravimetric soil moisture
!   (mass of water / mass of soil in %) assuming a constant density 
!   for soil
         pclay=beta(1,izob)*100.
         wp=0.0014*pclay*pclay+0.17*pclay
         smois_grav=(smois(i,1,j)/2.6)*100.
         if(smois_grav.gt.wp) then
           wetfactor=sqrt(1.0+1.21*(smois_grav-wp)**0.68)
         else
           wetfactor=1.0
         endif
!        wetfactor=1.0

! part5 - dustflux
! lower bound on ustar = 20 cm/s as in Shaw et al, but set upper
! bound to 100 cm/s

         ustar1=ust(i,j)*100.0
         if(ustar1.gt.100.0) ustar1=100.0
         ustart0=20.0
         ustart=ustart0*wetfactor
         if(ustar1.le.ustart) then
           dustflux=0.0
         else
           dustflux=1.0e-14*(ustar1**4)*(1.0-(ustart/ustar1))
         endif
         dustflux=dustflux*10.0
! units kg m-2 s-1
         ftot=0.0
         do ii=1,2
           ftot=ftot+dustflux*alphamask*delta(ii)
         enddo
! convert to ug m-2 s-1
         ftot=ftot*1.0e+09

!   apportion other inorganics only
         factaa = (dtstep/dz8w(i,k,j))*alt(i,k,j)
         factbb = factaa * ftot
         fractot = 1.0

         do imode = 1, n_dust_mode
!   loop over the cam-mam modes that have dust
         itype = itype_dust_mode(imode)
         isize = isize_dust_mode(imode)

         l_dust = lptr_dust_aer(isize,itype,iphase)
         if ((l_dust < p1st) .or. (l_dust > num_chem)) cycle
         l_numb = numptr_aer(isize,itype,iphase)
         if ((l_numb < p1st) .or. (l_numb > num_chem)) l_numb = -1

!   if (ivgtyp(i,j) .eq. 8) print*,'jdf',i,j,ustar1,ustart0,factaa,ftot
!        do 1810 n = 1, nsize_aer(itype)
         do n = 1, 8
!   calculate contribution of each mosaic-8bin bin to the current cam-mam mode
            if (sz_wght_dust_mode(n,imode) <= 0.0) cycle
            factcc = factbb * sz(n) * fractot * sz_wght_dust_mode(n,imode)

            chem(i,k,j,l_dust) = chem(i,k,j,l_dust) + factcc

            if (l_numb >= p1st) &
            chem(i,k,j,l_numb) = chem(i,k,j,l_numb) + factcc/mass1part_mos8bin(n)

        end do ! n

        end do ! imode


1820    continue
1830    continue

        return

   END subroutine cam_mam_mosaic_dust_emiss

!====================================================================================
!add dust emission scheme following CLM mechanism (DEAD)  --czhao  08/10/2011
!Reference: Mahowald et al. (2006); Zender et al. (2003)
!====================================================================================
  subroutine cam_mam_dust_DEADemis (               &
         dt,num_soil_layers,alt,u_phy,v_phy,  &
         chem,rho_phy,dz8w,smois,ust,u10,v10,                  &
         ivgtyp,isltyp,xland,g,e_dust,        &
         num_land,landf,num_soil,soilf,sh2o,snowc,xlai,        &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

  USE module_configure
  USE module_state_description
  USE module_model_constants, ONLY: mwdry
! USE module_data_mosaic_asect
  USE module_data_gocart_dust, ONLY: porosity
  USE module_data_cam_mam_asect, only:  ai_phase,   &
         lptr_dust_aer, ntype_aer, numptr_aer
  USE modal_aero_data, only:  &
         modeptr_accum, modeptr_coarse, modeptr_coardust, modeptr_finedust
  USE module_data_sorgam,  ONLY: rhosoil,soilfac

  IMPLICIT NONE

  INTEGER,      INTENT(IN   ) :: num_soil_layers,num_land,num_soil,     &
                                 ids,ide, jds,jde, kds,kde,               &
                                 ims,ime, jms,jme, kms,kme,               &
                                 its,ite, jts,jte, kts,kte

  INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
         INTENT(IN   ) ::   ivgtyp,isltyp

  REAL, DIMENSION(ims:ime,num_land,jms:jme),               &
         INTENT(IN) ::  landf

  REAL, DIMENSION(ims:ime,num_soil,jms:jme),               &
         INTENT(IN) ::  soilf

  REAL, DIMENSION(ims:ime,jms:jme),               &
         INTENT(IN) ::  snowc,xlai

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
        INTENT(INOUT ) ::                                   chem

  REAL,  DIMENSION( ims:ime , jms:jme ), INTENT(INOUT ) ::  e_dust

  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(IN) ::               smois,sh2o

  REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                           u10,                  &
                                                     v10,                  &
                                                     ust,                  &
                                                     xland
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                     dz8w,             &
                                              u_phy,v_phy,rho_phy

  REAL, INTENT(IN   ) :: dt,g

  !local variables
  REAL, PARAMETER  :: Tfactor=5.0*1.e-4 !factor for resolution sensitivity, highly tunable 
  REAL, PARAMETER  :: Sfactor=1.0       !source erodibility factor, place holder 
  INTEGER, PARAMETER  :: BARESOIL=19    !bare soil with sparse veg. from USGS 24-class data
  INTEGER, PARAMETER  :: WATERBODY=16    !water body from USGS 24-class data
  LOGICAL,PARAMETER :: prescribe_soil=.false.
  real baresoilf
  real f_v
  real,parameter  :: xlai_thresh=0.3    !threshhold for LAI
  REAL alpha !sandblasting mass efficiency (m^-1)
  REAL f_soil    !grid cell fraction of exposed bare soil 
  real mclay,mclay1 !mass fraction of clay particles 
  real f_clay !fraction of clay soil
  real index1
  real qflux  !total horizontally saltating mass flux (kg m^-1 s^-1)
  real term1,term2,term3
  real ustar   ! wind friction speed (m/s) accounting for the Owen effect
  real ustar_t ! threshold wind friction speed for saltation (m/s)
  real wind10m,wind10m_t
  REAL, PARAMETER  :: CSfactor=2.61    !saltation constant
  real totalemis
  !Accumulation mode:
  real jdustemis
  real, parameter :: dgvem_j1 = 0.832E-6 !m, mass median diameter
  real, parameter :: sgem_j1 = 2.1 !geometric standard deviation
  real, parameter :: emass_jfrac1 = 0.04 !mass fraction
  !Coarse mode:
  real cdustemis
  real, parameter :: dgvem_c1 = 4.82E-6 !m, mass median diameter
  real, parameter :: sgem_c1 = 1.9 !geometric standard deviation
  real, parameter :: emass_cfrac1 = 0.96 !mass fraction
  real accfrac,corfrac
  real jdustcon,cdustcon
  real convert_j,convert_c
  REAL :: sand_f(num_soil) ! sand fraction of 16-class soil type from Noah land surface model 
  REAL :: silt_f(num_soil) ! silt fraction of 16-class soil type from Noah land surface model 
  REAL :: clay_f(num_soil) ! clay fraction of 16-class soil type from Noah land surface model 
  integer isoil,i,j,k
  real gwet,gwet_t,bulk_dens !bulk density of soil 
  real,parameter :: liquid_dens = 1.e3  !kg/m3 
  real,parameter :: dosp = 75*10e-6  ! m, optimal saltation particles 
  real fwet !factor dependent on soil moisture
  real,parameter :: fzfactor=1.0 !factor dependent on surface roughness, set as a place holder for now
  real re_t_f,re_t

  integer, parameter :: max_dust_mode = 2
  integer :: ii,l,l_dust,l_numb
  integer :: imode, iphase, isize,itype
  integer :: p1st
  integer :: isize_dust_mode(max_dust_mode)
  integer :: itype_dust_mode(max_dust_mode)
  integer :: n_dust_mode

    p1st = PARAM_FIRST_SCALAR



!**********************************************************
!1. Initializing parameters 
!**********************************************************

   !-----------------------------------------
   !1.1 Checking land & soil information 
   !-----------------------------------------
    if (num_land.ne.24) then
      write(0,*) 'incorrect number of land category ',num_land,baresoil
      print*, 'incorrect number of land category ',num_land,baresoil
      stop
    endif

    if (num_soil.ne.16) then
      write(0,*) 'incorrect number of soil category ',num_soil
      print*, 'incorrect number of soil category ',num_soil
      stop
    endif

   !-----------------------------------------
   !1.2 Initialize soil categary percentage 
   !-----------------------------------------
    sand_f = (/92,82,58,17,10,43,58,10,32,52, 6,22,0,0,0,0/)
    silt_f = (/ 5,12,32,70,85,39,15,56,34, 6,47,20,0,0,0,0/)
    clay_f = (/ 3, 6,10,13, 5,18,27,34,34,42,47,58,0,0,0,0/)

   !-----------------------------------------
   !1.3 Initialize dust emission fraction of modes 
   !-----------------------------------------
    accfrac=emass_jfrac1
    corfrac=emass_cfrac1

!**********************************************************
!2. Dust emission calculation 
!**********************************************************

  k=kts   !surface only
  do j=jts,jte
  do i=its,ite

      !zero the emission array
      e_dust(i,j)=0.0
      totalemis = 0.0

      !baresoilf fraction
      if (prescribe_soil) then
      baresoilf = landf(i,baresoil,j)
      else
      f_v=xlai(i,j)/xlai_thresh
      f_v=max(0.0,f_v)
      f_v=min(1.0,f_v)
      baresoilf = (1.-snowc(i,j))*(1.-landf(i,waterbody,j))*(1.-f_v)
      endif

   !dust emission only over baresoil land 
    if(xland(i,j).lt.1.5.and.baresoilf.gt.0)then
   !if(landf(i,waterbody,j).lt.0.5.and.baresoilf.gt.0)then

      if (baresoilf.gt.0) then
        write(0,*) 'czhao check baresoil ',i,j,num_land,baresoilf,landf(i,baresoil,j),landf(i,waterbody,j)
        print*, 'czhao check baresoil ',i,j,num_land,baresoilf,landf(i,baresoil,j),landf(i,waterbody,j)
      endif

   !-----------------------------------------------------
   !2.1 grid cell fraction of exposed bare soil (f_soil) 
   !-----------------------------------------------------
      if (smois(i,1,j).gt.0) then
       f_soil = baresoilf*sh2o(i,1,j)/smois(i,1,j)
      else
       f_soil = baresoilf
      endif

   !-----------------------------------------------------
   !2.2 sandblasting mass efficiency (alpha) 
   !-----------------------------------------------------
      !2.2.1 clay percentage 
      f_clay=0.0
      do isoil=1,num_soil
       f_clay=f_clay+clay_f(isoil)*soilf(i,isoil,j)
      enddo

      !2.2.2 mass fraction of clay particles 
      if (f_clay.gt.20.and.f_clay.le.100) then
       mclay = 20*0.01
      else
       if (f_clay.lt.0) f_clay=0.0
       if (f_clay.gt.20) f_clay=20.0
       mclay = f_clay*0.01
      endif
      !2.2.3 alpha 
      index1= 13.4*mclay-6.0
      alpha = 100.0*exp(index1*alog(10.0))

   !-----------------------------------------------------
   !2.3 threshold wind friction speed for saltation (ustar_t)  
   !-----------------------------------------------------
      !2.3.1 soil moisture content (gwet)
      bulk_dens=(1.0-porosity(isltyp(i,j)))*2.7e3
      gwet=smois(i,1,j)*liquid_dens/bulk_dens

      !2.3.2 soil moisture dependent factor (fwet) 
      mclay1 = f_clay*0.01
      mclay1 = min(max(mclay1,0.0),1.0)
      if (mclay1==0) then
       gwet_t = 0.0
      else
       gwet_t=(0.17*mclay1+0.14*(mclay1**2))/mclay1
      endif

      if (gwet.le.gwet_t) then
       fwet = 1.0
      else
       fwet = (1+1.21*((100*(gwet-gwet_t))**0.68))**0.5
      endif

      !2.3.3 threshold friction Reynolds factor (re_t_f) 
      re_t = 0.38+1331*((100*dosp)**1.56)
      if (re_t.gt.10) then
      re_t_f = 0.0144*(1-0.0858*exp(-0.0617*(re_t-10)))**2
      else
       if (re_t.ge.0.03) then
        !re_t_f = 0.1291*0.1291/(-1+1.928*re_t)   !from CLM document
        re_t_f = 0.1291*0.1291/(-1+1.928*(re_t**0.0922)) !from CLM code
       endif
      endif

      !2.3.4 ustar_t 
      term1=re_t_f*rhosoil*g*dosp*(1+(6.e-7)/(rhosoil*g*(dosp**2.5)))
      ustar_t=fzfactor*(term1**0.5)*(rho_phy(i,k,j)**(-0.5))*fwet

   !-----------------------------------------------------
   !2.4 horizontally saltating mass flux (kg/m/s) (qflux) 
   !-----------------------------------------------------
      !2.4.1 friction velocity accounting for the Owen effect 
      wind10m=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
      ! we don't trust the u10,v10 values, if model layers are very thin near surface
      if(dz8w(i,kts,j).le.10.) wind10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
      wind10m_t=ustar_t*wind10m/ust(i,j)
      if (wind10m.lt.wind10m_t) then
       ustar = ust(i,j)
      else
       ustar = ust(i,j)+0.003*(wind10m-wind10m_t)**2
      endif

      write(0,*) 'czhao check ustar',i,j,ustar,ustar_t,wind10m,wind10m_t
      print*, 'czhao check ustar',i,j,ustar,ustar_t,wind10m,wind10m_t

      !2.4.2 friction velocity accounting for the Owen effect 
      if (ustar.gt.ustar_t) then
       term1 = csfactor*rho_phy(i,kts,j)*(ustar**3)/g
       term2 = 1-ustar_t/ustar
       term3 = (1+ustar_t/ustar)**2
       qflux = term1*term2*term3
      else
       qflux = 0.0
      endif

   !-----------------------------------------------------
   !2.5 vertical mass fluxes of dust (kg m-2 s-1) (totalemis) 
   !-----------------------------------------------------
      totalemis = tfactor*sfactor*f_soil*alpha*qflux

      write(0,*) 'czhao check dust emission',i,j,xlai(i,j),f_soil,alpha,qflux
      print*, 'czhao check dust emission',i,j,xlai(i,j),f_soil,alpha,qflux

      !2.5.1 unit convert from kg/m2/s to ug/m2/s
      totalemis = totalemis*1.e9

      !2.5.2 size distribution 
      jdustemis = totalemis*accfrac   ! accumulation mode
      cdustemis = totalemis*corfrac   ! coarse mode 

      e_dust(i,j)=cdustemis+jdustemis

      !2.5.3 convert emission to mass concentration
      !ug/m2/s to ug/kg-dry air (alt in m3/kg)
      jdustcon = (jdustemis*dt/dz8w(i,k,j))*alt(i,k,j)
      cdustcon = (cdustemis*dt/dz8w(i,k,j))*alt(i,k,j)

      !pre-processing for MAM aerosol scheme
#if (defined MODAL_AERO)
        ! 3 modes -- dust is in accum and coarse modes
        iphase=ai_phase
        n_dust_mode = 2
        itype_dust_mode(1) = modeptr_accum
        itype_dust_mode(2) = modeptr_coarse
        isize_dust_mode(:) = 1
#else
            write(msg,'(a)') 'DEAD dust emission is only available for MAM3 now'
            call wrf_error_fatal( msg )
#endif
      !2.5.4 update dust mass and number with emission 
        do imode = 1, n_dust_mode
         !   loop over the cam-mam modes that have dust
         itype = itype_dust_mode(imode)
         isize = isize_dust_mode(imode)

         l_dust = lptr_dust_aer(isize,itype,iphase)
         if ((l_dust < p1st) .or. (l_dust > num_chem)) cycle
         l_numb = numptr_aer(isize,itype,iphase)
         if ((l_numb < p1st) .or. (l_numb > num_chem)) l_numb = -1

          if (imode==1) then !accumulation mode
           convert_j = exp(4.5*log(sgem_j1)**2)/dgvem_j1**3
           if (l_dust >= p1st) & 
           chem(i,k,j,l_dust) = chem(i,k,j,l_dust) + jdustcon 
           if (l_numb >= p1st) &
           chem(i,k,j,l_numb) = chem(i,k,j,l_numb) + jdustcon*convert_j*soilfac 
          else !coarse mode
           convert_c = exp(4.5*log(sgem_c1)**2)/dgvem_c1**3
           if (l_dust >= p1st) &
           chem(i,k,j,l_dust) = chem(i,k,j,l_dust) + cdustcon
           if (l_numb >= p1st) &
           chem(i,k,j,l_numb) = chem(i,k,j,l_numb) + cdustcon*convert_c*soilfac
          endif
        end do ! imode

   endif  !xland 
  enddo !i
  enddo !j


end subroutine cam_mam_dust_DEADemis
!----------------------------------------------------------------------

     DOUBLE PRECISION FUNCTION FGAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      implicit none
      INTEGER I,N
      LOGICAL PARITY
      REAL (kind=8) :: X
      REAL                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO

      REAL, PARAMETER :: PI = 3.1415926535897932384626434
      REAL, PARAMETER :: SQRTPI = 0.9189385332046727417803297

      REAL, DIMENSION(7) :: C
      REAL, DIMENSION(8) :: P
      REAL, DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,                                 &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,    &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            RES=RES*Y
            Y=Y+ONE
          END DO
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
            SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE) RES=FACT/RES
  900 FGAMMA=REAL(RES,8)
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION FGAMMA

!----------------------------------------------------------------------
END MODULE module_cam_mam_addemiss
