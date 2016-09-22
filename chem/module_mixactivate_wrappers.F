!**********************************************************************************  
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of 
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! MOSAIC module: see module_mosaic_driver.F for references and terms of use
!**********************************************************************************  

!----------------------------------------------------------------------
! This module contains interface wrapper routines to couple the aerosol
! modules with mixactivate in the physics directory. Due to compiling
! dependencies, these cannot be placed in module_mixactivate.
!----------------------------------------------------------------------

MODULE module_mixactivate_wrappers

CONTAINS

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine mosaic_mixactivate (                        &
           id, ktau, dtstep, config_flags, idrydep_onoff,    &
           !rho_phy, t_phy, w, cldfra, cldfra_old,            &
           p_phy,rho_phy, t_phy, w, cldfra, cldfra_old,            &
           ddvel, z, dz8w, p_at_w, t_at_w, exch_h,           &
           !qv, qc, qi, qndrop3d, f_qc, f_qi, chem,           &
     !czhao
     qv,qc,qi,qs,f_ice_phy,f_rain_phy,qndrop3d,qnice3d,f_qv,f_qc,f_qi,f_qs,chem, &
           ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
        !czhao 
        nihf3d,niimm3d,nidep3d,nimey3d,  &
        rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4,  &
        tke_pbl,relhum3d,wsub3d,  &
           ids,ide, jds,jde, kds,kde,                        &
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

    USE module_configure, only: grid_config_rec_type
	use module_state_description, only:  num_chem
	use module_data_mosaic_asect
	use module_mixactivate, only:  icenucleat,mixactivate

! wrapper to call mixactivate for mosaic description of aerosol

	implicit none

!   subr arguments
	integer, intent(in) ::               &
         id, ktau,            &
         ids, ide, jds, jde, kds, kde,   &
         ims, ime, jms, jme, kms, kme,   &
         its, ite, jts, jte, kts, kte,   &
         idrydep_onoff

	real, intent(in) :: dtstep

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		rho_phy, t_phy, w, p_phy,  &
		z, dz8w, p_at_w, t_at_w, exch_h

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: cldfra, cldfra_old

	real, intent(in),   &
		dimension( its:ite, jts:jte, num_chem ) :: ddvel

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		qv, qc, qi,qs, f_ice_phy,f_rain_phy

    LOGICAL, intent(in) :: f_qv,f_qc, f_qi,f_qs

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		qndrop3d,qnice3d

      real, intent(inout), dimension(ims:ime,kms:kme,jms:jme) ::  &
        nihf3d,niimm3d,nidep3d,nimey3d
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) ::  &
       rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4  ! dust for contact freezing 

      REAL, intent(out), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
       relhum3d,wsub3d

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
         INTENT(IN ) :: tke_pbl

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource,&
	     ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  ! number conc of aerosols activated at supersat

	type(grid_config_rec_type), intent(in) :: config_flags
! local vars
	real qsrflx(ims:ime, jms:jme, num_chem) ! dry deposition flux of aerosol
	real sumhygro,sumvol
	integer i,j,k,l,m,n
	real hygro( its:ite, kts:kte, jts:jte, maxd_asize, maxd_atype ) ! bulk

! aerosol mass and number #/cm^3 or m2/cm^3, for 3 size ranges
  real,dimension(ims:ime,kms:kme,jms:jme,1:3) ::  &
  so4_num,soot_num,organic_num,dst_num,no3_num,nh4_num,oin_num,so4_sfc,soot_sfc,organic_sfc,dst_sfc,no3_sfc,nh4_sfc,oin_sfc
  integer lnum,lmass
  real scale0 
  real chem1d(1:num_chem) ! aerosol number and mass mixing ratios (#/m3 or kg/m3) 
  REAL nacon1d(4), rndst1d(4)




! calculate volume-weighted bulk hygroscopicity for each type and size
      do 100 j=jts,jte
      do 100 k=kts,kte
      do 100 i=its,ite
       do n=1,ntype_aer
       do m=1,nsize_aer(n)
	       sumhygro=0.
	       sumvol=0.
	       do l=1,ncomp_aer(n)
                  lmass=massptr_aer(l,m,n,ai_phase)
	          sumhygro = sumhygro+hygro_aer(l,n)*   &
                   chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	          sumvol = sumvol+chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	       end do ! comp
           hygro(i,k,j,m,n)=sumhygro/sumvol
	end do ! size
	end do ! type
  100 continue

!============================================================
!aerosol ice nucleation  -czhao 12/2010
!============================================================
    ! get aerosol number and mass for modal aerosols 
    do 200 j=jts,jte
    do 200 k=kts,kte
    do 200 i=its,ite

    so4_num(i,k,j,1:3)=0.0
    no3_num(i,k,j,1:3)=0.0
    nh4_num(i,k,j,1:3)=0.0
    oin_num(i,k,j,1:3)=0.0
    soot_num(i,k,j,1:3)=0.0
    organic_num(i,k,j,1:3)=0.0
    dst_num(i,k,j,1:3)=0.0
    so4_sfc(i,k,j,1:3)=0.0
    no3_sfc(i,k,j,1:3)=0.0
    nh4_sfc(i,k,j,1:3)=0.0
    oin_sfc(i,k,j,1:3)=0.0
    soot_sfc(i,k,j,1:3)=0.0
    organic_sfc(i,k,j,1:3)=0.0
    dst_sfc(i,k,j,1:3)=0.0

    nacon1d(1:4)=0. ! initialize 
    rndst1d(1:4)=1.e-8 ! 0.01 um to avoid divid by zero

    do n=1,ntype_aer
    do m=1,nsize_aer(n)
         lnum=numptr_aer(m,n,ai_phase)
         if(lnum>0)then
            scale0 = rho_phy(i,k,j) !air density (kg/m3)
            chem1d(lnum)=chem(i,k,j,lnum)*scale0 ! #/kg to #/m3
         endif
         do l=1,ncomp_aer(n)
            lmass=massptr_aer(l,m,n,ai_phase)
            scale0 = rho_phy(i,k,j)*1.e-9
            chem1d(lmass)=chem(i,k,j,lmass)*scale0 ! ug/kg to kg/m3
         enddo
    enddo ! size
    enddo ! type

    call get_sectional_aer(chem1d,num_chem, &
    so4_num(i,k,j,1:3),soot_num(i,k,j,1:3),organic_num(i,k,j,1:3),dst_num(i,k,j,1:3),no3_num(i,k,j,1:3),nh4_num(i,k,j,1:3),oin_num(i,k,j,1:3), & 
    so4_sfc(i,k,j,1:3),soot_sfc(i,k,j,1:3),organic_sfc(i,k,j,1:3),dst_sfc(i,k,j,1:3),no3_sfc(i,k,j,1:3),nh4_sfc(i,k,j,1:3),oin_sfc(i,k,j,1:3), &
    nacon1d,rndst1d,i,k,j)

    !dust for contact freezing
    rndst1(i,k,j)=rndst1d(1)
    rndst2(i,k,j)=rndst1d(2)
    rndst3(i,k,j)=rndst1d(3)
    rndst4(i,k,j)=rndst1d(4)
    nacon1(i,k,j)=nacon1d(1)
    nacon2(i,k,j)=nacon1d(2)
    nacon3(i,k,j)=nacon1d(3)
    nacon4(i,k,j)=nacon1d(4)

200 continue

    call icenucleat (qv,qc,qi,qs,qnice3d,t_phy,p_phy,w, &
           !czhao 
           nihf3d,niimm3d,nidep3d,nimey3d,  &
           ! sunny 
           rndst1, rndst2, rndst3, rndst4, nacon1, nacon2, nacon3, nacon4, & 
           so4_num,soot_num,organic_num,dst_num,no3_num,nh4_num,oin_num, &
           so4_sfc,soot_sfc,organic_sfc,dst_sfc,no3_sfc,nh4_sfc,oin_sfc, &
           tke_pbl,relhum3d,wsub3d,     &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old,id, ktau, dtstep, &
           f_qv,f_qc,f_qi,f_qs,f_ice_phy,f_rain_phy,    &
           config_flags%ra_sw_physics,config_flags%bl_pbl_physics,config_flags%wsubscheme,config_flags%ice_para,config_flags%mixcloud )

!-------------------------------------------------------------------------

!============================================================
!aerosol droplet activation and dry depostion
!============================================================
! check arguments of mixactivate for consistency between send, receive
! 06-nov-2005 rce - id & ktau added to arg list
      call mixactivate(  msectional, &
           !chem, num_chem, qv, qc, qi, qndrop3d,   &
           !t_phy, w, ddvel, idrydep_onoff,  &
           chem, num_chem, qv, qc, qi,qs, qndrop3d,   &
           t_phy,p_phy, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, dlo_sect, dhi_sect, sigmag_aer, dcen_sect,  &
           dens_aer, mw_aer,           &
           waterptr_aer, hygro,  ai_phase, cw_phase,                &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old, qsrflx, &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
	       id, ktau, dtstep, &
        f_qv,f_qc, f_qi,f_qs,f_ice_phy,f_rain_phy,config_flags%ra_sw_physics)

      end subroutine mosaic_mixactivate


!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine mosaic_mixactivate_init(                    &
           config_flags, chem, scalar,                       & 
           chem_in_opt,                                      & !BSINGH(12/04/13):  Added for SAPRC 8 bin vbs
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

      USE module_configure, only: grid_config_rec_type
      use module_state_description, only:  num_chem, num_scalar, p_qndrop
      use module_data_mosaic_asect

	implicit none

! subr arguments
      type(grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) ::               &
           ims, ime, jms, jme, kms, kme,   &
           its, ite, jts, jte, kts, kte
      INTEGER,      INTENT(IN   ) :: chem_in_opt !BSINGH(12/04/13):  Added for SAPRC 8 bin vbs
      real, intent(inout),   &
           dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
           chem

      real, intent(inout),   &
           dimension( ims:ime, kms:kme, jms:jme, 1:num_scalar ) :: &
           scalar

      integer :: i, j, k, m, n, l

      do j=jts,jte
         do k=kts,kte
            do i=its,ite
               scalar(i,k,j,p_qndrop) = 0.               
            end do
         end do
      end do

      if( cw_phase > 0 ) then   !sanity check in case using prognostic
         !drop number without aq. chemistry
         if(config_flags%chem_in_opt.ne.1) then !BSINGH(12/04/13):  Added for SAPRC 8 bin vbs
            
            do n=1,ntype_aer
               do m=1,nsize_aer(n)
                  chem(its:ite,kts:kte,jts:jte,numptr_aer(m,n,cw_phase)) = 0.
                  do l=1,ncomp_aer(n)
                     chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,cw_phase)) = 0.
                  end do              ! comp
               end do                 ! size
            end do                 ! type
         elseif (config_flags%chem_in_opt == 1) then !added by MS to transfer cloud borne species to interstitial species when using chem_in_opt=1
            
            do n=1,ntype_aer
               do m=1,nsize_aer(n)
                  chem(its:ite,kts:kte,jts:jte,numptr_aer(m,n,cw_phase)) = 0.
                  
                  
                  do l=1,ncomp_aer(n)
                     
                     if( ai_phase > 0 ) then
                        
                        ! add cloud borne aerosol mass to interstitial aerosol mass
                        
                        chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,ai_phase))= &
                             chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,ai_phase)) + &
                             chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,cw_phase))
                        
                     endif ! for ai_phase>0
                     
                     chem(its:ite,kts:kte,jts:jte,massptr_aer(l,m,n,cw_phase)) = 0.
                  end do              ! comp
               end do                 ! size
            end do                 ! type
            
         endif ! for chem_in_opt
      end if

      end subroutine mosaic_mixactivate_init


!----------------------------------------------------------------------
!----------------------------------------------------------------------


      subroutine sorgam_mixactivate (                        &
           id, ktau, dtstep, config_flags, idrydep_onoff,    &
           !rho_phy, t_phy, w, cldfra, cldfra_old,            &
           p_phy,rho_phy, t_phy, w, cldfra, cldfra_old,            &
           ddvel, z, dz8w, p_at_w, t_at_w, exch_h,           &
           !qv, qc, qi, qndrop3d, f_qc, f_qi, chem,           &
        !czhao
     qv,qc,qi,qs,f_ice_phy,f_rain_phy,qndrop3d,qnice3d,f_qv,f_qc,f_qi,f_qs,chem, &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
        !czhao 
        nihf3d,niimm3d,nidep3d,nimey3d,  &
        rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4,  &
        dgnuc,dgacc,dgcor,    &
        tke_pbl,relhum3d,wsub3d,   &
           ids,ide, jds,jde, kds,kde,                        &
           ims,ime, jms,jme, kms,kme,                        &
           its,ite, jts,jte, kts,kte                         )

    USE module_configure, only: grid_config_rec_type
	use module_state_description, only:  num_chem
	use module_data_sorgam
	use module_mixactivate, only: icenucleat, mixactivate

! wrapper to call mixactivate for sorgam description of aerosol

	implicit none

!   subr arguments
	integer, intent(in) ::                  &
		id, ktau,                       &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte,   &
                idrydep_onoff

	real, intent(in) :: dtstep

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		rho_phy, t_phy, w, p_phy,  &
		z, dz8w, p_at_w, t_at_w, exch_h

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: cldfra, cldfra_old

	real, intent(in),   &
		dimension( its:ite, jts:jte, num_chem ) :: ddvel

	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		qv, qc, qi,qs,f_ice_phy,f_rain_phy

    LOGICAL, intent(in) :: f_qv,f_qc, f_qi,f_qs

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme  ) :: &
		qndrop3d,qnice3d

	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) :: nsource, &
	     ccn1,ccn2,ccn3,ccn4,ccn5,ccn6  ! number conc of aerosols activated at supersat
      real, intent(inout), dimension(ims:ime,kms:kme,jms:jme) ::  &
        nihf3d,niimm3d,nidep3d,nimey3d
      real, intent(out), dimension(ims:ime,kms:kme,jms:jme) ::  &
       rndst1,rndst2,rndst3,rndst4,nacon1,nacon2,nacon3,nacon4  ! dust for contact freezing 
      real, intent(in), dimension(ims:ime,kms:kme,jms:jme) ::  &
       dgnuc,dgacc,dgcor 

      REAL, intent(out), DIMENSION( ims:ime, kms:kme, jms:jme ) :: &
       relhum3d,wsub3d

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
         INTENT(IN ) :: tke_pbl

	type(grid_config_rec_type), intent(in) :: config_flags

! local vars
	real qsrflx(ims:ime, jms:jme, num_chem) ! dry deposition flux of aerosol
	real sumhygro,sumvol
	integer i,j,k,l,m,n
	real hygro( its:ite, kts:kte, jts:jte,maxd_asize, maxd_atype )

! aerosol mass and number #/cm^3 or m2/cm^3
  real,dimension(ims:ime,kms:kme,jms:jme,1:3) ::  &
  so4_num,soot_num,organic_num,dst_num,no3_num,nh4_num,oin_num,so4_sfc,soot_sfc,organic_sfc,dst_sfc,no3_sfc,nh4_sfc,oin_sfc 
! real,dimension(1:3) :: &
! so4_num0,soot_num0,organic_num0,dst_num0,no3_num0,nh4_num0,oin_num0,so4_sfc0,soot_sfc0,organic_sfc0,dst_sfc0,no3_sfc0,nh4_sfc0,oin_sfc0 

  integer lnum,lmass
  real scale0
  real :: dgnum_m(maxd_asize, maxd_atype) ! median diameter (m) of number distrib of mode 
  real :: alogsig(maxd_asize, maxd_atype) ! natl log of geometric standard dev of aerosol
  real dpvolmean_m(maxd_asize, maxd_atype)  ! mean-volume diameter (m) of mode
  real dpsrfmean_m(maxd_asize, maxd_atype)  ! mean-surface diameter (m) of mode
  real chem1d(1:num_chem) ! aerosol number and mass mixing ratios (#/m3 or kg/m3) 
  REAL nacon1d(4), rndst1d(4)

! calculate volume-weighted bulk hygroscopicity for each type and size

      do 100 j=jts,jte
      do 100 k=kts,kte
      do 100 i=its,ite
       do n=1,ntype_aer
       do m=1,nsize_aer(n)
	       sumhygro=0
	       sumvol=0
	       do l=1,ncomp_aer(n)
	          sumhygro = sumhygro+hygro_aer(l,n)*   &
                   chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	          sumvol = sumvol+chem(i,k,j,massptr_aer(l,m,n,ai_phase))/dens_aer(l,n)
	       end do ! comp
               hygro(i,k,j,m,n)=sumhygro/sumvol
	end do ! size
	end do ! type
  100 continue

!============================================================
!aerosol ice nucleation  -czhao 12/2010
!============================================================
    ! get aerosol number and mass for modal aerosols 
    do 200 j=jts,jte
    do 200 k=kts,kte
    do 200 i=its,ite

    so4_num(i,k,j,1:3)=0.0
    no3_num(i,k,j,1:3)=0.0
    nh4_num(i,k,j,1:3)=0.0
    oin_num(i,k,j,1:3)=0.0
    soot_num(i,k,j,1:3)=0.0
    organic_num(i,k,j,1:3)=0.0
    dst_num(i,k,j,1:3)=0.0
    so4_sfc(i,k,j,1:3) = 0.0
    no3_sfc(i,k,j,1:3) = 0.0
    nh4_sfc(i,k,j,1:3) = 0.0
    oin_sfc(i,k,j,1:3) = 0.0
    soot_sfc(i,k,j,1:3) = 0.0
    organic_sfc(i,k,j,1:3) = 0.0
    dst_sfc(i,k,j,1:3) = 0.0

    nacon1d(1:4)=0. ! initialize 
    rndst1d(1:4)=1.e-8 ! 0.01 um to avoid divid by zero

    !--------------------------------------------------------
    !Get volume mean diameter (not volume medium diameter) 
    do n=1,ntype_aer
    do m=1,nsize_aer(n)
      alogsig(m,n)=alog(sigmag_aer(m,n))
      if (n==1.and.m==1) dgnum_m(m,n) = dgnuc(i,k,j)  ! Aitken mode
      if (n==1.and.m==2) dgnum_m(m,n) = dgacc(i,k,j)  ! Accumulaion mode 
      if (n==2.and.m==1) dgnum_m(m,n) = dgcor(i,k,j)  ! Coarse mode 
      dpvolmean_m(m,n) = dgnum_m(m,n)*exp(1.5*alogsig(m,n)*alogsig(m,n)) 
      dpsrfmean_m(m,n) = dgnum_m(m,n)*exp(alogsig(m,n)*alogsig(m,n) )
    enddo
    enddo

    do n=1,ntype_aer
    do m=1,nsize_aer(n)
         lnum=numptr_aer(m,n,ai_phase)
         if(lnum>0)then
            scale0 = rho_phy(i,k,j) !air density (kg/m3)
            chem1d(lnum)=chem(i,k,j,lnum)*scale0 ! #/kg to #/m3
         endif
         do l=1,ncomp_aer(n)
            lmass=massptr_aer(l,m,n,ai_phase)
            scale0 = rho_phy(i,k,j)*1.e-9
            chem1d(lmass)=chem(i,k,j,lmass)*scale0 ! ug/kg to kg/m3
         enddo
    enddo ! size
    enddo ! type

    call get_modal_aer(chem1d,num_chem,dpvolmean_m,dpsrfmean_m,dgnum_m, &
    so4_num(i,k,j,1:3),soot_num(i,k,j,1:3),organic_num(i,k,j,1:3),dst_num(i,k,j,1:3),no3_num(i,k,j,1:3),nh4_num(i,k,j,1:3),oin_num(i,k,j,1:3), &
    so4_sfc(i,k,j,1:3),soot_sfc(i,k,j,1:3),organic_sfc(i,k,j,1:3),dst_sfc(i,k,j,1:3),no3_sfc(i,k,j,1:3),nh4_sfc(i,k,j,1:3),oin_sfc(i,k,j,1:3), &
    nacon1d,rndst1d,i,k,j)

!   so4_num(i,k,j)=max(1.e-20,so4_num0)
!   soot_num(i,k,j)=max(1.e-20,soot_num0)
!   organic_num(i,k,j)=max(1.e-20,organic_num0)
!   dst_num(i,k,j)=max(1.e-20,dst_num0)
!   soot_sfc(i,k,j)=max(1.e-20,soot_sfc0)
!   organic_sfc(i,k,j)=max(1.e-20,organic_sfc0)
!   dst_sfc(i,k,j)=max(1.e-20,dst_sfc0)

!   so4_num(i,k,j)=so4_num0
!   soot_num(i,k,j)=soot_num0
!   organic_num(i,k,j)=organic_num0
!   dst_num(i,k,j)=dst_num0
!   soot_sfc(i,k,j)=soot_sfc0
!   organic_sfc(i,k,j)=organic_sfc0
!   dst_sfc(i,k,j)=dst_sfc0

    !dust for contact freezing
    rndst1(i,k,j)=rndst1d(1)
    rndst2(i,k,j)=rndst1d(2)
    rndst3(i,k,j)=rndst1d(3)
    rndst4(i,k,j)=rndst1d(4)
    nacon1(i,k,j)=nacon1d(1)
    nacon2(i,k,j)=nacon1d(2)
    nacon3(i,k,j)=nacon1d(3)
    nacon4(i,k,j)=nacon1d(4)

200 continue

    call icenucleat (qv,qc,qi,qs,qnice3d,t_phy,p_phy,w, &
           !czhao 
           nihf3d,niimm3d,nidep3d,nimey3d,  &
           ! sunny 
           rndst1, rndst2, rndst3, rndst4, nacon1, nacon2, nacon3, nacon4, &
           so4_num,soot_num,organic_num,dst_num,no3_num,nh4_num,oin_num, &
           so4_sfc,soot_sfc,organic_sfc,dst_sfc,no3_sfc,nh4_sfc,oin_sfc, &
           tke_pbl,relhum3d,wsub3d,    &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old,id, ktau, dtstep, &
           f_qv,f_qc,f_qi,f_qs,f_ice_phy,f_rain_phy,       &
           config_flags%ra_sw_physics,config_flags%bl_pbl_physics,config_flags%wsubscheme,config_flags%ice_para,config_flags%mixcloud)

!-------------------------------------------------------------------------

!============================================================
!aerosol droplet activation and dry depostion
!============================================================

! check arguments of mixactivate for consistency between send, receive
! 06-nov-2005 rce - id & ktau added to arg list
      call mixactivate(  msectional, &
           !chem, num_chem, qv, qc, qi, qndrop3d,   &
           !t_phy, w, ddvel, idrydep_onoff,  &
           chem, num_chem, qv, qc, qi,qs, qndrop3d,   &
           t_phy,p_phy, w, ddvel, idrydep_onoff,  &
           maxd_acomp, maxd_asize, maxd_atype, maxd_aphase,   &
           ncomp_aer, nsize_aer, ntype_aer, nphase_aer,  &
           numptr_aer, massptr_aer, dlo_sect, dhi_sect, sigmag_aer, dcen_sect,  &
           dens_aer, mw_aer,           &
           waterptr_aer, hygro,  ai_phase, cw_phase,                 &
           ids,ide, jds,jde, kds,kde,                            &
           ims,ime, jms,jme, kms,kme,                            &
           its,ite, jts,jte, kts,kte,                            &
           rho_phy, z, dz8w, p_at_w, t_at_w, exch_h,      &
           cldfra, cldfra_old, qsrflx,                      &
	       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
	       id, ktau, dtstep, &
      f_qv,f_qc,f_qi,f_qs,f_ice_phy,f_rain_phy,config_flags%ra_sw_physics  )

      end subroutine sorgam_mixactivate

!----------------------------------------------------------------------
  subroutine get_modal_aer(chem1d,num_chem,dpvolmean,dpsrfmean,dgnum_m, &
    so4_num0,soot_num0,organic_num0,dst_num0,no3_num0,nh4_num0,oin_num0, &
    so4_sfc0,soot_sfc0,organic_sfc0,dst_sfc0,no3_sfc0,nh4_sfc0,oin_sfc0, &
    nacon1d,rndst1d,iclm,kclm,jclm)

  USE module_data_sorgam !, only: dens_so4_aer,dens_oc_aer,dens_ec_aer,dens_dust_aer
  USE module_state_description, only:  param_first_scalar

  implicit none
!-----------------------------------------------------
! Input Arguments
!
  INTEGER, intent(in) ::  num_chem
  integer, intent(in) ::  iclm,jclm,kclm
  real,intent(in) :: chem1d(1:num_chem) ! aerosol number and mass mixing ratios (#/m3 or kg/m3) 
  real,intent(in) :: dpvolmean(maxd_asize, maxd_atype)  ! mean-volume diameter (m) of mode
  real,intent(in) :: dpsrfmean(maxd_asize, maxd_atype)  ! mean-surface diameter (m) of mode
  real,intent(in) :: dgnum_m(maxd_asize, maxd_atype) ! median diameter (m) of number distrib of mode 

!-----------------------------------------------------
! Output Arguments
!
  !3 size:total, and >0.1um, >0.5um
  real, intent(inout),dimension(1:3) :: so4_num0,soot_num0,organic_num0,dst_num0,no3_num0,nh4_num0,oin_num0, &
                            so4_sfc0,soot_sfc0,organic_sfc0,dst_sfc0,no3_sfc0,nh4_sfc0,oin_sfc0

  REAL, intent(inout), DIMENSION(1:4) ::  nacon1d, rndst1d

!-----------------------------------------------------
! Local Arguments
!
  integer isize,itype
  integer p1st
  integer n,l
  !Aiken
  real  mass_so4i , mass_oci ,  mass_bci  , mass_dusti,   &
        mass_aro1i, mass_aro2i, mass_alk1i, mass_ole1i,   &
        mass_ba1i , mass_ba2i,  mass_ba3i , mass_ba4i ,   &
        mass_pai,mass_no3i,mass_nh4i,mass_oini
  !Accumulation
  real  mass_so4j , mass_ocj ,  mass_bcj  , mass_dustj,   &
        mass_aro1j, mass_aro2j, mass_alk1j, mass_ole1j,   &
        mass_ba1j , mass_ba2j,  mass_ba3j , mass_ba4j ,   &
        mass_paj,mass_no3j,mass_nh4j,mass_oinj
  !Coarse      
  real  mass_soil,mass_anth
  real  vol_so4i   , vol_oci   ,  vol_bci    , vol_dusti, &
        vol_so4j   , vol_ocj   ,  vol_bcj    , vol_dustj, &
        vol_soil   , vol_no3i  ,  vol_no3j   , vol_nh4i, vol_nh4j,vol_oini,vol_oinj, vol_anth

  real  sfc_so4i   , sfc_oci   ,  sfc_bci    , sfc_dusti, &
        sfc_so4j   , sfc_ocj   ,  sfc_bcj    , sfc_dustj, &
        sfc_soil   , sfc_no3i  ,  sfc_no3j   , sfc_nh4i,  sfc_nh4j,sfc_oini,sfc_oinj, sfc_anth

  real  num_so4i   , num_oci   ,  num_bci    , num_dusti, &
        num_so4j   , num_ocj   ,  num_bcj    , num_dustj, &
        num_soil   , num_no3i, num_no3j,num_nh4i, num_nh4j,num_oini,num_oinj,num_anth

   real volume0,number0,surface0 
   real vol_so4   , vol_dust, vol_oc, vol_bc, vol_no3, vol_nh4 , vol_oin
   real sfc_so4   , sfc_dust, sfc_oc, sfc_bc, sfc_no3, sfc_nh4 , sfc_oin
   real num_so4   , num_dust, num_oc, num_bc, num_no3, num_nh4 , num_oin


  integer, parameter :: nbin_o=16
  real dcen_bin(nbin_o),dlo_bin(nbin_o),dhi_bin(nbin_o),dum  !in cm
  real dlo_um,dhi_um,duma,dgmin,dgnum_um
  real, dimension(1:nbin_o) :: xmas_secti,xmas_sectj,xmas_sectc
  real, dimension (1:3) :: IN_size
  integer :: in_bin, binsize


      IN_size(1) = 0.0
      IN_size(2) = 0.1*1.e-4  !0.1 um in cm
      IN_size(3) = 0.5*1.e-4  !0.5 um in cm

        dgmin=1.0e-03 ! in (um)
        duma=1.0
        dlo_um=0.0390625
        dhi_um=10.0

        !set bin size arrays for integration
        dlo_bin(1) = dlo_um*1.e-4  ! cm 
        dhi_bin(nbin_o) = dhi_um*1.e-4 ! cm 
        dum = alog(dhi_bin(nbin_o)/dlo_bin(1))/nbin_o
        do n = 2, nbin_o
          dlo_bin(n) = dlo_bin(1)*exp((n-1)*dum)
          dhi_bin(n-1) = dlo_bin(n)
        end do
        do n = 1, nbin_o 
          dcen_bin(n) = sqrt(dlo_bin(n)*dhi_bin(n))
        end do

        p1st = param_first_scalar

        ! Aitken mode...
        isize = 1 ; itype = 1
        l=lptr_so4_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_so4i= chem1d(l)
        l=lptr_no3_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_no3i= chem1d(l)
        l=lptr_nh4_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_nh4i= chem1d(l)
        l=lptr_dust_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_dusti= chem1d(l)
        l=lptr_orgaro1_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_aro1i= chem1d(l)
        l=lptr_orgaro2_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_aro2i= chem1d(l)
        l=lptr_orgalk_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_alk1i= chem1d(l)
        l=lptr_orgole_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ole1i= chem1d(l)
        l=lptr_orgba1_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba1i= chem1d(l)
        l=lptr_orgba2_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba2i= chem1d(l)
        l=lptr_orgba3_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba3i= chem1d(l)
        l=lptr_orgba4_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba4i= chem1d(l)
        l=lptr_orgpa_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_pai= chem1d(l)
        l=lptr_ec_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_bci= chem1d(l)
        l=lptr_p25_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_oini= chem1d(l)
        mass_oci=mass_aro1i+mass_aro2i+mass_alk1i+mass_ole1i+ &
                 mass_ba1i+mass_ba2i+mass_ba3i+mass_ba4i+mass_pai

        vol_so4i = mass_so4i/(dens_so4_aer*1.e3)   ! m3/m3
        vol_oci = mass_oci/(dens_oc_aer*1.e3)
        vol_bci = mass_bci/(dens_ec_aer*1.e3)
        vol_dusti = mass_dusti/(dens_dust_aer*1.e3)
        vol_no3i = mass_no3i/(dens_no3_aer*1.e3)   ! m3/m3
        vol_nh4i = mass_nh4i/(dens_nh4_aer*1.e3)   ! m3/m3
        vol_oini = mass_oini/(dens_oin_aer*1.e3)   ! m3/m3

        dgnum_um=amax1(dgmin,dgnum_m(isize,itype)*1.e6)
        call sect16bin(dgnum_um,sginin,duma,nbin_o,dlo_um,dhi_um,xmas_secti)

! Accumulation mode...
        isize = 2 ; itype = 1
        l=lptr_so4_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_so4j= chem1d(l)
        l=lptr_no3_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_no3j= chem1d(l)
        l=lptr_nh4_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_nh4j= chem1d(l)
        l=lptr_dust_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_dustj= chem1d(l)
        l=lptr_orgaro1_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_aro1j= chem1d(l)
        l=lptr_orgaro2_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_aro2j= chem1d(l)
        l=lptr_orgalk_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_alk1j= chem1d(l)
        l=lptr_orgole_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ole1j= chem1d(l)
        l=lptr_orgba1_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba1j= chem1d(l)
        l=lptr_orgba2_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba2j= chem1d(l)
        l=lptr_orgba3_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba3j= chem1d(l)
        l=lptr_orgba4_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_ba4j= chem1d(l)
        l=lptr_orgpa_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_paj= chem1d(l)
        l=lptr_ec_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_bcj= chem1d(l)
        l=lptr_p25_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_oinj= chem1d(l)

        mass_ocj=mass_aro1j+mass_aro2j+mass_alk1j+mass_ole1j+ &
                 mass_ba1j+mass_ba2j+mass_ba3j+mass_ba4j+mass_paj

        vol_so4j = mass_so4j/(dens_so4_aer*1.e3)  ! m3/m3
        vol_ocj = mass_ocj/(dens_oc_aer*1.e3)   ! m3/m3
        vol_bcj = mass_bcj/(dens_ec_aer*1.e3)   ! m3/m3
        vol_dustj = mass_dustj/(dens_dust_aer*1.e3)  !m3/m3
        vol_no3j = mass_no3j/(dens_no3_aer*1.e3)   ! m3/m3
        vol_nh4j = mass_nh4j/(dens_nh4_aer*1.e3)   ! m3/m3
        vol_oinj = mass_oinj/(dens_oin_aer*1.e3)   ! m3/m3

        dgnum_um=amax1(dgmin,dgnum_m(isize,itype)*1.e6)
        call sect16bin(dgnum_um,sginin,duma,nbin_o,dlo_um,dhi_um,xmas_sectj)


! Coarse mode...
        isize = 1 ; itype = 2   ! after  march-2008 ordering
        l=lptr_soil_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_soil= chem1d(l)
        l=lptr_anth_aer(isize,itype,ai_phase)
        if (l .ge. p1st)  mass_anth= chem1d(l)
        vol_soil = mass_soil/(dens_dust_aer*1.e3)   !m3/m3
        vol_anth = mass_soil/(dens_oin_aer*1.e3)   !m3/m3

        dgnum_um=amax1(dgmin,dgnum_m(isize,itype)*1.e6)
        call sect16bin(dgnum_um,sginin,duma,nbin_o,dlo_um,dhi_um,xmas_sectc)

      !for contact freezing
      nacon1d(1:4)=0.0
      rndst1d(1:4)=1.0

      do in_bin=1, 3 

      vol_so4=0.0
      vol_dust=0.0
      vol_oc=0.0
      vol_bc=0.0
      vol_no3=0.0
      vol_nh4=0.0
      vol_oin=0.0
      sfc_so4=0.0
      sfc_dust=0.0
      sfc_oc=0.0
      sfc_bc=0.0
      sfc_no3=0.0
      sfc_nh4=0.0
      sfc_oin=0.0
      num_so4=0.0
      num_dust=0.0
      num_oc=0.0
      num_bc=0.0
      num_no3=0.0
      num_nh4=0.0
      num_oin=0.0

      do isize = 1, nbin_o  ! 16 bins 

       if (dcen_bin(isize).ge.IN_size(in_bin)) then

          volume0 = vol_so4i*xmas_secti(isize)+vol_so4j*xmas_sectj(isize) 
          vol_so4 = vol_so4+volume0
          number0 = volume0*6/pirs/((dcen_bin(isize)*1.e-2)**3)  ! #/m3
          num_so4 = num_so4+number0
          surface0 = number0*pirs*((dcen_bin(isize)*1.e-2)**2)  ! m2/m3
          sfc_so4 = sfc_so4+surface0

          volume0 = vol_no3i*xmas_secti(isize)+vol_no3j*xmas_sectj(isize)
          vol_no3 = vol_no3+volume0
          number0 = volume0*6/pirs/((dcen_bin(isize)*1.e-2)**3)  ! #/m3
          num_no3 = num_no3+number0
          surface0 = number0*pirs*((dcen_bin(isize)*1.e-2)**2)  ! m2/m3
          sfc_no3 = sfc_no3+surface0

          volume0 = vol_nh4i*xmas_secti(isize)+vol_nh4j*xmas_sectj(isize) 
          vol_nh4 = vol_nh4+volume0
          number0 = volume0*6/pirs/((dcen_bin(isize)*1.e-2)**3)  ! #/m3
          num_nh4 = num_nh4+number0
          surface0 = number0*pirs*((dcen_bin(isize)*1.e-2)**2)  ! m2/m3
          sfc_nh4 = sfc_nh4+surface0

          volume0 = vol_oini*xmas_secti(isize)+vol_oinj*xmas_sectj(isize)+vol_anth*xmas_sectc(isize) 
          vol_oin = vol_oin+volume0
          number0 = volume0*6/pirs/((dcen_bin(isize)*1.e-2)**3)  ! #/m3
          num_oin = num_oin+number0
          surface0 = number0*pirs*((dcen_bin(isize)*1.e-2)**2)  ! m2/m3
          sfc_oin = sfc_oin+surface0

          volume0 = vol_oci*xmas_secti(isize)+vol_ocj*xmas_sectj(isize)
          vol_oc = vol_oc+volume0
          number0 = volume0*6/pirs/((dcen_bin(isize)*1.e-2)**3)  ! #/m3
          num_oc = num_oc+number0
          surface0 = number0*pirs*((dcen_bin(isize)*1.e-2)**2)  ! m2/m3
          sfc_oc = sfc_oc+surface0

          volume0 = vol_bci*xmas_secti(isize)+vol_bcj*xmas_sectj(isize) 
          vol_bc = vol_bc+volume0
          number0 = volume0*6/pirs/((dcen_bin(isize)*1.e-2)**3)  ! #/m3
          num_bc = num_bc+number0
          surface0 = number0*pirs*((dcen_bin(isize)*1.e-2)**2)  ! m2/m3
          sfc_bc = sfc_bc+surface0

          volume0 = vol_dusti*xmas_secti(isize)+vol_dustj*xmas_sectj(isize)+vol_soil*xmas_sectc(isize)
          vol_dust = vol_dust+volume0
          number0 = volume0*6/pirs/((dcen_bin(isize)*1.e-2)**3)  ! #/m3
          num_dust = num_dust+number0
          surface0 = number0*pirs*((dcen_bin(isize)*1.e-2)**2)  ! m2/m3
          sfc_dust = sfc_dust+surface0

          if (in_bin.eq.1) then !full size of dust 
          !--------------------------------------------------------------------------------
          ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
          ! estimate rndst and nacon for 4 dust bins here to pass to MG microphysics
          ! rndst: radius of 4 dust bins for contact freezing
          ! nacon: number in 4 dust bins for contact freezing 
          ! the estimate here has large uncertainty, can be improved in future
          ! For sectional aerosols: match size for dust 1-4 mode
          !--------------------------------------------------------------------------------
            binsize=int((isize+3)/4)
            nacon1d(binsize) = nacon1d(binsize)+number0  !#/m3
            rndst1d(binsize) = rndst1d(binsize)*(0.5*dcen_bin(isize)*1.e-2) !m
          !--------------------------------------------------------------------------------
          end if !in_bin

       endif !IN_size(in_bin)
      enddo !isize

        so4_num0(in_bin) = num_so4*1.e-6   ! #/cm3 
        no3_num0(in_bin) = num_no3*1.e-6   ! #/cm3 
        nh4_num0(in_bin) = num_nh4*1.e-6   ! #/cm3 
        oin_num0(in_bin) = num_oin*1.e-6   ! #/cm3 
        organic_num0(in_bin) = num_oc*1.e-6
        soot_num0(in_bin) = num_bc*1.e-6
        dst_num0(in_bin) = num_dust*1.e-6

        so4_sfc0(in_bin)=sfc_so4*1.e-6 !m2/cm3
        no3_sfc0(in_bin)=sfc_no3*1.e-6 !m2/cm3
        nh4_sfc0(in_bin)=sfc_nh4*1.e-6 !m2/cm3
        oin_sfc0(in_bin)=sfc_oin*1.e-6 !m2/cm3
        organic_sfc0(in_bin) = sfc_oc*1.e-6  ! m2/cm3 
        soot_sfc0(in_bin) = sfc_bc*1.e-6
        dst_sfc0(in_bin) = sfc_dust*1.e-6

        if (in_bin.eq.1) then !full size of dust
        rndst1d(1:4)=rndst1d(1:4)**0.25 !each 4-bin is integrated 
        endif

     enddo !in_bin=1, 3 

  return

  end subroutine get_modal_aer


!----------------------------------------------------------------------
  subroutine get_sectional_aer(chem1d,num_chem, &
    so4_num,soot_num,organic_num,dst_num,no3_num,nh4_num,oin_num, &
    so4_sfc,soot_sfc,organic_sfc,dst_sfc,no3_sfc,nh4_sfc,oin_sfc, &
    nacon1d,rndst1d,iclm,kclm,jclm)

  USE module_data_mosaic_asect
  USE module_data_mosaic_other
  USE module_state_description, only:  param_first_scalar

!-----------------------------------------------------
! Input Arguments
!
  INTEGER, intent(in) ::  num_chem
  INTEGER, intent(in) ::  iclm,kclm,jclm 
  real,intent(in) :: chem1d(1:num_chem) ! aerosol number and mass mixing ratios (#/m3 or kg/m3) 
!-----------------------------------------------------  
! Output Arguments
!
  !3 size:total, and >0.1um, >0.5um
  real, intent(inout), DIMENSION(1:3) :: so4_num,soot_num,organic_num,dst_num,no3_num,nh4_num,oin_num, &
                            so4_sfc,soot_sfc,organic_sfc,dst_sfc,no3_sfc,nh4_sfc,oin_sfc
  real, dimension (1:3) :: IN_size
  integer :: in_bin
  REAL, intent(inout), DIMENSION(1:4) ::  nacon1d, rndst1d

!-----------------------------------------------------
! Local Arguments
!
  integer isize,itype,iphase
  integer p1st
  integer l,l1,l2,l3,l4
  integer binsize
  real*8 pirs

   real  mass_so4  , mass_co3  , mass_ca   , mass_oc   , mass_bc   , &
         mass_aro1 , mass_aro2 , mass_alk1 , mass_ole1 , mass_api1 , &
         mass_api2 , mass_lim1 , mass_lim2 , mass_dust , mass_no3, mass_nh4,mass_oin
   real volume0,number0,surface0
   real  vol_so4   , vol_dust, vol_oc, vol_bc, vol_no3, vol_nh4 , vol_oin 
   real  sfc_so4   , sfc_dust, sfc_oc, sfc_bc, sfc_no3, sfc_nh4 , sfc_oin
   real  num_so4   , num_dust, num_oc, num_bc, num_no3, num_nh4 , num_oin

      IN_size(1) = 0.0
      IN_size(2) = 0.1*1.e-4  !0.1 um in cm
      IN_size(3) = 0.5*1.e-4  !0.5 um in cm

      pirs = 3.14159265358979324
      p1st = param_first_scalar
      itype = 1 ; iphase = 1

      !for contact freezing
      nacon1d(1:4)=0.0
      rndst1d(1:4)=1.0

      do in_bin=1, 3 

      vol_so4=0.0 
      vol_dust=0.0
      vol_oc=0.0
      vol_bc=0.0
      vol_no3=0.0
      vol_nh4=0.0
      vol_oin=0.0
      sfc_so4=0.0 
      sfc_dust=0.0
      sfc_oc=0.0
      sfc_bc=0.0
      sfc_no3=0.0
      sfc_nh4=0.0
      sfc_oin=0.0
      num_so4=0.0 
      num_dust=0.0
      num_oc=0.0
      num_bc=0.0
      num_no3=0.0
      num_nh4=0.0
      num_oin=0.0

      do isize = 1, nsize_aer(itype) 
       
       if (dcen_sect(isize,itype).ge.IN_size(in_bin)) then         
 
          mass_so4=0.0
          l=lptr_so4_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_so4= chem1d(l)
          mass_dust=0.0
          l=lptr_dust_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_dust=mass_dust+chem1d(l)
          l1=lptr_tr1r1_aer(isize,itype,iphase)
          if (l1 .ge. p1st) mass_dust=mass_dust+chem1d(l1) 
          l2=lptr_tr1r2_aer(isize,itype,iphase)
          if (l2 .ge. p1st) mass_dust=mass_dust+chem1d(l2) 
          l3=lptr_tr1r3_aer(isize,itype,iphase)
          if (l3 .ge. p1st) mass_dust=mass_dust+chem1d(l3) 
          l4=lptr_tr1r4_aer(isize,itype,iphase)
          if (l4 .ge. p1st) mass_dust=mass_dust+chem1d(l4) 
          mass_oc=0.0
          l=lptr_oc_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_oc= chem1d(l)
          mass_bc=0.0
          l=lptr_bc_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_bc= chem1d(l)
          mass_co3=0.0
          l=lptr_co3_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_co3= chem1d(l)
          mass_ca=0.0
          l=lptr_ca_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_ca= chem1d(l)
          mass_no3=0.0
          l=lptr_no3_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_no3= chem1d(l)
          mass_nh4=0.0
          l=lptr_nh4_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_nh4= chem1d(l)
          mass_oin=0.0
          l=lptr_oin_aer(isize,itype,iphase)
          if (l .ge. p1st)  mass_oin= chem1d(l)


          volume0 = mass_so4/(dens_so4_aer*1.e3)   ! m3/m3
          vol_so4 = vol_so4+volume0
          number0 = volume0*6/pirs/((dcen_sect(isize,itype)*1.e-2)**3)  ! #/m3
          num_so4 = num_so4+number0
          surface0 = number0*pirs*((dcen_sect(isize,itype)*1.e-2)**2)  ! m2/m3
          sfc_so4 = sfc_so4+surface0

          volume0 = mass_no3/(dens_no3_aer*1.e3)   ! m3/m3
          vol_no3 = vol_no3+volume0
          number0 = volume0*6/pirs/((dcen_sect(isize,itype)*1.e-2)**3)  ! #/m3
          num_no3 = num_no3+number0
          surface0 = number0*pirs*((dcen_sect(isize,itype)*1.e-2)**2)  ! m2/m3
          sfc_no3 = sfc_no3+surface0

          volume0 = mass_nh4/(dens_nh4_aer*1.e3)   ! m3/m3
          vol_nh4 = vol_nh4+volume0
          number0 = volume0*6/pirs/((dcen_sect(isize,itype)*1.e-2)**3)  ! #/m3
          num_nh4 = num_nh4+number0
          surface0 = number0*pirs*((dcen_sect(isize,itype)*1.e-2)**2)  ! m2/m3
          sfc_nh4 = sfc_nh4+surface0

          volume0 = mass_oin/(dens_oin_aer*1.e3)   ! m3/m3
          vol_oin = vol_oin+volume0
          number0 = volume0*6/pirs/((dcen_sect(isize,itype)*1.e-2)**3)  ! #/m3
          num_oin = num_oin+number0
          surface0 = number0*pirs*((dcen_sect(isize,itype)*1.e-2)**2)  ! m2/m3
          sfc_oin = sfc_oin+surface0

          volume0 = mass_oc/(dens_oc_aer*1.e3)
          vol_oc = vol_oc+volume0
          number0 = volume0*6/pirs/((dcen_sect(isize,itype)*1.e-2)**3)  ! #/m3
          num_oc = num_oc+number0
          surface0 = number0*pirs*((dcen_sect(isize,itype)*1.e-2)**2)  ! m2/m3
          sfc_oc = sfc_oc+surface0

          volume0 = mass_bc/(dens_bc_aer*1.e3)
          vol_bc = vol_bc+volume0
          number0 = volume0*6/pirs/((dcen_sect(isize,itype)*1.e-2)**3)  ! #/m3
          num_bc = num_bc+number0
          surface0 = number0*pirs*((dcen_sect(isize,itype)*1.e-2)**2)  ! m2/m3
          sfc_bc = sfc_bc+surface0

          volume0 = (mass_dust+mass_co3+mass_ca)/(dens_dust_aer*1.e3)
          vol_dust = vol_dust+volume0
          number0 = volume0*6/pirs/((dcen_sect(isize,itype)*1.e-2)**3)  ! #/m3
          num_dust = num_dust+number0
          surface0 = number0*pirs*((dcen_sect(isize,itype)*1.e-2)**2)  ! m2/m3
          sfc_dust = sfc_dust+surface0

          if (in_bin.eq.1) then !full size of dust 
          !--------------------------------------------------------------------------------
          ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
          ! estimate rndst and nacon for 4 dust bins here to pass to MG microphysics
          ! rndst: radius of 4 dust bins for contact freezing
          ! nacon: number in 4 dust bins for contact freezing 
          ! the estimate here has large uncertainty, can be improved in future
          ! For sectional aerosols: match size for dust 1-4 mode
          !--------------------------------------------------------------------------------
          if (nsize_aer(itype).eq.4) then
            nacon1d(isize) = number0  !#/m3
            rndst1d(isize) = 0.5*dcen_sect(isize,itype)*1.e-2 !radius in m 
          else if (nsize_aer(itype).eq.8) then
                 binsize=int((isize+1)/2)
                 nacon1d(binsize) = nacon1d(binsize)+number0  !#/m3
                 rndst1d(binsize) = rndst1d(binsize)*0.5*dcen_sect(isize,itype)*1.e-2 
          end if
          !--------------------------------------------------------------------------------
          end if !in_bin

       endif !IN_size(in_bin)
     enddo !isize

        so4_num(in_bin) = num_so4*1.e-6   ! #/cm3 
        no3_num(in_bin) = num_no3*1.e-6   ! #/cm3 
        nh4_num(in_bin) = num_nh4*1.e-6   ! #/cm3 
        oin_num(in_bin) = num_oin*1.e-6   ! #/cm3 
        organic_num(in_bin) = num_oc*1.e-6
        soot_num(in_bin) = num_bc*1.e-6
        dst_num(in_bin) = num_dust*1.e-6

        so4_sfc(in_bin)=sfc_so4*1.e-6 !m2/cm3
        no3_sfc(in_bin)=sfc_no3*1.e-6 !m2/cm3
        nh4_sfc(in_bin)=sfc_nh4*1.e-6 !m2/cm3
        oin_sfc(in_bin)=sfc_oin*1.e-6 !m2/cm3
        organic_sfc(in_bin) = sfc_oc*1.e-6  ! m2/cm3 
        soot_sfc(in_bin) = sfc_bc*1.e-6
        dst_sfc(in_bin) = sfc_dust*1.e-6
      
        if (in_bin.eq.1) then !full size of dust
        if (nsize_aer(itype).eq.8) then
        rndst1d(1:4)=sqrt(rndst1d(1:4))
        endif
        endif

     enddo !in_bin=1, 3 

!    write(0,*) 'czhao check mosaic size ', nsize_aer(itype),dcen_sect(:,1)
!    print*, 'czhao check mosaic size ', nsize_aer(itype),dcen_sect(:,1)

  return

  end subroutine get_sectional_aer

    subroutine sect16bin(dgnum_um,sigmag,duma,nbin,dlo_um,dhi_um,xmas_sect)
!
        implicit none
        REAL, DIMENSION(nbin), INTENT(OUT) :: xmas_sect
        integer n, nbin
        real &
          dgnum, dgnum_um, dhi, dhi_um, dlo, dlo_um,   &
          dstar, duma, dumfrac, dx,   & 
          sigmag,summas,   &
          sx, sxroot2, thi, tlo, vtot,   &
          x0, x3, xhi, xlo, xmtot
        real dlo_sect(nbin), dhi_sect(nbin)
        real pi
        parameter (pi = 3.1415926536)
!
        xmtot = duma
!   compute section boundaries
        dlo = dlo_um*1.0e-4
        dhi = dhi_um*1.0e-4
        xlo = log( dlo )
        xhi = log( dhi )
        dx = (xhi - xlo)/nbin
        do n = 1, nbin
            dlo_sect(n) = exp( xlo + dx*(n-1) )
            dhi_sect(n) = exp( xlo + dx*n )
        end do
!   compute modal "working" parameters including total num/vol/mass
        dgnum = dgnum_um*1.0e-4
        sx = alog( sigmag )
        x0 = alog( dgnum )
        x3 = x0 + 3.*sx*sx
        dstar = dgnum * exp(1.5*sx*sx)
!   compute number and mass for each section
        sxroot2 = sx * sqrt( 2.0 )
        summas = 0.
        do n = 1, nbin
            xlo = alog( dlo_sect(n) )
            xhi = alog( dhi_sect(n) )
            tlo = (xlo - x3)/sxroot2
            thi = (xhi - x3)/sxroot2
            if (tlo .le. 0.) then
                dumfrac = 0.5*( erfc_num_recipes(-thi) - erfc_num_recipes(-tlo) )
            else
                dumfrac = 0.5*( erfc_num_recipes(tlo) - erfc_num_recipes(thi) )
            end if
            xmas_sect(n) = xmtot*dumfrac
            summas = summas + xmas_sect(n)
        end do

      end subroutine  sect16bin

      real function erfc_num_recipes( x )
!   from press et al, numerical recipes, 1990, page 164
        implicit none
        real x
        double precision erfc_dbl, dum, t, z
        z = abs(x)
        t = 1.0/(1.0 + 0.5*z)
        dum =  ( -z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 +   &
          t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +   &
                                           t*(-1.13520398 +   &
          t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))
        erfc_dbl = t * exp(dum)
        if (x .lt. 0.0) erfc_dbl = 2.0d0 - erfc_dbl
        erfc_num_recipes = erfc_dbl
        return 
      end function erfc_num_recipes


END MODULE module_mixactivate_wrappers
