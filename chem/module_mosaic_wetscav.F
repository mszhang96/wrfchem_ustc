MODULE module_mosaic_wetscav

CONTAINS

!===========================================================================
!===========================================================================
   subroutine wetscav_cbmz_mosaic (id,ktau,dtstep,ktauc,config_flags,      &
               dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,        &
	       qlsink,ilsink,r2evap,i2evap,precr,preci,precs,precg,        &
	       gas_aqfrac, numgas_aqfrac,                                  &
               qv,qc,qi,qs,   &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )

!  wet removal by grid-resolved precipitation
!  scavenging of cloud-phase aerosols and gases by collection, freezing, ...
!  scavenging of interstitial-phase aerosols by impaction
!  scavenging of gas-phase gases by mass transfer and reaction

!----------------------------------------------------------------------
   USE module_configure, only: grid_config_rec_type
   USE module_state_description
   USE module_data_mosaic_asect

!----------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac
   REAL, INTENT(IN   ) :: dtstep,dtstepc
!
! all advected chemical species
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem

! fraction of gas species in cloud water
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac

!
!
! input from meteorology
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
         INTENT(IN   ) ::                                           &
                                                        alt,        &
                                                      t_phy,        &
                                                      p_phy,        &
                                                   t8w,p8w,         &
                            qlsink,ilsink,precr,preci,precs,precg, &
                            r2evap,i2evap, &
                                                    rho_phy,cldfra
        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                qv, qc, qi,qs

   call wetscav (id,ktau,dtstep,ktauc,config_flags,                      &
        dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,             &
        qlsink,ilsink,precr,preci,precs,precg,                   &
        r2evap,i2evap,                                                   &
        gas_aqfrac, numgas_aqfrac,                                       &
        ntype_aer, nsize_aer, ncomp_aer,                                 &
        massptr_aer, dens_aer, numptr_aer,hygro_aer,                     &
        maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase, cw_phase, &
        volumcen_sect, volumlo_sect, volumhi_sect,                       &
        waterptr_aer, dens_water_aer,                                    &
        scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
        qv,qc,qi,qs,   &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte                                        )

   end subroutine wetscav_cbmz_mosaic



!===========================================================================
!===========================================================================
   subroutine wetscav (id,ktau,dtstep,ktauc,config_flags,                &
        dtstepc,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,             &
        qlsink,ilsink,precr,preci,precs,precg,                    &
        r2evap,i2evap,                             &
        gas_aqfrac, numgas_aqfrac,                                       &
        ntype_aer, nsize_aer, ncomp_aer,                                 &
        massptr_aer, dens_aer, numptr_aer,hygro_aer,                     &
        maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase, cw_phase, &
        volumcen_sect, volumlo_sect, volumhi_sect,                       &
        waterptr_aer, dens_water_aer,                                    &
        scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
        qv,qc,qi,qs,   &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte                                        )

!  wet removal by grid-resolved precipitation
!  scavenging of cloud-phase aerosols and gases by collection, freezing, ...
!  scavenging of interstitial-phase aerosols by impaction
!  scavenging of gas-phase gases by mass transfer and reaction

!----------------------------------------------------------------------
!Updates and Notes by Chun Zhao,PNNL
!1. update the current scheme working without aerosol-cloud-interaction/cloud-born-aerosol, only change the in-cloud scavenging part.
!2. Now the below cloud scavenging only accounts for large-scale rain
!3. The convective updraft scavening is dealt during the convective transport
!4. adding the new wet-scavenging scheme following GEOS-CHEM for in-cloud (rainout) scavenging and below-cloud (washout) scavenging in the cases without online aerosol-cloud interaction/aerosol activation
!----------------------------------------------------------------------

!----------------------------------------------------------------------
   USE module_configure, only: grid_config_rec_type
   USE module_state_description
   USE module_model_constants, only: g, rhowater, mwdry
   USE module_dep_simple, only: is_aerosol

   IMPLICIT NONE

!======================================================================

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac
      REAL,      INTENT(IN   ) :: dtstep,dtstepc
!
! all advected chemical species
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem

! fraction of gas species in cloud water
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac

!
! input from meteorology
!
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
          INTENT(IN   ) ::                                          &
                                                        alt,        &
                                                      t_phy,        &
                                                      p_phy,        &
                                                    t8w,p8w,        &
                             qlsink,ilsink,precr,preci,precs,precg,        &
                             r2evap,i2evap,        &
                                             rho_phy,cldfra
   integer, intent(in) :: maxd_atype, maxd_asize, maxd_acomp, maxd_aphase
   integer, intent(in) :: ai_phase, cw_phase
   integer, intent(in) :: ntype_aer
   integer, intent(in) :: nsize_aer( maxd_atype ),   & ! number of size bins
      	  ncomp_aer( maxd_atype ),   & ! number of chemical components
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase), & ! index for mixing ratio
      	  waterptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol water
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) ! index for the number mixing ratio
   real, intent(in) :: dens_aer( maxd_acomp, maxd_atype ),   & ! material density (g/cm3)
	  dens_water_aer, &  ! water density (g/m3)
          hygro_aer( maxd_acomp, maxd_atype )
   real, intent(in) :: volumcen_sect( maxd_asize, maxd_atype ),   & ! single-particle volumes (cm3)
                        volumlo_sect( maxd_asize, maxd_atype ),   & ! at section center, lower boundary,
                        volumhi_sect( maxd_asize, maxd_atype )      ! upper boundary

   real, intent(in) :: dlndg_nimptblgrow
   integer, intent(in) :: nimptblgrow_mind, nimptblgrow_maxd
   real, intent(in) :: scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype), &
     	     scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)
        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                qv, qc, qi,qs


!
! LOCAL  VAR
!
     integer :: i,j,k,l,m,n
     integer :: lmasscw,lnumcw
     real scale,scale1,scale2
   logical :: isprx(ims:ime,kms:kme,jms:jme)

   real :: pdel(ims:ime,kms:kme,jms:jme)
   real :: delpf(kms:kme)
   real :: delpfhalf
   real :: dump
   real :: fac_pwght_ls(kms:kme)   !
   real :: fapincld, fapoutcld, faptot
   real :: fapmin_ls     !
   real :: fapx(ims:ime,kms:kme,jms:jme)
   real :: fracscav
   real :: pf_above, pf_below, pdel_fac
   real :: pf_ls(kms:kme)
   real :: pfoutcld
   real :: pfsmall_ls    ! l-s precip fluxes (kg/m2/s) smaller than this
                             ! are ignored (treated as zero)
   real :: pfsmall_min
   real :: pfx(ims:ime,kms:kme,jms:jme)
   real :: pfx_inrain(ims:ime,kms:kme,jms:jme)
   real :: sumfa, sumpf, sumpffa
   REAL :: dqdt( ims:ime, kms:kme, jms:jme, num_chem )
   logical dotend(num_chem)


   dotend(1:num_chem) = .false.
   do j=jts,jte
   do i=its,ite
   do k=kts,kte
     dqdt(i,k,j,1:num_chem)=0.0
   enddo
   end do
   end do

!******************************************
!1. in-cloud scavenging
!******************************************

   !scavenging cloud-phase aerosol assuming precip falls to surface
   incloud_select: SELECT CASE(config_flags%chem_opt)   

   !-----------------------------------------
   !1.1 in-cloud scavenging for cases with cloud-born aerosol and 
   !    aqueous phase chemistry
   !-----------------------------------------
   CASE (RADM2SORG_AQ,RACMSORG_AQ,CBMZSORG_AQ, & 
         CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
         CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, &
         CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
         CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE, CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE )

    if (config_flags%progn == 0 .or. config_flags%cldchem_onoff == 0) then 
      CALL wrf_error_fatal('Wet scavenging with cloud-born-aerosol must work with aerosol-indirect effect/online aerosol activation and cloud chemistry')
    endif

   !-----------------------------------------
   !1.1.1 in-cloud scavenging of aerosols 
   !-----------------------------------------
      do 100 j=jts,jte
        do k=kts,kte
         do i=its,ite
            pdel(i,k,j)=p8w(i,k,j)-p8w(i,k+1,j)
         end do
        end do

      do 100 k=kts,kte
      do 100 i=its,ite
         scale=1.0-dtstepc*qlsink(i,k,j)
         scale=max(0.0,min(1.0,scale)) ! make sure 0 <= scale <= 1
         if (qlsink(i,k,j) >  0.0) then
            pdel_fac = pdel(i,k,j)/g
            do n=1,ntype_aer
               do m=1,nsize_aer(n)
                  do l=1,ncomp_aer(n)
                     lmasscw=massptr_aer(l,m,n,cw_phase)
                     if (lmasscw < param_first_scalar) cycle
!czhao                     qsrflx(i,j,lmasscw)=qsrflx(i,j,lmasscw)+chem(i,k,j,lmasscw)*(scale-1.)*pdel_fac ! aerosol mass (ug/m2)
                     chem(i,k,j,lmasscw)=chem(i,k,j,lmasscw)*scale
                  end do ! comp
                  lnumcw=numptr_aer(m,n,cw_phase)
                  if (lnumcw < param_first_scalar) cycle
!czhao                  qsrflx(i,j,lnumcw)=qsrflx(i,j,lnumcw)+chem(i,k,j,lnumcw)*(scale-1.)*pdel_fac ! aerosol number (1/m2)
                  chem(i,k,j,lnumcw)=chem(i,k,j,lnumcw)*scale
               end do ! size
            end do ! type
         end if    ! qlsink > 0
  100 continue

   !-----------------------------------------
   !1.1.2 in-cloud scavenging of gases 
   !-----------------------------------------
      do 290 l = param_first_scalar, min( num_chem, numgas_aqfrac )
      if ( is_aerosol(l) ) goto 290
      do 270 j = jts,jte
      do 270 k = kts,kte
      do 270 i = its,ite
         fracscav = dtstepc*qlsink(i,k,j)*gas_aqfrac(i,k,j,l)
         if (fracscav > 0.0) then ! make sure fracscav > 0
            fracscav = min(1.0,fracscav) ! make sure fracscav <= 1
            scale = 1.0 - fracscav
            pdel_fac = pdel(i,k,j)/(g*mwdry)
!czhao            qsrflx(i,j,l) = qsrflx(i,j,l)+chem(i,k,j,l)*(scale-1.)*pdel_fac ! mmol/m2
            chem(i,k,j,l) = chem(i,k,j,l)*scale
         end if
 270   continue
 290   continue

   !-----------------------------------------
   !1.2 in-cloud scavenging for cases without cloud-born aerosol and
   !    aqueous phase chemistry
   !    in this case, we should still call cloud phase gas chemistry
   !    to get the scavenging coefficient (gas_aqfrac). 
   !    This should be done, but not yet, now gas_aqfrac is re-calculated for soluble gases
   !-----------------------------------------
   CASE (RADM2SORG,CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN)

    if (config_flags%progn == 1.or.config_flags%cldchem_onoff == 1) then
      CALL wrf_error_fatal('Wet scavenging without cloud-born-aerosol must work without aerosol-indirect effect/online aerosol activation and cloud chemistry')
    endif

    if (config_flags%mp_physics/=MORR_TWO_MOMENT) then 
     CALL wrf_error_fatal('In-cloud wet scavenging without cloud-born-aerosol is currently only coupled with morrison mp_physics scheme')
    endif

        do j=jts,jte
        do k=kts,kte
        do i=its,ite
            pdel(i,k,j)=p8w(i,k,j)-p8w(i,k+1,j)
        enddo
        enddo
        enddo

    call wetscav_incloud(id,ktau,dtstep,ktauc,config_flags,           &
        dtstepc,pdel,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,             &
        qlsink,ilsink,                    &
        r2evap,i2evap,                    &
        gas_aqfrac, numgas_aqfrac,                                       &
        ntype_aer, nsize_aer, ncomp_aer,                                 &
        massptr_aer, dens_aer, numptr_aer, hygro_aer,                    &
        maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase,        &
        volumcen_sect, volumlo_sect, volumhi_sect,                       &
        waterptr_aer, dens_water_aer,                                    &
        qv,qc,qi,qs,      &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte                                        )


   CASE DEFAULT

   END SELECT incloud_select


!******************************************
!2. below-cloud scavenging
!******************************************

! precip rates -- 1.0 kgwtr/m2/s = 1.0e-3 m3wtr/m2/s = 1.0e-3 m/s
!     7.06e-5 kg/m2/s = 7.06e-8 m/s = 0.01 inch/h
!     3.17e-5 kg/m2/s = 3.17e-8 m/s = 1 m/y = global annual average
!     1.00e-7 kg/m2/s = 1.00e-10 m/s is a very small precip rate!

   fapmin_ls = 1.0e-3
   pfsmall_ls = 1.0e-7

    do 5900 j=jts,jte
    do 5900 i=its,ite

      do k = kts,kte
        isprx(i,k,j) = .false.
        pfx(i,k,j) = 0.0
        pfx_inrain(i,k,j) = 0.0
        fapx(i,k,j) = 0.0
      enddo

!----------------------------------------
!2.1 compute l-s precip rates
!----------------------------------------
! pf_ls(k) = precip flux at center of level
      pf_below = 0.0
      do k = kte,kts,-1
         pf_above = pf_below
         !pf_below = precr(i,k,j) + preci(i,k,j) + precs(i,k,j) + precg(i,k,j) ! total precip rate
         pf_below = precr(i,k,j) + preci(i,k,j) + precs(i,k,j) !problem with precg values,not used for now
         if (pf_below < pfsmall_ls) pf_below = 0.0
         delpf(k) = pf_below - pf_above
         pf_ls(k) = 0.5*(pf_below + pf_above)
         if (pf_ls(k) < pfsmall_ls) pf_ls(k) = 0.0
      end do

! compute fac_pwght_ls which is an average of cloud fractional area in and
!    above the current level, weighted by precip production in each level
! basically this reflect the cloud area at levels where precip is produced
      do k = kte, kts,-1
         if (k == kte) then
!     compute change from (k-1/k) interface to level k center
            delpfhalf = 0.5*delpf(k)
            if (delpfhalf > 0.0) then
               fac_pwght_ls(k) = max( cldfra(i,k,j), fapmin_ls )
               sumpffa = delpfhalf * fac_pwght_ls(k)
               sumpf = delpfhalf
            else
               fac_pwght_ls(k) = fapmin_ls
               sumpffa = 0.0
               sumpf = 0.0
	    end if
         else
!     compute change from level (k+1) center to (k+1/k) interface
            delpfhalf = 0.5*delpf(k+1)
            if (delpfhalf > 0.0) then
               sumpffa = sumpffa + delpfhalf*max( cldfra(i,k+1,j), fapmin_ls )
               sumpf = sumpf + delpfhalf
               fac_pwght_ls(k) = max( (sumpffa/sumpf), fapmin_ls )
            else
               fac_pwght_ls(k) = fac_pwght_ls(k+1)
               sumpffa = max( (sumpffa + delpfhalf*fac_pwght_ls(k)), 0.0 )
               sumpf = max( (sumpf + delpfhalf), 0.0 )
	    end if
!     compute change from (k-1/k) interface to level k center
            delpfhalf = 0.5*delpf(k)
            if (delpfhalf > 0.0) then
               sumpffa = sumpffa + delpfhalf*max( cldfra(i,k,j), fapmin_ls )
               sumpf = sumpf + delpfhalf
               fac_pwght_ls(k) = max( (sumpffa/sumpf), fapmin_ls )
            else
!              here, fac_pwght_ls(k) is unchanged from its value computed just above
               sumpffa = max( (sumpffa + delpfhalf*fac_pwght_ls(k)), 0.0 )
               sumpf = max( (sumpf + delpfhalf), 0.0 )
	    end if
	 end if
      end do

! loop through levels
      do 4900 k = kte, kts,-1
!----------------------------------------------------------------------
! for each cloud/precip type (ls, dp, sh), compute
!   fapx = fractional area with precip
!   pfx = precip flux based on entire grid-cell area (kg/m2/s)
!   pfx_inrain = precip flux within the precip fractional area (kg/m2/s)
!
      sumpf = 0.0
      sumfa = 0.0

! l-s cloud
! assume  faptot = total (in + out-of-cloud) precip area
!                = 0.5*fac_pwght_ls(k)
! then  fapincld = in-cloud precip area
!                = min( faptot, cloud area )
!      fapoutcld = out-of-cloud precip area
!                = max( 0.0, faptot-fapincld )
! also  pfoutcld = out-of-cloud precip flux = fraction of total precip flux
!                = pf_ls(k)*(fapoutcld/faptot)
      if (pf_ls(k) > 0.0) then
         faptot= 0.5*fac_pwght_ls(k)
         fapincld = min( faptot, cldfra(i,k,j) )
         fapoutcld = max( 0.0, faptot-fapincld )
         pfoutcld = pf_ls(k)*(fapoutcld/faptot)

         if (pfoutcld >= pfsmall_ls) then
            isprx(i,k,j) = .true.
            pfx(i,k,j) = pfoutcld
            fapx(i,k,j) = fapoutcld
            if (faptot.le.0.or.faptot.ge.1.e10) then
             write(0,*) 'check faptot ',k,faptot
             print*, 'check faptot ',k,faptot
             stop
            endif
            pfx_inrain(i,k,j) = pf_ls(k)/faptot
            if (pfx_inrain(i,k,j).ge.1.e20) then
             write(0,*) 'check pfx_inrain ',i,j,k,pfx_inrain(i,k,j),pf_ls(k),faptot
             print*, 'check pfx_inrain ',i,j,k,pfx_inrain(i,k,j),pf_ls(k),faptot
             stop
            endif
            sumpf = sumpf + pfx(i,k,j)
            sumfa = sumfa + fapx(i,k,j)
         end if
      end if

4900 continue  ! "k = 1, pver"

5900 continue  ! "i = 1, ncol", j

!----------------------------------------
!2.2 do below cloud scavenging
!----------------------------------------
      call aerimpactscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte, num_chem,       &
           ntype_aer, nsize_aer, ncomp_aer, massptr_aer, dens_aer, numptr_aer,  &
           maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase,          &
           volumcen_sect, volumlo_sect, volumhi_sect, &
           waterptr_aer, dens_water_aer,                &
           scavimptblvol, scavimptblnum,      &
           nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow,     &
           dtstepc, t_phy, p_phy, pdel, chem,          &
           isprx, fapx, pfx, pfx_inrain,           &
           dqdt, dotend      )

      call gasrainscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte, num_chem,     &
                      config_flags,      &
                      dtstepc,    t_phy,      p_phy,      pdel,  chem,      &
                      isprx,      fapx,       pfx,        pfx_inrain, &
                      dqdt,       dotend      )

!----------------------------------------
!2.3 update chemistry for below-cloud scavenging 
!----------------------------------------
   do n=1,num_chem
      if(dotend(n))then
         do 6000 j=jts,jte
         do 6000 k=kts,kte
         do 6000 i=its,ite
            !czhao found dqdt(1,1,1,2) undefined
            chem(i,k,j,n) = chem(i,k,j,n) + dqdt(i,k,j,n)*dtstepc
 6000    continue
      end if
   end do


!===========================================================================
!additional option of large scale wet scavenging following GEOS-CHEM
!     call wetscav_largescale (
!       ids,ide, jds,jde, kds,kde,                                       &
!       ims,ime, jms,jme, kms,kme,                                       &
!       its,ite, jts,jte, kts,kte                                        )
!===========================================================================

   end subroutine wetscav

!===========================================================================
 subroutine wetscav_incloud (id,ktau,dtstep,ktauc,config_flags,           &
        dtstepc,pdel,alt,t_phy,p8w,t8w,p_phy,chem,rho_phy,cldfra,             &
        qlsink,ilsink,                    &
        r2evap,i2evap,                    &
        gas_aqfrac, numgas_aqfrac,                                       &
        ntype_aer, nsize_aer, ncomp_aer,                                 &
        massptr_aer, dens_aer, numptr_aer, hygro_aer,                    &
        maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase,        &
        volumcen_sect, volumlo_sect, volumhi_sect,                       &
        waterptr_aer, dens_water_aer,                                    &
        qv,qc,qi,qs,      &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte                                        )

!-----------------------------------------------------------------------
! Purpose:
! Does in-cloud scavenging of aerosols without cloud-born aerosol 
!
! Authors: Chun Zhao, PNNL 
!-----------------------------------------------------------------------
  USE module_configure, only: grid_config_rec_type
  USE module_model_constants, only: g,rhowater, mwdry
  USE module_state_description
  USE module_dep_simple, only: is_aerosol

! USE module_radiation_driver, only: cal_cldfra,cal_cldfra2


   implicit none

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    ::                                &
                                      ids,ide, jds,jde, kds,kde,    &
                                      ims,ime, jms,jme, kms,kme,    &
                                      its,ite, jts,jte, kts,kte,    &
                                      id, ktau, ktauc, numgas_aqfrac
      REAL,      INTENT(IN   ) :: dtstep,dtstepc
!
! all advected chemical species
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),          &
         INTENT(INOUT ) ::                                chem

! fraction of gas species in cloud water
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, numgas_aqfrac ),     &
         INTENT(IN ) ::                                   gas_aqfrac

!
! input from meteorology
!
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,        &
          INTENT(IN   ) ::                                          &
                                                        alt,        &
                                                      t_phy,        &
                                                      p_phy,        &
                                                    t8w,p8w,        &
                             qlsink,ilsink,        &
                             r2evap,i2evap,        &
                                             rho_phy,cldfra
   integer, intent(in) :: maxd_atype, maxd_asize, maxd_acomp, maxd_aphase
   integer, intent(in) :: ai_phase

   integer, intent(in) :: ntype_aer 
   integer, intent(in) :: nsize_aer( maxd_atype ),   & ! number of size bins
          ncomp_aer( maxd_atype ),   & ! number of chemical components
          massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase), & ! index for mixing ratio 
          waterptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol water
          numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) ! index for the number mixing ratio
   real, intent(in) :: dens_aer( maxd_acomp, maxd_atype ),   & ! material density (g/cm3)
          dens_water_aer, hygro_aer( maxd_acomp, maxd_atype) ! water density (g/m3)
   real, intent(in) :: volumcen_sect( maxd_asize, maxd_atype ),   & ! single-particle volumes (cm3)
                        volumlo_sect( maxd_asize, maxd_atype ),   & ! at section center, lower boundary,
                        volumhi_sect( maxd_asize, maxd_atype )      ! upper boundary

   real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: &
                qv, qc, qi,qs
   real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: &
                pdel 


! LOCAL  VAR
!
   real :: pdel_fac
   real :: qcld,lcldfra,icldfra
   real :: scale,scale1,scale2
   real :: fracscav

   integer :: i,j,k,l,m,n,lmod
   integer :: lmass,lnum

   real, parameter :: tmelt = 273.15     ! freezing T of fresh water (K)


!Prescribed scavenging ratios for cases without online aerosol-activation:
!--- stratifrom cloud water in water and mixed clouds (T>238.5K):
   REAL :: csr_strat_wat(3) ! for three modes
!--- stratifrom cloud ice in mixed clouds (T>238.5K):
   REAL :: csr_strat_mix(3) ! for three modes
!--- stratifrom cloud ice in ice clouds (T<238.5K):
   REAL :: csr_strat_ice(3) ! for three modes
!--- convective clouds: 
!  REAL :: csr_conv(3) ! is not used for now

   !aerosol in-cloud scavenging coefficent
   real :: lscav_coef(3)  !3 modes
   real :: iscav_coef(3)  !3 modes
   !gas in-cloud scavenging coefficent
   real :: gasscav_coef(numgas_aqfrac) 
   real :: totprep,sumclcp,sump
   real,dimension( ims:ime, kms:kme, jms:jme ) :: feffpr
   real, parameter :: minp = 1.e-20  
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem ) ::wdep_mass
   real :: sumhygro,sumvol,hygro,sol_f
   real,parameter :: chem_minv=1.e-20


   wdep_mass=0.0

!--- 1. Pre-set parameters
   
   !--- 1.1 Pre-set scavenging ratios following ECHAM5:
   !reference: based on the Thesis of Stier
   !since we don't have soluble and insoluble modes, all is assumed to be 
   !soluble for now
   !--- Stratiform:
   csr_strat_wat = (/0.25, 0.85, 0.99/)
   csr_strat_mix = (/0.40, 0.75, 0.75/)
   csr_strat_ice = (/0.10, 0.10, 0.10/)
   !--- Convective:
   !csr_conv      = (/0.60, 0.99, 0.99/) ! not used for now

   !--- 1.2 Calculate the effective grid-box fraction 
   !   affected by precipitation (feffpr):
   do j=jts,jte
   do i=its,ite
      totprep=0.0
      sumclcp=0.0
      sump=0.0
   do k=kte,kts,-1  ! from top to bottom
      !--- Preciptiation formation:
      totprep=qlsink(i,k,j)+ilsink(i,k,j)
      !--- Precipitaion-formation weighted averaging:
      sumclcp=sumclcp+cldfra(i,k,j)*totprep
      sump=sump+totprep
      IF(sump.GT.minp) THEN
         feffpr(i,k,j)=sumclcp/sump
      ELSE
         feffpr(i,k,j)=0.
      ENDIF
   enddo
   enddo
   enddo

 do 100 j=jts,jte
 do 100 k=kts,kte
 do 100 i=its,ite

   !--- 1.3 Calculate cloud fraction for liquid and ice 
      if(f_qi)then
        qcld=qc(i,k,j)+qi(i,k,j)
      else
        qcld=qc(i,k,j)
      endif
      if(qcld.lt.-1..or.qcld.gt.1.)then
        write(6,'(a,g12.2,a,3i5)')'qcld=',qcld,' for i,k,j=',i,k,j
        call wrf_error_fatal("1")
      endif
      if(qcld.gt.1.e-20)then
        lcldfra=cldfra(i,k,j)*qc(i,k,j)/qcld
        icldfra=1.0-lcldfra
      else
        lcldfra=0.
        icldfra=0.
      endif

   !--- 1.4 Assign scavenging ratios for different modes of aerosols and types of clouds 
      !--- water clouds:
      if (t_phy(i,k,j)> tmelt) then 
        lscav_coef(:) = csr_strat_wat
        iscav_coef(:) = csr_strat_wat

      !--- mixing clouds:
      else if (t_phy(i,k,j)>238.15) then 
        lscav_coef(:) = csr_strat_mix
        iscav_coef(:) = csr_strat_mix

      !--- cold clouds:
      else 
        lscav_coef(:) = csr_strat_wat
        iscav_coef(:) = csr_strat_ice
      endif


!--- 2. Wet scavenging 

   !--- 2.1 Wet scavenging of aerosol in liquid and ice cloud 
      if (qlsink(i,k,j) >  0.0 .or. ilsink(i,k,j) >  0.0) then

         pdel_fac = pdel(i,k,j)/(g*mwdry)
         do n=1,ntype_aer
          do m=1,nsize_aer(n)

            !calculate volume-weighted bulk hygroscopicity
            sumhygro=0.
            sumvol=0. 
            do l=1,ncomp_aer(n)
              lmass=massptr_aer(l,m,n,ai_phase)
              chem(i,k,j,lmass)=max(chem(i,k,j,lmass),chem_minv)
              sumhygro = sumhygro+hygro_aer(l,n)*chem(i,k,j,lmass)/dens_aer(l,n)
              sumvol = sumvol+chem(i,k,j,lmass)/dens_aer(l,n)
              if (sumvol.le.0) then 
                write(0,*) 'incorrect vol ',n,m,l,sumvol,chem(i,k,j,lmass)
                print*, 'incorrect vol ',n,m,l,sumvol,chem(i,k,j,lmass)
              endif
            enddo ! comp
            hygro=sumhygro/sumvol
            !calculate the soluble fraction of the internal-mixing aerosol based on hygroscopicity
            !this is not accurate, just for temporary application !!!!
            sol_f=min(max(hygro/1.0,0.0),1.0)
           
           if (ntype_aer.gt.1) then ! modal approach 
            if (n==1.and.m==1) lmod = 1 !Aitken mode
            if (n==1.and.m==2) lmod = 2 !Accumulation mode
            if (n==2.and.m==1) lmod = 3 !Coarse mode
           else   ! bin approach,now only support 1 type, i.e., ntype_aer=1 
            if (nsize_aer(1).le.4) then  ! 4-bin
             if (m==1) lmod = 1 !Aitken mode
             if (m==2) lmod = 2 !Accumulation mode
             if (m==3) lmod = 2 !Accumulation mode
             if (m==4) lmod = 3 !Coarse mode
            else    ! 8-bin
             if (m==1) lmod = 1 !Aitken mode
             if (m==2) lmod = 1 !Aitken mode
             if (m==3) lmod = 2 !Accumulation mode
             if (m==4) lmod = 2 !Accumulation mode
             if (m==5) lmod = 2 !Accumulation mode
             if (m==6) lmod = 3 !Coarse mode
             if (m==7) lmod = 3 !Coarse mode
             if (m==8) lmod = 3 !Coarse mode
            endif
           endif
 
           do l=1,ncomp_aer(n)

             scale1=max(dtstepc*qlsink(i,k,j)*lcldfra*lscav_coef(lmod)*sol_f,0.)
             if(scale1.gt.1.)then
              write(0,*)'qlsink,scale1=',qlsink(i,k,j),scale1,' i,k,j=',i,k,j
              print *,'qlsink,scale1=',qlsink(i,k,j),scale1,' i,k,j=',i,k,j
              scale1=1.
             endif

             !scale2=max(dtstepc*ilsink(i,k,j)*icldfra*iscav_coef(lmod)*sol_f,0.)
             scale2=max(dtstepc*ilsink(i,k,j)*icldfra*iscav_coef(lmod),0.)
             if(scale2.gt.1.)then
              write(0,*)'ilsink,scale2=',ilsink(i,k,j),scale2,' i,k,j=',i,k,j
              print *,'ilsink,scale2=',ilsink(i,k,j),scale2,' i,k,j=',i,k,j
              scale2=1.
             endif
             !turn off the ice-cloud scavenging for now.   
             !The auto-conversion rate should be correctly treated in microphysics in future.
             scale2=0.0

             scale=max(1.0-scale1-scale2,0.)
             if(scale.gt.1.)then
              write(0,*)'scale=',scale,' i,k,j=',i,k,j
              print *,'scale=',scale,' i,k,j=',i,k,j
              scale=1.
             endif
 
             lmass=massptr_aer(l,m,n,ai_phase)
             if (lmass < param_first_scalar) cycle
!            qsrflx(i,j,lmass)=qsrflx(i,j,lmass)+chem(i,k,j,lmass)*(scale-1.)*pdel_fac
             wdep_mass(i,k,j,lmass)=chem(i,k,j,lmass)*(1.0-scale)
             chem(i,k,j,lmass)=chem(i,k,j,lmass)*scale
           enddo ! comp
             lnum=numptr_aer(m,n,ai_phase)
             if (lnum < param_first_scalar) cycle
!            qsrflx(i,j,lnum)=qsrflx(i,j,lnum)+chem(i,k,j,lnum)*(scale-1.)*pdel_fac
             wdep_mass(i,k,j,lnum)=chem(i,k,j,lnum)*(1.0-scale)
             chem(i,k,j,lnum)=chem(i,k,j,lnum)*scale
          enddo ! size
         enddo ! type
      endif   ! qlsink > 0 or ilsink > 0

   !--- 2.2 Wet scavenging of gases in liquid cloud 

      gasscav_coef(1:numgas_aqfrac) = gas_aqfrac(i,k,j,1:numgas_aqfrac)
      !The coefficient should be calculated in cloud-chemistry,
      !but now,it's not available, i.e., gas_aqfrac==0.0, 
      !so the coefficient is calculated following GEOS-CHEM 
      call gascloudscav_coef(numgas_aqfrac,config_flags,dtstepc,  & 
                 t_phy(i,k,j),p_phy(i,k,j),gasscav_coef,   &
                 qlsink(i,k,j),ilsink(i,k,j),lcldfra,icldfra,feffpr(i,k,j), &
                 qc(i,k,j),qi(i,k,j),rho_phy(i,k,j),num_chem,chem(i,k,j,:), &
                             ims,ime,kms,kme,jms,jme,              &
                             its,ite,kts,kte,jts,jte)
             
      do 101 l = param_first_scalar, min(num_chem,numgas_aqfrac)
        if ( is_aerosol(l) ) goto 101
        fracscav = gasscav_coef(l)
        if (fracscav > 0.0) then
           fracscav = max( 0.0, min( 1.0, fracscav ) )
           scale = 1.0 - fracscav
           pdel_fac = pdel(i,k,j)/(g*mwdry)
!          qsrflx(i,j,l) = qsrflx(i,j,l)+chem(i,k,j,l)*(scale-1.)*pdel_fac
           wdep_mass(i,k,j,l)=chem(i,k,j,l)*(1.0-scale)
           chem(i,k,j,l) = chem(i,k,j,l)*scale
        endif
 101  continue

 100 continue

   !--- 2.3 Re-evaporation  
   !The evaporation rates (r2evap,i2evap) are outputed, but not treated for RAINOUT for now 


 end subroutine wetscav_incloud
!===========================================================================

subroutine gascloudscav_coef (numgas_aqfrac, config_flags,deltat,  &
                             t,pmid,gasscav_coef,  &
                             qlsink,ilsink,lcldfra,icldfra,feffpr,  &
                             qc,qi,airdens,num_chem,chem0,  &
                             ims,ime,kms,kme,jms,jme,        &
                             its,ite,kts,kte,jts,jte)

!-----------------------------------------------------------------------
!
! Purpose:
! Calculate the gas scavenging coefficient in cloud water 
!
! Authors: Chun Zhao, PNNL 
!
!-----------------------------------------------------------------------

  use module_configure, only:  grid_config_rec_type,   &
                param_first_scalar,   &
                p_so2, p_h2o2, p_sulf, p_msa,   &
                p_hno3, p_hcl, p_nh3,p_hcho

   implicit none

   integer, intent(in)  :: ims,ime            ! column dimension
   integer, intent(in)  :: kms,kme            ! level dimension
   integer, intent(in)  :: jms,jme            ! column dimension
   integer, intent(in)  :: its,ite            ! column identifier
   integer, intent(in)  :: kts,kte            ! level identifier
   integer, intent(in)  :: jts,jte            ! column identifier
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   integer, intent(in)  :: num_chem 
   REAL, DIMENSION(num_chem ),          &
       INTENT(INOUT ) ::                          chem0


   real, intent(in) :: deltat           ! model timestep

   real, intent(in) :: qlsink,ilsink,lcldfra,icldfra,feffpr,qc,qi,airdens


   real, intent(in) :: t    ! temperature
   real, intent(in) :: pmid ! pressure at model levels

   INTEGER,      INTENT(IN)    :: numgas_aqfrac
   ! fraction of gas species in cloud water
   REAL, DIMENSION(numgas_aqfrac ),INTENT(INOUT) :: gasscav_coef


!--------------------------Local Variables------------------------------
   integer :: i, j, k         ! x, y, z work index
   integer :: p1st

   integer, parameter :: ng = 8
   integer, parameter :: ig_so2   = 1
   integer, parameter :: ig_h2o2  = 2
   integer, parameter :: ig_h2so4 = 3
   integer, parameter :: ig_msa   = 4
   integer, parameter :: ig_hno3  = 5
   integer, parameter :: ig_hcl   = 6
   integer, parameter :: ig_nh3   = 7
   integer, parameter :: ig_hcho   = 8
   integer :: ig, lg, lg_ptr(ng)

   real :: rainrate,rainrate1 
   real :: rainfrac 

   REAL, PARAMETER    :: CONV_H2O2 = 8.27042925126d-1
   REAL, PARAMETER    :: CONV_NH3  = 5.69209978831d-1
   real :: c_h2o
   real :: qcv,qiv,l2g,i2g
   real :: c_tot,f_l,f_i 
   real :: so2loss


   gasscav_coef = 0.0
   rainrate=qlsink+ilsink

   qcv=(qc/(1.0*1.e3))*airdens  !kg/kg to m3/m3
   qiv=(qi/(0.9*1.e3))*airdens

   p1st = param_first_scalar

   lg_ptr(ig_so2  ) = p_so2
   lg_ptr(ig_h2o2 ) = p_h2o2
   lg_ptr(ig_h2so4) = p_sulf
   lg_ptr(ig_msa  ) = p_msa
   lg_ptr(ig_hno3 ) = p_hno3
   lg_ptr(ig_hcl  ) = p_hcl
   lg_ptr(ig_nh3  ) = p_nh3
   lg_ptr(ig_hcho ) = p_hcho

   !------------------------------------
   ! Coefficient begins here!
   !------------------------------------
   do ig = 1, ng
    if (lg_ptr(ig)==p_hno3) then !HNO3 (treat it as aerosol)
      rainfrac = GET_RAINFRAC(rainrate, feffpr, deltat)
      gasscav_coef(lg_ptr(ig))=min(rainfrac,0.9999)
!     if (rainfrac>0.5) then 
!     write(0,*)  'czhao check hno3 coef',rainfrac
!     print*,  'czhao check hno3 coef',rainfrac
!     endif
    else if (lg_ptr(ig)==p_h2o2) then ! H2O2 (liquid & ice phases) 
         ! Compute ice to gas ratio for H2O2 by co-condensation
         ! (Eq. 9, Jacob et al, 2000)
         c_h2o= E_ICE( t ) / (pmid*0.01)  ! hPa 

         IF ( c_h2o > 0d0 ) THEN
            I2G = ( qiv / C_H2O ) * CONV_H2O2
         ELSE
            I2G = 0d0
         ENDIF

         ! Compute liquid to gas ratio for H2O2, using
         ! the appropriate parameters for Henry's law
         ! (Eqs. 7, 8 and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( 8.3d4, -7.4d3, t, qcv, L2G )

         ! Fraction of H2O2 in liquid & ice phases
         ! (Eqs. 4, 5, 6, Jacob et al, 2000)
         C_TOT = 1d0 + L2G + I2G
         F_L   = L2G / C_TOT
         F_I   = I2G / C_TOT

         ! Compute the rate constant K.  The retention factor  
         ! for liquid H2O2 is 0.05 for 248 K < T < 268 K, and 
         ! 1.0 for T >= 268 K.  (Eq. 1, Jacob et al, 2000)
         IF ( t >= 268d0 ) THEN
            rainrate1 = rainrate * ( F_L + F_I )
         ELSE IF ( t > 248d0  .and. t < 268d0 ) THEN
            rainrate1 = rainrate * ( ( 5d-2 * F_L ) + F_I )
         ELSE
            rainrate1 = rainrate * F_I
         ENDIF

         ! Compute RAINFRAC, the fraction of rained-out H2O2
         ! (Eq. 10, Jacob et al, 2000)
         rainfrac = GET_RAINFRAC(rainrate1, feffpr, deltat)
         gasscav_coef(lg_ptr(ig))=min(rainfrac,0.9999)

    else if (lg_ptr(ig)==p_hcho) then ! CH2O (liquid phase only) 
         ! Compute liquid to gas ratio for CH2O, using
         ! the appropriate parameters for Henry's law
         ! (Eqs. 7, 8 and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( 3.0d3, -7.2d3, T, qcv, L2G )

         ! Fraction of CH2O in liquid phase 
         ! NOTE: CH2O does not exist in the ice phase!
         ! (Eqs. 4, 5, Jacob et al, 2000)
         C_TOT = 1d0 + L2G
         F_L   = L2G / C_TOT

         ! Compute the rate constant K.  The retention factor  
         ! for liquid CH2O is 0.02 for 248 K < T < 268 K, and 
         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
         IF ( t >= 268d0 ) THEN
            rainrate1 = rainrate * F_L
         ELSE IF ( t > 248d0 .and. t < 268d0 ) THEN
            rainrate1 = rainrate * ( 2d-2 * F_L )
         ELSE
            rainrate1 = 0d0
         ENDIF
         ! Compute RAINFRAC, the fraction of rained-out CH2O
         ! (Eq. 10, Jacob et al, 2000)
         rainfrac = GET_RAINFRAC(rainrate1, feffpr, deltat)
         gasscav_coef(lg_ptr(ig))=min(rainfrac,0.9999)
    else if (lg_ptr(ig)==p_so2) then ! SO2
         !==============================================================
         ! NOTE: SO2 and H2O2 are in [v/v] and here RAINFRAC contains 
         ! the amount of SO2 lost due to rainout normalized by the
         ! total SO2 -- so that in WETDEP routine mulitiplying SO2 in 
         ! [kg] will produce correct amount.  Need to verify this. 
         ! (rjp, 01/16/02)
         !==============================================================
         ! Treat SO2 as an aerosol
         RAINFRAC = GET_RAINFRAC(rainrate,feffpr,deltat)
         ! Update SO2 and H2O2
         IF ( chem0(p_so2) > 1.e-16 ) THEN

            ! Limit RAINFRAC 
            SO2LOSS      = MIN( chem0(p_so2), chem0(p_h2o2) )
            RAINFRAC     = SO2LOSS * RAINFRAC / chem0(p_so2)
            RAINFRAC     = MAX( RAINFRAC, 0d0 )

            ! Update saved H2O2 concentration
            chem0(p_h2o2) = chem0(p_h2o2) - ( chem0(p_so2) * RAINFRAC )
            chem0(p_h2o2) = MAX(chem0(p_h2o2),1.e-16)
         ELSE
            RAINFRAC = 0D0
         ENDIF
         gasscav_coef(lg_ptr(ig))=min(rainfrac,0.9999)
    else if (lg_ptr(ig)==p_sulf.or.lg_ptr(ig)==p_hcl) then !SO4, assume HCL the same as SO4
         RAINFRAC = GET_RAINFRAC(rainrate,feffpr,deltat)
         gasscav_coef(lg_ptr(ig))=min(rainfrac,0.9999)
    else if (lg_ptr(ig)==p_msa) then !MSA
         RAINFRAC = GET_RAINFRAC(rainrate,feffpr,deltat)
         gasscav_coef(lg_ptr(ig))=min(rainfrac,0.9999)
    else if (lg_ptr(ig)==p_nh3) then !NH3 (liquid & ice phases)
         ! Compute ice to gas ratio for NH3 by co-condensation
         ! (Eq. 9, Jacob et al, 2000)
         c_h2o= E_ICE( t ) / (pmid*0.01)  ! hPa 
         IF ( C_H2O > 0d0 ) THEN
            I2G = ( qiv / C_H2O) * CONV_NH3
         ELSE
            I2G = 0d0
         ENDIF
         ! Compute liquid to gas ratio for NH3, using
         ! the appropriate parameters for Henry's law
         ! (Seinfeld and Pandis, p343 eq. 6.8)
         ! PH    = 4.5  ! Assumed PH for typical cloud drop
         ! Hstar = 1.054d11 * (10.**(-PH)) == 3.3d6
         CALL COMPUTE_L2G( 3.3d6, -4.1d3, t, qcv, L2G )
         ! Fraction of NH3 in liquid & ice phases
         ! (Eqs. 4, 5, 6, Jacob et al, 2000)
         C_TOT = 1d0 + L2G + I2G
         F_L   = L2G / C_TOT
         F_I   = I2G / C_TOT
         ! Compute the rate constant K.  The retention factor  
         ! for liquid NH3 is 0.05 for 248 K < T < 268 K, and 
         ! 1.0 for T >= 268 K.  (Eq. 1, Jacob et al, 2000)
         IF ( t >= 268d0 ) THEN
            rainrate1 = rainrate * ( F_L + F_I )
         ELSE IF ( t > 248d0  .and. t < 268d0 ) THEN
            rainrate1 = rainrate * ( ( 5d-2 * F_L ) + F_I )
         ELSE
            rainrate1 = rainrate * F_I
         ENDIF
         ! Compute RAINFRAC, the fraction of rained-out NH3
         ! (Eq. 10, Jacob et al, 2000)
         RAINFRAC = GET_RAINFRAC(rainrate1,feffpr,deltat)
         gasscav_coef(lg_ptr(ig))=min(rainfrac,0.9999)
    endif

   enddo  !ig

   return
end subroutine gascloudscav_coef

!============================================================================
SUBROUTINE COMPUTE_L2G( Kstar298, H298_R, TK, H2OLIQ, L2G )
!
!******************************************************************************
!  Subroutine COMPUTE_L2G computes the ratio L2G = Cliq / Cgas, which is 
!  the mixing ratio of tracer in the liquid phase, divided by the mixing 
!  ratio of tracer in the gas phase.
!
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate 
!  values of Kstar298 and H298_R must be supplied for each tracer.  
!  (cf Jacob et al 2000, p. 3)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) Kstar298 (REAL*8) : Eff. Henry's law constant @ 298 K   [moles/atm]
!  (2 ) H298_R   (REAL*8) : Molar heat of formation @ 298 K / R [K]
!  (3 ) TK       (REAL*8) : Temperature at grid box (I,J,L)     [K]
!  (4 ) H2OLIQ   (REAL*8) : Liquid water content at (I,J,L)     [cm3 H2O/cm3 air]
!  Arguments as Output:
!  ============================================================================
!  (5 ) L2G      (REAL*8) : Cliq/Cgas ratio for given tracer  [unitless]
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Jacob et al, 2000
!******************************************************************************
!
      ! Arguments
      REAL*8,INTENT(IN) :: KStar298, H298_R
      REAL, INTENT(IN) :: TK, H2OLIQ
      REAL, INTENT(OUT) :: L2G

      ! Local variables
      REAL              :: Kstar

      ! R = universal gas constant [atm/moles/K]
      REAL, PARAMETER   :: R = 8.32d-2

      ! INV_T0 = 1/298 K
      REAL, PARAMETER   :: INV_T0 = 1d0 / 298d0

      !=================================================================
      ! COMPUTE_L2G begins here!
      !=================================================================

      ! Get Kstar, the effective Henry's law constant for temperature TK
      Kstar = Kstar298 * EXP( -H298_R * ( ( 1d0 / TK ) - INV_T0 ) )

      ! Use Henry's Law to get the ratio:
      ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
      L2G   = Kstar * H2OLIQ * R * TK

      ! Return to calling program
END SUBROUTINE COMPUTE_L2G

!=================================================================

FUNCTION E_ICE( TK ) RESULT( VALUE )
! 
!******************************************************************************
!  Subroutine E_ICE computes Eice(T), the saturation vapor pressure of ice
!  at a given Celsius temperature. 
!  
!  Arguments as Input:
!  (1 ) TK (REAL*8) : Ambient temperature [K] 
! 
!  References:
!  (1 ) Marti & Mauersberber (GRL '93) formulation of saturation 
!        vapor pressure of ice [Pa] is: log P = A/TK + B
!  NOTES:
!  (1 ) Now use the same analytic function as the Goddard CTM (bmy, 2/8/05)
!******************************************************************************
      ! Arguments as Input
      REAL, INTENT(IN) :: TK
      ! Return value
      REAL             :: VALUE
      ! Parameters
      REAL, PARAMETER  :: A = -2663.5d0
      REAL, PARAMETER  :: B =  12.537d0

      ! E_ICE begins here!

      ! Saturation vap press of Ice [Pa] -- divide by 100 for [hPa]
      VALUE = ( 10d0**( A/TK + B ) ) / 100d0

      ! Return to calling program
END FUNCTION E_ICE


!=================================================================
FUNCTION GET_RAINFRAC( K, F, DT ) RESULT( RAINFRAC )
!
!  Function GET_RAINFRAC computes the fraction of tracer lost to rainout 
!  according to Jacob et al 2000. 
!
!  Arguments as Input:
!******************************************************************************
!  (1 ) K  (REAL*8) : Rainout rate constant              [1/s]
!  (2 ) DT (REAL*8) : Timestep for rainout event         [s]
!  (3 ) F  (REAL*8) : Fraction of grid box precipitating [unitless]
!******************************************************************************
!
      ! Arguments
      REAL, INTENT(IN) :: K, F, DT
      ! Local variables
      REAL             :: RAINFRAC

      ! (Eq. 10, Jacob et al, 2000 ) 
      RAINFRAC = F * ( 1 - EXP( -K * DT ) )

      ! Return to RAINOUT
      END FUNCTION GET_RAINFRAC
!===========================================================================

!===========================================================================
subroutine aerimpactscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,num_chem,  &
     ntype_aer, nsize_aer, ncomp_aer, massptr_aer, dens_aer, numptr_aer, &
     maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase,           &
     volumcen_sect, volumlo_sect, volumhi_sect,                          &
     waterptr_aer, dens_water_aer,                                       &
     scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
     deltat,     t,          pmid,       pdel, chem,                     &
     isprx,      fapx,       pfx,        pfx_inrain,                     &
     dqdt,       dotend                                           )


!-----------------------------------------------------------------------
!
! Purpose:
! Does below cloud scavenging of aerosols through impaction-interception
!
! Removal rates for aersol number, (surface - still to do), and volume
! are computed for each mode using precalculated lookup tables.
! The tables account for variables in wet-dgnum, wet density,
! air temperature, and air pressure.
!
! Authors: R. Easter
!
!-----------------------------------------------------------------------
  USE module_model_constants, only: g,rhowater, mwdry
  USE module_state_description, only: param_first_scalar

   implicit none

!-----------------------------------------------------------------------
!
! Input arguments
!
! abbreviations & acronyms
!    TMR = tracer mixing ratio
!    l-s = large scale
!
   integer, intent(in)  :: num_chem           ! number of chemical species
   integer, intent(in)  :: ims,ime            ! column dimension
   integer, intent(in)  :: kms,kme            ! level dimension
   integer, intent(in)  :: jms,jme            ! column dimension
   integer, intent(in)  :: its,ite            ! column identifier
   integer, intent(in)  :: kts,kte            ! level identifier
   integer, intent(in)  :: jts,jte            ! column identifier
   real, intent(in) :: deltat           ! model timestep

   real, intent(in) :: t(ims:ime,kms:kme,jms:jme)    ! temperature
   real, intent(in) :: pmid(ims:ime,kms:kme,jms:jme) ! pressure at model levels
   real, intent(in) :: pdel(ims:ime,kms:kme,jms:jme) ! pressure thickness of levels
   real, intent(in) :: chem(ims:ime,kms:kme,jms:jme,num_chem) ! chem array

   logical, intent(in)  :: isprx(ims:ime,kms:kme,jms:jme) ! true if precip at i,k
   real, intent(in) :: fapx(ims:ime,kms:kme,jms:jme)    ! frac. area for precip
   real, intent(in) :: pfx(ims:ime,kms:kme,jms:jme)     ! grid-avg precip
                                                 ! flux (kg/m2/s)
   real, intent(in) :: pfx_inrain(ims:ime,kms:kme,jms:jme)  ! in-rain-area precip flux (kg/m2/s)

   real, intent(inout) :: dqdt(ims:ime,kms:kme,jms:jme,num_chem) ! TMR tendency array
   logical,  intent(inout) :: dotend(num_chem)     ! flag for doing scav
!czhao   real, intent(inout) :: qsrflx(ims:ime,jms:jme,num_chem) 
                        ! changes to column tracer burdens by wet scavenging over current timestep
                        ! this routine adds on the contribution from aerosol impaction and brownian
                        !    diffusion scavenging by precipitation
   integer, intent(in)  :: maxd_atype, maxd_asize, maxd_acomp, maxd_aphase
   integer, intent(in) :: ai_phase
   integer, intent(in) :: ntype_aer
   integer, intent(in) ::  nsize_aer( maxd_atype ),   & ! number of size bins
      	  ncomp_aer( maxd_atype ),   & ! number of chemical components
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase), & ! index for mixing ratio
      	  waterptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol water
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) ! index for the number mixing ratio
   real, intent(in) :: volumcen_sect( maxd_asize, maxd_atype ),   & ! single-particle volumes (cm3)
                        volumlo_sect( maxd_asize, maxd_atype ),   & ! at section center, lower boundary,
                        volumhi_sect( maxd_asize, maxd_atype )      ! upper boundary
   real, intent(in) :: dens_aer( maxd_acomp, maxd_atype ),   & ! material density (g/cm3)
                       dens_water_aer ! water density (g/m3)

   real, intent(in) :: dlndg_nimptblgrow
   integer, intent(in) :: nimptblgrow_mind, nimptblgrow_maxd
   real, intent(in) :: scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype), &
     	     scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)

!--------------------------Local Variables------------------------------

   integer :: i,j               ! Work index
   integer :: jgrow, jp       ! Work index
   integer :: k               ! Work index
   integer :: l, ll, m, n        ! Work index

   logical :: ispr_anywhere

   real :: dry_mass, dry_volu, dry_volu_1p
   real :: duma
   real :: dumfhi, dumflo
   real :: dumimpactamt0, dumimpactamt3, dumimpactamtsum
   real :: dumimpactratea0, dumimpactratea3
   real :: dumimpactrateb0, dumimpactrateb3
   real :: dumdgratio, dumrate, dumwetdens
   real :: dumlogdens, dumlogptot, dumlogtemp
   real :: dumnumb
   real :: dumscavratenum, dumscavratevol
   real :: pdel_dt_fac
   real :: pf_to_prmmh
   real :: scavimptbl1, scavimptbl2, scavimptbl3, scavimptbl4
   real :: xgrow
   real :: wet_mass, wet_volu
   real, save :: third = 0.33333333

!-----------------------------------------------------------------------
!

!  if (ncol .ne. -987654321) return

! precip rates -- 1.0 kgwtr/m2/s = 1.0e-3 m3wtr/m2/s = 1.0e-3 m/s
!                                                    = 1.0 mm/s = 3600 mm/h

   ispr_anywhere = .false.

   do 5900 i = its,ite
   do 5900 j = jts,jte

      do 4900 k = kte,kts,-1

! skip this level if no precip
      if ( isprx(i,k,j) ) then
        ispr_anywhere = .true.
      else
        goto 4900
      end if


	dumlogtemp = log( t(i,k,j) )
!   dumlogptot = log( pressure in dynes/cm2 )
	dumlogptot = log( 10.0*pmid(i,k,j) )

!   compute removal amounts for each aerosol mode
       do 3900 n=1,ntype_aer
       do 3900 m=1,nsize_aer(n)
          dry_volu = 0.
          dry_mass = 0
          do ll = 1, ncomp_aer(n)
             l = massptr_aer(ll,m,n,ai_phase)
             if (l < param_first_scalar) cycle
             duma = max( chem(i,k,j,l)*1.0e-6, 0.0 )
             dry_mass = dry_mass + duma                 ! g-AP/kg-air
             dry_volu = dry_volu + duma/dens_aer(ll,n)  ! cm3-AP/kg-air
          end do
!   If no aerosol is present at this size and type, there is nothing
!   to scavenge so go onto the next size bin. wig, 25-Oct-2005
          if( dry_volu < 1.0e-26 ) goto 3900

!   wet volume
          l = waterptr_aer(m,n)
          if (l >= param_first_scalar) then
             duma = max( chem(i,k,j,l)*1.0e-6, 1.e-20)
          else
             duma = 0.0
          end if
          wet_mass = dry_mass + duma
          wet_volu = dry_volu + duma/dens_water_aer

!   calc dry volume of 1 particle
          dumnumb = chem(i,k,j,numptr_aer(m,n,ai_phase))
          if (dry_volu > volumhi_sect(m,n)*dumnumb) then
             dry_volu_1p = volumhi_sect(m,n)
          else if (dry_volu < volumlo_sect(m,n)*dumnumb) then
             dry_volu_1p = volumlo_sect(m,n)
          else
             dry_volu_1p = dry_volu/dumnumb
          end if
!   interpolate table values using 
!      dumdgratio = [(actual-wet-diam)/(base-dry-diam)]
!                 = [(actual-wet-volu)/(base-dry-volu)]**1/3
          if (wet_volu <= 1.0e9*dry_volu) then
             duma = wet_volu/dry_volu
          else
             duma = 1.0e9
          end if
          dumdgratio = ((dry_volu_1p/volumcen_sect(m,n))*duma)**third

!   dumwetdens = wet aerosol density in g/cm3
	  dumwetdens = wet_mass/wet_volu
	  dumlogdens = log( dumwetdens )
	  dumimpactamt3 = 0
	  dumimpactamt0 = 0
!
!         compute impaction scavenging removal amount for volume
!

	  if ((dumdgratio .ge. 0.99) .and. (dumdgratio .le. 1.01)) then
	    scavimptbl1 = scavimptblvol(1,0,m,n)
	    scavimptbl2 = scavimptblvol(2,0,m,n)
	    scavimptbl3 = scavimptblvol(3,0,m,n)
	    scavimptbl4 = scavimptblvol(4,0,m,n)

	  else
	    xgrow = log( dumdgratio ) / dlndg_nimptblgrow
	    jgrow = int( xgrow )
	    if (xgrow .lt. 0.) jgrow = jgrow - 1
	    if (jgrow .lt. nimptblgrow_mind) then
		jgrow = nimptblgrow_mind
		xgrow = jgrow
	    else 
		jgrow = min( jgrow, nimptblgrow_maxd-1 )
	    end if

	    dumfhi = xgrow - jgrow
	    dumflo = 1. - dumfhi

	    scavimptbl1 = dumflo*scavimptblvol(1,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(1,jgrow+1,m,n)

	    scavimptbl2 = dumflo*scavimptblvol(2,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(2,jgrow+1,m,n)

	    scavimptbl3 = dumflo*scavimptblvol(3,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(3,jgrow+1,m,n)

	    scavimptbl4 = dumflo*scavimptblvol(4,jgrow,m,n) +   &
			  dumfhi*scavimptblvol(4,jgrow+1,m,n)
	  end if

!         apply temperature and pressure corrections
	  dumimpactratea3 = exp( scavimptbl1 + scavimptbl2*dumlogtemp +    &
      		scavimptbl3*dumlogptot + scavimptbl4*dumlogdens )

!   dumimpactratea3 = impaction scav rate (1/h) for precip = 1 mm/h
!   dumimpactrateb3 = impaction scav rate (1/s) for precip = pfx_inrain
!       (dumimpactratea3/3600) = impaction scav rate (1/s) for precip = 1 mm/h
!       (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
!       dumimpactrateb3 = (dumimpactratea3/3600) * (pfx_inrain*3600)
!   dumimpactamt3   = fraction of aerosol removed from entire grid cell
!                         in time deltat
	dumimpactamtsum = 0.0
	    dumimpactrateb3 = dumimpactratea3 * pfx_inrain(i,k,j)
	    dumimpactamt3 =  (1. - exp(-deltat*dumimpactrateb3)) * fapx(i,k,j)
	    dumimpactamtsum = dumimpactamtsum + dumimpactamt3
	dumimpactamt3 = min( dumimpactamtsum, 1.0 )

!   diagnostic output
!	dump = 10.0*pmid(i,k)
!	call calc_1_impact_rate( dgncur_awet(i,k,j,m,n),   &
!      		sigmag_amode(n), dumwetdens,   &
!      		t(i,k), dump,   &
!      		dumscavratenum, dumscavratevol, lun )
!
!	dumr = -9.99e35
!	if (dumscavratevol > 1.0e-35)   &
!		dumr = (dumimpactratea3/dumscavratevol) - 1.0
! 	write(lun,9440) nstep, lchnk, i, k, jp,   &
!		(dumtemp-273.16), dumpress*.001
! 	write(lun,9442) 'rhowet, dgnwet, dgratio, xgrow',   &
!		dumwetdens, dgncur_awet(i,k,n), dumdgratio, xgrow
!	write(lun,9442) 'exa&approx vol rt, relerr, amt',   &
!		dumscavratevol, dumimpactratea3, dumr, dumimpactamt3
!	write(lun,9442) 'pfx_inrain(1:3)               ',   &
!		(pfx_inrain(jp,i,k), jp=1,3)
!	write(lun,9442) 'fapx(1:3)                     ',   &
!		(fapx(jp,i,k), jp=1,3)
!9440	format( / 'ns,lc,i,k,jp,   T(C), p(mb)', i6, 2i4, 2i3, 2f7.1 )
!9442	format( a, 4(1pe11.3) )
!   end diagnostic output


!
!   compute impaction scavenging removal amount to number
!
	if (numptr_aer(m,n,ai_phase) < param_first_scalar) then
	    dumimpactamt0 = dumimpactamt3
	    goto 3700
	end if

!   interpolate table values using log of (actual-wet-size)/(base-dry-size)
	if ((dumdgratio .ge. 0.99) .and. (dumdgratio .le. 1.01)) then
	    scavimptbl1 = scavimptblnum(1,0,m,n)
	    scavimptbl2 = scavimptblnum(2,0,m,n)
	    scavimptbl3 = scavimptblnum(3,0,m,n)
	    scavimptbl4 = scavimptblnum(4,0,m,n)

	else
	    scavimptbl1 = dumflo*scavimptblnum(1,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(1,jgrow+1,m,n)

	    scavimptbl2 = dumflo*scavimptblnum(2,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(2,jgrow+1,m,n)

	    scavimptbl3 = dumflo*scavimptblnum(3,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(3,jgrow+1,m,n)

	    scavimptbl4 = dumflo*scavimptblnum(4,jgrow,m,n) +   &
			  dumfhi*scavimptblnum(4,jgrow+1,m,n)
	end if

!   apply temperature and pressure corrections
	dumimpactratea0 = exp( scavimptbl1 + scavimptbl2*dumlogtemp +    &
      		scavimptbl3*dumlogptot + scavimptbl4*dumlogdens )

	dumimpactamt0 = 0.0
	    dumimpactrateb0 = dumimpactratea0 * pfx_inrain(i,k,j)
	    dumimpactamt0 = dumimpactamt0 +   &
		(1. - exp( -deltat*dumimpactrateb0 )) * fapx(i,k,j)
	dumimpactamt0 = min( dumimpactamt0, 1.0 )

!   diagnostic output
!	dumr = -9.99e35
!	if (dumscavratenum > 1.0e-35)   &
!		dumr = (dumimpactratea0/dumscavratenum) - 1.0
!	write(lun,9442) 'exa&approx num rt, relerr, amt',   &
!		dumscavratenum, dumimpactratea0, dumr, dumimpactamt0
!   end diagnostic output


3700	continue

!
!   compute tendencies
!
	pdel_dt_fac = deltat*pdel(i,k,j)/g
	dumrate = -dumimpactamt3/(deltat*(1.0 + 1.0e-8))
	dumrate = min(0.0,max(-1.0/deltat,dumrate)) ! make sure -1 <= dumrate*deltat <= 0
	do ll = 1, ncomp_aer(n)
	    l = massptr_aer(ll,m,n,ai_phase)
	    if (l < param_first_scalar) cycle
	    dqdt(i,k,j,l) = chem(i,k,j,l)*dumrate
!czhao	    qsrflx(i,j,l) = qsrflx(i,j,l) + dqdt(i,k,j,l)*pdel_dt_fac ! aerosol mass (ug/m2)
	end do
	l = waterptr_aer(m,n)
	if (l >= param_first_scalar) then 
	    dqdt(i,k,j,l) = chem(i,k,j,l)*dumrate
!czhao	    qsrflx(i,j,l) = qsrflx(i,j,l) + dqdt(i,k,j,l)*pdel_dt_fac ! aerosol water mass (ug/m2)
	end if
	l = numptr_aer(m,n,ai_phase)
	if (l >= param_first_scalar) then 
	    dumrate = -dumimpactamt0/(deltat*(1.0 + 1.0e-8))
	    dqdt(i,k,j,l) = chem(i,k,j,l)*dumrate
!czhao	    qsrflx(i,j,l) = qsrflx(i,j,l) + dqdt(i,k,j,l)*pdel_dt_fac ! aerosol number (1/m2)
	end if



3900 continue  ! "n = 1, ntot_amode"

4900 continue  ! "k = kts,kte"

5900 continue  ! "i = its,ite, j=jts,jte"


! set dotend's
   if ( ispr_anywhere ) then
      do n=1,ntype_aer
      do m=1,nsize_aer(n)
         do ll = 1, ncomp_aer(n)
            if (massptr_aer(ll,m,n,ai_phase) >= param_first_scalar) &
               dotend(massptr_aer(ll,m,n,ai_phase)) = .true.
         end do
         if (waterptr_aer(m,n) >= param_first_scalar) &
            dotend(waterptr_aer(m,n)) = .true.
         if (numptr_aer(m,n,ai_phase) >= param_first_scalar) &
            dotend(numptr_aer(m,n,ai_phase)) = .true.
      end do
      end do
   end if


   return
end subroutine aerimpactscav



!===========================================================================
!===========================================================================
  subroutine initwet( ntype_aer, nsize_aer, ncomp_aer, massptr_aer,      &
       dens_aer, numptr_aer, maxd_acomp, maxd_asize,maxd_atype,          &
       maxd_aphase, dcen_sect, sigmag_aer, waterptr_aer, dens_water_aer, &
       scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, &
       dlndg_nimptblgrow)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes lookup table for aerosol impaction/interception scavenging rates
!
! Authors: R. Easter
!
!-----------------------------------------------------------------------
  implicit none

   integer, intent(in) :: maxd_atype, maxd_asize, maxd_acomp, maxd_aphase
   integer, intent(in) :: ntype_aer
   integer, intent(in) ::  nsize_aer( maxd_atype ) ! number of size bins
   integer, intent(in) :: ncomp_aer( maxd_atype ) ! number of chemical components
   integer, intent(in) :: massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase) ! index for mixing ratio
   integer, intent(in) :: waterptr_aer( maxd_asize, maxd_atype ) ! index for aerosol water
   integer, intent(in) :: numptr_aer( maxd_asize, maxd_atype, maxd_aphase ) ! index for the number mixing ratio
   real, intent(in) :: dens_aer( maxd_acomp, maxd_atype ) ! material density (g/cm3)
   real, intent(in) :: dens_water_aer ! water density (g/m3)
   real, intent(in) :: sigmag_aer(maxd_asize, maxd_atype)
   real, intent(in) :: dcen_sect( maxd_asize, maxd_atype ) 
       ! "base" volume-mean diameter (cm) == section center diameter
       ! scav rates for aerosol number and volume are pre-calculated for a range of volume-mean diameters
       !    as well as aerosol densities and air pressure and temperature
       ! for specific time/location/aerosol properties, these rates are interpolated, and
       !    the size interpolation uses [(actual volume-mean diameter)/(base volume-mean diameter)]

   real, intent(out) :: dlndg_nimptblgrow
   integer, intent(in) :: nimptblgrow_mind, nimptblgrow_maxd
   real, intent(out) :: scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)
   real, intent(out) :: scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)

!   local variables
	integer nnfit_maxd
	parameter (nnfit_maxd=27)

	integer i, j, m, n, jgrow, jdens, jpress, jtemp, ll, mode, nnfit
        integer lunerr

	real dg0, dg0_base, press, rhodryaero, rhowetaero, rmserr, &
     	scavratenum, scavratevol, sigmag,                &
     	temp, wetdiaratio, wetvolratio
	real aafitnum(4), xxfitnum(4,nnfit_maxd), yyfitnum(nnfit_maxd)
	real aafitvol(4), xxfitvol(4,nnfit_maxd), yyfitvol(nnfit_maxd)


	lunerr = 6
	dlndg_nimptblgrow = log( 1.25d00 )

        do 7900 n=1,ntype_aer
        do 7900 m=1,nsize_aer(n)

	sigmag = sigmag_aer(m,n)
	dg0_base = dcen_sect(m,n)*exp( -1.5*((log(sigmag))**2) )

!   for setting up the lookup table, use the dry density of the first 
!   chemical component of the aerosol type (which currently will be so4)
	rhodryaero = dens_aer(1,n)

	do 7800 jgrow = nimptblgrow_mind, nimptblgrow_maxd

	wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
	dg0 = dg0_base*wetdiaratio

	wetvolratio = exp( jgrow*dlndg_nimptblgrow*3. )
	rhowetaero = 1.0 + (rhodryaero-1.0)/wetvolratio
	rhowetaero = min( rhowetaero, rhodryaero )

!
!   compute impaction scavenging rates at 9 temp-press pairs and save
!
	nnfit = 0

	do 5900 jtemp = -1, 1
	temp = 273.16 + 25.*jtemp

	do 5900 jpress = -1, 1
	press = 0.75e6 + 0.25e6*jpress

	do 5900 jdens = 0, 2
	rhowetaero = rhodryaero**(0.5*jdens)

	call calc_1_impact_rate( &
     		dg0, sigmag, rhowetaero, temp, press, &
     		scavratenum, scavratevol, lunerr ) 

	nnfit = nnfit + 1
	if (nnfit .gt. nnfit_maxd) then
	    call wrf_error_fatal('*** subr. aerosol_impact_setup -- nnfit too big')
	end if

	xxfitnum(1,nnfit) = 1.
	xxfitnum(2,nnfit) = log( temp )
	xxfitnum(3,nnfit) = log( press )
	xxfitnum(4,nnfit) = log( rhowetaero )
	yyfitnum(nnfit) = log( scavratenum )

	xxfitvol(1,nnfit) = 1.
	xxfitvol(2,nnfit) = log( temp )
	xxfitvol(3,nnfit) = log( press )
	xxfitvol(4,nnfit) = log( rhowetaero )
	yyfitvol(nnfit) = log( scavratevol )

5900	continue

!
!   do linear regression
!	log(scavrate) = a1 + a2*log(temp) + a3*log(press) + a4*log(wetdens)
!
	call mlinft( xxfitnum, yyfitnum, aafitnum, nnfit, 4, 4, rmserr )
	call mlinft( xxfitvol, yyfitvol, aafitvol, nnfit, 4, 4, rmserr )

	do i = 1, 4
	    scavimptblnum(i,jgrow,m,n) = aafitnum(i)
	    scavimptblvol(i,jgrow,m,n) = aafitvol(i)
	end do


7800	continue
7900	continue

        return
      end subroutine initwet



!===========================================================================
!===========================================================================
	subroutine calc_1_impact_rate(             &
     		dg0, sigmag, rhoaero, temp, press, &
     		scavratenum, scavratevol, lunerr )
!
!   routine computes a single impaction scavenging rate
!	for precipitation rate of 1 mm/h
!
!   dg0 = geometric mean diameter of aerosol number size distrib. (cm)
!   sigmag = geometric standard deviation of size distrib.
!   rhoaero = density of aerosol particles (g/cm^3)
!   temp = temperature (K)
!   press = pressure (dyne/cm^2)
!   scavratenum = number scavenging rate (1/h)
!   scavratevol = volume or mass scavenging rate (1/h)
!   lunerr = logical unit for error message
!
	implicit none

!   subr. parameters
	integer lunerr
	real dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol

!   local variables
	integer nrainsvmax
	parameter (nrainsvmax=50)
	real rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
     		vfallrainsv(nrainsvmax)

	integer naerosvmax
	parameter (naerosvmax=51)
	real aaerosv(naerosvmax), &
     	ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

	integer i, ja, jr, na, nr
	real a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
     	real anumsum, avolsum, boltzmann, cair, chi
     	real d, dr, dum, dumfuchs, dx
     	real ebrown, eimpact, eintercept, etotal, freepath, gravity
     	real pi, precip, precipmmhr, precipsum
     	real r, rainsweepout, reynolds, rhi, rhoair, rhowater, rlo, rnumsum
     	real scavsumnum, scavsumnumbb
     	real scavsumvol, scavsumvolbb
     	real schmidt, sqrtreynolds, sstar, stokes, sx
     	real taurelax, vfall, vfallstp
     	real x, xg0, xg3, xhi, xlo, xmuwaterair


	rlo = .005
	rhi = .250
	dr = 0.005
	nr = 1 + nint( (rhi-rlo)/dr )
	if (nr .gt. nrainsvmax) then
	    call wrf_error_fatal('*** subr. calc_1_impact_rate -- nr > nrainsvmax')
	end if

	precipmmhr = 1.0
	precip = precipmmhr/36000. ! cm/s

	ag0 = dg0/2.
	sx = log( sigmag )
	xg0 = log( ag0 )
	xg3 = xg0 + 3.*sx*sx

	xlo = xg3 - 4.*sx
	xhi = xg3 + 4.*sx
	dx = 0.2*sx

	dx = max( 0.2*sx, 0.01 )
	xlo = xg3 - max( 4.*sx, 2.*dx )
	xhi = xg3 + max( 4.*sx, 2.*dx )

	na = 1 + nint( (xhi-xlo)/dx )
	if (na .gt. naerosvmax) then
	    call wrf_error_fatal('*** subr. calc_1_impact_rate -- na > naerosvmax')
	end if

!   air molar density
	cair = press/(8.31436e7*temp)
!   air mass density
	rhoair = 28.966*cair
!   molecular freepath
	freepath = 2.8052e-10/cair
!   boltzmann constant
	boltzmann = 1.3807e-16
!   water density
	rhowater = 1.0
!   gravity
	gravity = 980.616
!   air dynamic viscosity
	airdynvisc = 1.8325e-4 * (416.16/(temp+120.)) *    &
                                        ((temp/296.16)**1.5)
!   air kinemaic viscosity
	airkinvisc = airdynvisc/rhoair
!   ratio of water viscosity to air viscosity (from Slinn)
	xmuwaterair = 60.0

	pi = 3.1415926536

!
!   compute rain drop number concentrations
!	rrainsv = raindrop radius (cm)
!	xnumrainsv = raindrop number concentration (#/cm^3)
!		(number in the bin, not number density)
!	vfallrainsv = fall velocity (cm/s)
!
	precipsum = 0.
	do i = 1, nr
	    r = rlo + (i-1)*dr
	    rrainsv(i) = r
	    xnumrainsv(i) = exp( -r/2.7e-2 )

	    d = 2.*r
	    if (d .le. 0.007) then
		vfallstp = 2.88e5 * d**2.
	    else if (d .le. 0.025) then
		vfallstp = 2.8008e4 * d**1.528
	    else if (d .le. 0.1) then
		vfallstp = 4104.9 * d**1.008
	    else if (d .le. 0.25) then
		vfallstp = 1812.1 * d**0.638
	    else
		vfallstp = 1069.8 * d**0.235
	    end if

	    vfall = vfallstp * sqrt(1.204e-3/rhoair)
	    vfallrainsv(i) = vfall
	    precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
	end do
	precipsum = precipsum*pi*1.333333

	rnumsum = 0.
	do i = 1, nr
	    xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
	    rnumsum = rnumsum + xnumrainsv(i)
	end do

!
!   compute aerosol concentrations
!	aaerosv = particle radius (cm)
!	fnumaerosv = fraction of total number in the bin (--)
!	fvolaerosv = fraction of total volume in the bin (--)
!
	anumsum = 0.
	avolsum = 0.
	do i = 1, na
	    x = xlo + (i-1)*dx
	    a = exp( x )
	    aaerosv(i) = a
	    dum = (x - xg0)/sx
	    ynumaerosv(i) = exp( -0.5*dum*dum )
	    yvolaerosv(i) = ynumaerosv(i)*1.3333*pi*a*a*a
	    anumsum = anumsum + ynumaerosv(i)
	    avolsum = avolsum + yvolaerosv(i)
	end do

	do i = 1, na
	    ynumaerosv(i) = ynumaerosv(i)/anumsum
	    yvolaerosv(i) = yvolaerosv(i)/avolsum
	end do


!
!   compute scavenging
!
	scavsumnum = 0.
	scavsumvol = 0.
!
!   outer loop for rain drop radius
!
	do 5900 jr = 1, nr

	r = rrainsv(jr)
	vfall = vfallrainsv(jr)

	reynolds = r * vfall / airkinvisc
	sqrtreynolds = sqrt( reynolds )

!
!   inner loop for aerosol particle radius
!
	scavsumnumbb = 0.
	scavsumvolbb = 0.

	do 5500 ja = 1, na

	a = aaerosv(ja)

	chi = a/r

	dum = freepath/a
	dumfuchs = 1. + 1.246*dum + 0.42*dum*exp(-0.87/dum)
	taurelax = 2.*rhoaero*a*a*dumfuchs/(9.*rhoair*airkinvisc)

	aeromass = 4.*pi*a*a*a*rhoaero/3.
	aerodiffus = boltzmann*temp*taurelax/aeromass

	schmidt = airkinvisc/aerodiffus
	stokes = vfall*taurelax/r

	ebrown = 4.*(1. + 0.4*sqrtreynolds*(schmidt**0.3333333)) /  &
     			(reynolds*schmidt)

	dum = (1. + 2.*xmuwaterair*chi) /         &
     			(1. + xmuwaterair/sqrtreynolds)
	eintercept = 4.*chi*(chi + dum)

	dum = log( 1. + reynolds )
	sstar = (1.2 + dum/12.) / (1. + dum)
	eimpact = 0.
	if (stokes .gt. sstar) then
	    dum = stokes - sstar
	    eimpact = (dum/(dum+0.6666667)) ** 1.5
	end if

	etotal = ebrown + eintercept + eimpact
	etotal = min( etotal, 1.0 )

	rainsweepout = xnumrainsv(jr)*4.*pi*r*r*vfall

	scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
	scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

5500	continue

	scavsumnum = scavsumnum + scavsumnumbb
	scavsumvol = scavsumvol + scavsumvolbb
5900	continue

	scavratenum = scavsumnum*3600.
	scavratevol = scavsumvol*3600.



  return
  end subroutine calc_1_impact_rate



!===========================================================================
!===========================================================================
subroutine gasrainscav (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,num_chem,  &
                      config_flags,   &
                      deltat,     t,          pmid,       pdel,       &
                      chem,                                  &
                      isprx,      fapx,       pfx,        pfx_inrain, &
                      dqdt,       dotend      )


!-----------------------------------------------------------------------
!
! Purpose:
! Does below cloud scavenging of gases by rain
!
! Currently does
!    Irreversible uptake of h2so4 and msa
!    Reactive uptake of so2 and h2o2.   This is assumed to be rate limited
!    by the mass transfer to rain (and not by aqueous reaction)
!
! Authors: R. Easter
!
!-----------------------------------------------------------------------
  USE module_model_constants, only: g,rhowater, mwdry
  use module_configure, only:  grid_config_rec_type,   &
		param_first_scalar,   &
		p_so2, p_h2o2, p_sulf,p_h2so4, p_msa,   &
		p_hno3, p_hcl, p_nh3

   implicit none

!-----------------------------------------------------------------------
!
! Input arguments
!
! abbreviations & acronyms
!    TMR = tracer mixing ratio
!    l-s = large scale
!    dp-cnv = deep convective
!    sh-cnv = shallow convective
!
   integer, intent(in)  :: num_chem           ! number of chemical species
   integer, intent(in)  :: ims,ime            ! column dimension
   integer, intent(in)  :: kms,kme            ! level dimension
   integer, intent(in)  :: jms,jme            ! column dimension
   integer, intent(in)  :: its,ite            ! column identifier
   integer, intent(in)  :: kts,kte            ! level identifier
   integer, intent(in)  :: jts,jte            ! column identifier
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   real, intent(in) :: deltat           ! model timestep

   real, intent(in) :: t(ims:ime,kms:kme,jms:jme)    ! temperature
   real, intent(in) :: pmid(ims:ime,kms:kme,jms:jme) ! pressure at model levels
   real, intent(in) :: pdel(ims:ime,kms:kme,jms:jme) ! pressure thickness of levels
   real, intent(in) :: chem(ims:ime,kms:kme,jms:jme,num_chem) ! TMR array including chemistry

   logical, intent(in)  :: isprx(ims:ime,kms:kme,jms:jme) ! true if precip at i,k
   real, intent(in) :: fapx(ims:ime,kms:kme,jms:jme)    ! frac. area for precip
   real, intent(in) :: pfx(ims:ime,kms:kme,jms:jme)     ! grid-avg precip
                                                 ! flux (kg/m2/s)
   real, intent(in) :: pfx_inrain(ims:ime,kms:kme,jms:jme)  !  precip flux (kg/m2/s)

   real, intent(inout) :: dqdt(ims:ime,kms:kme,jms:jme,num_chem) ! TMR tendency array
   logical,  intent(inout) :: dotend(num_chem)     ! flag for doing scav
!  real, intent(inout) :: qsrflx(ims:ime,jms:jme,num_chem)
                        ! changes to column tracer burdens by wet scavenging over current timestep
                        ! this routine adds on the contribution from gas scavenging
                        !    by mass transfer to rain

!--------------------------Local Variables------------------------------

   integer :: i, j, k         ! x, y, z work index
   integer :: jp              ! precip type index
   integer :: p1st
   logical :: ispr_anywhere

   integer, parameter :: ng = 7
   integer, parameter :: ig_so2   = 1
   integer, parameter :: ig_h2o2  = 2
   integer, parameter :: ig_h2so4 = 3
   integer, parameter :: ig_msa   = 4
   integer, parameter :: ig_hno3  = 5
   integer, parameter :: ig_hcl   = 6
   integer, parameter :: ig_nh3   = 7
   integer :: ig, lg, lg_ptr(ng)

   real :: amtscav(ng), amtscav_sub(ng)
   real :: fracgas(ng)
   real :: fracscav(ng), fracscav_sub(ng)
   real :: deltatinv
   real :: dum, dumamt, dumprecipmmh, dumpress, dumtemp
   real :: pdel_fac
   real :: r_gc(ng)
   real :: scavrate_hno3
   real :: scavrate(ng), scavrate_factor(ng)

!-----------------------------------------------------------------------
!

!  if (ncol .ne. -987654321) return

! precip rates -- 1.0 kgwtr/m2/s = 1.0e-3 m3wtr/m2/s = 1.0e-3 m/s
!                                                    = 1.0 mm/s = 3600 mm/h

   ispr_anywhere = .false.
   deltatinv = 1.0/(deltat*(1.0d0 + 1.0d-15))

   p1st = param_first_scalar

   lg_ptr(ig_so2  ) = p_so2  
   lg_ptr(ig_h2o2 ) = p_h2o2 
   lg_ptr(ig_h2so4) = p_sulf
   lg_ptr(ig_h2so4) = p_h2so4
   lg_ptr(ig_msa  ) = p_msa  
   lg_ptr(ig_hno3 ) = p_hno3 
   lg_ptr(ig_hcl  ) = p_hcl  
   lg_ptr(ig_nh3  ) = p_nh3  

   scavrate_factor(ig_so2  ) = 1.08
   scavrate_factor(ig_h2o2 ) = 1.38
   scavrate_factor(ig_h2so4) = 0.80
   scavrate_factor(ig_msa  ) = 0.80
   scavrate_factor(ig_hno3 ) = 1.00
   scavrate_factor(ig_hcl  ) = 1.15
   scavrate_factor(ig_nh3  ) = 1.59


   do 5900 j = jts,jte
   do 5900 i = its,ite

      do 4900 k = kts,kte

! skip this level if no precip
      if ( isprx(i,k,j)  ) then
        ispr_anywhere = .true.
      else
        goto 4900
      end if


! skip this level if below freezing
	dumtemp = t(i,k,j)
	if (dumtemp .le. 273.16) goto 4900
	dumpress = 10.0*pmid(i,k,j)

	do ig = 1, ng
	    fracscav(ig) = 0.0
	    fracgas(ig) = 1.0
	    lg = lg_ptr(ig)
	    if (lg .ge. p1st) then
		r_gc(ig) = max( chem(i,k,j,lg), 0.0 )
! activate this after gas_aqfrac is added to arguments
!		if (lg .le. numgas_aqfrac)   &
!			fracgas(ig) = gas_aqfrac(lg)
	    else
		r_gc(ig) = 0.0
	    end if
	end do
		
	if ( .not. isprx(i,k,j) ) goto 3600

!   precip rate in mm/h over rainy portion of the subarea
	dumprecipmmh = pfx_inrain(i,k,j)*3600.0

!   rain scavenging rate for hno3 (power law fit to schwarz and levine,
!   with temperature and pressure adjustments) -- units are (1/s)
	scavrate_hno3 = 6.262e-5*(dumprecipmmh**0.7366)   &
      		* ((dumtemp/298.0)**1.12)   &
      		* ((1.013e6/dumpress)**.75)

	do ig = 1, ng
	    scavrate(ig) = scavrate_hno3*scavrate_factor(ig)
	    fracscav_sub(ig) = (1. - exp(-scavrate(ig)*deltat))   &
			*fracgas(ig)*fapx(i,k,j)
	    amtscav_sub(ig) = r_gc(ig)*min( fracscav_sub(ig), 1.0 )
	end do

!   for so2 & h2o2, assume aqueous oxidation is fast, so reactive 
!   uptake is limited by the smaller of the two mass transfer rates
	dumamt = min( amtscav_sub(ig_so2), amtscav_sub(ig_h2o2) )
	fracscav_sub(ig_so2 ) = dumamt/max( r_gc(ig_so2 ), 1.0e-30 )
	fracscav_sub(ig_h2o2) = dumamt/max( r_gc(ig_h2o2), 1.0e-30 )
	amtscav_sub(ig_so2 ) = r_gc(ig_so2 )*min( fracscav_sub(ig_so2 ), 1.0 )
	amtscav_sub(ig_h2o2) = r_gc(ig_h2o2)*min( fracscav_sub(ig_h2o2), 1.0 )

!   for nh3, limit uptake by uptake of all acid gases combined
	dumamt = 2.0*amtscav_sub(ig_so2)   &
		+ 2.0*amtscav_sub(ig_h2so4) + amtscav_sub(ig_msa)   &
		+ amtscav_sub(ig_hno3) + amtscav_sub(ig_hcl)
	dumamt = min( dumamt, amtscav_sub(ig_nh3) )
	fracscav_sub(ig_nh3) = dumamt/max( r_gc(ig_nh3), 1.0e-30 )
	amtscav_sub(ig_nh3 ) = r_gc(ig_nh3 )*min( fracscav_sub(ig_nh3 ), 1.0 )

	do ig = 1, ng
	    fracscav(ig) = fracscav(ig) + fracscav_sub(ig)
	end do

!   diagnostic output
! 	write(lun,9440) nstep, lchnk, i, k, jp,   &
!		(dumtemp-273.16), dumpress*.001
!	write(lun,9442) 'pfx, pfx_inrain, fapx               ',   &
!		pfx(jp,i,k), pfx_inrain(jp,i,k), fapx(jp,i,k)
!	write(lun,9442) 'scavrate_so2, h2o2, msa, h2so4      ',   &
!		scavrate(ig_so2), scavrate(ig_h2o2),   &
!		scavrate(ig_msa), scavrate(ig_h2so4)
!	write(lun,9442) 'rso2gc, rso2g, rh2o2gc, rh2o2g      ',   &
!		r_gc(ig_so2), r_gc(ig_so2)*fracgas(ig)so2),   &
!		r_gc(ig_h2o2), r_gc(ig_h2o2)*fracgas(ig)h2o2),
!	write(lun,9442) 'amtscav_sub so2, h2o2               ',   &
!		amtscav_sub(ig_so2), amtscav_sub(ig_h2o2)
!	write(lun,9442) 'fracscav_sub so2, h2o2, msa, h2so4  ',   &
!		fracscav_sub(ig_so2), fracscav_sub(ig_h2o2),   &
!		fracscav_sub(ig_msa), fracscav_sub(ig_h2so4)
!9440	format( / 'ns,lc,i,k,jp,   T(C), p(mb)', i6, 2i4, 2i3, 2f7.1 )
!9442	format( a, 4(1pe11.3) )
!   end diagnostic output


3600	continue

!
!   compute tendencies
!
	pdel_fac = (pdel(i,k,j)/(g*mwdry))

	do ig = 1, ng
	    fracscav(ig) = max(0.0,min(1.0,fracscav(ig))) ! make sure 0 <= fracscav <= 1
	    amtscav(ig)  = fracscav(ig)*r_gc(ig)
	    lg = lg_ptr(ig)
	    if (lg .ge. p1st) then
		dqdt(i,k,j,lg) = -deltatinv*amtscav(ig)  
!czhao		qsrflx(i,j,lg) = qsrflx(i,j,lg) - pdel_fac*amtscav(ig) ! mmol/m2
	    end if
	end do

4900 continue  ! "k = 1, pver"

5900 continue  ! "i = 1, ncol"


! set dotend's
   if ( ispr_anywhere ) then
       do ig = 1, ng
           if (lg_ptr(ig) .ge. p1st) dotend(lg_ptr(ig)) = .true.
       end do
   end if


   return
end subroutine gasrainscav



!===========================================================================
!===========================================================================
	subroutine mlinft( x, y, a, n, m, mmaxd, rmserr )
!
!   fits y = a(1)*x(1) + a(2)*x(2) + ... + a(m)*x(m)
!
!	x - array containing x values
!	    x(i,k) is parameter i, observation k
!	y - array containing y values
!	    y(k) is observation
!	a - array !ontaining the regression coefficients
!	n - number of observations
!	m - number of parameters
!	mmaxd - first dimension of the x array
!	rmserr - root mean square residual
!	    rmserr = sqrt( avg-sq-err )
!	    avg-sq-err = (sum of residuals squared)/(number of values)
!	    residual = y - (a1*x1 + a2*x2 + ... + am*xm)
!
	implicit none

!   subr. parameters
	integer n, m, mmaxd
	real x(mmaxd,n), y(n), a(mmaxd), rmserr

!   local variables
	integer i, j, jflag, k
	real aa(10,10), bb(10), errsq, resid, ydum

	if (n .le. 1) then
	    a(1) = 1.e30
	    rmserr = 0.
	    return
	end if

	do 2900 i = 1, m
	    do 2100 j = 1, m
		aa(i,j) = 0.0
2100	    continue
	    bb(i) = 0.0

	    do 2500 k = 1, n
		do 2300 j = 1, m
		    aa(i,j) = aa(i,j) + x(i,k)*x(j,k)
2300		continue
		bb(i) = bb(i) + x(i,k)*y(k)
2500	    continue

2900	continue

!	do 4100 i = 1, m
!	    write(13,9300) i, (aa(i,j), j=1,m), bb(i)
!4100	continue
!9300	format( i5, 5f15.2 )

!	subr linsolv( a, x, b, n, m1, m2, jflag )
	call linsolv( aa, a, bb, m, 10, 10, jflag )


	errsq = 0.
	do 3300 k = 1, n
	    ydum = 0.0
	    do 3100 i = 1, m
		ydum = ydum + a(i)*x(i,k)
3100	    continue
	    resid = ydum - y(k)
	    errsq = errsq + resid*resid
3300	continue
	rmserr = sqrt( errsq/n )

	return
	end subroutine mlinft



!===========================================================================
!===========================================================================
	subroutine linsolv( a, x, b, n, m1, m2, jflag )
!
!   solves linear eqn system a*x = b using gaussian-elimination
!	with partial pivoting
!
!    n = order of the system
!    m1, m2 = fortran dimensions of a array
!    jflag = completion flag
!	1 - system solved successfully
!	0 - system is singular or close to it, and could not be solved.
!		computation was halted to avoid overflow or divide by zero.
!
!   *** note ***   rsmall should be defined as close to but somewhat larger
!	than the smallest single precision real on the computer.
!
!   initial coding on 29-aug-86 by r.c. easter
!   change on 4-feb-89 by r.c.easter -- added jflag to parameter list
!	and checking for singularity
!
	implicit none

!   subr. parameters
	integer n, m1, m2, jflag
	real a(m1,m2), x(n), b(n)

!   local variables
	integer i, imax, iup, j, k
	real amax, asmall, dmy, rsmall
	parameter (rsmall = 1.0e-16)

	jflag = 0

!
!   reduce coef. matrix to upper triangular form
!
	do 1900 k = 1, n
!
!   find pivot element, and
!   move pivot row into row k if necessary
!
	    imax = k
	    amax = abs( a(imax,k) )
	    do 1200 i = k+1, n
		if (abs(a(i,k)) .gt. amax) then
		    imax = i
		    amax = abs(a(i,k))
		end if
1200	    continue
	    if (amax .eq. 0.) return

	    if (imax .ne. k) then
		do 1400 j = k, n
		    dmy = a(imax,j)
		    a(imax,j) = a(k,j)
		    a(k,j) = dmy
1400		continue
		dmy = b(imax)
		b(imax) = b(k)
		b(k) = dmy
	    end if

!
!   reduce
!
	    asmall = abs(a(k,k))
	    do 1700 i = k+1, n
		if (a(i,k) .ne. 0.0) then
		    if (asmall .le. abs(rsmall*a(i,k))) return
		    dmy = a(i,k)/a(k,k)
		    a(i,k) = 0.0
		    do 1600 j = k+1, n
			a(i,j) = a(i,j) - dmy*a(k,j)
1600		    continue
		    b(i) = b(i) - dmy*b(k)
		end if
1700	    continue

1900	continue

!
!   backsolve
!
	do 2900 iup = 1, n
	    i = n + 1 - iup
	    dmy = b(i)
	    do 2500 j = i+1, n
		dmy = dmy - a(i,j)*x(j)
2500	    continue
	    asmall = abs(a(i,i))
	    if (abs(a(i,i)) .le. abs(rsmall*dmy)) return
	    x(i) = dmy/a(i,i)
2900	continue

	jflag = 1

	return
	end subroutine linsolv


!===========================================================================
!subroutine wetscav_largescale (ims,ime,kms,kme,jms,jme,its,ite,kts,kte,jts,jte,num_chem,  &
!    ntype_aer, nsize_aer, ncomp_aer, massptr_aer, dens_aer, numptr_aer, &
!    maxd_acomp, maxd_asize,maxd_atype, maxd_aphase, ai_phase,           &
!    volumcen_sect, volumlo_sect, volumhi_sect,                          &
!    waterptr_aer, dens_water_aer,                                       &
!    scavimptblvol, scavimptblnum, nimptblgrow_mind, nimptblgrow_maxd, dlndg_nimptblgrow, &
!    deltat,     t,          pmid,       pdel, chem,                     &
!    isprx,      fapx,       pfx,        pfx_inrain,                     &
!    dqdt,       dotend                                      )


!-----------------------------------------------------------------------
! Purpose:
! Does in-cloud and below-cloud scavenging of aerosols following GEOS-Chem 
!
! Authors: Chun Zhao, PNNL
!
! References:
! (1) Liu,H., D.J. Jacob, I. Bey and R.M. Yantosca, "Constraints from 210Pb 
!        and 7Be on wet deposition and transport in a global three-dimensional
!        chemical tracer model driven by assimilated meteorological fields", 
!        JGR, Vol 106, pp 12109-12128, 2001.
! (2) D.J. Jacob, H. Liu, C. Mari, and R. M. Yantosca, "Harvard wet 
!        deposition scheme for GMI", Harvard Atmospheric Chemistry Modeling 
!        Group, March 2000.
!-----------------------------------------------------------------------
!******************************************************************************
!  Subroutine computes the downward mass flux of tracer due to washout 
!  and rainout of aerosols and soluble tracers in a column.  The timestep is 
!  the dynamic timestep.
!
!  The precip fields through the bottom of each level are indexed as follows:
!
!       Layer          GISS-CTM II         GEOS-CTM
!
!      ------------------------------------------------- Top of Atm.
!        LM            PSSW4(I,J,LM-1)   PDOWN(LM,I,J)
!                          |                  |
!      ====================V==================V========= Max Extent 
!        LM-1          PSSW4(I,J,LM)     PDOWN(LM-1,I,J)  of Clouds
!                          |                  |
!      --------------------V------------------V---------
!                         ...                ...             
!
!      -------------------------------------------------
!        4             PSSW4(I,J,3)      PDOWN(4,I,J)
!                          |                  |
!      --------------------V------------------V----------
!        3             PSSW4(I,J,2)      PDOWN(3,I,J)
!                          |                  |
!      --------------------V------------------V--------- Cloud base
!        2             PSSW4(I,J,1)      PDOWN(2,I,J) 
!                          |                  |
!      -  -  -  -  -  -  - V -  -   -   -   - V -  -  - 
!        1                               PDOWN(1,I,J) 
!                                             |
!      =======================================V========= Ground
!
!  From the diagram, we have the following for layer L:
!     
!  GISS-CTM:
!  (a) Precip coming in  thru top    of layer L = PSSW4(I,J,L  )
!  (b) Precip going  out thru bottom of layer L = PSSW4(I,J,L-1)
!
!  GEOS-CHEM
!  (a) Precip coming in  thru top    of layer L = PDOWN(L+1,I,J)
!  (b) Precip going  out thru bottom of layer L = PDOWN(L,  I,J) 
!
!  Thus: Precip coming in:  PSSW4(I,J,L  ) is analogous to PDOWN(L+1,I,J).
!        Precip going  out: PSSW4(I,J,L-1) is analogous to PDOWN(L,I,J  ).
!
!******************************************************************************




! USE module_model_constants, only: g,rhowater, mwdry
! USE module_state_description, only: param_first_scalar

!  implicit none

!  INTEGER,      INTENT(IN   )    ::                                &
!                                     ids,ide, jds,jde, kds,kde,    &
!                                     ims,ime, jms,jme, kms,kme,    &
!                                     its,ite, jts,jte, kts,kte,    &


!  integer, intent(in)  :: num_chem           ! number of chemical species
!  real, intent(in) :: deltat           ! model timestep

!  !=================================================================
!  ! WETDEP begins here!
!  !
!  ! (1)  I n i t i a l i z e   V a r i a b l e s
!  !=================================================================


!  !=================================================================
!  ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
!  !
!  ! Process rainout / washout by columns.
!  !=================================================================

!  DO J = jts,jte 
!  DO I = its,ite

!        ! Zero FTOP
!        FTOP = 0d0

!        ! Zero accumulator array
!        DO L  = 1, LLPAR
!        DO NN = 1, NSOL
!           DSTT(NN,L,I,J) = 0d0
!        ENDDO
!        ENDDO

!  !==============================================================
!  ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR) 
!  !
!  ! Assume that rainout is happening in the top layer if 
!  ! QQ(LLPAR,I,J) > 0.  In other words, if any precipitation 
!  ! forms in grid box (I,J,LLPAR), assume that all of it falls 
!  ! down to lower levels.
!  !
!  ! Soluble gases/aerosols are incorporated into the raindrops 
!  ! and are completely removed from grid box (I,J,LLPAR).  There 
!  ! is no evaporation and "resuspension" of aerosols during a 
!  ! rainout event.
!  !
!  ! For large-scale (a.k.a. stratiform) precipitation, the first 
!  ! order rate constant for rainout in the grid box (I,J,L=LLPAR) 
!  ! (cf. Eq. 12, Jacob et al, 2000) is given by:
!  !                        Q        
!  !    K_RAIN = K_MIN + -------    [units: s^-1]
!  !                      L + W    
!  ! and the areal fraction of grid box (I,J,L=LLPAR) that 
!  ! is actually experiencing large-scale precipitation 
!  ! (cf. Eq. 11, Jacob et al, 2000) is given by: 
!  !                  Q               
!  !    F' =  -------------------   [unitless]
!  !           K_RAIN * ( L + W )    
!  ! Where:
!  !
!  !    K_MIN  = minimum value for K_RAIN         
!  !           = 1.0e-4 [s^-1]
!  !
!  !    L + W  = condensed water content in cloud 
!  !           = 1.5e-6 [cm3 H2O/cm3 air]
!  !
!  !    Q = QQ = rate of precipitation formation 
!  !             [ cm3 H2O / cm3 air / s ]
!  !
!  ! For convective precipitation, K_RAIN = 5.0e-3 [s^-1], and the
!  ! expression for F' (cf. Eq. 13, Jacob et al, 2000) becomes:
!  !                                  { DT        }
!  !                    FMAX * Q * MIN{ --- , 1.0 }
!  !                                  { TAU       }
!  !  F' = ------------------------------------------------------
!  !               { DT        }
!  !        Q * MIN{ --- , 1.0 }  +  FMAX * K_RAIN * ( L + W )
!  !               { TAU       } 
!  ! Where:
!  !    Q = QQ = rate of precipitation formation 
!  !             [cm3 H2O/cm3 air/s]
!  !    FMAX   = maximum value for F' 
!  !           = 0.3
!  !    DT     = dynamic time step from the CTM [s]
!  !    TAU    = duration of rainout event 
!  !           = 1800 s (30 min)
!  !    L + W  = condensed water content in cloud 
!  !           = 2.0e-6 [cm3 H2O/cm3 air]
!  ! K_RAIN and F' are needed to compute the fraction of tracer
!  ! in grid box (I,J,L=LLPAR) lost to rainout.  This is done in 
!  ! module routine RAINOUT.
!  !==============================================================

!        ! Zero variables for this level
!        ALPHA       = 0d0
!        ALPHA2      = 0d0
!        F           = 0d0
!        F_PRIME     = 0d0
!        GAINED      = 0d0
!        K_RAIN      = 0d0
!        LOST        = 0d0
!        Q           = 0d0
!        QDOWN       = 0d0
!        MASS_NOWASH = 0d0
!        MASS_WASH   = 0d0
!        RAINFRAC    = 0d0
!        WASHFRAC    = 0d0
!        WETLOSS     = 0d0

!        ! Start at the top of the atmosphere
!        L = LLPAR

!        ! If precip forms at (I,J,L), assume it all rains out
!        IF ( QQ(L,I,J) > 0d0 ) THEN

!           ! Q is the new precip that is forming within grid box (I,J,L)
!           Q = QQ(L,I,J)

!           ! Compute K_RAIN and F' for either large-scale or convective
!           ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000) 
!           IF ( LS ) THEN
!              K_RAIN  = LS_K_RAIN( Q )
!              F_PRIME = LS_F_PRIME( Q, K_RAIN )
!           ELSE
!              K_RAIN  = 1.5d-3
!              F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
!           ENDIF

!           ! Set F = F', since there is no FTOP at L = LLPAR
!           F = F_PRIME

!           ! Only compute rainout if F > 0. 
!           ! This helps to eliminate unnecessary CPU cycles.
!           IF ( F > 0d0 ) THEN

!              ! ND16 diagnostic...save LS and Conv fractions
!              IF ( ND16 > 0 .and. L <= LD16 ) THEN
!                 AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                 CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!              ENDIF

!              ! ND17 diagnostic...increment counter
!              IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                 CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
!              ENDIF

!              ! Loop over soluble tracers and/or aerosol tracers    
!              DO NN = 1, NSOL
!                 N = IDWETD(NN)

!                 ! Call subroutine RAINOUT to compute the fraction
!                 ! of tracer lost to rainout in grid box (I,J,L=LLPAR)
!                 CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )

!                 ! WETLOSS is the amount of soluble tracer 
!                 ! lost to rainout in grid box (I,J,L=LLPAR)
!                 WETLOSS = STT(I,J,L,N) * RAINFRAC

!                 ! Remove rainout losses in grid box (I,J,L=LLPAR) from STT
!                 STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

!                 ! DSTT is an accumulator array for rained-out tracers.  
!                 ! The tracers in DSTT are in the liquid phase and will 
!                 ! precipitate to the levels below until a washout occurs.
!                 ! Initialize DSTT at (I,J,L=LLPAR) with WETLOSS.
!                 DSTT(NN,L,I,J) = WETLOSS

!                 ! ND17 diagnostic...LS and conv rainout fractions [unitless]
!                 IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                    AD17(I,J,L,NN,IDX) =
!    &                    AD17(I,J,L,NN,IDX) + RAINFRAC / F
!                 ENDIF

!                 ! ND39 diag - save rainout losses in [kg/s]
!                 ! Add LGTMM in condition for AD39 (ccc, 11/18/09)
!                 IF ( ( ND39 > 0 .or. LGTMM ) .and. L <= LD39 ) THEN
!                    AD39(I,J,L,NN) = AD39(I,J,L,NN) + WETLOSS / DT
!                 ENDIF

!                 ! Negative tracer...call subroutine SAFETY
!                 IF ( STT(I,J,L,N) < 0d0 ) THEN
!                    CALL SAFETY( I, J, L, N, 3,
!    &                            LS,             PDOWN(L,I,J),
!    &                            QQ(L,I,J),      ALPHA,
!    &                            ALPHA2,         RAINFRAC,
!    &                            WASHFRAC,       MASS_WASH,
!    &                            MASS_NOWASH,    WETLOSS,
!    &                            GAINED,         LOST,
!    &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
!                 ENDIF
!              ENDDO
!           ENDIF

!           ! Save FTOP for the next lower level 
!           FTOP = F
!        ENDIF

!        !==============================================================
!        ! (4)  R a i n o u t   i n   t h e   M i d d l e   L e v e l s
!        ! 


!end subroutine wetscav_largescale 

END MODULE module_mosaic_wetscav
