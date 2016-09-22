MODULE module_mosaic_settling

!---------------------------------------------------------------------------------
! Purpose:
! calculate gravitational settling for MOSAIC aerosols 
!
! Author: Chun Zhao 
! v1 : 2016/05/30
!---------------------------------------------------------------------------------


CONTAINS

SUBROUTINE mosaic_settling_driver(dt,config_flags,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
         ids,ide, jds,jde, kds,kde,                 &
         ims,ime, jms,jme, kms,kme,            &
         its,ite, jts,jte, kts,kte                   )

  USE module_configure
  USE module_state_description
  USE module_model_constants, ONLY: mwdry
  IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) ::                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::  t_phy,p_phy,dz8w,p8w,rho_phy
  REAL, INTENT(IN   ) :: dt,dx,g

!  INTEGER, PARAMETER :: settling_option=1

  settling_scheme: SELECT CASE(config_flags%settling_option)
   
  CASE (0)
      CALL wrf_debug(15,'no gravitational settling scheme called for mosaic aerosol')

  CASE (1)
      CALL wrf_debug(15,'gocart gravitational settling scheme called for mosaic aerosol')
       call gocart_settling_scheme(dt,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
         ids,ide, jds,jde, kds,kde,                  &
         ims,ime, jms,jme, kms,kme,             &
         its,ite, jts,jte, kts,kte            )

  CASE (2)
      ! call cam_settling()

  CASE DEFAULT
  END SELECT settling_scheme


END SUBROUTINE mosaic_settling_driver


SUBROUTINE gocart_settling_scheme(dt,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,p_phy,dx,g, &
         ids,ide, jds,jde, kds,kde,                 &
         ims,ime, jms,jme, kms,kme,            &
         its,ite, jts,jte, kts,kte                   )

  USE module_configure
  USE module_state_description
  USE module_data_mosaic_therm, only: nbin_a
  USE module_data_mosaic_asect, only: ncomp_aer,nsize_aer,massptr_aer,hygro_aer,dens_aer,dcen_sect
  USE module_data_gocart_dust, only: dyn_visc
  USE module_model_constants, ONLY: mwdry
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) ::                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist  
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
    REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),        &
          INTENT(IN   ) ::  t_phy,p_phy,dz8w,p8w,rho_phy
    REAL, INTENT(IN   ) :: dt,dx,g

!   REAL, DIMENSION( ims:ime, kms:kme, jms:jme,1:nbin_a),                 &
!        INTENT(INOUT ) :: buget_settling_mass 

 integer :: kkk,nmx,i,j,k,kk,lmx,l
! real*8, DIMENSION (kms:kme) :: tmp,airden,airmas,p_mid,delz,rh
  real*8, DIMENSION (kms:kme) :: tmp,airmas,p_mid,delz,rh
!
  real*8 conver,converi

  REAL*8, DIMENSION(kms:kme,1:nbin_a) :: & 
  totmass0,totmass,bstl_mass,avghygro,avgden,avgreff
  REAL*8 :: sumhygro,sumvol,lmass
  INTEGER :: ITYPE, IPHASE, ISIZE

    conver=1.e-9  !ug/kg to kg/kg
    converi=1.e9
    lmx=kte-kts+1

    do j=jts,jte
    do i=its,ite

      totmass(:,:) = 0.0
      bstl_mass(:,:)=0.  !ug/kg
      avghygro(:,:)=0.
      itype=1
      iphase=1
      do k=kts,kte
       do isize = 1, nbin_a
       !total mass and volume-weighted hygroscopicity for each type and size
         sumhygro=0.
         sumvol=0.
         do l=1,ncomp_aer(itype)
          lmass=massptr_aer(l,isize,itype,iphase)
          sumhygro = sumhygro+hygro_aer(l,itype)*   &
              chem(i,k,j,lmass)/(dens_aer(l,itype)*1.e12) !dens from g/cm3 to ug/m3
          sumvol = sumvol+chem(i,k,j,lmass)/(dens_aer(l,itype)*1.e12) !m3/kg-air
          totmass(k,isize) = totmass(k,isize)+chem(i,k,j,lmass)*conver !ug/kg to kg/kg
         end do ! comp
         avghygro(k,isize)=sumhygro/sumvol
         avgden(k,isize)=totmass(k,isize)/sumvol ! kg/m3 
         avgreff(k,isize)=dcen_sect(isize,itype)*0.01 !cm to m
        end do ! size
      enddo !k 

      do k=kts,kte
        p_mid(k)=.01*p_phy(i,k,j) ! Pa to hPa
        delz(k)=dz8w(i,k,j)
        airmas(k)=-(p8w(i,k+1,j)-p8w(i,k,j))*dx*dx/g  !kg
!       airden(k)=rho_phy(i,k,j)
        tmp(k)=t_phy(i,k,j)
        rh(k) = .95
        rh(k) = MIN(0.95, moist(i,k,j,p_qv) / &
           (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
           (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
        rh(k)=max(1.0D-1,rh(k))
      enddo

      totmass0=totmass
      call settling(kts,kte,kms,kme, &
                  nbin_a,g,dyn_visc, &
                  totmass, tmp, p_mid, delz, airmas, &
                  avghygro,avgden,avgreff, dt, bstl_mass, rh)

      do k=kts,kte
       do isize = 1, nbin_a
         !scale each aerosol component for total mass 
         do l=1,ncomp_aer(itype)
          lmass=massptr_aer(l,isize,itype,iphase)   
          chem(i,k,j,lmass) = chem(i,k,j,lmass)*totmass(k,isize)/totmass0(k,isize) !ug/kg
         end do ! comp
       end do ! isize 

       if (i==20.and.j==20.and.(k==kts.or.k==3.or.k==10.or.k==kte)) then 
       write(0,*) 'check grav settling budget for size 1,7,8',k,bstl_mass(k,1),bstl_mass(k,nbin_a-1),bstl_mass(k,nbin_a),totmass(k,nbin_a)/totmass0(k,nbin_a)
       print*, 'check grav settling budget for size 1,7,8',k,bstl_mass(k,1),bstl_mass(k,nbin_a-1),bstl_mass(k,nbin_a),totmass(k,nbin_a)/totmass0(k,nbin_a)
       endif

      end do ! k

    enddo !i
    enddo !j

END SUBROUTINE gocart_settling_scheme

subroutine settling(kts,kte,kms,kme, &
                    nmx,g0,dyn_visc, &
                    tc, tmp, p_mid, delz, airmas, &
                    hygro,den,reff, dt, bstl, rh)

! ****************************************************************************
! *                                                                          *
! *  Calculate the loss by settling, using an implicit method                *
! *                                                                          *
! *  Input variables:                
! *
! *    SIGE(k)         - sigma coordinate of the vertical edges              *
! *    PS(i,j)         - Surface pressure (mb)                               *
! *    TMP(i,j,k)      - Air temperature  (K)                                *
! *    CT(i,j)         - Surface exchange coeff for moisture
! *                                                                          *
! **************************************************************************** 

  USE module_data_mosaic_asect, only: hygro_na_aer, hygro_dust_aer 


  IMPLICIT  NONE

  INTEGER, INTENT(IN) :: kts,kte,kms,kme, nmx
  INTEGER :: ntdt
  REAL, INTENT(IN) :: dt,g0,dyn_visc
  REAL*8,    INTENT(IN) :: tmp(kms:kme), delz(kms:kme), p_mid(kms:kme),  &
                         airmas(kms:kme), rh(kms:kme), &
                         den(kms:kme,1:nmx), reff(kms:kme,1:nmx),hygro(kms:kme,1:nmx)
  REAL*8, INTENT(INOUT) :: tc(kms:kme,1:nmx)
  REAL*8, INTENT(OUT)   :: bstl(kms:kme,1:nmx)

  REAL*8    :: tc1(kms:kme,1:nmx), dt_settl(1:nmx), rcm(kms:kme,1:nmx), rho(kms:kme,1:nmx)
  INTEGER :: ndt_settl(1:nmx)
  REAL*8    :: dzmin, vsettl(kms:kme), dtmax, pres, rhb, rwet(kms:kme,1:nmx), ratio_r(kms:kme,1:nmx)
  REAL*8    :: c_stokes, free_path, c_cun, viscosity, vd_cor(kms:kme)
  REAL*8    :: growth_fac_dust,growth_fac_seas,growth_fac(kms:kme,1:nmx)
  INTEGER :: k, n, i, j, l, l2,k2,isize

  ! for sea-salt:
  REAL*8, PARAMETER :: c1=0.7674, c2=3.079, c3=2.573E-11, c4=-1.424

  ! for OMP:
  REAL*8 :: rwet_priv_dust(kms:kme,1:nmx), rho_priv_dust(kms:kme,1:nmx)
  REAL*8 :: rwet_priv_seas(kms:kme,1:nmx), rho_priv_seas(kms:kme,1:nmx)
  REAL*8 :: rwet_priv(kms:kme,1:nmx), rho_priv(kms:kme,1:nmx)

  ! executable statements

! dzmin = MINVAL(delz(:))

  DO n = 1,nmx

    growth_fac_dust = 1.0
    growth_fac_seas = 3.0
    !internal-mixed aerosol growth factor is scaled between dust and seasalt
    do k=kts,kte
      growth_fac(k,n) = growth_fac_dust + & 
      ((hygro(k,n)-hygro_dust_aer)/(hygro_na_aer-hygro_dust_aer))* &
      (growth_fac_seas-growth_fac_dust)
    enddo

      ! Settling velocity (m/s) for each tracer (Stokes Law)
      ! DEN       density                        (kg/m3)
      ! REFF      effective radius               (m)
      ! dyn_visc  dynamic viscosity              (kg/m/s)
      ! g0        gravity                        (m/s2)
      ! 3.0       corresponds to a growth of a factor 3 of radius with 100% RH
      ! 0.5       upper limit with temp correction

      tc1(:,n) = tc(:,n)
    do k=kts,kte
      vsettl(k) = 2.0/9.0 * g0 * den(k,n) * (growth_fac(k,n)*reff(k,n))**2 / &
              (0.5*dyn_visc)
    enddo

     ! Determine the maximum time-step satisying the CFL condition:
     ! dt <= (dz)_min / v_settl
     ntdt=INT(dt)
     dtmax = minval(delz(kts:kte)/vsettl(kts:kte))
     ndt_settl(n) = MAX(1, INT(ntdt/dtmax) )
     ! limit maximum number of iterations
     IF (ndt_settl(n) > 15) ndt_settl(n) = 15
     dt_settl(n) = REAL(ntdt) / REAL(ndt_settl(n))

     ! Particles radius in centimeters
     rcm(:,n) = reff(:,n)*100.0
!cz    IF (idust.eq.1)then
     rwet(:,n) = reff(:,n)
!         ratio_r(k) = 1.0
     rho(:,n) = den(:,n)
!cz     endif
  END DO !n

  ! Solve the bidiagonal matrix (l,l)

!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( i,   j,   l,   l2, n,   k,   rhb, rwet_priv, ratio_r, c_stokes)&
!$OMP PRIVATE( free_path, c_cun, viscosity, rho_priv, vd_cor )

  DO isize = 1,nmx
!cz       IF (idust.eq.1) THEN
!         rwet_priv(k) = rwet(k)
!         rho_priv(k)  = rho(k)
!cz       END IF

   DO n = 1,ndt_settl(isize)

      ! Solve each vertical layer successively (layer l)
      DO k = kte,kts,-1
        k2 = kte - k + kts 

        rhb = MIN(9.9D-1, rh(k))
        ! Aerosol growth with relative humidity (Gerber, 1985)
        rwet_priv_dust(k,isize) = rwet(k,isize)
        rho_priv_dust(k,isize)  = rho(k,isize)

        rwet_priv_seas(k,isize) = 0.01*(c1*rcm(k,isize)**c2/(c3*rcm(k,isize)**c4 - &
                 LOG10(rhb)) + rcm(k,isize)**3)**0.33

        rwet_priv(k,isize) = rwet_priv_dust(k,isize) + &
        ((hygro(k,isize)-hygro_dust_aer)/(hygro_na_aer-hygro_dust_aer))* &
        (rwet_priv_seas(k,isize)-rwet_priv_dust(k,isize))
 
        ratio_r(k,isize) = (reff(k,isize)/rwet_priv(k,isize))**3.0

        rho_priv(k,isize) = ratio_r(k,isize)*den(k,isize) + (1.0 - ratio_r(k,isize))*1000.0

        ! Dynamic viscosity
        c_stokes = 1.458E-6 * tmp(k)**1.5/(tmp(k) + 110.4)

        ! Mean free path as a function of pressure (mb) and 
        ! temperature (K)
        ! order of p_mid is top->sfc
        free_path = 1.1E-3/p_mid(k)/SQRT(tmp(k))
!!!     free_path = 1.1E-3/p_edge(i,j,l2)/SQRT(tmp(i,j,l))

        ! Slip Correction Factor
        c_cun = 1.0+ free_path/rwet_priv(k,isize)* &
           (1.257 + 0.4*EXP(-1.1*rwet_priv(k,isize)/free_path))

        ! Corrected dynamic viscosity (kg/m/s)
        viscosity = c_stokes / c_cun

        ! Settling velocity
!cz       IF (iseas.eq.1) THEN
!cz           rho_priv(k) = ratio_r(k)*den(k) + (1.0 - ratio_r(k))*1000.0
!cz       END IF
        vd_cor(k) = 2.0/9.0*g0*rho_priv(k,isize)*rwet_priv(k,isize)**2/viscosity

        ! Update mixing ratio in an implicit method -czhao
        ! Order of delz is top->sfc
        IF (k == kte) THEN !only remove
         tc(k,isize) = tc(k,isize) / &
            (1.0 + dt_settl(isize)*vd_cor(k)/delz(k))
        ELSE IF (k == kts) THEN ! only receive, the removal of the surface layer will be done with dry deposition part, avoid double counting
         tc(k,isize) = tc(k,isize)+tc(k+1,isize)*dt_settl(isize)*vd_cor(k+1)/delz(k+1)
        ELSE 
         tc(k,isize) = 1.0/(1.0+dt_settl(isize)*vd_cor(k)/delz(k)) &
         *(tc(k,isize)+tc(k+1,isize)*dt_settl(isize)*vd_cor(k+1)/delz(k+1))
        END IF

        END DO  !k
     END DO  !n
  END DO  !isize

!$OMP END PARALLEL DO

  DO n = 1,nmx
   DO k = kts,kte 
     bstl(k,n) = 0.0
     IF (tc(k,n) < 1.0D-32) tc(k,n) = 1.0D-32
     bstl(k,n)=(tc(k,n)-tc1(k,n))*1.e9 !kg/kg to ug/kg
   END DO
  END DO

END SUBROUTINE settling


END MODULE module_mosaic_settling
