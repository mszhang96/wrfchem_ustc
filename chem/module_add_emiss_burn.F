Module module_add_emiss_burn
CONTAINS
       subroutine add_emis_burn(id,dtstep,ktau,dz8w,rho_phy,chem,       &
            julday,gmt,xlat,xlong,t_phy,p_phy,                          &
            ebu,chem_opt,tracer_opt,biomass_burn_opt,                   &
            do_plumerisefire,   &
            num_c,ids,ide, jds,jde, kds,kde,                            &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )
  USE module_configure, only: grid_config_rec_type
  USE module_state_description
  IMPLICIT NONE


!  TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,julday,chem_opt,biomass_burn_opt,     &
                                  num_c,ids,ide, jds,jde, kds,kde,      &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte,tracer_opt
   INTEGER,      INTENT(IN   ) ::                                          &
                                  ktau
   LOGICAL,      INTENT(IN   ) ::  do_plumerisefire 
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_c ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_ebu ),                 &
         INTENT(IN    ) ::                           ebu
!
!
!
   REAL,  DIMENSION( ims:ime ,  jms:jme )         ,               &
          INTENT(IN   ) ::                                                 &
                                                      xlat,xlong
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::                                                 &
                                                      t_phy,               &
                                                      p_phy,               &
                                                      dz8w,                &
                                                    rho_phy


   !czhao distribute the biomass burning following AEROCOM up to 6000m
   INTEGER,PARAMETER :: firelevel=6 
   REAL  total_emis(firelevel),total_air(firelevel),bbemis_frac(firelevel)
   REAL  fireheight(firelevel+1)
   INTEGER :: nk,emisbot(firelevel),emistop(firelevel)

   REAL,      INTENT(IN   ) ::                                          &
                             dtstep,gmt
    integer ::imonth1,idate1,iyear1,itime1
    integer :: i,j,k
    real :: time,conv_rho
    integer :: iweek,idays
    real :: tign,timeq,r_q,r_antro
    real, dimension(7) :: week_CYCLE
    !                     dia da semana:  DOM   SEG   TER   QUA   QUI   SEX   SAB
    !                            iweek=   1     2     3     4     5     6     7
    !- dados cetesb/campinas/2005
    data (week_CYCLE(iweek),iweek=1,7) /0.67, 1.1, 1.1, 1.1, 1.1, 1.1, 0.83/ !total = 7
    real, parameter :: bx_bburn  = 18.041288 * 3600., & !- peak at 18 UTC
                  cx        =  2.184936 * 3600., &
                  rinti     =  2.1813936e-8    , &
                  ax        = 2000.6038        , &
                  bx_antro  = 15.041288 * 3600.    !- peak em 15 UTC
    !itime1 : initial time of simulation (hour*100)
    ! time  : time elapsed in seconds
    ! r_q : gaussian function in 1/sec

    !-------------biomass burning diurnal cycle --------------------
    !number of days of simulation
    itime1=0
    time=0.
    idays = int(( float(itime1)/100. + time/3600.)/24.+.00001)
    tign  = real(idays)*24.*3600.
    ! Modulacao da queimada media durante o ciclo diurno(unidade: 1/s)
    ! com a int( r_q dt) (0 - 24h)= 1.
    timeq= ( time + float(itime1)*0.01*3600. - tign )
    timeq=gmt*3600.+float(ktau)*dtstep
    timeq=mod(timeq,86400.)
    r_q  = rinti*( ax * exp( -(timeq-bx_bburn)**2/(2.*cx**2) ) + 100. -  &
           5.6712963e-4*( timeq ))


    !------------- anthropogenic diurnal cycle (industrial,residencial, ...)
    ! weekly cycle
    ! week day
    iweek= int(((float(julday)/7. - &
           int(julday/7))*7.)) + 1
    if(iweek.gt.7) iweek = iweek-7
    !- weekly + diurnal cycle
    r_q=r_q*86400.
    r_q=1. ! no diurnal cycle
!
!tracer_opt has to come in as zero, if chem_opt is not =0!
      temiss_select:  SELECT CASE(tracer_opt)
         CASE (TRACER_SMOKE)
!
! for smoke only
!
          do j=jts,jte
          do i=its,ite
          do k=kts,kte
             conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
             chem(i,k,j,p_smoke)  =  chem(i,k,j,p_smoke)+ebu(i,k,j,p_ebu_co)*conv_rho
          enddo
          enddo
          enddo
         CASE (TRACER_TEST2)
          do j=jts,jte
          do i=its,ite
          do k=kts,kte
             conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
             chem(i,k,j,p_tr17_3)  =  chem(i,k,j,p_tr17_3)+ebu(i,k,j,p_ebu_co)*conv_rho
             chem(i,k,j,p_tr17_4)  =  chem(i,k,j,p_tr17_4)+ebu(i,k,j,p_ebu_co)*conv_rho
          enddo
          enddo
          enddo
         CASE DEFAULT
             call wrf_debug(15,'nothing done with burn emissions for tracers here')
      END SELECT temiss_select

      ! no biomass burning
      IF (biomass_burn_opt .eq. 0) THEN
          call wrf_debug(15,'fire emission, return from biomass_add')
          return
      ENDIF

      emiss_select:  SELECT CASE(chem_opt)
      CASE (RACMPM_KPP)
          do j=jts,jte
          do i=its,ite
           do k=kts,kte
        conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
        chem(i,k,j,p_so2)  =  chem(i,k,j,p_so2)                        &
                         +ebu(i,k,j,p_ebu_so2)*conv_rho
        chem(i,k,j,p_sulf)  =  chem(i,k,j,p_sulf)                        &
                         +ebu(i,k,j,p_ebu_sulf)*conv_rho
        chem(i,k,j,p_csl)  =  chem(i,k,j,p_csl)                        &
                         +ebu(i,k,j,p_ebu_csl)*conv_rho
        chem(i,k,j,p_iso)  = chem(i,k,j,p_iso)                         &
                         +ebu(i,k,j,p_ebu_iso)*conv_rho
        chem(i,k,j,p_no)   = chem(i,k,j,p_no)                          &
                         +ebu(i,k,j,p_ebu_no)*conv_rho
        chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)                         &
                         +ebu(i,k,j,p_ebu_no2)*conv_rho
        chem(i,k,j,p_ald)  = chem(i,k,j,p_ald)                         &
                         +ebu(i,k,j,p_ebu_ald)*conv_rho
        chem(i,k,j,p_hcho) = chem(i,k,j,p_hcho)                        &
                         +ebu(i,k,j,p_ebu_hcho)*conv_rho
        chem(i,k,j,p_ora2)  = chem(i,k,j,p_ora2)                       &
                         +ebu(i,k,j,p_ebu_ora2)*conv_rho
        chem(i,k,j,p_hc3)  = chem(i,k,j,p_hc3)                         &
                         +ebu(i,k,j,p_ebu_hc3)*conv_rho
        chem(i,k,j,p_hc5)  = chem(i,k,j,p_hc5)                         &
                         +ebu(i,k,j,p_ebu_hc5)*conv_rho
        chem(i,k,j,p_hc8)  = chem(i,k,j,p_hc8)                         &
                         +ebu(i,k,j,p_ebu_hc8)*conv_rho
        chem(i,k,j,p_eth)  = chem(i,k,j,p_eth)                         &
                         +ebu(i,k,j,p_ebu_eth)*conv_rho
        chem(i,k,j,p_co)  = chem(i,k,j,p_co)                           &
                         +ebu(i,k,j,p_ebu_co)*conv_rho
        chem(i,k,j,p_olt)  = chem(i,k,j,p_olt)                         &
                         +ebu(i,k,j,p_ebu_olt)*conv_rho
        chem(i,k,j,p_oli)  = chem(i,k,j,p_oli)                         &
                         +ebu(i,k,j,p_ebu_oli)*conv_rho
        chem(i,k,j,p_tol)  = chem(i,k,j,p_tol)                         &
                         +ebu(i,k,j,p_ebu_tol)*conv_rho
        chem(i,k,j,p_xyl)  = chem(i,k,j,p_xyl)                         &
                         +ebu(i,k,j,p_ebu_xyl)*conv_rho
        chem(i,k,j,p_ket)  =  chem(i,k,j,p_ket)                        &
                         +ebu(i,k,j,p_ebu_ket)*conv_rho
        chem(i,k,j,p_pm_25)  =  chem(i,k,j,p_pm_25)                        &
                         +r_q*ebu(i,k,j,p_ebu_pm25)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_pm_10)  =  chem(i,k,j,p_pm_10)                        &
                         +r_q*ebu(i,k,j,p_ebu_pm10)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        enddo
        enddo
        enddo
      CASE (RADM2,RACM_KPP,RACM_MIM_KPP,RADM2SORG,RACMSORG_KPP, RADM2SORG_KPP, RACM_ESRLSORG_KPP, RACM_SOA_VBS_KPP,  &
            RADM2SORG_AQ, RADM2SORG_AQCHEM, RACMSORG_AQ, RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP)
        do j=jts,jte
        do i=its,ite

          if (kts.ne.1) then
           print*, 'kts not start from the first layer'
           stop
          endif

          fireheight(1)=0
          fireheight(2)=100
          fireheight(3)=500
          fireheight(4)=1000
          fireheight(5)=2000
          fireheight(6)=3000
          fireheight(7)=6000
          emisbot=1       
          emistop=1    
         
          do nk=1,firelevel 
           if (nk.gt.1) emisbot(nk)=emistop(nk-1)+1
           do k=kts,kte
            if (sum(dz8w(i,kts:k,j)).gt.fireheight(nk+1)) then 
             if (k.eq.1) then 
             emistop(nk)=k
             else 
             emistop(nk)=k-1   
             endif
             if (emistop(nk).lt.emisbot(nk)) emistop(nk)=emisbot(nk)
            goto  9999 
            endif
           enddo !k
9999      continue
          enddo !nk
         
          !fire emission fraction is region dependent
          if (abs(xlat(i,j)).lt.30) then !tropical 
          bbemis_frac(1)=0.2
          bbemis_frac(2)=0.4
          bbemis_frac(3)=0.4
          bbemis_frac(4)=0.0
          bbemis_frac(5)=0.0
          bbemis_frac(6)=0.0
          endif 
          if (abs(xlat(i,j)).ge.30.and.abs(xlat(i,j)).le.60) then  !temperate
          bbemis_frac(1)=0.2
          bbemis_frac(2)=0.2
          bbemis_frac(3)=0.2
          bbemis_frac(4)=0.4
          bbemis_frac(5)=0.0
          bbemis_frac(6)=0.0
          endif
          if (xlat(i,j).ge.60.and.xlong(i,j).ge.0) then !East Boreal
          bbemis_frac(1)=0.1
          bbemis_frac(2)=0.1
          bbemis_frac(3)=0.2
          bbemis_frac(4)=0.2
          bbemis_frac(5)=0.4
          bbemis_frac(6)=0.0
          endif
          if (xlat(i,j).ge.60.and.xlong(i,j).lt.0) then !West Boreal
          bbemis_frac(1)=0.1
          bbemis_frac(2)=0.1
          bbemis_frac(3)=0.1
          bbemis_frac(4)=0.1
          bbemis_frac(5)=0.2
          bbemis_frac(6)=0.4
          endif

          !for gas phase
          total_air=0.0
          total_emis=0.0
          do nk=1,firelevel
          do k=emisbot(nk),emistop(nk)
            !rho_phy kg/m3, convert to kg/m2
            total_air(nk)=total_air(nk)+rho_phy(i,k,j)*dz8w(i,k,j) !air per unit area
          enddo
          enddo
          !kg/m2 to mole/m2
          total_air=total_air*1000.0/29

           if (i==28.and.j==25) then 
           do nk=1,firelevel
           write(0,*) 'czhao check fire distribution ',nk,emisbot(nk),emistop(nk) 
           print*, 'czhao check fire distribution ',nk,emisbot(nk),emistop(nk) 
           enddo

           write(0,*) 'czhao check fire emissions ',ebu(i,kts:kte,j,p_ebu_co)
           print*, 'czhao check fire emissions ',ebu(i,kts:kte,j,p_ebu_co)

           endif


          do nk=1,firelevel
          do k=emisbot(nk),emistop(nk)
          !conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/60./dz8w(i,k,j)
          !czhao change the formula to be easy-understanding
!czhao         conv_rho=r_q*dtstep*29.0/(rho_phy(i,k,j)*1000.0*dz8w(i,k,j)*3600.0)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_so2))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_sulf))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_csl))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_csl)=chem(i,k,j,p_csl)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_iso))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_iso)=chem(i,k,j,p_iso)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_no))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_no)=chem(i,k,j,p_no)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_ald))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_ald)=chem(i,k,j,p_ald)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_hcho))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_hcho)=chem(i,k,j,p_hcho)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_ora2))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_ora2)=chem(i,k,j,p_ora2)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_hc3))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_hc3)=chem(i,k,j,p_hc3)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_hc5))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_hc5)=chem(i,k,j,p_hc5)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_hc8))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_hc8)=chem(i,k,j,p_hc8)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_eth))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_eth)=chem(i,k,j,p_eth)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_co))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_co)=chem(i,k,j,p_co)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_olt))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_olt)=chem(i,k,j,p_olt)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_oli))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_oli)=chem(i,k,j,p_oli)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_tol))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_tol)=chem(i,k,j,p_tol)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_xyl))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_xyl)=chem(i,k,j,p_xyl)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_ket))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_ket)=chem(i,k,j,p_ket)+1.e6*total_emis(nk)/total_air(nk)

          enddo !k
          enddo !nk

        enddo !j
        enddo !i
    CASE (CBMZ, CBMZ_BB, CBMZ_BB_KPP, CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN, &
          CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
          CBMZSORG, CBMZSORG_AQ, CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, &
          CBMZ_SORGAM_KPP,CBMZ_MOSAIC_KPP, &
          CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
          CBMZ_MOSAIC_DMS_4BIN_AQ_SOURCE, CBMZ_MOSAIC_DMS_8BIN_AQ_SOURCE)
          !czhao add biomass burning emission for CBMZ chemistry 
          do j=jts,jte
          do i=its,ite

          if (kts.ne.1) then
           print*, 'kts not start from the first layer'
           stop
          endif

          fireheight(1)=0
          fireheight(2)=100
          fireheight(3)=500
          fireheight(4)=1000
          fireheight(5)=2000
          fireheight(6)=3000
          fireheight(7)=6000
          emisbot=1
          emistop=1

          do nk=1,firelevel
           if (nk.gt.1) emisbot(nk)=emistop(nk-1)+1
           do k=kts,kte
            if (sum(dz8w(i,kts:k,j)).gt.fireheight(nk+1)) then
             if (k.eq.1) then
             emistop(nk)=k
             else
             emistop(nk)=k-1
             endif
             if (emistop(nk).lt.emisbot(nk)) emistop(nk)=emisbot(nk)
            goto  999
            endif
           enddo !k
999      continue
          enddo !nk

          !fire emission fraction is region dependent
          if (abs(xlat(i,j)).lt.30) then  !tropical 
          bbemis_frac(1)=0.2
          bbemis_frac(2)=0.4
          bbemis_frac(3)=0.4
          bbemis_frac(4)=0.0
          bbemis_frac(5)=0.0
          bbemis_frac(6)=0.0
          endif
          if (abs(xlat(i,j)).ge.30.and.abs(xlat(i,j)).le.60) then !temperate
          bbemis_frac(1)=0.2
          bbemis_frac(2)=0.2
          bbemis_frac(3)=0.2
          bbemis_frac(4)=0.4
          bbemis_frac(5)=0.0
          bbemis_frac(6)=0.0
          endif
          if (xlat(i,j).ge.60.and.xlong(i,j).ge.0) then !East Boreal
          bbemis_frac(1)=0.1
          bbemis_frac(2)=0.1
          bbemis_frac(3)=0.2
          bbemis_frac(4)=0.2
          bbemis_frac(5)=0.4
          bbemis_frac(6)=0.0
          endif
          if (xlat(i,j).ge.60.and.xlong(i,j).lt.0) then !West Boreal
          bbemis_frac(1)=0.1
          bbemis_frac(2)=0.1
          bbemis_frac(3)=0.1
          bbemis_frac(4)=0.1
          bbemis_frac(5)=0.2
          bbemis_frac(6)=0.4
          endif

          !for gas phase
          total_air=0.0
          total_emis=0.0
          do nk=1,firelevel
          do k=emisbot(nk),emistop(nk)
            !rho_phy kg/m3, convert to kg/m2
            total_air(nk)=total_air(nk)+rho_phy(i,k,j)*dz8w(i,k,j) !air per unit area
          enddo
          enddo
          !kg/m2 to mole/m2
          total_air=total_air*1000.0/29

           if (i==28.and.j==25) then
           do nk=1,firelevel
           write(0,*) 'czhao check fire distribution ',nk,emisbot(nk),emistop(nk)
           print*, 'czhao check fire distribution ',nk,emisbot(nk),emistop(nk)
           enddo

           write(0,*) 'czhao check fire emissions ',ebu(i,kts:kte,j,p_ebu_co)
           print*, 'czhao check fire emissions ',ebu(i,kts:kte,j,p_ebu_co)

           endif

          do nk=1,firelevel
          do k=emisbot(nk),emistop(nk)
        !mole/km2/hr to ppm
!czhao        conv_rho=r_q*dtstep*29.0/(rho_phy(i,k,j)*1000.0*dz8w(i,k,j)*3600.0)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_so2))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_sulf))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_csl))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_csl)=chem(i,k,j,p_csl)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_iso))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_iso)=chem(i,k,j,p_iso)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_no))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_no)=chem(i,k,j,p_no)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_ald))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_ald)=chem(i,k,j,p_ald)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_hcho))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_hcho)=chem(i,k,j,p_hcho)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_ora2))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_ora2)=chem(i,k,j,p_ora2)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_eth))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_eth)=chem(i,k,j,p_eth)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_co))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_co)=chem(i,k,j,p_co)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_olt))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_olt)=chem(i,k,j,p_olt)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_oli))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_oli)=chem(i,k,j,p_oli)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*sum(ebu(i,kts:kte,j,p_ebu_tol))*(dtstep/3600.)/1.e6
           chem(i,k,j,p_tol)=chem(i,k,j,p_tol)+1.e6*total_emis(nk)/total_air(nk)

           total_emis(nk)=bbemis_frac(nk)*(0.4*sum(ebu(i,kts:kte,j,p_ebu_ald))   &
                                          +2.9*sum(ebu(i,kts:kte,j,p_ebu_hc3))   &
                                          +4.8*sum(ebu(i,kts:kte,j,p_ebu_hc5))   &
                                          +7.9*sum(ebu(i,kts:kte,j,p_ebu_hc8))   &
                                          +0.9*sum(ebu(i,kts:kte,j,p_ebu_ket))   &
                                          +2.8*sum(ebu(i,kts:kte,j,p_ebu_oli))   &
                                          +1.8*sum(ebu(i,kts:kte,j,p_ebu_olt))   &
                                          +1.0*sum(ebu(i,kts:kte,j,p_ebu_ora2)) )   &
                                         *(dtstep/3600.)/1.e6
           chem(i,k,j,p_par)=chem(i,k,j,p_par)+1.e6*total_emis(nk)/total_air(nk)


          enddo !k
          enddo !nk

        enddo !j
        enddo !i

      CASE (GOCART_SIMPLE)
          do j=jts,jte
          do i=its,ite
          do k=kts,kte
        conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
        chem(i,k,j,p_so2)  =  chem(i,k,j,p_so2)                        &
                         +ebu(i,k,j,p_ebu_so2)*conv_rho
        chem(i,k,j,p_sulf)  =  chem(i,k,j,p_sulf)                        &
                         +ebu(i,k,j,p_ebu_sulf)*conv_rho
        chem(i,k,j,p_dms)  =  chem(i,k,j,p_dms)                        &
                         +ebu(i,k,j,p_ebu_dms)*conv_rho
        chem(i,k,j,p_oc1)  =  chem(i,k,j,p_oc1)                        &
                         +r_q*ebu(i,k,j,p_ebu_oc)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_bc1)  =  chem(i,k,j,p_bc1)                        &
                         +r_q*ebu(i,k,j,p_ebu_bc)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_p25)  =  chem(i,k,j,p_p25)                        &
                         +r_q*ebu(i,k,j,p_ebu_pm25)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_p10)  =  chem(i,k,j,p_p10)                        &
                         +r_q*ebu(i,k,j,p_ebu_pm10)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        enddo
        enddo
        enddo
      CASE (GOCARTRACM_KPP,GOCARTRADM2_KPP,GOCARTRADM2)
          do j=jts,jte
          do i=its,ite
           do k=kts,kte
        conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
        chem(i,k,j,p_so2)  =  chem(i,k,j,p_so2)                        &
                         +ebu(i,k,j,p_ebu_so2)*conv_rho
        chem(i,k,j,p_sulf)  =  chem(i,k,j,p_sulf)                        &
                         +ebu(i,k,j,p_ebu_sulf)*conv_rho
        chem(i,k,j,p_dms)  =  chem(i,k,j,p_dms)                        &
                         +ebu(i,k,j,p_ebu_dms)*conv_rho
        chem(i,k,j,p_csl)  =  chem(i,k,j,p_csl)                        &
                         +ebu(i,k,j,p_ebu_csl)*conv_rho
        chem(i,k,j,p_iso)  = chem(i,k,j,p_iso)                         &
                         +ebu(i,k,j,p_ebu_iso)*conv_rho
        chem(i,k,j,p_no)   = chem(i,k,j,p_no)                          &
                         +ebu(i,k,j,p_ebu_no)*conv_rho
        chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)                         &
                         +ebu(i,k,j,p_ebu_no2)*conv_rho
        chem(i,k,j,p_ald)  = chem(i,k,j,p_ald)                         &
                         +ebu(i,k,j,p_ebu_ald)*conv_rho
        chem(i,k,j,p_hcho) = chem(i,k,j,p_hcho)                        &
                         +ebu(i,k,j,p_ebu_hcho)*conv_rho
        chem(i,k,j,p_ora2)  = chem(i,k,j,p_ora2)                       &
                         +ebu(i,k,j,p_ebu_ora2)*conv_rho
        chem(i,k,j,p_hc3)  = chem(i,k,j,p_hc3)                         &
                         +ebu(i,k,j,p_ebu_hc3)*conv_rho
        chem(i,k,j,p_hc5)  = chem(i,k,j,p_hc5)                         &
                         +ebu(i,k,j,p_ebu_hc5)*conv_rho
        chem(i,k,j,p_hc8)  = chem(i,k,j,p_hc8)                         &
                         +ebu(i,k,j,p_ebu_hc8)*conv_rho
        chem(i,k,j,p_eth)  = chem(i,k,j,p_eth)                         &
                         +ebu(i,k,j,p_ebu_eth)*conv_rho
        chem(i,k,j,p_co)  = chem(i,k,j,p_co)                           &
                         +ebu(i,k,j,p_ebu_co)*conv_rho
        chem(i,k,j,p_olt)  = chem(i,k,j,p_olt)                         &
                         +ebu(i,k,j,p_ebu_olt)*conv_rho
        chem(i,k,j,p_oli)  = chem(i,k,j,p_oli)                         &
                         +ebu(i,k,j,p_ebu_oli)*conv_rho
        chem(i,k,j,p_tol)  = chem(i,k,j,p_tol)                         &
                         +ebu(i,k,j,p_ebu_tol)*conv_rho
        chem(i,k,j,p_xyl)  = chem(i,k,j,p_xyl)                         &
                         +ebu(i,k,j,p_ebu_xyl)*conv_rho
        chem(i,k,j,p_ket)  =  chem(i,k,j,p_ket)                        &
                         +ebu(i,k,j,p_ebu_ket)*conv_rho
        chem(i,k,j,p_oc1)  =  chem(i,k,j,p_oc1)                        &
                         +r_q*ebu(i,k,j,p_ebu_oc)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_bc1)  =  chem(i,k,j,p_bc1)                        &
                         +r_q*ebu(i,k,j,p_ebu_bc)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_p25)  =  chem(i,k,j,p_p25)                        &
                         +r_q*ebu(i,k,j,p_ebu_pm25)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_p10)  =  chem(i,k,j,p_p10)                        &
                         +r_q*ebu(i,k,j,p_ebu_pm10)/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        enddo
        enddo
        enddo

      CASE (SAPRC99_KPP,SAPRC99_MOSAIC_4BIN_VBS2_KPP,SAPRC99_MOSAIC_4BIN_VBS9_KPP, &
           SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_KPP)!czhao add biomass burning FINN speciation for SAPRC99, gas phase only

        do k=kts,5
        if (sum(ebu(its:ite,k,jts:jte,p_ebu_co)).gt.0) then 
        write(0,*) 'czhao check plume fire in module_add_emiss_burn ',k,sum(ebu(its:ite,k,jts:jte,p_ebu_co))
        print*, 'czhao check plume fire in module_add_emiss_burn ',k,sum(ebu(its:ite,k,jts:jte,p_ebu_co))
        endif 
        enddo

        do j=jts,jte
        do i=its,ite
        do k=kts,kte
!       conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
        !mole/km2/hr to ppmV
        conv_rho=8.047e-6/rho_phy(i,k,j)*dtstep/dz8w(i,k,j)
        chem(i,k,j,p_no)  =  chem(i,k,j,p_no)+ebu(i,k,j,p_ebu_no)*conv_rho
        chem(i,k,j,p_no2)  =  chem(i,k,j,p_no2)+ebu(i,k,j,p_ebu_no2)*conv_rho
        chem(i,k,j,p_co)  =  chem(i,k,j,p_co)+ebu(i,k,j,p_ebu_co)*conv_rho
        chem(i,k,j,p_co2)  =  chem(i,k,j,p_co2)+ebu(i,k,j,p_ebu_co2)*conv_rho
        chem(i,k,j,p_ch4)  =  chem(i,k,j,p_ch4)+ebu(i,k,j,p_ebu_ch4)*conv_rho
        chem(i,k,j,p_c2h6)  =  chem(i,k,j,p_c2h6)+ebu(i,k,j,p_ebu_c2h6)*conv_rho
        chem(i,k,j,p_acet)  =  chem(i,k,j,p_acet)+ebu(i,k,j,p_ebu_acet)*conv_rho
        chem(i,k,j,p_alk3)  =  chem(i,k,j,p_alk3)+ebu(i,k,j,p_ebu_alk3)*conv_rho
        chem(i,k,j,p_alk4)  =  chem(i,k,j,p_alk4)+ebu(i,k,j,p_ebu_alk4)*conv_rho
        chem(i,k,j,p_alk5)  =  chem(i,k,j,p_alk5)+ebu(i,k,j,p_ebu_alk5)*conv_rho
        chem(i,k,j,p_aro1)  =  chem(i,k,j,p_aro1)+ebu(i,k,j,p_ebu_aro1)*conv_rho
        chem(i,k,j,p_aro2)  =  chem(i,k,j,p_aro2)+ebu(i,k,j,p_ebu_aro2)*conv_rho
        chem(i,k,j,p_bald)  =  chem(i,k,j,p_bald)+ebu(i,k,j,p_ebu_bald)*conv_rho
        chem(i,k,j,p_ccho)  =  chem(i,k,j,p_ccho)+ebu(i,k,j,p_ebu_ccho)*conv_rho
        chem(i,k,j,p_cco_oh)  =  chem(i,k,j,p_cco_oh)+ebu(i,k,j,p_ebu_cco_oh)*conv_rho
        chem(i,k,j,p_ethene)  =  chem(i,k,j,p_ethene)+ebu(i,k,j,p_ebu_ethene)*conv_rho
        chem(i,k,j,p_hcho)  =  chem(i,k,j,p_hcho)+ebu(i,k,j,p_ebu_hcho)*conv_rho
        chem(i,k,j,p_hcooh)  =  chem(i,k,j,p_hcooh)+ebu(i,k,j,p_ebu_hcooh)*conv_rho
        chem(i,k,j,p_hono)  =  chem(i,k,j,p_hono)+ebu(i,k,j,p_ebu_hono)*conv_rho
        chem(i,k,j,p_isoprene)  =  chem(i,k,j,p_isoprene)+ebu(i,k,j,p_ebu_isoprene)*conv_rho
        chem(i,k,j,p_mek)  =  chem(i,k,j,p_mek)+ebu(i,k,j,p_ebu_mek)*conv_rho
        chem(i,k,j,p_meoh)  =  chem(i,k,j,p_meoh)+ebu(i,k,j,p_ebu_meoh)*conv_rho
        chem(i,k,j,p_methacro)  =  chem(i,k,j,p_methacro)+ebu(i,k,j,p_ebu_methacro)*conv_rho
        chem(i,k,j,p_mgly)  =  chem(i,k,j,p_mgly)+ebu(i,k,j,p_ebu_mgly)*conv_rho
        chem(i,k,j,p_mvk)  =  chem(i,k,j,p_mvk)+ebu(i,k,j,p_ebu_mvk)*conv_rho
        chem(i,k,j,p_ole1)  =  chem(i,k,j,p_ole1)+ebu(i,k,j,p_ebu_ole1)*conv_rho
        chem(i,k,j,p_ole2)  =  chem(i,k,j,p_ole2)+ebu(i,k,j,p_ebu_ole2)*conv_rho
        chem(i,k,j,p_phen)  =  chem(i,k,j,p_phen)+ebu(i,k,j,p_ebu_phen)*conv_rho
        chem(i,k,j,p_prod2)  =  chem(i,k,j,p_prod2)+ebu(i,k,j,p_ebu_prod2)*conv_rho
        chem(i,k,j,p_rcho)  =  chem(i,k,j,p_rcho)+ebu(i,k,j,p_ebu_rcho)*conv_rho
        chem(i,k,j,p_rno3)  =  chem(i,k,j,p_rno3)+ebu(i,k,j,p_ebu_rno3)*conv_rho
        chem(i,k,j,p_so2)  =  chem(i,k,j,p_so2)+ebu(i,k,j,p_ebu_so2)*conv_rho
        chem(i,k,j,p_nh3)  =  chem(i,k,j,p_nh3)+ebu(i,k,j,p_ebu_nh3)*conv_rho

       enddo
       enddo
       enddo

!     CASE (SAPRC99_KPP,SAPRC99_MOSAIC_4BIN_VBS2_KPP,SAPRC99_MOSAIC_4BIN_VBS9_KPP, &
!          SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_KPP)!BSINGH(12/03/2013): Added SAPRC 8 bin and non-aq on (04/07/2014)
!         do j=jts,jte
!         do i=its,ite
!          do k=kts,kte
!       conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/60./dz8w(i,k,j)
!       chem(i,k,j,p_csl)  =  chem(i,k,j,p_csl)                        &
!                        +ebu(i,k,j,p_ebu_csl)*conv_rho
!       chem(i,k,j,p_iso)  = chem(i,k,j,p_iso)                         &
!                        +ebu(i,k,j,p_ebu_iso)*conv_rho
!       chem(i,k,j,p_no)   = chem(i,k,j,p_no)                          &
!                        +ebu(i,k,j,p_ebu_no)*conv_rho
!       chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)                         &
!                        +ebu(i,k,j,p_ebu_no2)*conv_rho
!       chem(i,k,j,p_ald)  = chem(i,k,j,p_ald)                         &
!                        +ebu(i,k,j,p_ebu_ald)*conv_rho
!       chem(i,k,j,p_hcho) = chem(i,k,j,p_hcho)                        &
!                        +ebu(i,k,j,p_ebu_hcho)*conv_rho
!       chem(i,k,j,p_ora2)  = chem(i,k,j,p_ora2)                       &
!                        +ebu(i,k,j,p_ebu_ora2)*conv_rho
!       chem(i,k,j,p_hc3)  = chem(i,k,j,p_hc3)                         &
!                        +ebu(i,k,j,p_ebu_hc3)*conv_rho
!       chem(i,k,j,p_hc5)  = chem(i,k,j,p_hc5)                         &
!                        +ebu(i,k,j,p_ebu_hc5)*conv_rho
!       chem(i,k,j,p_hc8)  = chem(i,k,j,p_hc8)                         &
!                        +ebu(i,k,j,p_ebu_hc8)*conv_rho
!       chem(i,k,j,p_eth)  = chem(i,k,j,p_eth)                         &
!                        +ebu(i,k,j,p_ebu_eth)*conv_rho
!       chem(i,k,j,p_co)  = chem(i,k,j,p_co)                           &
!                        +ebu(i,k,j,p_ebu_co)*conv_rho
!       chem(i,k,j,p_olt)  = chem(i,k,j,p_olt)                         &
!                        +ebu(i,k,j,p_ebu_olt)*conv_rho
!       chem(i,k,j,p_oli)  = chem(i,k,j,p_oli)                         &
!                        +ebu(i,k,j,p_ebu_oli)*conv_rho
!       chem(i,k,j,p_tol)  = chem(i,k,j,p_tol)                         &
!                        +ebu(i,k,j,p_ebu_tol)*conv_rho
!       chem(i,k,j,p_xyl)  = chem(i,k,j,p_xyl)                         &
!                        +ebu(i,k,j,p_ebu_xyl)*conv_rho
!       chem(i,k,j,p_ket)  =  chem(i,k,j,p_ket)                        &
!                        +ebu(i,k,j,p_ebu_ket)*conv_rho
!       enddo
!       enddo
!       enddo
      CASE (MOZART_KPP,MOZCART_KPP,MOZART_MOSAIC_4BIN_VBS0_KPP )
        if( biomass_burn_opt == BIOMASSB_MOZC .or. biomass_burn_opt == BIOMASSB_MOZ ) then
          do j=jts,jte
            do k=kts,kte
              do i=its,ite
                conv_rho = (r_q*4.828e-4*dtstep)/(rho_phy(i,k,j)*60.*dz8w(i,k,j))
                chem(i,k,j,p_co)  = chem(i,k,j,p_co) + ebu(i,k,j,p_ebu_co)*conv_rho
                chem(i,k,j,p_no)  = chem(i,k,j,p_no) + ebu(i,k,j,p_ebu_no)*conv_rho
                chem(i,k,j,p_no2) = chem(i,k,j,p_no2) + ebu(i,k,j,p_ebu_no2)*conv_rho
                chem(i,k,j,p_bigalk) =  chem(i,k,j,p_bigalk) + ebu(i,k,j,p_ebu_bigalk)*conv_rho
                chem(i,k,j,p_bigene) =  chem(i,k,j,p_bigene) + ebu(i,k,j,p_ebu_bigene)*conv_rho
                chem(i,k,j,p_c2h4)   = chem(i,k,j,p_c2h4) + ebu(i,k,j,p_ebu_c2h4)*conv_rho
                chem(i,k,j,p_c2h5oh) = chem(i,k,j,p_c2h5oh) + ebu(i,k,j,p_ebu_c2h5oh)*conv_rho
                chem(i,k,j,p_c2h6) = chem(i,k,j,p_c2h6) + ebu(i,k,j,p_ebu_c2h6)*conv_rho
                chem(i,k,j,p_c3h6) = chem(i,k,j,p_c3h6) + ebu(i,k,j,p_ebu_c3h6)*conv_rho
                chem(i,k,j,p_c3h8) = chem(i,k,j,p_c3h8) + ebu(i,k,j,p_ebu_c3h8)*conv_rho
                chem(i,k,j,p_hcho) = chem(i,k,j,p_hcho) +ebu(i,k,j,p_ebu_ch2o)*conv_rho
                chem(i,k,j,p_ald)  = chem(i,k,j,p_ald) +ebu(i,k,j,p_ebu_ch3cho)*conv_rho
                chem(i,k,j,p_acetol) = chem(i,k,j,p_acetol) +ebu(i,k,j,p_ebu_acetol)*conv_rho
                chem(i,k,j,p_isopr)  = chem(i,k,j,p_isopr) +ebu(i,k,j,p_ebu_isop)*conv_rho
                chem(i,k,j,p_macr) = chem(i,k,j,p_macr) +ebu(i,k,j,p_ebu_macr)*conv_rho
                chem(i,k,j,p_mvk)  = chem(i,k,j,p_mvk) +ebu(i,k,j,p_ebu_mvk)*conv_rho
                chem(i,k,j,p_acet)  = chem(i,k,j,p_acet) + ebu(i,k,j,p_ebu_ch3coch3)*conv_rho
                chem(i,k,j,p_ch3oh) = chem(i,k,j,p_ch3oh) + ebu(i,k,j,p_ebu_ch3oh)*conv_rho
                chem(i,k,j,p_mek)  = chem(i,k,j,p_mek) + ebu(i,k,j,p_ebu_mek)*conv_rho
                chem(i,k,j,p_so2)  =  chem(i,k,j,p_so2) +ebu(i,k,j,p_ebu_so2)*conv_rho
                chem(i,k,j,p_tol)  = chem(i,k,j,p_tol) +ebu(i,k,j,p_ebu_toluene)*conv_rho
                chem(i,k,j,p_nh3)  = chem(i,k,j,p_nh3) + ebu(i,k,j,p_ebu_nh3)*conv_rho
                chem(i,k,j,p_open) = chem(i,k,j,p_open) + ebu(i,k,j,p_ebu_open)*conv_rho
                chem(i,k,j,p_c10h16) = chem(i,k,j,p_c10h16) + ebu(i,k,j,p_ebu_c10h16)*conv_rho
                chem(i,k,j,p_cres)   = chem(i,k,j,p_cres) + ebu(i,k,j,p_ebu_cres)*conv_rho
                chem(i,k,j,p_glyald) = chem(i,k,j,p_glyald) + ebu(i,k,j,p_ebu_glyald)*conv_rho
                chem(i,k,j,p_gly)    = chem(i,k,j,p_gly) + ebu(i,k,j,p_ebu_gly)*conv_rho
              enddo
            enddo
          enddo
          if( biomass_burn_opt == BIOMASSB_MOZC ) then
            do j=jts,jte
              do k=kts,kte
                do i=its,ite
                  conv_rho = (r_q*dtstep)/(rho_phy(i,k,j)*dz8w(i,k,j))
                  chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + conv_rho*ebu(i,k,j,p_ebu_oc)
                  chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + conv_rho*ebu(i,k,j,p_ebu_bc)
                  chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + conv_rho*ebu(i,k,j,p_ebu_pm10)
                  chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + conv_rho*ebu(i,k,j,p_ebu_pm25)
                enddo
              enddo
            enddo
          endif
        endif
      CASE (CHEM_TRACE2)
          do j=jts,jte
          do i=its,ite
          do k=kts+1,kte-1
        conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/60./dz8w(i,k,j)
        chem(i,k,j,p_tracer_1)  =  chem(i,k,j,p_tracer_1)                        &
                         +ebu(i,k,j,p_ebu_co)*conv_rho
        enddo
        k=kts
        conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
        chem(i,k,j,p_tracer_1)  =  chem(i,k,j,p_tracer_1)                        &
                         +ebu(i,k,j,p_ebu_co)*conv_rho
        enddo
        enddo
      CASE (GHG_TRACER)
        if( biomass_burn_opt == BIOMASSB_GHG ) then
          do j=jts,jte
          do k=kts,kte
          do i=its,ite
             conv_rho=r_q*4.828e-4/rho_phy(i,k,j)*dtstep/60./dz8w(i,k,j)
             chem(i,k,j,p_co_bbu)  = chem(i,k,j,p_co_bbu)   +ebu(i,k,j,p_ebu_co)*conv_rho
             chem(i,k,j,p_co2_bbu)   = chem(i,k,j,p_co2_bbu)  +ebu(i,k,j,p_ebu_co2)*conv_rho
             chem(i,k,j,p_ch4_bbu) =  chem(i,k,j,p_ch4_bbu) + ebu(i,k,j,p_ebu_ch4)*conv_rho
          enddo
          enddo
          enddo
        endif
    CASE DEFAULT
       call wrf_debug(15,'nothing done with burn emissions for chem array')
    END SELECT emiss_select


    END subroutine add_emis_burn

END Module module_add_emiss_burn
