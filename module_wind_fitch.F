!WRF:MODEL_LAYER:PHYSICS

MODULE module_wind_fitch
!
!Represents kinetic energy extracted by wind turbines and turbulence
! (TKE) they produce at model levels within the rotor area.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NOTICE
!! The following paper should be cited whenever presenting results using this scheme
!! (using either the original version or any modified versions of the scheme):
!! Fitch, A. C. et al. 2012: Local and Mesoscale Impacts of Wind Farms as Parameterized in a
!! Mesoscale NWP Model. Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-11-00352.1
!!
!! Anna C. Fitch, National Center for Atmospheric Research (formerly University of Bergen)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! History of changes:
!
! WRFV3.5.1: 
! WRFV3.6:   Modified by Pedro A. Jimenez to include:
!             - Initialize the wind turbines in this module.
!             - Introduce z_at_walls to avoid instabilities due to neglecting
!                the perturbation of the geopotential height.
!             - User friendly interface to introduce the technical characteritics of
!                the wind turbines.
!             - Only uses one set of turbine coefficients using the wind speed at hub height
!             - Two standing coefficients.
!             - Calculates the power produced by the wind turbines.
!
! References:
!
! Fitch, A. C. et al. 2012: Local and Mesoscale Impacts of Wind Farms as Parameterized in a
!    Mesoscale NWP Model. Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-11-00352.1
! Fitch, A. C. et al. 2013: Mesoscale Influences of Wind Farms Throughout a Diurnal Cycle.
!    Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-12-00185.1
! Fitch, A. C. et al. 2013: Parameterization of Wind Farms in Climate Models.
!    Journal of Climate, doi:http://dx.doi.org/10.1175/JCLI-D-12-00376.1
! Jimenez, P.A., J. Navarro, A.M. Palomares and J. Dudhia:  Mesoscale modeling of offshore wind turbines
!    wakes at the wind farm resolving scale: a composite-based analysis with the WRF model over Horns Rev. 
!    Wind Energy, (In Press.).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  USE module_driver_constants, ONLY : max_domains
  USE module_model_constants, ONLY :  piconst
!
  USE module_llxy
  USE module_dm, ONLY : wrf_dm_min_real
  USE module_configure, ONLY : grid_config_rec_type

  IMPLICIT NONE

  INTEGER, PARAMETER :: MAXVALS  = 100   
  INTEGER, PARAMETER :: MAXVALS2 = 100     
!
  INTEGER           :: nt
  INTEGER, DIMENSION(:), ALLOCATABLE :: NKIND
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ival,jval
  REAL, DIMENSION(:), ALLOCATABLE :: hubheight,diameter,stc,stc2,cutin,cutout,npower
!
  REAL :: turbws(maxvals,maxvals2),turbtc(maxvals,maxvals2),turbpw(maxvals,maxvals2)
  REAL :: correction_factor
!
CONTAINS

 SUBROUTINE tg_interpolate(id,z_at_w,u,v,rho                 &
       &,dx,dz,dt                            &
       &,du,dv                               & !,upwp_tg,vpwp_tg                                     &
       &,windfarm_opt                      &
       &,ids,ide,jds,jde,kds,kde                 &
       &,ims,ime,jms,jme,kms,kme                 &
       &,its,ite,jts,jte,kts,kte                 &
       &,speed_84,u_84,v_84,upwp_TG,vpwp_TG,U_140,V_140,UW_TG140,UW_TG28,VW_TG140,VW_TG28)!,u_140,v_140,speed_140,uw_tg140,uw_tg28,vw_tg140,vw_tg28)
  INTEGER, INTENT(IN) :: id,windfarm_opt
  INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde
  REAL, INTENT(IN) :: dx,dt
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) ::dz,z_at_w,rho,u,v,upwp_TG,vpwp_TG
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: du,dv
  REAL, DIMENSION(ims:ime,jms:jme),INTENT(INOUT)::speed_84,U_84,V_84,U_140,V_140,UW_TG140,UW_TG28,VW_TG140,VW_TG28
  INTEGER i_TG,j_TG,k_TG
! Local
!
  REAL blade_l_point1,blade_u_point1,zheightl1,zheightl128,zheightl1140,zheightu1,z1tg,z2tg,z1tg28,z2tg28,z1tg140,z2tg140
  REAL     speed84,u84,v84,u140,v140,uwtg140,uwtg28,vwtg140,vwtg28  !TG2 Added airdensity
  INTEGER  itf,jtf,ktf
  INTEGER  i1,j1,k1,n1
  INTEGER  k_turbine_bot1, k_turbine_top1

  LOGICAL :: kfound1
!
!
  REAL :: speed1,speed2
  INTEGER :: kbot1,ktop1,kt1,kbot1_28,ktop1_28,ktop1_140,kbot1_140

  itf=MIN0(ite,ide-1)
  jtf=MIN0(jte,jde-1)
  ktf=MIN0(kte,kde-1)


  k_TG=1
  speed_84=0.
 ! speed_140=0.
  u_84=0.
  v_84=0.
  u_140=0.
  v_140=0.
  uw_tg28=0.
  uw_tg140=0.
  vw_tg28=0.
  vw_tg140=0.

!  DO i_TG=1,999
!         DO j_TG=1,999
!                ival1(k_TG)=i_TG
!                jval1(k_TG)=j_TG
!                k_TG=k_TG+1


!   ENDDO
!  ENDDO

DO j1=JTS,JTF
        DO i1=ITS,ITF
          IF ( windfarm_opt .eq. 1 ) THEN

        k_turbine_bot1=0      !bottom level
        k_turbine_top1=-1     !top
         if (i1.ne.-9999.and.j1.ne.-9999) then
                 IF (( its .LE. i1 .AND. i1 .LE. itf ) .AND. &
                     ( jts .LE. j1 .AND. j1 .LE. jtf )  ) THEN
!
          blade_l_point1=28!=hubheight(kt1)-diameter(kt1)/2. ! height of lower blade tip above ground (m)
          blade_u_point1=140!hubheight(kt1)+diameter(kt1)/2. ! height of upper blade tip above ground (m)
!*!
          kfound1 = .false.
          zheightl1=0.0
!*          ! find vertical levels cut by turbine bladeS
!               print*,"TG kts and ktf",kts,ktf
                  DO k1=kts,ktf
                    IF(.NOT. kfound1) THEN
                      zheightu1 = zheightl1 + dz(i1,k1,j1) ! increment height
                        IF(blade_l_point1 .GE. zheightl1 .AND. blade_l_point1 .LE.zheightu1) THEN
                        k_turbine_bot1=k1 ! lower blade tip cuts this level
                       ENDIF
!*
                      IF(blade_u_point1 .GE. zheightl1 .AND. blade_u_point1 .LE.zheightu1) THEN
                         k_turbine_top1=k1 ! upper blade tip cuts this level
                        kfound1 = .TRUE.
                       ENDIF
!*
                  zheightl1 = zheightu1
                  ENDIF
               ENDDO
!         print*,"TG k_turbine_bot1 k_turbine_top1 kfound1",k_turbine_bot1,k_turbine_top1,kfound1
         CALL wrf_debug(100,'in phys/module_wind_fitch.F 2nd do loop')
!*

       IF ( kfound1 ) THEN
       kfound1 = .false.

          zheightl1=0.
          ! find vertical levels (half levels) within the hub height
          DO k1=kts,ktf
            IF(.NOT. kfound1) THEN
              z2tg = zheightl1 + 0.5*dz(i1,k1,j1)
!              print*,"TG1 k1 and z2tg",k1,z2tg
                IF(84 .GE. z2tg ) THEN
                kbot1=k1
                !print*,"TG2 k1 or kbot1",kbot1
                ELSE
                ktop1=k1
                kfound1 = .TRUE.
                !print*,"TG3 k1 or ktop1",ktop1
                ENDIF
!
              if (.NOT. kfound1) z1tg=z2tg
              zheightl1 = z2tg + 0.5*dz(i1,k1,j1)
            ENDIF
          ENDDO
!
       zheightl128=0.
       kfound1=.false.
        DO k1=kts,ktf
            IF(.NOT. kfound1) THEN
               z2tg28 = zheightl128 + 0.5*dz(i1,k1,j1)
               !print*,"TG1 k1, dz, z2tg",k1,dz(i1,k1,j1),z2tg28
                IF(28 .GE. z2tg28 ) THEN
                kbot1_28=k1
              !  print*,"TG2 k1 or kbot1",kbot1_28
                ELSE
                ktop1_28=k1
                kfound1 = .TRUE.
               ! print*,"TG3 k1 or ktop1",ktop1_28
                ENDIF
              !print*,"TG kfound1", kfound1
              if (.NOT. kfound1) z1tg28=z2tg28
              !print*,"TG3 z1tg28",z1tg28
              zheightl128 = z2tg28 + 0.5*dz(i1,k1,j1)
            ENDIF
          ENDDO


      zheightl1140=0.
      kfound1=.false.
        DO k1=kts,ktf
            IF(.NOT. kfound1) THEN
              z2tg140 = zheightl1140 + 0.5*dz(i1,k1,j1)
             ! print*,"TG1 k1 and z2tg",k1,z2tg140
                IF(140 .GE. z2tg140 ) THEN
                kbot1_140=k1
                !print*,"TG2 k1 or kbot1",kbot1_28
                ELSE
                ktop1_140=k1
                kfound1 = .TRUE.
                !print*,"TG3 k1 or ktop1",ktop1_28
                ENDIF
!
              if (.NOT. kfound1) z1tg140=z2tg140
              zheightl1140 = z2tg140 + 0.5*dz(i1,k1,j1)
            ENDIF
          ENDDO

         speed1=0.
          speed2=0.
!          print*,"TG4 ktop1 ",ktop1
          if (ktop1.eq.1) then
           speed84=sqrt(u(i1,1,j1)**2.+v(i1,1,j1)**2.)*84/z1tg
           u84=u(i1,1,j1)*84/z1tg
           v84=v(i1,1,j1)*84/z1tg

          else
           speed1=sqrt(u(i1,kbot1,j1)**2.+v(i1,kbot1,j1)**2.)
           speed2=sqrt(u(i1,ktop1,j1)**2.+v(i1,ktop1,j1)**2.)
           speed84=speed1+((speed2-speed1)/(z2tg-z1tg))*(hubheight(kt1)-z1tg)
           u84=u(i1,kbot1,j1)+((u(i1,ktop1,j1)-u(i1,kbot1,j1))/(z2tg-z1tg))*(84-z1tg)
           v84=v(i1,kbot1,j1)+((v(i1,ktop1,j1)-v(i1,kbot1,j1))/(z2tg-z1tg))*(84-z1tg)



!TG 31/12/2015 starts
          speed_84(i1,j1)=speed84  ! TG Added to O/P over all domain
          u_84(i1,j1)=u84
          v_84(i1,j1)=v84
      !    print*, "uat84 and vat84",u_84(i1,j1),v_84(i1,j1),"at i1 and j1",i1,j1
         CALL wrf_debug(100,'in phys/module_wind_fitch.F computed u84')
         endif



       ! print*,"TG1 ktop1_28 kbot1_28",ktop1_28,kbot1_28
        if (ktop1_28.eq.1) then
         !  u28=u(i1,1,j1)*28/z1tg28
          ! v28=v(i1,1,j1)*28/z1tg28
           uwtg28=upwp_TG(i1,1,j1)*28/z1tg28
           vwtg28=vpwp_TG(i1,1,j1)*28/z1tg28
       !    print*,"TG uwtg28 upwp_TG z2tg28 ktop1_28",uwtg28,upwp_TG(i1,1,j1),z2tg28,ktop1_28
          else
          ! u28=u(i1,kbot1_28,j1)+((u(i1,ktop1_28,j1)-u(i1,kbot1_28,j1))/(z2tg28-z1tg28))*(28-z1tg28)
          ! v28=v(i1,kbot1_28,j1)+((v(i1,ktop1_28,j1)-v(i1,kbot1_28,j1))/(z2tg28-z1tg28))*(28-z1tg28)
           uwtg28=upwp_TG(i1,kbot1_28,j1)+((upwp_TG(i1,ktop1_28,j1)-upwp_TG(i1,kbot1_28,j1))/(z2tg28-z1tg28))*(28-z1tg28)
           vwtg28=vpwp_TG(i1,kbot1_28,j1)+((vpwp_TG(i1,ktop1_28,j1)-vpwp_TG(i1,kbot1_28,j1))/(z2tg28-z1tg28))*(28-z1tg28)

       !   v_28(i1,j1)=v28
          uw_tg28(i1,j1)=uwtg28
          vw_tg28(i1,j1)=vwtg28
      !    print*, "uat84 and vat84",u_84(i1,j1),v_84(i1,j1),"at i1 and j1",i1,j1
     !    CALL wrf_debug(100,'in phys/module_wind_fitch.F computed u84')
         endif


  if (ktop1_140.eq.1) then
           u140=u(i1,1,j1)*140/z1tg140
           v140=v(i1,1,j1)*140/z1tg140
           uwtg140=upwp_TG(i1,1,j1)*140/z1tg140
           vwtg140=vpwp_TG(i1,1,j1)*140/z1tg140

          else
           u140=u(i1,kbot1_140,j1)+((u(i1,ktop1_140,j1)-u(i1,kbot1_140,j1))/(z2tg140-z1tg140))*(140-z1tg140)
           v140=v(i1,kbot1_140,j1)+((v(i1,ktop1_140,j1)-v(i1,kbot1_140,j1))/(z2tg140-z1tg140))*(140-z1tg140)
           uwtg140=upwp_TG(i1,kbot1_140,j1)+((upwp_TG(i1,ktop1_140,j1)-upwp_TG(i1,kbot1_140,j1))/(z2tg140-z1tg140))*(140-z1tg140)
           vwtg140=vpwp_TG(i1,kbot1_140,j1)+((vpwp_TG(i1,ktop1_140,j1)-vpwp_TG(i1,kbot1_140,j1))/(z2tg140-z1tg140))*(140-z1tg140)

!TG 31/12/2015 starts
          u_140(i1,j1)=u140
          v_140(i1,j1)=v140
          uw_tg140(i1,j1)=uwtg140
          vw_tg140(i1,j1)=vwtg140
      !    print*, "uat84 and vat84",u_84(i1,j1),v_84(i1,j1),"at i1 and
      !    j1",i1,j1
     !    CALL wrf_debug(100,'in phys/module_wind_fitch.F computed u84')
         endif




ENDIF
ENDIF
       ENDIF
      ENDIF

        ENDDO
ENDDO
 END SUBROUTINE tg_interpolate


  SUBROUTINE  dragforce(                      &
       & id                                      &
       &,z_at_w,u,v,rho                &
       &,dx,dz,dt,qke                            &
       &,du,dv                                   &
       &,windfarm_opt,power,power_2                      &
       &,ids,ide,jds,jde,kds,kde                 &
       &,ims,ime,jms,jme,kms,kme                 &
       &,its,ite,jts,jte,kts,kte                 &
       &,speedhub1,airdensity1,momdef)  
!
!
!
  INTEGER, INTENT(IN) :: id,windfarm_opt 
  INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde
  REAL, INTENT(IN) :: dx,dt
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz,rho,u,v,z_at_w
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: du,dv,qke
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) ::power,speedhub1,airdensity1,power_2
  REAL,DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(INOUT) ::momdef
   REAL,DIMENSION(ims:ime,kms:kme,jms:jme)::momdef_u,momdef_v
!
! Local
!
  REAL     blade_l_point,blade_u_point,zheightl,zheightu,z1,z2,tarea
  REAL     speed,tkecof,powcof,thrcof,wfdensity,airdensity
  INTEGER  itf,jtf,ktf
  INTEGER  i,j,k,n
  INTEGER  k_turbine_bot, k_turbine_top

  LOGICAL :: kfound
!
! ... PAJ: more variables ...
!
  REAL :: speedhub,speed1,speed2
  real :: power1,power2,area,ec
  INTEGER :: kbot,ktop,kt

  itf=MIN0(ite,ide-1)
  jtf=MIN0(jte,jde-1)
  ktf=MIN0(kte,kde-1)

    wfdensity = 1.0/(dx*dx)   !  per turbine, so numerator is 1i
 !TG Nov.2020
    speedhub1=0.
    airdensity1=0.
    momdef=0.
    momdef_u=0.
    momdef_v=0.

    power=0.

    DO kt = 1,nt  
      IF ( windfarm_opt .eq. 1 ) THEN
!
! vertical layers cut by turbine blades
!
        k_turbine_bot=0      !bottom level
        k_turbine_top=-1     !top level
        i = ival(kt,id)
        j = jval(kt,id)
!
         if (i.ne.-9999.and.j.ne.-9999) then
        IF (( its .LE. i .AND. i .LE. itf ) .AND. &
            ( jts .LE. j .AND. j .LE. jtf )  ) THEN
!
          blade_l_point=hubheight(kt)-diameter(kt)/2. ! height of lower blade tip above ground (m)
          blade_u_point=hubheight(kt)+diameter(kt)/2. ! height of upper blade tip above ground (m)
!
          kfound = .false.
          zheightl=0.0
          ! find vertical levels cut by turbine blades
          DO k=kts,ktf
            IF(.NOT. kfound) THEN
              zheightu = zheightl + dz(i,k,j) ! increment height

              IF(blade_l_point .GE. zheightl .AND. blade_l_point .LE. zheightu) THEN
                k_turbine_bot=k ! lower blade tip cuts this level
              ENDIF

              IF(blade_u_point .GE. zheightl .AND. blade_u_point .LE. zheightu) THEN
                k_turbine_top=k ! upper blade tip cuts this level
                kfound = .TRUE.
              ENDIF

              zheightl = zheightu
            ENDIF
          ENDDO
          IF ( kfound ) THEN
!
! ... PAJ: Changes introduced to compute only one set of turbine coefficients ...
!          First computes the wind speed at the hub height.
!
          kfound = .false.
          zheightl=0.
          ! find vertical levels (half levels) within the hub height
          DO k=kts,ktf
            IF(.NOT. kfound) THEN
              z2 = zheightl + 0.5*dz(i,k,j) 
!
              IF(hubheight(kt) .GE. z2 ) THEN
                kbot=k 
              ELSE
                ktop=k
                kfound = .TRUE.
              ENDIF
!
              if (.NOT. kfound) z1=z2
              zheightl = z2 + 0.5*dz(i,k,j)
            ENDIF
          ENDDO
!
          speed1=0.
          speed2=0.
          if (ktop.eq.1) then
           speedhub=sqrt(u(i,1,j)**2.+v(i,1,j)**2.)*hubheight(kt)/z1
          else
           speed1=sqrt(u(i,kbot,j)**2.+v(i,kbot,j)**2.)
           speed2=sqrt(u(i,ktop,j)**2.+v(i,ktop,j)**2.)
           speedhub=speed1+((speed2-speed1)/(z2-z1))*(hubheight(kt)-z1)
   
          !TG 31/12/2015 starts
	  speedhub1(i,j)=speedhub  ! TG Added to O/P over all domain
	!TG 31/12/2015 end
       endif
        !TG adding air density
         airdensity=rho(i,kbot,j)+((rho(i,ktop,j)-rho(i,kbot,j))/(z2-z1))*(hubheight(kt)-z1) !TG16_11_15 added airdensity linear interpolation
          airdensity1(i,j)=airdensity
! ... calculate TKE, power and thrust coeffs
!
              CALL dragcof(tkecof,powcof,thrcof,               &
                           speedhub,cutin(kt),cutout(kt),   &
                           npower(kt),diameter(kt),stc(kt),stc2(kt),nkind(kt),airdensity)
!
! ... PAJ: Computation of power generated by the wind turbine ...
!
          area=piconst/4.*diameter(kt)**2.          ! area swept by turbine blades


!!!!!!!!!!!!!!TG Commenting out origial code for power calcuation
!          power1=0.5*airdensity*speedhub**3.*area*powcof !TG removed 1.23 to airdensity
!          power(i,j)=power1+power(i,j)
!          power2=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!! TG adding Power curve 3.75 MW turbine
   IF(speedhub .LE. cutin(k) .OR. speedhub .GE. cutout(kt)) THEN  !TG
     power1=0.
     !print*,"TG cutin speed is:",cutin(k),"cutout speed is",cutout(kt) ,"speedhub :",speedhub  !TG 26/11/16
   ELSEIF (speedhub .GE. cutin(k) .AND. speedhub .LT. 12.0) THEN   !TG 25/11/16
        power1=0.5*airdensity*speedhub**3.*area*powcof
   ELSE   !TG
   power1=3075000. !TG
   ENDIF    !TG
        power(i,j)=power1+power(i,j)
        power_2(i,j)=power2+power_2(i,j)
        power2=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            DO k=k_turbine_bot,k_turbine_top ! loop over turbine blade levels
              z1=z_at_w(i,k,j)-blade_l_point-z_at_w(i,1,j)  ! distance between k level and lower blade tip
              z2=z_at_w(i,k+1,j)-blade_l_point-z_at_w(i,1,j) ! distance between k+1 level and lower blade tip
              IF(z1 .LT. 0.) z1=0.0 ! k level lower than lower blade tip
              IF(z2 .GT. diameter(kt)) z2=diameter(kt) ! k+1 level higher than turbine upper blade tip
              CALL turbine_area(z1,z2,diameter(kt),wfdensity,tarea)
!
              speed=sqrt(u(i,k,j)**2.+v(i,k,j)**2.)
              power2=power2+0.5*powcof*1.23*(speed**3.)*tarea/wfdensity
            ENDDO
!
! ... PAJ: Computes the tendencies of TKE and momentum ...
!
            DO k=k_turbine_bot,k_turbine_top ! loop over turbine blade levels
              z1=z_at_w(i,k,j)-blade_l_point-z_at_w(i,1,j)  ! distance between k lev and lower blade tip
              z2=z_at_w(i,k+1,j)-blade_l_point-z_at_w(i,1,j) !distance between k+1 lev and lower blade tip
              IF(z1 .LT. 0.) z1=0.0 ! k level lower than lower blade tip
              IF(z2 .GT. diameter(kt)) z2=diameter(kt) ! k+1 level higher than turbine upper blade tip
!
              CALL turbine_area(z1,z2,diameter(kt),wfdensity,tarea)
!
              speed=sqrt(u(i,k,j)**2.+v(i,k,j)**2.)
!`
! ... PAJ: normalization introduced to conserve energy ...
!
              if (power1.eq.0.or.power2.eq.0) then
              ec=1.
              else
              ec=power1/power2
              endif
!
              ! output TKE
              qke(i,k,j) = qke(i,k,j)+speed**3.*tarea*tkecof*dt/dz(i,k,j)*ec
              ! output u tendency
              !!! TG added momentum
                momdef_u(i,k,j)=momdef_u(i,k,j)+0.5*u(i,k,j)*thrcof*speed*tarea/dz(i,k,j)*ec
                momdef_v(i,k,j)=momdef_v(i,k,j)+.5*v(i,k,j)*thrcof*speed*tarea/dz(i,k,j)*ec
                momdef(i,k,j)=sqrt((momdef_u(i,k,j)*momdef_u(i,k,j))+(momdef_v(i,k,j)*momdef_v(i,k,j)))
              !!!!!!!!!!

	      du(i,k,j) = du(i,k,j)-.5*u(i,k,j)*thrcof*speed*tarea/dz(i,k,j)*ec
              ! output v tendency
              dv(i,k,j) = dv(i,k,j)-.5*v(i,k,j)*thrcof*speed*tarea/dz(i,k,j)*ec
               ENDDO
	      !!!! TG Removed Momdef from other layers
                 DO k=kms,kme
            		if(k.LE.k_turbine_bot.OR.k.GE.k_turbine_top) then
            		momdef(i,k,j)=0.
            		momdef_u(i,k,j)=0.
            		momdef_v(i,k,j)=0.
                	endif
		 ENDDO
          ENDIF
        ENDIF
        endif
      ENDIF
    ENDDO

  END SUBROUTINE dragforce

! This subroutine calculates area of turbine between two vertical levels
! Input variables : 
!            z1 = distance between k level and lower blade tip
!            z2 = distance between k+1 level and lower blade tip
!            wfdensity = wind farm density in m^-2
!     tdiameter = turbine diameter
! Output variable :
!         tarea = area of turbine between two levels * wfdensity
  SUBROUTINE turbine_area(z1,z2,tdiameter,wfdensity,tarea)

  REAL, INTENT(IN) ::tdiameter,wfdensity
  REAL, INTENT(INOUT) ::z1,z2
  REAL, INTENT(OUT):: tarea
  REAL r,zc1,zc2

  r=tdiameter/2.              !r = turbine radius
  z1=r-z1                   !distance of kth level from turbine center 
  z2=r-z2                   !distance of k+1 th level from turbine center
  zc1=abs(z1)
  zc2=abs(z2)
  !turbine area between z1 and z2
  IF(z1 .GT. 0. .AND. z2 .GT. 0.) THEN
     tarea=zc1*sqrt(r*r-zc1*zc1)+r*r*asin(zc1/r)- &
     (zc2*sqrt(r*r-zc2*zc2)+r*r*asin(zc2/r))
  ELSE IF(z1 .LT. 0. .AND. z2 .LT. 0.) THEN
     tarea=zc2*sqrt(r*r-zc2*zc2)+r*r*asin(zc2/r)- &
     (zc1*sqrt(r*r-zc1*zc1)+r*r*asin(zc1/r))
  ELSE
     tarea=zc2*sqrt(r*r-zc2*zc2)+r*r*asin(zc2/r)+ &
     zc1*sqrt(r*r-zc1*zc1)+r*r*asin(zc1/r)
  ENDIF
  tarea=tarea*wfdensity      !turbine area * wind farm density 

  END SUBROUTINE turbine_area


  SUBROUTINE dragcof(tkecof,powcof,thrcof,speed,cispeed,cospeed, &
                     tpower,tdiameter,stdthrcoef,stdthrcoef2,nkind,airdensity)


  REAL, INTENT(IN):: speed, cispeed, cospeed, tpower,tdiameter,stdthrcoef,stdthrcoef2,airdensity !TG  
  REAL, INTENT(OUT):: tkecof,powcof,thrcof
  REAL :: power,area,mspeed,hspeed
!
! ... PAJ ...
!
   INTEGER :: nkind,k,nu,nb
   LOGICAL :: vfound
   REAL :: fac1,fac2

  area=piconst/4.*tdiameter**2.          ! area swept by turbine blades

      vfound=.false.
      DO k=1,maxvals2
            IF(.NOT. vfound) THEN
              IF(turbws(nkind,k).GT.speed) THEN
                nu=k 
                nb=k-1
                vfound=.true.
              ENDIF
            ENDIF
      ENDDO
!
!!!!!!!!!!!!!!!!!! TG commenting out original code
!  IF (speed .LE. cispeed) THEN
!     thrcof = stdthrcoef
!  ELSE
!    IF (speed .GE. cospeed) THEN
!     thrcof = stdthrcoef2
!     ELSE
!     thrcof = turbtc(nkind,nb)+(turbtc(nkind,nu)-turbtc(nkind,nb))/(turbws(nkind,nu)-turbws(nkind,nb))*(speed-turbws(nkind,nb))
!    ENDIF
!  ENDIF
!!
!! ... power coeficient ...
!!
!  IF(speed .LE. cispeed .OR. speed .GE. cospeed) THEN
!     power=0.
!     powcof=0.
!  ELSE
!      fac1=1000./(0.5*1.23*turbws(nkind,nb)**3.*area)
!      fac2=1000./(0.5*1.23*turbws(nkind,nu)**3.*area)
!      power = turbpw(nkind,nb)+(turbpw(nkind,nu)-turbpw(nkind,nb))/(turbws(nkind,nu)-turbws(nkind,nb)) &
!                               *(speed-turbws(nkind,nb))
!     powcof = turbpw(nkind,nb)*fac1+(turbpw(nkind,nu)*fac2-turbpw(nkind,nb)*fac1)/(turbws(nkind,nu)-turbws(nkind,nb)) &
!                                     *(speed-turbws(nkind,nb))
!  ENDIF
!!
!  ! tke coefficient calculation 
!
!  tkecof=thrcof-powcof
!  tkecof = correction_factor * tkecof
!  IF(tkecof .LT. 0.) tkecof=0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TG commenting ends


!!!!******** TG adding new code
powcof=0.
thrcof=0.
tkecof=0.
!!!!!!!!!!!!!! TG adds air density depedent power curves
IF(airdensity .GT. 0. .AND. airdensity .LT. 0.95) THEN

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=(-7.24254284e6) + 1.11115488e7*speed + (-7.21755964e6)*speed**2 + &
!TG*    2.60343126e6*speed**3 + (-578362.544)*speed**4 + 83445.1807*speed**5 + (-7997.9659)*speed**6 + &
!TG*    511.91011*speed**7 + (-21.5926809)*speed**8 + 0.575797509*speed**9 + &
!TG*    (-0.00879450146)*speed**10 + 5.86229169e-5*speed**11

!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof= (( -0.042412 )*speed**4+  2.756276 *speed**3+   (-49.195333 )*speed**2 +    403.595389*speed+(  -763.930369 ))/ (speed**4+ (  -25.753501 )*speed**3+   242.673623  *speed**2+ &
(-870.990195)*speed+ (1783.558339 ) )
        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.056
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                thrcof=( ( -1.508647)*speed**4+   106.806300*speed**3+  (-2086.758428)*speed**2+14402.562299*speed+ (-7227.247840))/( speed**5 +(-37.027925)*speed**4+  672.633049*speed**3+ &
                ( -6923.593249)*speed**2+   35471.100766*speed+  (-37723.233104))  !TG Added new data
        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.056
        ELSE
                thrcof = 0.056
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
        tkecof=thrcof-powcof
       tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

ELSE IF (airdensity .GE. 0.95 .AND. airdensity .LT. 0.975) THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.068768 )*speed**4+  4.228612 *speed**3+   (    -76.810044 )*speed**2 +    596.522220*speed+( -1109.141114 ))/ (speed**4+ (-24.625814 )*speed**3+212.025380  *speed**2+ &
(-643.628030)*speed+ (1520.519641 ) )
        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.055
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof=( ( -1.479477)*speed**4+   103.446978*speed**3+  (-1979.484341)*speed**2+13499.262713*speed+ (-9185.740624))/( speed**5 +(-35.462349)*speed**4+  616.201059*speed**3+ &
        ( -6175.989471)*speed**2+   31321.474273*speed+  (-34726.520430)) !TG Added thr curve

        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.055
        ELSE
                thrcof = 0.055
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
       tkecof=thrcof-powcof
       tkecof = correction_factor * tkecof  
       IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

!BIN3 0.975 to 1.00
ELSE IF (airdensity .GE. 0.975 .AND. airdensity .LT. 1.000) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=830400.945 + 1.20113511e6*speed + (-2.09450677e6)*speed**2 + 1.13096578e6*speed**3 + &
!TG*    (-316290.849)*speed**4 + 53011.1622*speed**5 + (-5632.85507)*speed**6 + 388.219026*speed**7 + &
!TG*    (-17.3041778)*speed**8 + 0.481517325*speed**9 + (-0.00760800041)*speed**10 + &
!TG*    5.21362168e-5*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.081523 )*speed**4+  4.923965 *speed**3+   (-89.499173 )*speed**2 +    682.861879*speed+( -1273.432153 ))/ (speed**4+ (  -23.820955 )*speed**3+   190.939884  *speed**2+ &
(-482.958685)*speed+ (1212.033297 ) )

        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.054
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof=( ( -0.850984)*speed**4+    66.351438*speed**3+  (-1251.800362)*speed**2+ 8061.429917*speed+ (-3290.162373))/( speed**5 +(-34.990669)*speed**4+  560.396626*speed**3+ &
        ( -5010.743396)*speed**2+   22747.107885*speed+  (-23287.039024)) !TG Added thrcof

        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.054
        ELSE
                thrcof = 0.054
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
        tkecof = correction_factor * tkecof 
	IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

!BIN4 1 to 1.025 ******************
ELSE IF (airdensity .GE. 1.0 .AND. airdensity .LT. 1.025) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
!TG*    power=4.91317077e6 + (-3.7990088e6)*speed + 484887.782*speed**2 + 390862.284*speed**3 + &
!TG*    (-184738.002)*speed**4 + 37748.3903*speed**5 + (-4447.48299)*speed**6 + 326.249002*speed**7 + &
!TG*    (-15.1559336)*speed**8 + 0.434287103*speed**9 + (-0.0070134649)*speed**10 + &
!TG*    4.88842654e-5*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.086546 )*speed**4+  5.167414 *speed**3+   (    -93.263281 )*speed**2 +    697.571940*speed+( -1286.730088 ))/ (speed**4+ (  -23.708951 )*speed**3+   187.522690  *speed**2+ &
(-469.118198)*speed+ (1188.051099 ) )


        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.

  ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.053
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof= ( ( -0.635447)*speed**4+    52.428431*speed**3+  ( -967.367261)*speed**2+ 5749.439702*speed+ (  414.791234))/( speed**5 +(-35.432109)*speed**4+  556.073972*speed**3+&
         ( -4731.968085)*speed**2+   19977.564001*speed+  (-18458.852395))


        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.053
        ELSE
                thrcof = 0.053
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
	 tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF
             !BIN5 1.025 to 1.05
ELSE IF (airdensity .GE. 1.025 .AND. airdensity .LT. 1.05) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*    power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
!TG*    power=8.7817508e6 + (-8.54064203e6)*speed + 2.93227788e6*speed**2 + (-311315.659)*speed**3 + &
!TG*    (-60049.7455)*speed**4 + 23312.774*speed**5 + (-3330.15936)*speed**6 + 268.11726*speed**7 + &
!TG*    (-13.1536951)*speed**8 + 0.390635887*speed**9 + (-0.00646996258)*speed**10 + &
!TG*    4.59538544e-5*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.105759 )*speed**4+  6.241584 *speed**3+   (   -113.499208 )*speed**2 +    841.536275*speed+( -1548.239609 ))/ (speed**4+ (  -22.616362 )*speed**3+   158.635166  *speed**2+ &
(-247.387193)*speed+ ( 851.429834 ) )


        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.051
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof= ( ( -0.497413)*speed**4+    42.457877*speed**3+  ( -755.972783)*speed**2+ 3837.041206*speed+ ( 4978.754667))/( speed**5 +(-36.678836)*speed**4+  578.697071*speed**3+&
        ( -4788.720224)*speed**2+   18872.871130*speed+  (-14626.435592))


        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.051
        ELSE
                thrcof = 0.051
        ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
       tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

!BIN6 1.05 to 1.075
ELSE IF (airdensity .GE. 1.05 .AND. airdensity .LT. 1.075) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=1.45483587e7 + (-1.58900179e7)*speed + 6.90773917e6*speed**2 + &
!TG*    (-1.51713492e6)*speed**3 + 168411.923*speed**4 + (-5184.49417)*speed**5 + &
!TG*    (-929.887197)*speed**6 + 130.818827*speed**7 + (-7.89757119)*speed**8 + &
!TG*    0.26172869*speed**9 + (-0.00463898502)*speed**10 + 3.44990675e-5*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.099784 )*speed**4+  5.840045 *speed**3+   (   -104.589876 )*speed**2 +    763.717357*speed+( -1389.683201 ))/ (speed**4+ (  -22.928764 )*speed**3+   168.460215  *speed**2+ &
(-343.129152)*speed+ ( 977.703606 ) )

        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.

 ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.05
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof= ( ( -0.448738)*speed**4+    39.857156*speed**3+  ( -706.972596)*speed**2+ 3643.814744*speed+ ( 3492.182204))/( speed**5 +(-35.676367)*speed**4+  546.879096*speed**3+&
                 ( -4414.170362)*speed**2+   17148.151536*speed+  (-13884.936469))

        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.05
        ELSE
                thrcof = 0.05
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
       tkecof=thrcof-powcof
       tkecof = correction_factor * tkecof         
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

!BIN6 ends***************
!BIN7 1.075 to 1.1
ELSE IF (airdensity .GE. 1.075 .AND. airdensity .LT. 1.1) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=1.99230964e7 + (-2.27668628e7)*speed + 1.06398358e7*speed**2 + &
!TG*    (-2.65186784e6)*speed**3 + 383735.996*speed**4 + (-32063.6112)*speed**5 + &
!TG*    1334.41748*speed**6 + 1.34025813*speed**7 + (-2.94414226)*speed**8 + &
!TG*    0.140356798*speed**9 + (-0.00291692225)*speed**10 + 2.37387317e-5*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.104839 )*speed**4+  6.085602 *speed**3+   (   -108.484820 )*speed**2 +    786.979709*speed+( -1430.211537 ))/ (speed**4+ (  -22.311274 )*speed**3+   154.311284  *speed**2+ &
(-241.170641)*speed+ ( 752.676803 ) )
        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.049
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof= ( ( -0.401960)*speed**4+    36.844081*speed**3+  ( -648.289285)*speed**2+ 3259.361016*speed+ ( 3325.829954))/( speed**5 +(-35.482687)*speed**4+  537.581530*speed**3+&
                 ( -4267.695165)*speed**2+   16265.783042*speed+  (-13270.764728)) !BIN7

        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.049
        ELSE
                thrcof = 0.049
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         
        tkecof=thrcof-powcof
       tkecof = correction_factor * tkecof
	 IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

        !BIN8 Starts 1.1to1.125
ELSE IF (airdensity .GE. 1.1 .AND. airdensity .LT. 1.125) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=2.55723142e7 + (-2.99698721e7)*speed + 1.45370428e7*speed**2 + &
!TG*    (-3.83391749e6)*speed**3 + 607629.519*speed**4 + (-59976.2756)*speed**5 + &
!TG*    3683.72604*speed**6 + (-132.923081)*speed**7 + 2.19059513*speed**8 + 0.0145645648*speed**9 + &
!TG*    (-0.00113224947)*speed**10 + 1.25869195e-5*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.104808 )*speed**4+  6.036608 *speed**3+   (   -106.649920 )*speed**2 +    764.012886*speed+( -1377.660136 ))/ (speed**4+ (  -22.319410 )*speed**3+   155.401889  *speed**2+ &
(-262.762453)*speed+ ( 771.020709 ) )

        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.048
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof= ( ( -0.313841)*speed**4+    31.333599*speed**3+  ( -540.374020)*speed**2+ 2450.245871*speed+ ( 4519.732033))/( speed**5 +(-35.497233)*speed**4+  532.418684*speed**3+&
                ( -4135.297515)*speed**2+   15203.597905*speed+  (-11576.773578))

        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.048
        ELSE
                thrcof = 0.048
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
	tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

        !BIN9 1.125 to 1.15
ELSE IF (airdensity .GE. 1.125 .AND. airdensity .LT. 1.15) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=3.08339865e7 + (-3.66915198e7)*speed + 1.81798213e7*speed**2 +(-4.94012906e6)*speed**3 + &
!TG*    817303.909*speed**4 + (-86121.1353)*speed**5 + 5883.69004*speed**6 + (-258.575809)*speed**7 + &
!TG*    6.99187689*speed**8 + (-0.102934255)*speed**9 + 0.000532769915*speed**10 + &
!TG*    2.19619732e-6*speed**11



!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.114170 )*speed**4+  6.532063 *speed**3+   (   -115.454545 )*speed**2 +    823.936243*speed+( -1485.811246 ))/ (speed**4+ (  -21.548118 )*speed**3+   136.769392  *speed**2+ &
(-124.839032)*speed+ ( 515.265320 ) )
        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.048
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof= ( ( -0.214338)*speed**4+    25.370478*speed**3+  ( -427.251065)*speed**2+ 1617.109827*speed+ ( 5841.306009))/( speed**5 +(-35.550784)*speed**4+  528.155431*speed**3+&
                 ( -4011.101315)*speed**2+   14176.452233*speed+  (-9834.849555 ))



        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.048
        ELSE
                thrcof = 0.048

        ENDIF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
          tkecof = correction_factor * tkecof
	IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF


        !BIN10 1.15 to 1.175
ELSE IF (airdensity .GE. 1.15 .AND. airdensity .LT. 1.175) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=3.55532748e7 + (-4.28710602e7)*speed + 2.16178087e7*speed**2 + (-6.01309e6)*speed**3 + &
!TG*    1.02644101e6*speed**4 + (-112940.709)*speed**5 + 8203.83998*speed**6 + &
!TG*    (-394.731005)*speed**7 + 12.3329518*speed**8 + (-0.236999754)*speed**9 + &
!TG*    0.00247945963*speed**10 + (-1.02401599e-5)*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.105509 )*speed**4+  6.004156 *speed**3+   (   -104.749993 )*speed**2 +    739.651502*speed+( -1325.008044 ))/ (speed**4+ (  -21.889466 )*speed**3+   147.187847  *speed**2+ &
(-216.718963)*speed+ ( 615.168283 ) )

        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
       !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.047
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof=( ( -0.202778)*speed**4+    24.509352*speed**3+  ( -407.695219)*speed**2+ 1499.137990*speed+ ( 5664.197392))/( speed**5 +(-35.161842)*speed**4+  516.644662*speed**3+&
                 ( -3880.671060)*speed**2+   13561.696507*speed+  (-9290.093827 ))


        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.047
        ELSE
                thrcof = 0.047
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
	 tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

        !BIN11 1.175 to 1.2
ELSE IF (airdensity .GE. 1.175 .AND. airdensity .LT. 1.2) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=4.01527821e7 + (-4.88878756e7)*speed + 2.49635106e7*speed**2 + &
!TG*    (-7.05698001e6)*speed**3 + 1.22987027e6*speed**4 + (-139020.696)*speed**5 + &
!TG*    10459.033*speed**6 + (-526.998749)*speed**7 + 17.5179753*speed**8 + (-0.367047919)*speed**9 + &
!TG*    0.0043662491*speed**10 + (-2.22834004e-5)*speed**11


!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.110001 )*speed**4+  6.223378 *speed**3+   (   -108.328625 )*speed**2 +    762.570300*speed+( -1364.983328 ))/ (speed**4+ (  -21.345502 )*speed**3+   135.290751  *speed**2+ &
(-134.037990)*speed+ ( 452.783735 ) )


        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
  IF(speed .LT. cispeed) THEN
                thrcof = 0.046
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof=( ( -0.148595)*speed**4+    21.320495*speed**3+  ( -348.875589)*speed**2+ 1118.638954*speed+ ( 5905.258239))/( speed**5 +(-34.964506)*speed**4+  508.341978*speed**3+&
                 ( -3756.656457)*speed**2+   12832.147311*speed+  (-8505.836192 ))


        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.046
        ELSE
                thrcof = 0.046
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
	 tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF

        !BIN 12 1.2 to 1.225
ELSE IF (airdensity .GE. 1.2 .AND. airdensity .LT. 1.225) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER  
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.107056 )*speed**4+  6.021184 *speed**3+   (   -103.871874 )*speed**2 +    723.874960*speed+( -1287.432451 ))/ (speed**4+ (  -21.467733 )*speed**3+   139.460182  *speed**2+&
(-176.674249)*speed+ ( 502.845473 ) )

        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.045
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof= ( ( -0.117564)*speed**4+    19.458595*speed**3+  ( -315.861546)*speed**2+  937.998670*speed+ ( 5698.964698))/( speed**5 +(-34.765599)*speed**4+  500.770054*speed**3+&
                 ( -3654.098257)*speed**2+   12297.098389*speed+  (-8180.646230 ))


        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.045
        ELSE
                thrcof = 0.045
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
	 tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF
   !BIN13 1.225 to 1.25
ELSE IF (airdensity .GE. 1.225 .AND. airdensity .LT. 1.25) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=4.8033481e7 + (-5.9349306e7)*speed + 3.08644627e7*speed**2 + (-8.92329508e6)*speed**3 + &
!TG*    1.59819883e6*speed**4 + (-186791.145)*speed**5 + 14633.5279*speed**6 + &
!TG*    (-774.180041)*speed**7 + 27.2920068*speed**8 + (-0.614137868)*speed**9 + &
!TG*    0.00797710716*speed**10 + (-4.54846137e-5)*speed**11


!Tg*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.106748 )*speed**4+  5.980493 *speed**3+   (   -102.747785 )*speed**2 +    716.006829*speed+( -1277.515608 ))/ (speed**4+ (  -20.976760 )*speed**3+   129.682706  *speed**2+ &
(-107.711862)*speed+ ( 322.263709 ) )

        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.044
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof=( ( -0.112184)*speed**4+    19.374107*speed**3+  ( -317.980122)*speed**2+ 1043.491904*speed+ ( 4841.270979))/( speed**5 +(-34.275795)*speed**4+  487.295294*speed**3+&
       ( -3520.969312)*speed**2+   11827.539906*speed+  (-8261.959983 ))
        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.044
        ELSE
    thrcof = 0.044
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
	 tkecof = correction_factor * tkecof
         IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF
        !BIN 14 1.25 to 1.275

ELSE IF (airdensity .GE. 1.25 .AND. airdensity .LT. 1.5) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TURBINE POWER      !!!!!!!!!!
!TG*    IF(speed .LT. cispeed) THEN
!TG*            power=0.
!TG*    ELSE IF(speed .GE. cispeed .AND. speed .LT. 12.0) THEN
        ! for testing, will comment out existing power curve for another directly below
        !power = 2.25E6*(tanh((speed - 9.1569)/2.8016447) + 1)
!TG*    power=5.14153879e7 + (-6.38873012e7)*speed + 3.34528644e7*speed**2 + &
!TG*    (-9.75122737e6)*speed**3 + 1.7634223e6*speed**4 + (-208447.624)*speed**5 + &
!Tg*    16544.6999*speed**6 + (-888.369783)*speed**7 + 31.844598*speed**8 + (-0.730094797)*speed**9 + &
!TG*    0.00968327825*speed**10 + (-5.65161915e-5)*speed**11



!TG*    ELSE IF(speed .GE. 12 .AND. speed .LE. 25.0) THEN
!TG*            power=3075000.
!TG*    ELSE
!TG*            power=0.
!TG*    ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! POWER  COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                powcof=0.
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
                powcof=(( -0.101708 )*speed**4+  5.669153 *speed**3+   (    -96.438986 )*speed**2 +    665.863297*speed+( -1180.588577 ))/ (speed**4+ (  -21.216195 )*speed**3+   136.904453  *speed**2+ &
(-172.518200)*speed+ ( 412.312728 ) )

        ELSE IF(speed .GT. 25.0) THEN
                powcof=0.
        ELSE
                powcof=0.
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! THRUST COEFFICIENT !!!!!!!!!!
        IF(speed .LT. cispeed) THEN
                thrcof = 0.044
        ELSE IF(speed .GE. cispeed .AND. speed .LE. 25.0) THEN
        thrcof=( ( -0.020511)*speed**4+    14.412174*speed**3+  ( -231.163559)*speed**2+  477.082458*speed+ ( 5594.877159))/( speed**5 +(-34.188057)*speed**4+  481.232789*speed**3+&
        ( -3407.747627)*speed**2+   11066.799360*speed+  (-7160.727565 ))


        ELSE IF(speed .GT. 25.0) THEN
                thrcof = 0.044
        ELSE
                thrcof = 0.044
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!! TKE COEFFICIENT !!!!!!!!!!!!!
         tkecof=thrcof-powcof
        tkecof = correction_factor * tkecof 
	IF(tkecof .LT. 0.) THEN
        !   OPEN(file=filename,unit=133,POSITION='APPEND')
        !   WRITE(133,*) "Congratulations! You just invented the perpetual motion machine!"
        !   CLOSE(133)
    tkecof=0.
        ENDIF
ELSE
        powcof=0.
        thrcof=0.
        tkecof=0.
ENDIF  

END SUBROUTINE dragcof


!
  SUBROUTINE init_module_wind_fitch(id,config_flags,xlong,xlat,windfarm_initialized,&
                                            ims,ime,jms,jme,its,ite,jts,jte,ids,ide,jds,jde)
!
  IMPLICIT NONE
!
   integer ims,ime,jms,jme,ids,ide,jds,jde
   integer its,ite,jts,jte
   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: xlong,xlat
   TYPE (grid_config_rec_type) :: config_flags
   TYPE (PROJ_INFO) :: ts_proj
   logical :: windfarm_initialized
! 
   CHARACTER*256 num,input,message_wind
   real lat,lon,ts_rx,ts_ry
   REAL :: known_lat, known_lon
   INTEGER i,j,nval,k,id

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
!

   correction_factor = config_flags%windfarm_tke_factor
      IF ( wrf_dm_on_monitor() ) THEN
!
! ... PAJ: Opens the file with the location of the wind turbines ...
!
        if ( config_flags%windfarm_ij .eq. 1 ) then
          open(70,file='windturbines-ij.txt',form='formatted',status='old')
        else
          open(70,file='windturbines.txt',form='formatted',status='old')
        end if
!
! ... PAJ: Counts the turbines ...
!
       nt=0
 10    read(70,*,end=100) 
       nt=nt+1
       goto 10
!
 100   continue
       rewind (70)
     END IF
!
     CALL wrf_dm_bcast_integer(nt,1)
!
! ... PAJ: Initializes the configuration of the wind farm(s) ...
!
     if (.not. windfarm_initialized) then
       allocate (nkind(nt),ival(nt,max_domains),jval(nt,max_domains))
       allocate (hubheight(nt),stc(nt),stc2(nt),cutin(nt),cutout(nt),diameter(nt),npower(nt))
       ival=-9999
       jval=-9999
       windfarm_initialized=.true.
     endif
!
     IF ( wrf_dm_on_monitor() ) THEN
     do k=1,nt
       if ( config_flags%windfarm_ij .eq. 1 ) then
         read(70,*) ival(k,id), jval(k,id), nkind(k)
         write(message_wind,*)'WINDFARM Turbine #',k,': I, J = ',ival(k,id), jval(k,id),'; Type = ',nkind(k)
         CALL wrf_message(message_wind)

       else

         read(70,*)lat,lon,nkind(k)
         write(message_wind,*)'WINDFARM Turbine #',k,': Lat, lon = ',lat,lon,'; Type = ',nkind(k)
         CALL wrf_message(message_wind)

         CALL map_init(ts_proj)

         known_lat = xlat(its,jts)
         known_lon = xlong(its,jts)

      ! Mercator
      IF (config_flags%map_proj == PROJ_MERC) THEN
         CALL map_set(PROJ_MERC, ts_proj,               &
                      truelat1 = config_flags%truelat1, &
                      lat1     = known_lat,             &
                      lon1     = known_lon,             &
                      knowni   = REAL(its),             &
                      knownj   = REAL(jts),             &
                      dx       = config_flags%dx)

      ! Lambert conformal
      ELSE IF (config_flags%map_proj == PROJ_LC) THEN
         CALL map_set(PROJ_LC, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      truelat2 = config_flags%truelat2,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = REAL(its),              &
                      knownj   = REAL(jts),              &
                      dx       = config_flags%dx)
!      ! Polar stereographic
      ELSE IF (config_flags%map_proj == PROJ_PS) THEN
         CALL map_set(PROJ_PS, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = REAL(its),              &
                      knownj   = REAL(jts),              &
                      dx       = config_flags%dx)
!#if (EM_CORE == 1)
!      ! Cassini (global ARW)
!      ELSE IF (config_flags%map_proj == PROJ_CASSINI) THEN
!         CALL map_set(PROJ_CASSINI, ts_proj,                            &
!                      latinc   = grid%dy*360.0/(2.0*EARTH_RADIUS_M*PI), &
!                      loninc   = grid%dx*360.0/(2.0*EARTH_RADIUS_M*PI), &
!                      lat1     = known_lat,                             &
!                      lon1     = known_lon,                             &
!                      lat0     = config_flags%pole_lat,                 &
!                      lon0     = config_flags%pole_lon,                 &
!                      knowni   = 1.,                                    &
!                      knownj   = 1.,                                    &
!                      stdlon   = config_flags%stand_lon)
!#endif
!
!      ! Rotated latitude-longitude
!      ELSE IF (config_flags%map_proj == PROJ_ROTLL) THEN
!         CALL map_set(PROJ_ROTLL, ts_proj,                      &
!! I have no idea how this should work for NMM nested domains
!                      ixdim    = grid%e_we-1,                   &
!                      jydim    = grid%e_sn-1,                   &
!                      phi      = real(grid%e_sn-2)*grid%dy/2.0, &
!                      lambda   = real(grid%e_we-2)*grid%dx,     &
!                      lat1     = config_flags%cen_lat,          &
!                      lon1     = config_flags%cen_lon,          &
!                      latinc   = grid%dy,                       &
!                      loninc   = grid%dx,                       &
!                      stagger  = HH)
!
      END IF
!
         CALL latlon_to_ij(ts_proj, lat, lon, ts_rx, ts_ry)
!
          ival(k,id)=nint(ts_rx)
          jval(k,id)=nint(ts_ry)
          if (ival(k,id).lt.ids.and.ival(k,id).gt.ide) then
            ival(k,id)=-9999
            jval(k,id)=-9999
          endif
!
!         write(73,*) k,id,ival(k,id),jval(k,id)
          write(message_wind,*)'WINDFARM Turbine #',k,': Lat, lon = ',lat,lon, &
                               ', (i,j) = (',ival(k,id),',',jval(k,id),'); Type = ',nkind(k)
          CALL wrf_debug(0,message_wind)
!
     end if
!
     enddo
      close(70)
!
! ... PAJ: Read the tables with the turbine's characteristics ...
!
         turbws=0.
         turbtc=0.
         turbpw=0.
         DO i=1,nt
          write(num,*) nkind(i)
          num=adjustl(num)
          input="wind-turbine-"//trim(num)//".tbl"
          OPEN(file=TRIM(input),unit=19,FORM='FORMATTED',STATUS='OLD')
          READ (19,*,ERR=132)nval
          READ(19,*,ERR=132)hubheight(i),diameter(i),stc(i),npower(i)
            DO k=1,nval
              READ(19,*,ERR=132)turbws(nkind(i),k),turbtc(nkind(i),k),turbpw(nkind(i),k)
            ENDDO
          cutin(i)  = turbws(nkind(i),1)
          cutout(i) = turbws(nkind(i),nval)
          stc2(i) = turbtc(nkind(i),nval)
          close (19)
         ENDDO

 132   continue
!
! ... ...
!
      endif

        CALL wrf_dm_bcast_integer(ival,nt*max_domains)
        CALL wrf_dm_bcast_integer(jval,nt*max_domains)
        CALL wrf_dm_bcast_real(hubheight,nt)
        CALL wrf_dm_bcast_real(diameter,nt)
        CALL wrf_dm_bcast_real(stc,nt)
        CALL wrf_dm_bcast_real(npower,nt)
        CALL wrf_dm_bcast_real(cutin,nt)
        CALL wrf_dm_bcast_real(cutout,nt)
        CALL wrf_dm_bcast_integer(nkind,nt) 
        CALL wrf_dm_bcast_real(turbws,maxvals*maxvals2) 
        CALL wrf_dm_bcast_real(turbtc,maxvals*maxvals2) 
        CALL wrf_dm_bcast_real(turbpw,maxvals*maxvals2) 

  END SUBROUTINE init_module_wind_fitch
  
END MODULE module_wind_fitch
