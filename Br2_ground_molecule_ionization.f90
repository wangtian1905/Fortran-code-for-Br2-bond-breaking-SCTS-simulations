!  Github_Br2_SCTS_codes_full_for_Nat_Comm_final_clean.f90 
!
!  FUNCTIONS:
!  Github_Br2_SCTS_codes_full_for_Nat_Comm_final_clean - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Github_Br2_SCTS_codes_full_for_Nat_Comm_final_clean
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

 
module lsr_para
implicit none
!-----------------------------------------------------------------------------------------------
real*8,parameter::pi=3.141592653589793d0
real*8,parameter::c=2.99792458E8 !<light speed><ISU:m/s> 
real*8,parameter::w0=4.1341374457476439368E16 !<unit conversion factor for circular frequency>
real*8,parameter::I0=35.092d0 !<unit conversion factor for intensity><*10^15W/cm^2><PW/cm^2>
real*8,parameter::t0=0.0242 !<unit conversion factor for time><fs>
real*8,parameter::L_L=800.0E-9 !<wavelength><ISU:m>
real*8,parameter::L_W=(2.0*pi*c/L_L)/w0 !<circular frequency><a.u.>
real*8,parameter::L_T=(2*pi)/L_W !<duration time of one laser cycle><a.u.>
real*8,parameter::L_I=0.05 !<laser intensity><ISU:PW/cm^2>
real*8,parameter::L_PC=4 !<cycle number of the plateau of the laser fields>
real*8,parameter::L_TC=6    !<cycle number of the slope of the laser fields>
real*8,parameter::L_FC=0   !<cycle number of free-fly of the electron>
real*8,parameter::ceph=0.0  !<carrier envelope phase of laser fields><radian>   
real*8,parameter::T_PD=L_PC*L_T !<duration time of plateau of the laser fields><a.u.>
real*8,parameter::T_TD=L_TC*L_T !<duration time of slope of the laser fields><a.u.>
real*8,parameter::T_FD=L_FC*L_T !<duration time of free-fly of the electron><a.u.>
!-----------------------------------------------------------------------------------------------
end module  

module sim_para
use lsr_para
implicit none
!-----------------------------------------------------------------------------------------------
integer,parameter::gn=500000000!500000000

real*8,parameter::r_impct=0.5 !<the assummed impact crosssection radius of tunneled electron projecting on bonding electron of Br2+ ion><a.u.>
real*8,parameter::ion_fbd_TC=4
real*8,parameter::ion_fbd_tt=ion_fbd_TC*L_T
!real*8,parameter::d0=4.31!<bond length of Br2 molecule><a.u.>
!real*8,parameter::algn_thta_d=0.0 !<molecular alignment angle: theta><degree> 
!real*8,parameter::algn_thta_r=(algn_thta_d/180.d0)*pi !<molecular alignment angle: theta><radian>
real*8,parameter::tnl_t_lb=0.0 !<lowerbound of electron tunneling time><o.c.>
real*8,parameter::tnl_t_ub=4.0!4.0 !<upperbound of electron tunneling time><o.c.>
real*8,parameter::tnl_t_lb2=0.0 !<lowerbound of electron tunneling time><o.c.>
real*8,parameter::tnl_t_ub2=4.0!4.0!<upperbound of electron tunneling time><o.c.>
real*8,parameter::ttpk_wdw=0.0*L_T
real*8,parameter::ttcs_wdw=0.0*L_T!<less than (1/16)*L_T>
real*8,parameter::E_jdg=0.1 !<electric field strength judgment to remove zero field strength points> 
real*8,parameter::dlt_px=0.0!<px momentum shift for compensation of few-cycle induced momentum drift in turn-down field>

!pk_wdw,cs_wdw 
real*8,parameter::tnl_tt_lb=tnl_t_lb*L_T !<lowerbound of electron tunneling time><a.u.>
real*8,parameter::tnl_tt_ub=tnl_t_ub*L_T !<upperbound of electron tunneling time><a.u.>
real*8,parameter::tnl_tt_lb2=tnl_t_lb2*L_T !<lowerbound of electron tunneling time><a.u.>
real*8,parameter::tnl_tt_ub2=tnl_t_ub2*L_T !<upperbound of electron tunneling time><a.u.>
!real*8,parameter::tol=0.000000001!<error control tolerance for RKF45 algorithm!>
!real*8,parameter::pper_l=-1.0 !<lowerbound of the transverse momentum of tunneling electron><a.u.>
!real*8,parameter::pper_u=1.0 !<upperbound of the transverse momentum of tunneling electron><a.u.>
!real*8,parameter::r_jdg=5.0 !<position ionization-criterion for judging ionization status of tunneled electrons><a.u.>
!-----------------------------------------------------------------------------------------------
end module

module Sym_ctrl
implicit none
real*8,parameter::r_crt1=20!<radius criterion for near nucleus region><a.u.>
real*8,parameter::h1=0.5  !< time step for near nucleus region(extremely small)><a.u.>!0.02！latest old value:0.05
real*8,parameter::r_crt2=200!<radius criterion for far nucleus region><a.u.>!200
real*8,parameter::h2=2.0!< time step for far nucleus region(moderate small)><a.u.>(h2=0.005)
real*8,parameter::h3=8.0!< time step for too far nulceus region(not very small)><a.u.>(h3=0.05)!0.1
end module
 
module omp_para 
use sim_para
!-----------------------------------------------------------------------------------------------
integer,parameter::n_threads=10!<the total number of parallel threads><OMP>
integer,parameter::sn=gn/n_threads !<the size of each sub-ensemble in each thread><OMP>
!----------------------------------------------------------------------------------------------- 
end module

module MO_para
use lsr_para
implicit none
real*8,parameter::nn_Z1=0.5!0.5!0.5!1.0 !<nucleus charge for the first cucleus potential(here is for molecular potential)><a.u.> 
real*8,parameter::nn_Z2=0.5!0.5 !<nucleus charge for the second cucleus potential(here is for molecular potential)><a.u.>
real*8,parameter::n_Z=1.0!1.0 !<nucleus charge for the molecular ion when electron has flied far away from it(only for asymptotic momentum)><a.u.>
real*8,parameter::ne1=0.0 !<soft-core parameter for the first nucleus><a.u.>
real*8,parameter::ne2=0.0 !<soft-core parameter for the second nucleus><a.u.>
real*8,parameter::Ip1=0.457!0.434!<the ionization potential of Br2 pi_u orbital:0.457;for Br 4p orbital: 0.434><a.u.>
real*8,parameter::Ip2=0.457!0.434!0.38655!0.38655!<the ionization potential of pi_g orbital><a.u.>
real*8,parameter::prb_orb1=0.5!probability of the first trajectory group<contribution>
real*8,parameter::orb_prb_chk1=0.0d0
real*8,parameter::orb_prb_chk2=prb_orb1
real*8,parameter::orb_prb_chk3=1.0d0
real*8,parameter::R_M1=1.1!1.1!1.05!1.0!1.0125!1.0125!<radius of GSZ potential for upper multielectron ion core><a.u.>
real*8,parameter::R_M2=1.1!1.1!1.05!1.0!1.0125!1.0125!<radius of GSZ potential for lower multielectron ion core><a.u.>
real*8,parameter::apha_dg=90.0 !<alignment angle of molecule><degree>!old varbl:thta_dg
real*8,parameter::apha=(apha_dg/180)*pi !<alignment angle of molecule><radian>
real*8,parameter::n_R_Atm=4.3!0.0!4.3!4.3!<internucleus distance of ground Br2 molecule; for Br case, set zero!><a.u.>
real*8,parameter::n_R_Mlc=4.3!0.0<internucleus distance of ground Br2 molecule; for Br case, set zero!><a.u.>
real*8,parameter::apha_N=0.0!can be set within small range: 0.0~2.0!<polarizability of atom for considering stark shifted potential in initial state preparation><a.u.><This value represents the polarizability difference,not real polarizability!>
real*8,parameter::apha_I=0.0!<polarizability of ion for considering stark shifted potential in initial state preparation><a.u.>
real*8,parameter::apha_II_Atm=40!<atomic polarizability contribution for BBS ionization dynamics><atomic Br+:12;molecular ground Br2+:40;molecular stretched Br2+: 80><a.u.>
real*8,parameter::apha_II_Mlc=40!<molecular polarizability contribution for BBS ionization dynamics>atomic Br+:12;molecular ground Br2+:40;molecular stretched Br2+: 80><a.u.>
integer,parameter::p_per_div_n=2000
real*8,parameter::p_per_lb=-2.5
real*8,parameter::p_per_ub=2.5
real*8,parameter::p_per_dlt=(p_per_ub-p_per_lb)/p_per_div_n
end module

module intrfr_para
!-----------------------------------------------------------------------------------------------
integer,parameter::n_x=1000!1000!<resolution control:bin size on x direction><previous value:1000>
integer,parameter::n_y=1000!1000!<resolution control:bin size on y direction><previous value:1000>
integer,parameter::mn_x=1!<detail resolution control:bin size on x direction>
integer,parameter::mn_y=1!<detail resolution control:bin size on y direction>
real*8,parameter::x_lb=-2.0!<lowerbound of the momentum space on x direction>
real*8,parameter::x_ub=2.0!<upperbound of the momentum space on x direction>
real*8,parameter::y_lb=-2.0!<lowerbound of the momentum space on y direction>
real*8,parameter::y_ub=2.0!<upperbound of the momentum space on y direction> 
integer,parameter::len_X=n_x*mn_x!<momemtum space matrix:length of the matrix Spec on x direction>
integer,parameter::len_Y=n_y*mn_y!<momemtum space matrix:length of the matrix Spec on y direction>
real*8,parameter::width_x=(x_ub-x_lb)/n_x
real*8,parameter::width_y=(y_ub-y_lb)/n_y
real*8,parameter::mwidth_x=width_x/mn_x
real*8,parameter::mwidth_y=width_y/mn_y
!-------------check initial transverse velocity----------
integer,parameter::py_n_v=200
real*8,parameter::py_lb=-1.0
real*8,parameter::py_ub=1.0
!-----------------------------------------------------------------------------------------------
end module

program main
use omp_lib
use lsr_para
use sim_para
use MO_para
use omp_para
use intrfr_para
use Sym_ctrl

!use array_pub
implicit none
!-----------------------------------------------------------------------------------------------
real*8::h,rr_crt
real*8::dice_tnlt
real*8::px,py,pz,rx,ry,rz,rr,t,tx,wght,MO_wght,tnlt0,pper,pper_max,Ex0,E00,tnl_extpt0,S_ph,S_phx,S_phx0,S_phx1,fnl_C_ph,MO_ph,Spec_max
real*8::E_e,elec_egy0,elec_rr,Ip_ph
integer::chcka,chckb,chckb_ion,rjct
real*8::chcka_t0,chcka_dt,chckb_px0,chckb_t0,chckb_dt
integer::nstp_fr
real*8::t_fr
!real*8::y(1:yn) 
real*8::f_t0,f_t1
real*8::asym_px,asym_py,asym_pz,px_tst,py_tst
!integer::n_stp,n_stpx
integer::i,j,k,ii,jj
integer::thrd_num
integer::ion_jdg
integer::good_trj
real*8,external::random_grab,tnl_extpt,tnl_wght,E,f,H_ph,elec_egy,MO_dnsty
real*8::px0,pxx,MO_intgrl
real*8::Spec(1:len_Y,1:len_X),Spec_wght(1:len_Y,1:len_X)
real*8::M_cos(1:len_Y,1:len_X),M_sin(1:len_Y,1:len_X),Spec_grdX(1:len_X),Spec_grdY(1:len_Y)
real*8::p_per_grid(1:p_per_div_n),p_per_wght1(1:p_per_div_n),p_per_wght2(1:p_per_div_n)
real*8::p_per1,p_per2,pper_wght1,pper_wght2
real*8::p_per_vlu
integer::p_per_i 
real*8::rn1,rn2,rn1_min,rn2_min
real*8::dice_orb,Ip
integer::orb_labl
real*8::n_Z1,n_Z2,n_R,apha_II
 real*8::orb_dis_lbl
!-----------------------------------------------------------------------------------------------
!real*8::M_cos(1:len_Y,1:len_X),M_sin(1:len_Y,1:len_X),Spec_grdX(1:len_X),Spec_grdY(1:len_Y)
!-----------------------------------------------------------------------------------------------
write(*,*)"checkpoint1:start"
open(unit=10,file="Ini_state_tnl.txt")
open(unit=11,file="Laser_field.txt")
open(unit=12,file="OpenMP_test.txt")
open(unit=13,file="Spec_gridX.txt")    
open(unit=14,file="Spec_gridY.txt")
open(unit=15,file="Spec_Z.txt")
open(unit=16,file="Spec_Z_norm.txt")
open(unit=17,file="MO_gridX.txt")
open(unit=18,file="MO_gridY.txt")
open(unit=19,file="MO_pxy.txt")
open(unit=20,file="MO_per.txt")
open(unit=21,file="tran_vlcty_test.txt")
open(unit=22,file="fnl_stat_asymp_test.txt")
open(unit=23,file="fnl_stat_free_evlvd.txt")
open(unit=24,file="Start.txt")
open(unit=25,file="p_per_dis_test.txt")
open(unit=26,file="orb_labl_dis_test.txt")
open(unit=27,file="orb_labl_dis_debug.txt")
open(unit=28,file="chk_pnt_debug.txt")  
!open(unit=29,file="Dis_tran_MO1.txt")
!open(unit=30,file="Dis_tran_MO2.txt")
write(*,*)"checkpoint2"
!Matrices initilization
do i=1,len_Y  
do j=1,len_X
    M_cos(i,j)=0.d0
    M_sin(i,j)=0.d0
    Spec_wght(i,j)=0.d0
end do
end do

do i=1,p_per_div_n
!read(29,*) p_per_grid(i),p_per_wght1(i)!p_per_grid is useless!
!read(30,*) p_per_grid(i),p_per_wght2(i)!p_per_grid is useless!
end do


write(24,*) 'Program starts!'

!$OMP PARALLEL NUM_THREADS(n_threads)  DEFAULT(private) SHARED(M_cos,M_sin,p_per_wght1,p_per_wght2)
write(*,*)"checkpoint3:OpenMP"
thrd_num = OMP_GET_THREAD_NUM()
write(*,*) "OMP_thread_num:",thrd_num
i=0
do while(i<sn)!i=1,sn
good_trj=0
dice_tnlt=random_grab(-1.d0,1.d0)
if(dice_tnlt<0.0)then
tnlt0=random_grab(tnl_tt_lb,tnl_tt_ub)
else if(dice_tnlt>0.0)then
tnlt0=random_grab(tnl_tt_lb2,tnl_tt_ub2)
end if
!call tnlt0_grab(ttpk_wdw,ttcs_wdw,tnl_t_lb,tnl_t_ub,tnlt0)
!pper=random_grab(pper_l,pper_u)
 
Ex0=E(tnlt0)
E00=sqrt(L_I/I0)
if(abs(Ex0)<(E_jdg*E00))then
    good_trj=1
    cycle
end if
 !tnl_extpt0=tnl_extpt(Ex0) 
dice_orb=random_grab(0.d0,1.0d0)
if((dice_orb>orb_prb_chk1).and.(dice_orb<orb_prb_chk2))then
    orb_labl=1!orbital label for sigms_u MO!
else if((dice_orb>orb_prb_chk2).and.(dice_orb<orb_prb_chk3))then
    orb_labl=2!orbital label for sigms_g MO!
end if

pper_max=3.0*(abs(Ex0)/sqrt(2*Ip))**0.5
pper=random_grab(-2.d0,2.d0) !pper=random_grab((-pper_max),pper_max)
px=0.0
py=pper
rx=tnl_extpt(Ex0,Ip)
!if(Ex0>0.0)then
!rx=tnl_extpt0
!else if(Ex0<0.0)then
!    rx=-tnl_extpt(Ip,(-Ex0))
!end if
ry=0.0

p_per_i=int((pper-p_per_lb)/p_per_dlt)

if(orb_labl==1)then
    Ip=Ip1
    n_Z1=nn_Z1
    n_Z2=nn_Z2
    n_R=n_R_Mlc
    apha_II=apha_II_Mlc
else if(orb_labl==2)then
    Ip=Ip2
    n_Z1=nn_Z1
    n_Z2=nn_Z2
    n_R=n_R_Atm
    apha_II=apha_II_Atm
end if
MO_ph=0.d0
wght=tnl_wght(Ip,abs(Ex0),pper)
!write(24,*) py,wght

S_ph=0.0
S_phx=0.0
f_t0=0.0
f_t1=0.0
!MO_ph=0.0
!write(25,500)py,tnl_wght(Ip,abs(Ex0),pper),pper_wght1,wght
!write(10,*)p_per_i
!write(10,500)px,py,rx,ry,(tnlt0/L_T),(p_per2-p_per1),Ex0,wght

!ion_jdg=0

!********check initial transverse velocity(1)*********
!py_0(thrd_num*sn+i)=py
!*******************************
!********time eveolution (include initial preparation and main loop!)*********
!********initial preparation********


t=tnlt0
 S_phx=0.d0 
S_phx0=H_ph(px,py,rx,ry,n_Z1,n_Z2,n_R,apha_II)
ion_jdg=0
h=h1

rn1_min=sqrt(rx**2+(ry-0.5*n_R)**2)
rn2_min=sqrt(rx**2+(ry+0.5*n_R)**2)
!********************************
!********main loop begin!********
do while(t<(T_PD+T_TD))
   
   
call sym_adp_algrthm(px,py,rx,ry,t,h,n_Z1,n_Z2,n_R,apha_II)
    !step-adaptive-Symplectic algrithm end*** 
         t=t+h 
    ! Do integral for clasical action***
   S_phx1=H_ph(px,py,rx,ry,n_Z1,n_Z2,n_R,apha_II)
   S_phx=S_phx+0.5*(S_phx0+S_phx1)*h
   S_phx0=S_phx1
   elec_rr=sqrt(rx**2+ry**2)
   rn1=sqrt(rx**2+(ry-0.5*n_R)**2)
   rn2=sqrt(rx**2+(ry+0.5*n_R)**2)
   
   if(rn1<=rn1_min)then
       rn1_min=rn1
   end if
   
   if(rn2<=rn2_min)then
       rn2_min=rn2
   end if

end do

if ((rn1_min<R_M1).or.(rn2_min<R_M2))then
    good_trj=1 
end if

if(good_trj==0)then
    i=i+1

!write(19,500)px,py,rx,ry 
 call ion_jdg_egy(px,py,rx,ry,ion_jdg)
if(ion_jdg==1)then
call fnl_Coulomb_phase(px,py,rx,ry,fnl_C_ph)
!end if
!end do 
!write(23,500)px,py,rx,ry 
!Ip_ph=(1/4)*(Ip+apha_N-apha_I)*Ex0**2
Ip_ph=Ip+(1/2)*(apha_N-apha_I)*Ex0**2
S_ph=S_phx+fnl_C_ph+Ip_ph*tnlt0+MO_ph
!S_ph=mod(S_phx,(2*pi))
pz=0.0
rz=0.0
call asympt_momx(px,py,pz,rx,ry,rz,asym_px,asym_py,asym_pz)
asym_px=asym_px+dlt_px
!asym_px=px
!asym_py=py
!write(22,500)asym_px,asym_py,asym_pz
call sta_ana_intrfr_2D(asym_px,asym_py,wght,S_ph,n_x,n_y,mn_x,mn_y,x_lb,x_ub,y_lb,y_ub,len_X,len_Y,M_cos,M_sin)

if(mod(i,5000)==0)then
   write(12,*)thrd_num,i 
end if
!call sta_ana_wght_2D(px,py,wght,n_x,n_y,mn_x,mn_y,x_lb,x_ub,y_lb,y_ub,len_X,len_Y,Spec_wght)
!write(19,500)px,py!,px_tst,py_tst
end if
end if
end do
!$OMP END PARALLEL
!-----------------------------------------------------------------------------------------------
do i=1,len_Y
do j=1,len_X
    Spec(i,j)=(abs(M_cos(i,j)))**2.0+(abs(M_sin(i,j)))**2.0
end do
end do
!-----------------------------------------------------------------------------------------------
Spec_max=Spec(1,1)

do i=1,len_Y
do j=1,len_X
    if(Spec(i,j)>Spec_max)then
        Spec_max=Spec(i,j)
    end if
end do
end do
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------


!do i=1,len_Y
!Spec(i,:)=Spec(i,:)+Spec((len_Y-i+1),:)
   !write(18,1300) Spec_wght(i,:)
!end do

do i=1,len_Y
   write(15,1300) Spec(i,:)
   !write(18,1300) Spec_wght(i,:)
end do
!----------------------------------------------------------------------------------------------- 
do i=1,len_X

    Spec_grdX(i)=x_lb+i*mwidth_x-0.5*mwidth_x
    
end do
                                                                                       !ÏòÍø¸ñÊý×éX_sta,Y_sta¸³Öµ
do i=1,len_Y

    Spec_grdY(i)=y_lb+i*mwidth_y-0.5*mwidth_y

end do
!-----------------------------------------------------------------------------------------------
do i=1,len_X
   write(13,*) Spec_grdX(i)
end do

do i=1,len_Y
   write(14,*) Spec_grdY(i)
end do
!-----------------------------------------------------------------------------------------------
!*******************************
500 format(10(f30.15))
700 format(3(f30.15))
1100 format(2500(E25.15E3,1x))!format(2500(f30.15,1x))
1300 format(2500(f30.15,1x))
end program main

!*****************************************************************************
!                            Random grabing function
!*****************************************************************************
function random_grab(lbound,ubound)
implicit none
real*8::random_grab
real*8::lbound,ubound
real*8::length
real*8::x
length=ubound-lbound
call random_number(x)
random_grab=lbound+x*length
return
end function random_grab 

!*****************************************************************************
!                  Electric vector function of Laser fields
!*****************************************************************************
real(kind=8) function f(t)
use lsr_para
implicit none
real*8::t
if(t<=T_PD) then
    f=1.0
else if((t>T_PD).and.(t<(T_PD+T_TD))) then
    f=(cos(0.5*pi*(t-T_PD)/T_TD))**2.0
else 
    f=0.0
end if
return
end function f

real(kind=8) function E(t)
use lsr_para
implicit none
real*8::t
real*8,external::f
real*8::E0
E0=sqrt(L_I/I0)
E=E0*f(t)*sin(L_W*t+ceph)
return
    end function E
!*****************************************************************************
!                  Subroutine for grabing electron initial tunneling time in all slopes of the laser waveform
!*****************************************************************************
subroutine tnlt0_grab(pk_wdw,cs_wdw,tc_lb,tc_ub,tnlt0)
use lsr_para
implicit none
real*8::pk_wdw,cs_wdw,tnlt0
real*8::tt_rg,tt0,dt0
real*8::tc_lb,tc_ub,T_TC,nc
real*8::t_wdw
real*8,external::random_grab

T_TC=tc_ub-tc_lb

t_wdw=L_T/4-(pk_wdw+cs_wdw)/2
tt_rg=T_TC*4*t_wdw
tt0=random_grab(0.d0,tt_rg)
nc=int(tt0/(4*t_wdw))
dt0=tt0-nc*4*t_wdw
if((dt0>0).and.(dt0<t_wdw))then
    tnlt0=nc*L_T+dt0+cs_wdw/2
else if((dt0>t_wdw).and.(dt0<(2*t_wdw)))then
    tnlt0=nc*L_T+dt0+cs_wdw/2+pk_wdw
else if((dt0>(2*t_wdw)).and.(dt0<(3*t_wdw)))then
    tnlt0=nc*L_T+dt0+cs_wdw/2+pk_wdw+cs_wdw
else if((dt0>(3*t_wdw)).and.(dt0<(4*t_wdw)))then
    tnlt0=nc*L_T+dt0+cs_wdw/2+2*pk_wdw+cs_wdw
end if
end subroutine tnlt0_grab
!*****************************************************************************
!                  Position function of tunneling exit point
!*****************************************************************************
real(kind=8) function tnl_extpt(E,Ip)
use  MO_para
implicit none
real*8::Ipp,beta,re,E,Ip

Ipp=Ip+0.5*(apha_N-apha_I)*E**2
beta=1-sqrt(2*Ipp)/2
re=(Ipp+sqrt(Ipp**2-4*beta*abs(E)))/(2*abs(E))
if(E>0.0)then
    tnl_extpt=-re
else if(E<0.0)then
    tnl_extpt=re
end if
 
return
end function tnl_extpt

!***************************************************************************** 
!Weight function of initial transverse momentum originates from tunneling ionization
!*****************************************************************************
real(kind=8) function tnl_wght(Ip,E0,p_per)
implicit none
real*8::Ip,E0,p_per,k
real*8::W0,W1
real*8::pi=3.141592653589793d0

!pi=3.141592653589793d0
W0=exp((-2*(2*Ip)**1.5)/(3*E0))!((4*(2*Ip)**2)/abs(E0))*exp(-(2*k**3)/(3*E0))!((((2*Ip)**2)/abs(E0))**((2.d0/((2*Ip)**0.5))-1))*exp((-2*(2*Ip)**(3/2))/(3*abs(E0)))!((((2*Ip)**2)/abs(E0))**((2.d0/((2*Ip)**0.5))-1))*exp((-2*(2*abs(Ip))**(3/2))/(3*abs(E0)))
W1=exp((-p_per**2)*sqrt(2*Ip)/E0)!(k/(abs(E0)*pi))*exp(-k*(p_per**2)/E0)!(sqrt(2*Ip)/abs(E0))*exp(-sqrt(2*Ip)*p_per**2/abs(E0))!abs(p_per)*(((2*abs(Ip))**(0.5))/(abs(E0)*pi))*exp((-(p_per**2)*((2*abs(Ip))**0.5))/abs(E0))


!W1=abs(p_per)*(((2*abs(Ip))**(0.5))/(abs(E0)*pi))*exp((-(p_per**2)*((2*abs(Ip))**0.5))/abs(E0))!from old CTMC program! for test
tnl_wght=W0*W1

return
end function tnl_wght


!***************************************************************************** 
!                   Probability density of molecular orbital in momentum sopace
!*****************************************************************************
!real(kind=8) function MO_dnsty(px,py)! Molecular orbital in momentum space(1pi_g)
!use lsr_para
!use MO_para
!implicit none
!real*8::px,py,ppx,ppy,pp
!real*8::kn,Spz,p
!kn=sqrt(2*Ip)
!Spz=(1+kn*R0+(1/5)*(kn*R0)**2-(2/15)*(kn*R0)**3-(1/15)*(kn*R0)**4)*exp(-kn*R0)
!ppx=px*cos(thta)-py*sin(thta)
!ppy=px*sin(thta)+py*cos(thta)
!pp=sqrt(px**2+py**2)
!MO_dnsty=((2**5)*pp*(kn**3.5)*cos(atan2(ppy,ppx))*cos(ppx*R0/2)/(pi*sqrt(2.d0)*(pp**2+kn**2)**3)/sqrt(2*(1-Spz)))**2 

!return

!end function MO_dnsty
!*****************************************************************************
!             potential function for integrating classical action(H_ph)  
!*****************************************************************************
real(kind=8) function H_ph(px,py,rx,ry,n_Z1,n_Z2,n_R,apha_II) ! ,n_Z1,n_Z2,n_R,apha_II
use sim_para
use MO_para
implicit none
real*8::px,py,rx,ry,t,n_Z1,n_Z2,n_R,apha_II
real*8::V_r,grd_rx,grd_ry
real*8::rr1,rr2
real*8,external::partialt_u1x,partialt_u2x,partialt_l1x,partialt_l2x,omega_u,omega_l,partialt1,partialt2
t=0.d0
V_r=-n_Z1/sqrt((rx-0.5*n_R*cos(apha))**2+(ry-0.5*n_R*sin(apha))**2+ne1**2)-n_Z2/sqrt((rx+0.5*n_R*cos(apha))**2+(ry+0.5*n_R*sin(apha))**2+ne2**2)
grd_rx=partialt1(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
grd_ry=partialt2(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
H_ph=0.5*(px**2.0+py**2.0)+V_r-(rx*grd_rx+ry*grd_ry)

!end if
return
end function H_ph
!*****************************************************************************
!             Subroutine of ionization status judgement by electron energy
!*****************************************************************************
subroutine ion_jdg_egy(px,py,rx,ry,jdg)
use sim_para
use MO_para
implicit none
real*8::px,py,rx,ry
real*8::tot_egy
integer::jdg
!real*8,external::H_ph
tot_egy=0.5*(px**2.0+py**2.0)-n_Z/sqrt(rx**2+ry**2+ne1**2)!n_Z1/sqrt(rx**2+(ry-0.5*n_R)**2+ne1**2)-n_Z2/sqrt(rx**2+(ry+0.5*n_R)**2+ne2**2)
if(tot_egy>0)then
    jdg=1
else if(tot_egy<0)then
    jdg=0
end if
    end subroutine ion_jdg_egy

!*****************************************************************************
!             Subroutine of ionization status judgement by electron energy
!*****************************************************************************
real(kind=8) function elec_egy(px,py,rx,ry,n_Z1,n_Z2,n_R)
use sim_para
use MO_para
implicit none
real*8::px,py,rx,ry,n_Z1,n_Z2,n_R
elec_egy=0.5*(px**2.0+py**2.0)-n_Z1/sqrt(rx**2+(ry-0.5*n_R)**2+ne1**2)-n_Z2/sqrt(rx**2+(ry+0.5*n_R)**2+ne2**2)
return
    end function elec_egy
!*****************************************************************************
!             Subroutine for calculating cross product
!*****************************************************************************
subroutine crs_pdt(a1,a2,a3,b1,b2,b3,pdt_x,pdt_y,pdt_z)
implicit none
real*8::a1,a2,a3,b1,b2,b3,pdt_x,pdt_y,pdt_z
pdt_x=a2*b3-a3*b2
pdt_y=a3*b1-a1*b3
pdt_z=a1*b2-a2*b1
end subroutine crs_pdt

!*****************************************************************************
!             Subroutine for calculating asymptotic momentum
!*****************************************************************************
subroutine asympt_momx(qx,qy,qz,rx,ry,rz,asym_px,asym_py,asym_pz)
use sim_para
use MO_para
implicit none
real*8::qx,qy,qz,q_abs,rx,ry,rz,r_abs,asym_px,asym_py,asym_pz,p_abs
real*8::Lx,Ly,Lz,L_squr!angular momentum(three components and absolute value)
real*8::qLx,qLy,qLz!cross-production between momentum q and angularmomentum L
real*8::ax,ay,az!Runge-Lenz vector(three components)
real*8::Lax,Lay,Laz!cross-production between angular momentum L and vector A
real*8::Deno,coeff1,coeff2!two coefficients for calculating the final momentum
r_abs=(rx**2+ry**2+rz**2)**0.5
q_abs=(qx**2+qy**2+qz**2)**0.5
p_abs=(q_abs**2-(2*n_Z/r_abs))**0.5
call crs_pdt(rx,ry,rz,qx,qy,qz,Lx,Ly,Lz)
L_squr=Lx**2+Ly**2+Lz**2
call crs_pdt(qx,qy,qz,Lx,Ly,Lz,qLx,qLy,qLz)
ax=qLx-rx/r_abs
ay=qLy-ry/r_abs
az=qLz-rz/r_abs
call crs_pdt(Lx,Ly,Lz,ax,ay,az,Lax,Lay,Laz)
Deno=1+(p_abs**2)*L_squr
coeff1=(p_abs**2)/Deno
coeff2=p_abs/Deno
asym_px=coeff1*Lax-coeff2*ax
asym_py=coeff1*Lay-coeff2*ay
asym_pz=coeff1*Laz-coeff2*az
return
end subroutine asympt_momx
!*****************************************************************************
!             function for arsinh function
!*****************************************************************************
real(kind=8) function arsinh(x)
implicit none
real*8::x
arsinh=dlog(x+sqrt(x**2+1))
return
end function arsinh
!*****************************************************************************
!             Subroutine for calculating final Coulomb phase
!*****************************************************************************
subroutine fnl_Coulomb_phase(px,py,rx,ry,C_ph)
use sim_para
use MO_para
implicit none
real*8::px,py,rx,ry,C_ph
real*8::E,L,b,g,V_r
real*8,external::arsinh
V_r=-n_Z/sqrt(rx**2+ry**2)
E=0.5*(px**2+py**2)+V_r
L=rx*py-ry*px
b=1/(2*E)
g=sqrt(1+2*E*L**2)
C_ph=-n_Z*sqrt(b)*(dlog(g)+arsinh((rx*px+ry*py)/(g*sqrt(b))))
return
end subroutine fnl_Coulomb_phase
!************************************************************************************************
! sta_ana_intrfr_2D subroutine for doing the interference of all the trajectories on a 2D momentum space
!************************************************************************************************
subroutine sta_ana_intrfr_2D(px,py,wght,S_ph,n_x,n_y,mn_x,mn_y,x_lb,x_ub,y_lb,y_ub,len_X,len_Y,M_cos,M_sin)!(px,py,wght,S_ph,M_cos,M_sin)
implicit none
!use intrfr_para
real*8::px,py,x_lb,x_ub,y_lb,y_ub,wght,S_ph
integer::n_x,n_y,mn_x,mn_y
integer::len_X,len_Y
real*8::width_x,width_y,mwidth_x,mwidth_y
real*8::M_cos(1:len_Y,1:len_X),M_sin(1:len_Y,1:len_X)
integer::A_x,A_y !<array label>
integer::Int_X,Int_Y
integer::ci,cj 

width_x=(x_ub-x_lb)/n_x
width_y=(y_ub-y_lb)/n_y
mwidth_x=width_x/mn_x
mwidth_y=width_y/mn_y

if(((px<x_ub).and.(px>x_lb)).and.((py<y_ub).and.(py>y_lb))) then

    Int_X=int((px-x_lb)/width_x)
    Int_Y=int((py-y_lb)/width_y)
do ci=1,mn_x       
   do cj=1,mn_y 
    A_x=Int_X*mn_x+ci
    A_y=Int_Y*mn_y+cj
    
    M_cos(A_y,A_x)=M_cos(A_y,A_x)+sqrt(wght)*cos(S_ph)
    M_sin(A_y,A_x)=M_sin(A_y,A_x)+sqrt(wght)*sin(S_ph)
    end do
end do
end if
end subroutine sta_ana_intrfr_2D
!************************************************************************************************
! sta_ana_wght_2D subroutine for adding the weights of all the trajectories on a 2D momentum space
!************************************************************************************************
subroutine sta_ana_wght_2D(px,py,wght,n_x,n_y,mn_x,mn_y,x_lb,x_ub,y_lb,y_ub,len_X,len_Y,Spec)!(px,py,wght,S_ph,M_cos,M_sin)
implicit none
!use intrfr_para
real*8::px,py,x_lb,x_ub,y_lb,y_ub,wght
integer::n_x,n_y,mn_x,mn_y
integer::len_X,len_Y
real*8::width_x,width_y,mwidth_x,mwidth_y
real*8::Spec(1:len_Y,1:len_X)
integer::A_x,A_y !<array label>
integer::Int_X,Int_Y
integer::ci,cj 

width_x=(x_ub-x_lb)/n_x
width_y=(y_ub-y_lb)/n_y
mwidth_x=width_x/mn_x
mwidth_y=width_y/mn_y

if(((px<x_ub).and.(px>x_lb)).and.((py<y_ub).and.(py>y_lb))) then

    Int_X=int((px-x_lb)/width_x)
    Int_Y=int((py-y_lb)/width_y)
do ci=1,mn_x       
   do cj=1,mn_y 
    A_x=Int_X*mn_x+ci
    A_y=Int_Y*mn_y+cj
    
    Spec(A_y,A_x)=Spec(A_y,A_x)+wght
    
    end do
end do
end if
    end subroutine sta_ana_wght_2D


!************************************************************************************************

!*****************************************************************************
!   subroutine of Symplectic algorithm(2D case)(with laser field)
!*****************************************************************************
subroutine sym_algrthm(px,py,rx,ry,t,h,n_Z1,n_Z2,n_R,apha_II)
use MO_para
implicit none
real*8::a,b,c1,c2,c3,c4,d1,d2,d3,d4
real*8::px,py,rx,ry
real*8::px1,py1,rx1,ry1
real*8::px2,py2,rx2,ry2
real*8::px3,py3,rx3,ry3
real*8::t,t1,t2,t3,h
real*8::n_Z1,n_Z2,n_R,apha_II
real*8,external::partialt_u1x,partialt_u2x,partialt_l1x,partialt_l2x,partialt1,partialt2
real*8::rr1,rr2
rr1=sqrt(rx**2+(ry-0.5*n_R)**2)
rr2=sqrt(rx**2+(ry+0.5*n_R)**2)

a=1.d0/(2.d0-2.d0**(1.d0/3.d0))
b=1.d0-2.d0*a
c1=0.d0
c2=-a
c3=-b
c4=-a
d1=0.5d0*a
d2=0.5d0*(a+b)
d3=0.5d0*(a+b)
d4=0.5d0*a

if((rr1<R_M1).and.(rr2>R_M2))then

px1=px+c1*h*partialt1(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
rx1=rx+d1*h*px1
py1=py+c1*h*partialt2(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
ry1=ry+d1*h*py1

t1=t+h*d1

px2=px1+c2*h*partialt1(rx1,ry1,t1,n_Z1,n_Z2,n_R,apha_II)
rx2=rx1+d2*h*px2
py2=py1+c2*h*partialt2(rx1,ry1,t1,n_Z1,n_Z2,n_R,apha_II)
ry2=ry1+d2*h*py2
 
t2=t1+h*d2

px3=px2+c3*h*partialt1(rx2,ry2,t2,n_Z1,n_Z2,n_R,apha_II)
rx3=rx2+d3*h*px3
py3=py2+c3*h*partialt2(rx2,ry2,t2,n_Z1,n_Z2,n_R,apha_II)
ry3=ry2+d3*h*py3
 
t3=t2+h*d3

px=px3+c4*h*partialt1(rx3,ry3,t3,n_Z1,n_Z2,n_R,apha_II)
rx=rx3+d4*h*px
py=py3+c4*h*partialt2(rx3,ry3,t3,n_Z1,n_Z2,n_R,apha_II)
ry=ry3+d4*h*py

else if((rr2<R_M2).and.(rr1>R_M1))then

px1=px+c1*h*partialt1(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
rx1=rx+d1*h*px1
py1=py+c1*h*partialt2(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
ry1=ry+d1*h*py1

t1=t+h*d1

px2=px1+c2*h*partialt1(rx1,ry1,t1,n_Z1,n_Z2,n_R,apha_II)
rx2=rx1+d2*h*px2
py2=py1+c2*h*partialt2(rx1,ry1,t1,n_Z1,n_Z2,n_R,apha_II)
ry2=ry1+d2*h*py2
 
t2=t1+h*d2

px3=px2+c3*h*partialt1(rx2,ry2,t2,n_Z1,n_Z2,n_R,apha_II)
rx3=rx2+d3*h*px3
py3=py2+c3*h*partialt2(rx2,ry2,t2,n_Z1,n_Z2,n_R,apha_II)
ry3=ry2+d3*h*py3
 
t3=t2+h*d3

px=px3+c4*h*partialt1(rx3,ry3,t3,n_Z1,n_Z2,n_R,apha_II)
rx=rx3+d4*h*px
py=py3+c4*h*partialt2(rx3,ry3,t3,n_Z1,n_Z2,n_R,apha_II)
ry=ry3+d4*h*py

else
    
px1=px+c1*h*partialt1(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
rx1=rx+d1*h*px1
py1=py+c1*h*partialt2(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
ry1=ry+d1*h*py1

t1=t+h*d1

px2=px1+c2*h*partialt1(rx1,ry1,t1,n_Z1,n_Z2,n_R,apha_II)
rx2=rx1+d2*h*px2
py2=py1+c2*h*partialt2(rx1,ry1,t1,n_Z1,n_Z2,n_R,apha_II)
ry2=ry1+d2*h*py2
 
t2=t1+h*d2

px3=px2+c3*h*partialt1(rx2,ry2,t2,n_Z1,n_Z2,n_R,apha_II)
rx3=rx2+d3*h*px3
py3=py2+c3*h*partialt2(rx2,ry2,t2,n_Z1,n_Z2,n_R,apha_II)
ry3=ry2+d3*h*py3
 
t3=t2+h*d3

px=px3+c4*h*partialt1(rx3,ry3,t3,n_Z1,n_Z2,n_R,apha_II)
rx=rx3+d4*h*px
py=py3+c4*h*partialt2(rx3,ry3,t3,n_Z1,n_Z2,n_R,apha_II)
ry=ry3+d4*h*py
end if
return
end subroutine sym_algrthm
!*****************************************************************************
!             Partial derivative function of symplectic algorithm(with laser field)
!*****************************************************************************
real(kind=8) function partialt1(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
use MO_para
implicit none
real*8::rx,ry,t
real*8::n_Z1,n_Z2,n_R,apha_II
real*8,external::E 

partialt1=n_Z1*(rx**2+(ry-0.5*n_R)**2)**(-1.5)*rx+n_Z2*(rx**2+(ry+0.5*n_R)**2)**(-1.5)*rx+E(t)-apha_II*E(t)*(rx*(-3)*(rx**2+ry**2)**(-2.5)*rx+(rx**2+ry**2)**(-1.5))
return
end function partialt1

real(kind=8) function partialt2(rx,ry,t,n_Z1,n_Z2,n_R,apha_II)
use MO_para 
implicit none
real*8::rx,ry,t
real*8::n_Z1,n_Z2,n_R,apha_II
real*8,external::E

partialt2=n_Z1*(rx**2+(ry-0.5*n_R)**2)**(-1.5)*(ry-0.5*n_R)+n_Z2*(rx**2+(ry+0.5*n_R)**2)**(-1.5)*(ry+0.5*n_R)-apha_II*E(t)*rx*(-3)*(rx**2+ry**2)**(-2.5)*ry
return
    end function partialt2
!************************************************************************************************

!*****************************************************************************
!   subroutine of step adaptive Symplectic algorithm(2D case)(with laser field)
!*****************************************************************************
subroutine sym_adp_algrthm(px,py,rx,ry,t,h,n_Z1,n_Z2,n_R,apha_II)!sym_algrthm(px,py,rx,ry,t,h,n_Z1,n_Z2,n_R,apha_II)
use Sym_ctrl
implicit none
real*8::px,py,rx,ry,t,h
real*8::n_Z1,n_Z2,n_R,apha_II
real*8::rr_crt

rr_crt=(rx**2+ry**2)**0.5
    !step-adaptive-Symplectic algrithm begin*** 
    if(rr_crt<r_crt1)then
    call sym_algrthm(px,py,rx,ry,t,h1,n_Z1,n_Z2,n_R,apha_II)
           h=h1
    else if((rr_crt>r_crt1).and.(rr_crt<r_crt2))then
    call sym_algrthm(px,py,rx,ry,t,h2,n_Z1,n_Z2,n_R,apha_II)
           h=h2
    else if(rr_crt>r_crt2)then
        call sym_algrthm(px,py,rx,ry,t,h3,n_Z1,n_Z2,n_R,apha_II)
           h=h3
    end if
    !step-adaptive-Symplectic algrithm end*** 
return
end subroutine sym_adp_algrthm

subroutine sta_ana_1D_wghtx(A_exp,W_exp,n_v,mn_v,v_lb,v_ub,len_spec,spec)
implicit none
integer::n_v 
integer::mn_v 
integer::len_spec
integer::i,ci,cj,V_i 
real*8::width 
real*8::mwidth 
real*8::A_exp,W_exp
real*8::v_lb,v_ub 
real*8::spec(1:len_spec)

 width=(v_ub-v_lb)/n_v   
 mwidth=width/mn_v
 
   if(((A_exp-v_lb)>=0).and.((A_exp-v_ub)<=0)) then
       V_i=(int((A_exp-v_lb)/width))*mn_v
       
        do i=1,mn_v                                                   
        spec(V_i+i)=spec(V_i+i)+W_exp
        end do
   end if

return
    end subroutine sta_ana_1D_wghtx
!---------------------------------------------------------------------------------------------------------------
!***************************Functions for GSZ potential (include multielectron screening effect)******************

!Deleted!! Go to find test2 project!!


subroutine sta_ana_1D_wght_lbl_dis(lbl_exp,len_dis,lbl_dis)
implicit none
integer::lbl_exp,len_dis,V_i
integer::lbl_dis(1:len_dis)  
lbl_dis(int(lbl_exp))=lbl_dis(int(lbl_exp))+1  
return
end subroutine sta_ana_1D_wght_lbl_dis
!---------------------------------------------------------------------------------------------------------------
