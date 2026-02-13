
    !  SCTS_spec_ND_plot_fan_area_analysis.f90 
!
!  FUNCTIONS:
!  SCTS_spec_ND_plot_fan_area_analysis - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: SCTS_spec_ND_plot_fan_area_analysis
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************


    program SCTS_spec_ND_plot

    implicit none
real*8,parameter::pi=3.14159265358979
integer,parameter::n_x=1000!<resolution control:bin size on x direction><previous value:1000>
integer,parameter::n_y=1000!<resolution control:bin size on y direction><previous value:1000>
integer,parameter::mn_x=1!<detail resolution control:bin size on x direction>
integer,parameter::mn_y=1!<detail resolution control:bin size on y direction>
integer,parameter::len_X=n_x*mn_x!<momemtum space matrix:length of the matrix Spec on x direction>
integer,parameter::len_Y=n_y*mn_y!<momemtum space matrix:length of the matrix Spec on y direction>
real*8,parameter::px_intgrl_lb=-1.5
real*8,parameter::px_intgrl_ub=1.5
real*8,parameter::py_intgrl_lb=-1.5
real*8,parameter::py_intgrl_ub=1.5
real*8,parameter::bkgrd_fltr=0.d0!1.0d-8!1.0d-8!0.d0!1.0d-20
real*8,parameter::infi_fltr=1.d0
integer,parameter::diff_r_spec_n=50!100!30
integer,parameter::diff_r_spec_n2=50!100!30
integer,parameter::diff_r_spec_mn=1
real*8,parameter:: diff_r_spec_lb=0.0
real*8,parameter:: diff_r_spec_ub=1.0
real*8,parameter::fan_angl_dgre=30!30
real*8,parameter::fan_angl=fan_angl_dgre/180*pi
integer,parameter::lgth_diff_r=diff_r_spec_n*diff_r_spec_mn
integer,parameter::lgth_diff_r2=diff_r_spec_n2*diff_r_spec_mn
real*8,parameter::width=(diff_r_spec_ub-diff_r_spec_lb)/diff_r_spec_n
real*8,parameter::mwidth=width/diff_r_spec_mn
real*8,parameter::width2=(diff_r_spec_ub-diff_r_spec_lb)/diff_r_spec_n2
real*8,parameter::mwidth2=width2/diff_r_spec_mn
    integer::i,ii,j,jj,k
    real*8::sum1,sum2
    real*8::Spec1(1:len_Y,1:len_X),Spec2(1:len_Y,1:len_X),ND_Spec(1:len_Y,1:len_X),ND_Spec_in_fan(1:len_Y,1:len_X),ND_Spec_out_fan(1:len_Y,1:len_X)
    real*8::Spec_gridX(1:len_X),Spec_gridY(1:len_Y)
    real*8::px,py
    real*8::diff_r_vlu,diff_r_spec,diff_r_vlu_max,diff_r_vlu_intgl
    real*8::diff_r_dis(1:lgth_diff_r),diff_r_dis_in_fan(1:lgth_diff_r),diff_r_dis_out_fan(1:lgth_diff_r)
    real*8::diff_in_fan_vlu_max,diff_out_fan_vlu_max
    real*8::theta
open(unit=10,file="Spec1.txt")
open(unit=11,file="Spec2.txt")
open(unit=12,file="ND_spec.txt")
open(unit=13,file="Spec1_ND.txt")
open(unit=14,file="Spec2_ND.txt")
open(unit=15,file="sum_test.txt")
open(unit=16,file="Spec_gridX.txt")
open(unit=17,file="Spec_gridY.txt")
open(unit=18,file="Sum_test.txt")
open(unit=19,file="diff_r_dis_test.txt")
open(unit=20,file="ND_spec_in_fan.txt")
open(unit=21,file="ND_spec_out_fan.txt")
open(unit=22,file="diff_r_dis_test_in_fan.txt")
open(unit=23,file="diff_r_dis_test_out_fan.txt")
do ii=1,len_Y
   do jj=1,len_Y
       Spec1(ii,jj)=0.d0
       Spec2(ii,jj)=0.d0
   end do
end do
sum1=0.d0
sum2=0.d0
do i=1,len_Y
   read(10,1300) Spec1(i,:)
   read(11,1300) Spec2(i,:)
end do
do i=1,len_X
    read(16,*) Spec_gridX(i)
end do

do i=1,len_Y
    read(17,*) Spec_gridY(i)
end do

do i=1,len_Y
   do j=1,len_X
       px=Spec_gridX(j)
       py=Spec_gridY(i)
       if((px<px_intgrl_ub).and.(px>px_intgrl_lb).and.(py<py_intgrl_ub).and.(py>py_intgrl_lb))then
       if(Spec1(i,j)>infi_fltr)then
           Spec1(i,j)=0.d0
       end if
       if(Spec2(i,j)>infi_fltr)then
       Spec2(i,j)=0.d0
       end if
       sum1=sum1+Spec1(i,j)
       sum2=sum2+Spec2(i,j)
       end if
   end do
end do
!write(18,*)sum1,sum2
do i=1,len_Y
   do j=1,len_X
       Spec1(i,j)=Spec1(i,j)/sum1
       Spec2(i,j)=Spec2(i,j)/sum2
   end do
end do

do i=1,len_Y
   do j=1,len_X
   if((Spec1(i,j)>bkgrd_fltr).and.(Spec2(i,j)>bkgrd_fltr))then
ND_Spec(i,j)=(Spec2(i,j)-Spec1(i,j))/(Spec2(i,j)+Spec1(i,j))
   else 
ND_Spec(i,j)=0.d0
end if
   end do
end do

do i=1,len_Y
   write(12,1300) ND_Spec(i,:)
   write(13,1300) Spec1(i,:)
   write(14,1300) Spec2(i,:)
end do

do i=1,len_Y
do j=1,len_X
    diff_r_vlu=sqrt(Spec_gridY(i)**2+Spec_gridX(j)**2)
    call sta_ana_1D_wghtx(diff_r_vlu,ND_Spec(i,j),diff_r_spec_n,diff_r_spec_mn,diff_r_spec_lb,diff_r_spec_ub,lgth_diff_r,diff_r_dis)
end do
end do

diff_r_vlu_max=diff_r_dis(1)

do i=1,lgth_diff_r
    if(diff_r_dis(i)>=diff_r_vlu_max)then
        diff_r_vlu_max=diff_r_dis(i)
    end if
end do

diff_r_vlu_intgl=0.d0

do i=1,lgth_diff_r
   diff_r_vlu_intgl=diff_r_vlu_intgl+diff_r_dis(i)
end do

do i=1,lgth_diff_r
    diff_r_spec=diff_r_spec_lb+i*mwidth+0.5*mwidth
    write(19,*)diff_r_spec,diff_r_dis(i),(diff_r_dis(i)/diff_r_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_intgl)
end do

write(15,*)sum1,sum2

do i=1,len_Y
   do j=1,len_X
theta=atan2(Spec_gridY(i),Spec_gridX(j))

if(((theta>fan_angl).and.(theta<(pi-fan_angl))).or.((theta>(fan_angl-pi)).and.(theta<(-fan_angl))))then
            ND_Spec_in_fan(i,j)=ND_Spec(i,j)
            ND_Spec_out_fan(i,j)=0.d0
        else
            ND_Spec_in_fan(i,j)=0.d0!pZ_sta_in_fan(i,j)=0.d0
            ND_Spec_out_fan(i,j)=ND_Spec(i,j)
        end if
   end do
end do

do i=1,len_Y
do j=1,len_X
    diff_r_vlu=sqrt(Spec_gridY(i)**2+Spec_gridX(j)**2)
    call sta_ana_1D_wghtx(diff_r_vlu,ND_Spec_in_fan(i,j),diff_r_spec_n,diff_r_spec_mn,diff_r_spec_lb,diff_r_spec_ub,lgth_diff_r,diff_r_dis_in_fan)
    call sta_ana_1D_wghtx(diff_r_vlu,ND_Spec_out_fan(i,j),diff_r_spec_n2,diff_r_spec_mn,diff_r_spec_lb,diff_r_spec_ub,lgth_diff_r,diff_r_dis_out_fan)
end do
end do

do i=1,len_Y
   write(20,1300) ND_Spec_in_fan(i,:)
   write(21,1300) ND_Spec_out_fan(i,:)
end do

diff_in_fan_vlu_max=diff_r_dis_in_fan(1)!diff_out_fan_vlu_max
diff_out_fan_vlu_max=diff_r_dis_out_fan(1)

do i=1,lgth_diff_r
    if(diff_r_dis_in_fan(i)>diff_in_fan_vlu_max)then
        diff_in_fan_vlu_max=diff_r_dis_in_fan(i)
    end if
    
    if(diff_r_dis_out_fan(i)>diff_out_fan_vlu_max)then
        diff_out_fan_vlu_max=diff_r_dis_out_fan(i)
    end if
end do
do i=1,lgth_diff_r
    diff_r_spec=diff_r_spec_lb+i*mwidth+0.5*mwidth
    write(22,*)diff_r_spec,diff_r_dis_in_fan(i),(diff_r_dis_in_fan(i)/diff_in_fan_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_intgl)
    ! write(23,*)diff_r_spec,diff_r_dis_out_fan(i),(diff_r_dis_out_fan(i)/diff_out_fan_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_intgl)
end do

do i=1,lgth_diff_r2
    diff_r_spec=diff_r_spec_lb+i*mwidth2+0.5*mwidth2
   ! write(22,*)diff_r_spec,diff_r_dis_in_fan(i),(diff_r_dis_in_fan(i)/diff_in_fan_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_intgl)
     write(23,*)diff_r_spec,diff_r_dis_out_fan(i),(diff_r_dis_out_fan(i)/diff_out_fan_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_max)!,(diff_r_dis(i)/diff_r_vlu_intgl)
end do
500 format(3(f30.15))
700 format(6(f30.15))
900 format(14(f35.20))
1100 format(2500(f25.15,1x))
1300 format(2500(f30.15,1x))

end program SCTS_spec_ND_plot

!-----------------------------------------------------------------------------------------------------------------

subroutine sta_ana_1D_wghtx(A_exp,W_exp,n_v,mn_v,v_lb,v_ub,len_spec,spec)
implicit none
integer::n_v !<调整生成图像的物理分辨率n_v><统计结果数组的维数>
integer::mn_v !<调整生成图像的细分物理分辨率mn_v>
integer::len_spec
integer::i,ci,cj,V_i !<循环变量和数组标号>
real*8::width !<图像像素点的宽度>
real*8::mwidth !<图像细分像素点的宽度>
real*8::A_exp,W_exp
real*8::v_lb,v_ub !<调整图像绘制的区域和范围>
real*8::spec(1:len_spec)

 width=(v_ub-v_lb)/n_v
 mwidth=width/mn_v
 
   if(((A_exp-v_lb)>=0).and.((A_exp-v_ub)<=0)) then
       V_i=(int((A_exp-v_lb)/width))*mn_v
       
        do i=1,mn_v                                                   !统计计算实验数据分布，数组C_sta做记录
        spec(V_i+i)=spec(V_i+i)+W_exp
        end do
   end if

return
end subroutine sta_ana_1D_wghtx
!--------------------------------------------------------------------------------------------------------------- 

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
