module parameters
    integer,parameter::stringlength=100
	integer,parameter::tmplength=1000
    integer,parameter::stringnumber=200000
	integer,parameter::satsum=150! 1-40--GPS  41-100--BDS  101-150--GAL      
    real*8,parameter::am=6378137.0d0
    real*8,parameter::e2=0.00669437999013d0
    real*8,parameter::pi=3.1415926535898d0
    real*8,parameter::vc=2.99792458d8
    real*8,parameter::miu=3.986005d14
    real*8,parameter::omegae=7.2921151467d-5
    integer,parameter::sitenumber=100
    integer,parameter::systemnum=50
    integer,parameter::typenumber=50
    integer,parameter::epochguess=5000
    integer,parameter::pathnum=100
    integer,parameter::obsguess=100
    integer,parameter::numsatmax=200
contains
function get_lines(a,b)
    implicit none
    integer::get_lines,a,b,flag
    if(mod(a,b)==0)then
        flag=b
    else
        flag=0
    endif
    get_lines=(a-mod(a,b)-flag)/b+1
end function

function get_satpos(a)
	implicit none
	character(len=3)::a
	integer::get_satpos,b
	read(a(2:3),'(I2)')get_satpos
	if(a(1:1)=="G")then
		b=0
	else if(a(1:1)=="C")then
		b=40
	else if(a(1:1)=="E")then
		b=100
	else 
		b=-get_satpos
	endif
	get_satpos=get_satpos+b
end function

function get_satprn(a)
	implicit none
	integer::a
	character(len=3)::get_satprn,b
	if(a<10.and.a>0)then
		write(b,'(a1,i1)')"0",a
		get_satprn="G"//b
	else if(a<=40)then
		write(b,'(i2)')a
		get_satprn="G"//b
	else if(a<50)then
		write(b,'(a1,i1)')"0",a-40
		get_satprn="C"//b
	else if(a<=100)then
		write(b,'(i2)')a-40
		get_satprn="C"//b
	else if(a<110)then
		write(b,'(a1,i1)')"0",a-100
		get_satprn="E"//b
	else if(a<=150)then
		write(b,'(i2)')a-100
		get_satprn="E"//b
	endif
end function
end module
