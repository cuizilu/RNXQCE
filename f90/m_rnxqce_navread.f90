module m_rnxqce_navread
use timetrans
use parameters
implicit none
type delta_utc
    real*8::a0
    real*8::a1
    integer*8::t
    integer*8::w
endtype

type timecorrection
    character*4::name
    real*8::a0
    real*8::a1
    integer::t
    integer::w
endtype

type ioncorrection
    character*4::name
    real*8::a1
    real*8::a2
    real*8::a3
    real*8::a4
    character*1::mark
    integer::sv
endtype

type T_rnxnhd
    real*8::version
    real*8::ion_alpha(4)
    real*8::ion_beta(4)
    type(delta_utc)::dutc
    type(timecorrection),allocatable::corr(:)
    type(ioncorrection),allocatable::ion(:)
    integer::corrnum
    integer::ionnum
    integer::leapsec
    integer::headlinerow
endtype

type T_rnxndt!G、E、J三个系统的参数混合在前，其余系统的依次后排
    character*3::name
    integer*4::prn
    type(calendarTime)::cal
    type(gpsTime)::toc
    real*8::clkbias
    real*8::clkdrift
    real*8::clkdriftrate
    real*8::iode
    real*8::iodnav
    real*8::crs
    real*8::deltan
    real*8::m0
    real*8::cuc
    real*8::e
    real*8::cus
    real*8::sqrta
    real*8::toe
    real*8::cic
    real*8::omega
    real*8::cis
    real*8::i0
    real*8::crc
    real*8::omega1
    real*8::omegadot
    real*8::idot
    real*8::codesonl2channel
    real*8::datasources
    real*8::weeks
    real*8::l2pdataflag
    real*8::svaccuracy
    real*8::sisa
    real*8::svhealth
    real*8::tgd
    real*8::bgda
    real*8::iodc
    real*8::bgdb
    real*8::transtimeofmsg
    real*8::fitinterval
    real*8::spare1
    real*8::spare2
    real*8::spare3
    real*8::spare4
    real*8::posx!R系统
    real*8::posy
    real*8::posz
    real*8::xdot
    real*8::ydot
    real*8::zdot
    real*8::xacceleration
    real*8::yacceleration
    real*8::zacceleration
    real*8::health
    real*8::frequency
    real*8::ageoper
    real*8::aode!C系统
    real*8::sath1
    real*8::tgd1
    real*8::tgd2
    real*8::aodc
    integer::num
    real*8::accuracy!S系统
    real*8::iodn
    real*8::iodec!I系统
endtype

contains
subroutine read_nav_head(path,hdr)
implicit none
	type(T_rnxnhd)::hdr
	integer::i,nh,n,pp,flag,corrnum,ionnum,corrflag,ionflag,ios,jump
	character(len=tmplength)::tmp,tmp_ion(1000),tmp_sys(1000)
	character(len=stringlength)::path
	real*8 ww,qq
	nh=0
	ionnum=0
	corrnum=0
	flag=0
	open(1234,file=path,blank="null",status='old',iostat=ios)
    if(ios>0)then
        print*,"undifined file, please check your input"
    else
        do while(.true.)
            read(1234,'(a)',iostat=jump)tmp
            if(jump/=0.or.tmp(61:80)=='END OF HEADER')then
			    exit
			endif
			if(tmp(61:73)=="RINEX VERSION")then
			    if(tmp(6:6)/='.'.and.tmp(7:7)/='.'.and.tmp(8:8)/='.'.and.tmp(9:9)/='.')then
			    	flag=1
			    	read(tmp(6:9),'(i4)')pp
			    else
			    	flag=2
			    	read(tmp(6:9),'(f4.2)')ww
			    endif
			    if(flag==1)then
			        hdr%version=pp
			    elseif(flag==2)then
			        hdr%version=ww
			    endif
			endif
			if(tmp(61:69)=="ION ALPHA")then
			    read(tmp(4:14),'(f11.4)')hdr%ion_alpha(1)
			    read(tmp(16:26),'(f11.4)')hdr%ion_alpha(2)
			    read(tmp(28:38),'(f11.4)')hdr%ion_alpha(3)
			    read(tmp(40:50),'(f11.4)')hdr%ion_alpha(4)
			endif
			if(tmp(61:68)=="ION BETA")then
			    read(tmp(4:14),'(f11.4)')hdr%ion_beta(1)
			    read(tmp(16:26),'(f11.4)')hdr%ion_beta(2)
			    read(tmp(28:38),'(f11.4)')hdr%ion_beta(3)
			    read(tmp(40:50),'(f11.4)')hdr%ion_beta(4)
			endif
			if(tmp(61:76)=="TIME SYSTEM CORR")then
			    corrnum=corrnum+1
				tmp_sys(corrnum)=tmp
			endif
			if(tmp(61:76)=="IONOSPHERIC CORR")then
			    ionnum=ionnum+1
				tmp_ion(ionnum)=tmp
			endif
			if(tmp(61:72)=="LEAP SECONDS")then
			    read(tmp(5:10),'(i6)')hdr%leapsec
			endif
			if(tmp(61:69)=="DELTA-UTC")then
			    read(tmp(4:22),'(f19.12)')hdr%dutc%a0
			    read(tmp(23:41),'(f19.12)')hdr%dutc%a1
			    read(tmp(44:50),'(i7)')hdr%dutc%t
			    read(tmp(55:59),'(i5)')hdr%dutc%w
			endif
		enddo
	endif
	hdr%corrnum=corrnum
	hdr%ionnum=ionnum
	if(corrnum/=0)then
    	allocate(hdr%corr(corrnum))
    endif
    if(ionnum/=0)then
        allocate(hdr%ion(ionnum))
    endif
	do i=1,corrnum
		read(tmp_sys(i),'(a4,1x,d17.10,d16.9,i7,i5,1x,a5,1x,i2,1x)')hdr%corr(i)%name,hdr%corr(i)%a0&
		,hdr%corr(i)%a1,hdr%corr(i)%t,hdr%corr(i)%w
	enddo
	do i=1,ionnum
		read(tmp_ion(i),'(a4,1x,4d12.4,1x,a1,1x,i2)')hdr%ion(i)%name,hdr%ion(i)%a1,hdr%ion(i)%a2,hdr%ion(i)%a3,hdr%ion(i)%a4&
		,hdr%ion(i)%mark,hdr%ion(i)%sv
	enddo
end subroutine

subroutine getnavsum(navsum,hdr)
implicit none 
	type(T_rnxnhd)hdr
	integer::navsum(satsum),jump,ix,prn
	character(len=tmplength)::tmp
	character(len=3)::name
	navsum(:)=0
	do while(.true.)
		read(1234,'(a)',iostat=jump)tmp
        if(jump/=0)then
            exit
        endif
		if(hdr%version<3)then
			if(tmp(1:3)/="")then
				read(tmp(1:3),'(i2)')prn
	            if(prn>=10)then
	                write(name,'(a,i2)')"G",prn
	            else
	                write(name,'(a,i1)')"G0",prn
	            endif
				navsum(get_satpos(name))=navsum(get_satpos(name))+1
			endif
		else
			if(tmp(1:3)/="")then
				name=tmp(1:3)
				navsum(get_satpos(name))=navsum(get_satpos(name))+1
			endif
		endif
	enddo
	close(1234)
end subroutine 

function get_max(a)
implicit none
	integer::i,a(satsum),get_max
	get_max=0
	do i=1,satsum
		if(a(i)>get_max)then
			get_max=a(i)
		endif
	enddo
end function

subroutine read_nav_data(path,rnxnhd,rnxndt,navsum)
implicit none
	type(T_rnxnhd)rnxnhd
	type(T_rnxndt),allocatable::rnxndt(:,:)
	type(calendarTime)::cal
	integer::navsum(satsum),ios,jump,navnum(satsum),ix,prn
	character(len=tmplength)::tmp
	character(len=stringlength)::path
	character(len=3)::name
	call getnavsum(navsum,rnxnhd)
	navnum(:)=0
	allocate(rnxndt(satsum,get_max(navsum)))
	open(1234,file=path,blank="null",status='old',iostat=ios)
    if(ios>0)then
        print*,"undifined file, please check your input"
    else
        do while(.true.)
            read(1234,'(a)',iostat=jump)tmp
            if(jump/=0.or.tmp(61:80)=='END OF HEADER')then
                exit
            endif
		enddo
	endif
    do while(.true.)
        read(1234,'(a)',iostat=jump)tmp
        if(jump/=0)then
            exit
        endif
		if(rnxnhd%version<3)then
            read(tmp(1:3),'(i2)')prn
            if(prn>=10)then
                write(name,'(a,i2)')"G",prn
            else
                write(name,'(a,i1)')"G0",prn
            endif
			ix=get_satpos(name)
			navnum(ix)=navnum(ix)+1
			read(tmp,'(i2,1x,i2.2,4(1x,i2),f5.1,3d19.12)')rnxndt(ix,navnum(ix))%prn,cal%year,cal%month,cal%day,cal%hour&
			,cal%minute,cal%second,rnxndt(ix,navnum(ix))%clkbias,rnxndt(ix,navnum(ix))%clkdrift,rnxndt(ix,navnum(ix))%clkdriftrate
	        if(rnxndt(ix,navnum(ix))%prn>=10)then
	            write(rnxndt(ix,navnum(ix))%name,'(a,i2)')"G",rnxndt(ix,navnum(ix))%prn
	        else
	            write(rnxndt(ix,navnum(ix))%name,'(a,i1)')"G0",rnxndt(ix,navnum(ix))%prn
	        endif
	        if(cal%year<30)then
	        	cal%year=cal%year+2000
	        else
	            cal%year=cal%year+1900
	        endif 
	        call caltogps(cal,rnxndt(ix,navnum(ix))%toc)
			read(1234,'(a)',iostat=jump)tmp
	        read(tmp,'(3x,4d19.12)')rnxndt(ix,navnum(ix))%iode,rnxndt(ix,navnum(ix))%crs,rnxndt(ix,navnum(ix))%deltan,&
			rnxndt(ix,navnum(ix))%m0
			read(1234,'(a)',iostat=jump)tmp
	        read(tmp,'(3x,4d19.12)')rnxndt(ix,navnum(ix))%cuc,rnxndt(ix,navnum(ix))%e,rnxndt(ix,navnum(ix))%cus,&
			rnxndt(ix,navnum(ix))%sqrta
			read(1234,'(a)',iostat=jump)tmp
	        read(tmp,'(3x,4d19.12)')rnxndt(ix,navnum(ix))%toe,rnxndt(ix,navnum(ix))%cic,rnxndt(ix,navnum(ix))%omega,&
			rnxndt(ix,navnum(ix))%cis
			read(1234,'(a)',iostat=jump)tmp
	        read(tmp,'(3x,4d19.12)')rnxndt(ix,navnum(ix))%i0,rnxndt(ix,navnum(ix))%crc,rnxndt(ix,navnum(ix))%omega1&
			,rnxndt(ix,navnum(ix))%omegadot
			read(1234,'(a)',iostat=jump)tmp
	        read(tmp,'(3x,4d19.12)')rnxndt(ix,navnum(ix))%idot,rnxndt(ix,navnum(ix))%codesonl2channel,rnxndt(ix,navnum(ix))%weeks&
			,rnxndt(ix,navnum(ix))%l2pdataflag
			read(1234,'(a)',iostat=jump)tmp
	        read(tmp,'(3x,4d19.12)')rnxndt(ix,navnum(ix))%svaccuracy,rnxndt(ix,navnum(ix))%svhealth,rnxndt(ix,navnum(ix))%tgd&
			,rnxndt(ix,navnum(ix))%iodc
			read(1234,'(a)',iostat=jump)tmp
	        read(tmp,'(3x,4d19.12)')rnxndt(ix,navnum(ix))%transtimeofmsg,rnxndt(ix,navnum(ix))%fitinterval,&
			rnxndt(ix,navnum(ix))%spare1,rnxndt(ix,navnum(ix))%spare2
	        rnxndt(ix,navnum(ix))%cal=cal
		else
			ix=get_satpos(tmp(1:3))!index('GCE',tmp(1:1))    
			navnum(ix)=navnum(ix)+1
			if(ix/=0)then
            	read(tmp,'(a3,1x,i4,4(1x,i2.2),1x,f2.0,3d19.12)')rnxndt(ix,navnum(ix))%name,cal%year,cal%month,cal%day,cal%hour&
				,cal%minute,cal%second,rnxndt(ix,navnum(ix))%clkbias,rnxndt(ix,navnum(ix))%clkdrift,rnxndt(ix,navnum(ix))%clkdriftrate 
            	rnxndt(ix,navnum(ix))%cal=cal
            	call caltogps(cal,rnxndt(ix,navnum(ix))%toc)
            	if(rnxndt(ix,navnum(ix))%name(1:1)=="G")then
					read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%iode,rnxndt(ix,navnum(ix))%crs,rnxndt(ix,navnum(ix))%deltan&
					,rnxndt(ix,navnum(ix))%m0
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%cuc,rnxndt(ix,navnum(ix))%e,rnxndt(ix,navnum(ix))%cus&
					,rnxndt(ix,navnum(ix))%sqrta
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%toe,rnxndt(ix,navnum(ix))%cic,rnxndt(ix,navnum(ix))%omega&
					,rnxndt(ix,navnum(ix))%cis
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%i0,rnxndt(ix,navnum(ix))%crc,rnxndt(ix,navnum(ix))%omega1&
					,rnxndt(ix,navnum(ix))%omegadot
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%idot,rnxndt(ix,navnum(ix))%codesonl2channel&
					,rnxndt(ix,navnum(ix))%weeks,rnxndt(ix,navnum(ix))%l2pdataflag
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%svaccuracy,rnxndt(ix,navnum(ix))%svhealth&
					,rnxndt(ix,navnum(ix))%tgd,rnxndt(ix,navnum(ix))%iodc
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%transtimeofmsg,rnxndt(ix,navnum(ix))%fitinterval&
					,rnxndt(ix,navnum(ix))%spare1,rnxndt(ix,navnum(ix))%spare2
            	elseif(rnxndt(ix,navnum(ix))%name(1:1)=="E")then
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%iodnav,rnxndt(ix,navnum(ix))%crs,rnxndt(ix,navnum(ix))%deltan&
					,rnxndt(ix,navnum(ix))%m0
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%cuc,rnxndt(ix,navnum(ix))%e,rnxndt(ix,navnum(ix))%cus&
					,rnxndt(ix,navnum(ix))%sqrta
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%toe,rnxndt(ix,navnum(ix))%cic,rnxndt(ix,navnum(ix))%omega&
					,rnxndt(ix,navnum(ix))%cis
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%i0,rnxndt(ix,navnum(ix))%crc,rnxndt(ix,navnum(ix))%omega1&
					,rnxndt(ix,navnum(ix))%omegadot
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%idot,rnxndt(ix,navnum(ix))%datasources,rnxndt(ix,navnum(ix))%weeks&
					,rnxndt(ix,navnum(ix))%spare1
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%sisa,rnxndt(ix,navnum(ix))%svhealth,rnxndt(ix,navnum(ix))%bgda&
					,rnxndt(ix,navnum(ix))%bgdb
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%transtimeofmsg,rnxndt(ix,navnum(ix))%spare2&
					,rnxndt(ix,navnum(ix))%spare3,rnxndt(ix,navnum(ix))%spare4
            	elseif(rnxndt(ix,navnum(ix))%name(1:1)=="R")then
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%posx,rnxndt(ix,navnum(ix))%xdot,rnxndt(ix,navnum(ix))%xacceleration&
					,rnxndt(ix,navnum(ix))%health
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%posy,rnxndt(ix,navnum(ix))%ydot,rnxndt(ix,navnum(ix))%yacceleration&
					,rnxndt(ix,navnum(ix))%frequency
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%posz,rnxndt(ix,navnum(ix))%zdot,rnxndt(ix,navnum(ix))%zacceleration&
					,rnxndt(ix,navnum(ix))%ageoper
            	elseif(rnxndt(ix,navnum(ix))%name(1:1)=="J")then
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%iode,rnxndt(ix,navnum(ix))%crs,rnxndt(ix,navnum(ix))%deltan&
					,rnxndt(ix,navnum(ix))%m0
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%cuc,rnxndt(ix,navnum(ix))%e,rnxndt(ix,navnum(ix))%cus&
					,rnxndt(ix,navnum(ix))%sqrta
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%toe,rnxndt(ix,navnum(ix))%cic,rnxndt(ix,navnum(ix))%omega&
					,rnxndt(ix,navnum(ix))%cis
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%i0,rnxndt(ix,navnum(ix))%crc,rnxndt(ix,navnum(ix))%omega1&
					,rnxndt(ix,navnum(ix))%omegadot
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%idot,rnxndt(ix,navnum(ix))%codesonl2channel&
					,rnxndt(ix,navnum(ix))%weeks,rnxndt(ix,navnum(ix))%l2pdataflag
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%svaccuracy,rnxndt(ix,navnum(ix))%svhealth,rnxndt(ix,navnum(ix))%tgd&
					,rnxndt(ix,navnum(ix))%iodc
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%transtimeofmsg,rnxndt(ix,navnum(ix))%fitinterval&
					,rnxndt(ix,navnum(ix))%spare1,rnxndt(ix,navnum(ix))%spare2
            	elseif(rnxndt(ix,navnum(ix))%name(1:1)=="C")then
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%aode,rnxndt(ix,navnum(ix))%crs,rnxndt(ix,navnum(ix))%deltan&
					,rnxndt(ix,navnum(ix))%m0
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%cuc,rnxndt(ix,navnum(ix))%e,rnxndt(ix,navnum(ix))%cus&
					,rnxndt(ix,navnum(ix))%sqrta
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%toe,rnxndt(ix,navnum(ix))%cic,rnxndt(ix,navnum(ix))%omega&
					,rnxndt(ix,navnum(ix))%cis
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%i0,rnxndt(ix,navnum(ix))%crc,rnxndt(ix,navnum(ix))%omega1&
					,rnxndt(ix,navnum(ix))%omegadot
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%idot,rnxndt(ix,navnum(ix))%spare1,rnxndt(ix,navnum(ix))%weeks&
					,rnxndt(ix,navnum(ix))%spare2
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%svaccuracy,rnxndt(ix,navnum(ix))%sath1,rnxndt(ix,navnum(ix))%tgd1&
					,rnxndt(ix,navnum(ix))%tgd2
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%transtimeofmsg,rnxndt(ix,navnum(ix))%aodc&
					,rnxndt(ix,navnum(ix))%spare3,rnxndt(ix,navnum(ix))%spare4
            	elseif(rnxndt(ix,navnum(ix))%name(1:1)=="S")then
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%posx,rnxndt(ix,navnum(ix))%xdot,rnxndt(ix,navnum(ix))%xacceleration&
					,rnxndt(ix,navnum(ix))%health
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%posy,rnxndt(ix,navnum(ix))%ydot,rnxndt(ix,navnum(ix))%yacceleration&
					,rnxndt(ix,navnum(ix))%accuracy
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%posz,rnxndt(ix,navnum(ix))%zdot,rnxndt(ix,navnum(ix))%zacceleration&
					,rnxndt(ix,navnum(ix))%iodn
            	elseif(rnxndt(ix,navnum(ix))%name(1:1)=="I")then
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%iodec,rnxndt(ix,navnum(ix))%crs,rnxndt(ix,navnum(ix))%deltan&
					,rnxndt(ix,navnum(ix))%m0
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%cuc,rnxndt(ix,navnum(ix))%e,rnxndt(ix,navnum(ix))%cus&
					,rnxndt(ix,navnum(ix))%sqrta
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%toe,rnxndt(ix,navnum(ix))%cic,rnxndt(ix,navnum(ix))%omega&
					,rnxndt(ix,navnum(ix))%cis
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%i0,rnxndt(ix,navnum(ix))%crc,rnxndt(ix,navnum(ix))%omega1&
					,rnxndt(ix,navnum(ix))%omegadot
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%idot,rnxndt(ix,navnum(ix))%spare1,rnxndt(ix,navnum(ix))%weeks&
					,rnxndt(ix,navnum(ix))%spare1
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%accuracy,rnxndt(ix,navnum(ix))%health,rnxndt(ix,navnum(ix))%tgd&
					,rnxndt(ix,navnum(ix))%spare1
            	    read(1234,'(a)',iostat=jump)tmp
            	    read(tmp,'(4x,4d19.12)')rnxndt(ix,navnum(ix))%transtimeofmsg,rnxndt(ix,navnum(ix))%spare1&
					,rnxndt(ix,navnum(ix))%spare1,rnxndt(ix,navnum(ix))%spare1
            	endif
			endif
		endif
    enddo
end subroutine
end module
