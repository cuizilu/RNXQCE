module vertical_angle
    use calculatesiteposition
    use coordinate
    use m_rnxqce_rinexread
    use m_rnxqce_navread
	use m_rnxqce_control
    use clockbias
    use earthrot
    use timetrans
	contains
	subroutine calculate_vertical_angel(obshdr,navhdr,obsrec,navrec,epochsum,navsum,vz,flag_all,CTRL)
	implicit none
    type(T_rnxohd)::obshdr
	type(T_rnxqce_ctrl)::CTRL
	type(T_rnxodt),allocatable::obsrec(:)
    type(T_rnxnhd)::navhdr
	type(T_rnxndt),allocatable::navrec(:,:)
    type(gpsTime)::toc
    type(geodetic)::geo,geo11
	type(cartesian)::x2,sitepos
	real,allocatable::vz(:,:)
	integer,allocatable::flag_all(:,:)
	integer::i,j,k,l,ix,navstart(satsum),navsum(satsum),epochsum,flag,flag1,test
	real*8::start,finish,deltat,hh,rh,tao,tr,ts,vts
	allocate(vz(epochsum,satsum))
!	test=1707
	vz(:,:)=0
	if(obshdr%version<3)then                                                                                                                              
    	do j=1,obshdr%obstypenumber
        	if(obshdr%obstype(j)=="C1")then
            	flag=j
        	endif
        enddo
    else
        flag=1
    endif
    navstart(:)=1
    do i=1,epochsum
        if(obsrec(i)%epochflag==0)then
        	!print*,"epoch",i
			!print*,obsrec(i)%numsat
        	do j=1,satsum
				if(flag_all(i,j)==9)then
					cycle
				endif
				ix=get_satpos(obsrec(i)%satobsdata(j)%prn)!index('GCE',obsrec(i)%satobsdata(j)%prn(1:1))
        	    call caltogps(obsrec(i)%cal,toc)
				do k=1,navsum(ix)!navstart(ix),navsum(ix)
					if((toc%tow%sn-navrec(ix,k)%toc%tow%sn)<7200)then
						navstart(ix)=k
						exit
					endif
				enddo
        	    flag1=0
        	    if(obsrec(i)%satobsdata(j)%obs(flag)/=0)then
        	        tao=obsrec(i)%satobsdata(j)%obs(flag)/vc
        	    else
        	        do while(obsrec(i)%satobsdata(j)%obs(flag+flag1)==0)
						flag1=flag1+1
        	        enddo
        	        tao=obsrec(i)%satobsdata(j)%obs(flag+flag1)/vc
        	    endif
        	    tr=obsrec(i)%epoch%wn*24*60*60+obsrec(i)%epoch%tow%sn+obsrec(i)%epoch%tow%tos
        	    ts=tr-tao
        	    call cpu_time(start)
        	    do l=navstart(ix),navsum(ix)!navstart(ix),navsum(ix)
        	        call CPU_TIME(finish)
        	        if(navrec(ix,l)%name==obsrec(i)%satobsdata(j)%prn.and.navrec(ix,l)%toc%tow%sn<=toc%tow%sn&
					.and.(toc%tow%sn-navrec(ix,l)%toc%tow%sn)<14400)then
						!if(i==test)then
						!	print*,navrec(ix,l)%cal,obsrec(i)%satobsdata(j)%prn
						!endif
        	            ts=ts-navrec(ix,l)%toc%wn*24*60*60!-(navrec(xingliwei(i))%toc%tow%tos+navrec(xingliwei(i))%toc%tow%sn)
        	            tr=tr-navrec(ix,l)%toc%wn*24*60*60!-(navrec(xingliwei(i))%toc%tow%tos+navrec(xingliwei(i))%toc%tow%sn)
        	            call calclock(navrec(ix,l)%clkbias,navrec(ix,l)%clkdrift,navrec(ix,l)%clkdriftrate,ts,vts)
        	            call calsitepos(navrec(ix,l),ts,deltat,sitepos)
        	            vts=vts+deltat
        	            ts=ts-vts
        	            call calsitepos(navrec(ix,l),ts,deltat,sitepos)
        	            call calearth(ts,tr,sitepos,x2)
        	            call cartogeo(obshdr%approxpos,geo11)
        	            rh=sqrt((x2%x-obshdr%approxpos%x)**2+(x2%y-obshdr%approxpos%y)**2+(x2%z-obshdr%approxpos%z)**2)
        	            hh=pi/2d0-acos((cos(geo11%b)*(cos(geo11%l)*(x2%x-obshdr%approxpos%x)+sin(geo11%l)&
						*(x2%y-obshdr%approxpos%y))+sin(geo11%b)*(x2%z-obshdr%approxpos%z))/rh)
        	            vz(i,get_satpos(obsrec(i)%satobsdata(j)%prn))=hh/pi*180
						exit
        	        else if(navrec(ix,l)%toc%tow%sn>toc%tow%sn)then
        	            exit
        	        endif
        	    enddo
				if(vz(i,get_satpos(obsrec(i)%satobsdata(j)%prn))<CTRL%cutoff.and.flag_all(i,get_satpos&
				(obsrec(i)%satobsdata(j)%prn))/=9)then
					flag_all(i,get_satpos(obsrec(i)%satobsdata(j)%prn))=6
				endif
        	   ! if(abs(obsrec(i)%satobsdata(j)%vertical_angle)>=10.0)then
        	   !     test1=test1+1
        	   ! endif
        	   ! if(abs(obsrec(i)%satobsdata(j)%vertical_angle)<10.and.abs(obsrec(i)%satobsdata(j)%vertical_angle)>0)then
        	   !     test2=test2+1
        	   ! endif
        	enddo
        endif
    enddo
	!print*,vz(test,:),obsrec(test)%cal
	end subroutine

	subroutine vzprint(vz,epochsum)
	implicit none
		real,allocatable::vz(:,:)
		integer::epochsum,i,j
		write(321123,"(a)")"vertical angel"
		write(321123,"(a)")"format:epoch prn vz prn vz prn vz..."
		do i=1,epochsum
			write(321123,"(i4,1x,$)")i
			do j=1,satsum
				if(vz(i,j)/=0)then
					write(321123,"(1x,a3,f5.1,$)")get_satprn(j),vz(i,j)
				endif
			enddo
			write(321123,"(/)")
		enddo
	end subroutine
	end module
