module m_qc_findlack
	use parameters
	use m_rnxqce_rinexread
	use m_rnxqce_control

	type lack_result
		character(len=1)::systemname
		integer::obstypenumber
		character*3,allocatable::obstype(:)
		integer,allocatable::obstypelack(:)
		integer,allocatable::obstypesum(:)
	endtype
	contains
	subroutine qc_findlack(rnxohd,rnxodt,epochsum,CTRL,flag_all,lack)
	implicit none
		type(T_rnxohd)::rnxohd
		type(T_rnxodt),allocatable::rnxodt(:)
		type(T_rnxqce_ctrl)::CTRL
		type(lack_result),allocatable::lack(:)
		integer::test,epochsum,i,j,k,ix,lack_flag1(3),lack_flag2(3),lacksum(3),lacknum(3),avernum=0,aversum=0,l1num(3),l2num(3)
		character(len=1)::obstype_flag(2,3)
		character(len=3)::obstype(2,3)
		integer*4,allocatable::flag_all(:,:)
		allocate(lack(rnxohd%systemnumber))
		lack_flag1(:)=0
		lack_flag2(:)=0
		lacksum(:)=0
		lacknum(:)=0
		l1num(:)=0
		l2num(:)=0
		do i=1,rnxohd%systemnumber
			lack(i)%systemname=rnxohd%system(i)%systemname
			lack(i)%obstypenumber=rnxohd%system(i)%obstypenumber
			allocate(lack(i)%obstype(lack(i)%obstypenumber))
			allocate(lack(i)%obstypelack(lack(i)%obstypenumber))
			allocate(lack(i)%obstypesum(lack(i)%obstypenumber))
			do j=1,lack(i)%obstypenumber
				lack(i)%obstype(j)=rnxohd%system(i)%obstype(j)
				lack(i)%obstypelack(j)=0
				lack(i)%obstypesum(j)=0
			enddo
		enddo
		do i=1,epochsum
			if(rnxodt(i)%epochflag==0)then
				do j=1,satsum
					if(flag_all(i,j)==0)then
						ix=index(rnxohd%systemsumstring,rnxodt(i)%satobsdata(j)%prn(1:1)) 
						do k=1,lack(ix)%obstypenumber
							if(rnxodt(i)%satobsdata(j)%obs(k)==0)then
								lack(ix)%obstypelack(k)=lack(ix)%obstypelack(k)+1
								!if(k==14.and.j>40.and.j<100)then
								!print*,rnxodt(i)%cal,j,i
								!endif
							endif
							lack(ix)%obstypesum(k)=lack(ix)%obstypesum(k)+1
						enddo
					endif
				enddo
			endif
		enddo
		do i=1,rnxohd%systemnumber
			ix=index('GCE',lack(i)%systemname)
			if(ix==0)then
				cycle
			endif
			do j=1,lack(i)%obstypenumber
				if(CTRL%sysfreq(ix)(2:2)==lack(i)%obstype(j)(2:2).and.lack(i)%obstype(j)(1:1)=="L".and.&
				(lack(i)%obstypelack(lack_flag1(ix))>lack(i)%obstypelack(j).or.lack_flag1(ix)==0))then
					lack_flag1(ix)=j
					obstype_flag(1,ix)=lack(i)%obstype(j)(3:3)
					obstype(1,ix)=lack(i)%obstype(j)
				endif
				if(CTRL%sysfreq(ix)(3:3)==lack(i)%obstype(j)(2:2).and.lack(i)%obstype(j)(1:1)=="L".and.&
				(lack(i)%obstypelack(lack_flag2(ix))>lack(i)%obstypelack(j).or.lack_flag2(ix)==0))then
					lack_flag2(ix)=j
					obstype_flag(2,ix)=lack(i)%obstype(j)(3:3)
					obstype(2,ix)=lack(i)%obstype(j)
				endif
			enddo
		enddo
		do i=1,epochsum
			if(rnxodt(i)%epochflag==0)then
				do j=1,satsum
					if(flag_all(i,j)==0)then
						ix=index('GCE',rnxodt(i)%satobsdata(j)%prn(1:1))
						if(rnxodt(i)%satobsdata(j)%obs(lack_flag1(ix))==0.or.rnxodt(i)%satobsdata(j)%&
						obs(lack_flag2(ix))==0)then
							flag_all(i,j)=5
							lacknum(ix)=lacknum(ix)+1
							avernum=avernum+1
						endif
						if(rnxodt(i)%satobsdata(j)%obs(lack_flag1(ix))==0)then
							l1num(ix)=l1num(ix)+1
						endif
						if(rnxodt(i)%satobsdata(j)%obs(lack_flag2(ix))==0)then
							l2num(ix)=l2num(ix)+1
						endif
						lacksum(ix)=lacksum(ix)+1
						aversum=aversum+1
					endif
				enddo
			endif
		enddo
		!print*,flag_all(2694,49)
		do i=1,3
			if(lacksum(i)/=0)then
				if(obstype(1,i)(2:2)>"0".and.obstype(1,i)(2:2)<"9".and.obstype(2,i)(2:2)>"0".and.obstype(2,i)(2:2)<"9")then
					!print*,obstype(1,i)
				else
					write(*,'(a,a2,1x,a)')"#ERROR!There is some error in ","GCE"(i:i),"system"
				endif
			endif
		enddo
		write(321123,'(a)')"A means average, and is the average of all system"
		write(*,*)"完整率"
		write(321123,*)"完整率"
		do i=1,3
			if(lacksum(i)/=0)then
				write(*,'(a2,1x,f6.2,a)')"GCE"(i:i),(lacksum(i)-lacknum(i))*100.0/lacksum(i),"%"
				write(321123,'(a2,1x,f6.2,a)')"GCE"(i:i),(lacksum(i)-lacknum(i))*100.0/lacksum(i),"%"
				write(321123,'(a,i6)')"actual observation",lacksum(i)-lacknum(i)
				write(321123,'(a,i6)')"expect observation",lacksum(i)
			endif
		enddo
		write(*,'(a2,1x,f6.2,a)')"A",(aversum-avernum)*100.0/aversum,"%"
		write(321123,'(a2,1x,f6.2,a)')"A",(aversum-avernum)*100.0/aversum,"%"
		write(321123,'(a,i6)')"actual observation",aversum-avernum
		write(321123,'(a,i6)')"expect observation",aversum
		write(*,*)"缺失率"
		write(321123,*)"缺失率"
		do i=1,3
			if(lacksum(i)/=0)then
				write(*,'(a2,1x,a4,f6.2,a1,1x,a4,f6.2,a1)')"GCE"(i:i),obstype(1,i),l1num(i)*100.0/lacksum(i),"%",obstype(2,i)&
				,l2num(i)*100.0/lacksum(i),"%" 
				write(321123,'(a2,a4,5x,f6.2,a1,1x,a4,5x,f6.2,a1)')"GCE"(i:i),obstype(1,i),l1num(i)*100.0/lacksum(i),"%",obstype(2,i)&
				,l2num(i)*100.0/lacksum(i),"%" 
				write(321123,'(a2,a4,i6,a1,i6,a4,i6,a1,i6)')"GCE"(i:i),obstype(1,i),l1num(i),"/",lacksum(i),obstype(2,i)&
				,l2num(i),"/",lacksum(i) 
				if(l1num(i)*100.0/lacksum(i)>30.or.l2num(i)*100.0/lacksum(i)>30)then
					write(*,'(a,a2,1x,a)')"#WARNING!There may be some error in ","GCE"(i:i),"system"
				endif
			endif
		enddo
		!test=2
		!print*,lack(test)%systemname,lack(test)%obstype(test),lack(test)%obstypelack(test),lack(test)%obstypesum(test)
	end subroutine
end module

