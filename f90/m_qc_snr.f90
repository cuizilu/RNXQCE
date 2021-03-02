module m_qc_snr
	use parameters
	use m_rnxqce_rinexread
	use m_rnxqce_control
	use m_qc_findlack,only:lack_result

	type snr_result
		character(len=3)::snrtype(2,3)
		real::snrsum(2,3)
		integer::snrnum(2,3)
	endtype
	contains
	subroutine qc_snr(rnxohd,rnxodt,epochsum,CTRL,flag_all,lack)
	implicit none
		type(T_rnxohd)::rnxohd
		type(T_rnxodt),allocatable::rnxodt(:)
		type(T_rnxqce_ctrl)::CTRL
		type(snr_result)::snr
		integer*4,allocatable::flag_all(:,:)
		integer::epochsum,i,ix,j,snr_flag1(3),snr_flag2(3),avernum1=0,avernum2=0
		real::aversum1=0,aversum2=0
		type(lack_result),allocatable::lack(:)
		snr%snrsum(:,:)=0
		snr%snrnum(:,:)=0
		snr_flag1(:)=0
		snr_flag2(:)=0
		do i=1,rnxohd%systemnumber
		    ix=index('GCE',rnxohd%system(i)%systemname)
			if(ix==0)then
				cycle
			endif
		    do j=1,rnxohd%system(i)%obstypenumber
		        if(CTRL%sysfreq(ix)(2:2)==rnxohd%system(i)%obstype(j)(2:2).and.rnxohd%system(i)%obstype(j)(1:1)=="S".and.&
		        (lack(i)%obstypelack(snr_flag1(ix))>lack(i)%obstypelack(j).or.snr_flag1(ix)==0))then
		            snr_flag1(ix)=j
					snr%snrtype(1,ix)=rnxohd%system(i)%obstype(j)
		        endif
		        if(CTRL%sysfreq(ix)(3:3)==rnxohd%system(i)%obstype(j)(2:2).and.rnxohd%system(i)%obstype(j)(1:1)=="S".and.&
		        (lack(i)%obstypelack(snr_flag2(ix))>lack(i)%obstypelack(j).or.snr_flag2(ix)==0))then
		            snr_flag2(ix)=j
					snr%snrtype(2,ix)=rnxohd%system(i)%obstype(j)
		        endif
		    enddo
		enddo
		do i=1,epochsum
			do j=1,satsum
				if(flag_all(i,j)==0)then
					ix=index('GCE',rnxodt(i)%satobsdata(j)%prn(1:1))
					snr%snrsum(1,ix)=rnxodt(i)%satobsdata(j)%obs(snr_flag1(ix))+snr%snrsum(1,ix)
					snr%snrnum(1,ix)=snr%snrnum(1,ix)+1
					aversum1=aversum1+rnxodt(i)%satobsdata(j)%obs(snr_flag1(ix))
					avernum1=avernum1+1
					snr%snrsum(2,ix)=rnxodt(i)%satobsdata(j)%obs(snr_flag2(ix))+snr%snrsum(2,ix)
					snr%snrnum(2,ix)=snr%snrnum(2,ix)+1
					aversum2=aversum2+rnxodt(i)%satobsdata(j)%obs(snr_flag2(ix))
					avernum2=avernum2+1
				endif
			enddo
		enddo
		!print*,snr%snrsum(1,1)/snr%snrnum(1,1),snr%snrsum(2,1)/snr%snrnum(2,1)
		write(*,*)"SNR    1      2"
		write(321123,*)"SNR       1         2"
		do i=1,3
			if(snr%snrnum(1,i)/=0)then
				write(*,'(1x,a1,1x,f6.2,1x,f6.2)')"GCE"(i:i),snr%snrsum(1,i)/snr%snrnum(1,i),snr%snrsum(2,i)/snr%snrnum(2,i)
				write(321123,'(1x,a1,1x,f9.5,1x,f9.5)')"GCE"(i:i),snr%snrsum(1,i)/snr%snrnum(1,i),snr%snrsum(2,i)/snr%snrnum(2,i)
			endif
		enddo
		write(*,'(1x,a1,1x,f6.2,1x,f6.2)')"A",aversum1/avernum1,aversum2/avernum2
		write(321123,'(1x,a1,1x,f9.5,1x,f9.5)')"A",aversum1/avernum1,aversum2/avernum2
	end subroutine
end module
