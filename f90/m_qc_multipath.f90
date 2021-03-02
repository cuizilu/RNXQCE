module m_qc_multipath
	use parameters
	use m_rnxqce_rinexread
	use m_rnxqce_control
	use m_qc_findlack,only:lack_result
	use m_rnxqce_frequency

	type mp_result
		real*8::mp1=0
		real*8::mp2=0
	end type

	type mp_sum
    	integer::num=0
        real*8::pn=0
        real*8::ex=0
        real*8::ex2=0
        real*8::dx=0
    endtype

	contains
	subroutine qc_mp(rnxohd,rnxodt,epochsum,CTRL,flag_all,lack)
		implicit none
		type(T_rnxohd)::rnxohd
		type(T_rnxodt),allocatable::rnxodt(:)
		type(T_rnxqce_ctrl)::CTRL
		type(mp_result),allocatable::mp(:,:)
		type(mp_result)::old
		type(mp_sum)::mp1sum(satsum),mp2sum(satsum)
		integer*4,allocatable::flag_all(:,:)
		integer::mp_flag(4,3),epochsum,i,j,ix,mp1_num(satsum),mp2_num(satsum),sys_num1(3),sys_num2(3),avernum1=0,avernum2=0
		real*8::f1,f2,alpha,mp1_sum(satsum),mp2_sum(satsum),sys_sum1(3),sys_sum2(3),aversum1=0,aversum2=0  
		character(len=3)::obstype1(3),obstype2(3)
		type(lack_result),allocatable::lack(:)
		allocate(mp(epochsum,satsum))
		mp_flag(:,:)=0
		mp1_num(:)=0
		mp1_sum(:)=0
		mp2_num(:)=0
		mp2_sum(:)=0
		sys_num1(:)=0
		sys_num2(:)=0
		sys_sum1(:)=0
		sys_sum2(:)=0
        do i=1,rnxohd%systemnumber
            ix=index('GCE',rnxohd%system(i)%systemname)
            if(ix==0)then
                cycle
            endif
            do j=1,rnxohd%system(i)%obstypenumber!1,2-->phase observation 3,4-->code observation
                if(CTRL%sysfreq(ix)(2:2)==rnxohd%system(i)%obstype(j)(2:2).and.rnxohd%system(i)%obstype(j)(1:1)=="L".and.&
                (lack(i)%obstypelack(mp_flag(1,ix))>lack(i)%obstypelack(j).or.mp_flag(1,ix)==0))then
                    mp_flag(1,ix)=j
                    obstype1(ix)=rnxohd%system(i)%obstype(j)
                endif
                if(CTRL%sysfreq(ix)(3:3)==rnxohd%system(i)%obstype(j)(2:2).and.rnxohd%system(i)%obstype(j)(1:1)=="L".and.&
                (lack(i)%obstypelack(mp_flag(2,ix))>lack(i)%obstypelack(j).or.mp_flag(2,ix)==0))then
                    mp_flag(2,ix)=j
                    obstype2(ix)=rnxohd%system(i)%obstype(j)
                endif
                if(CTRL%sysfreq(ix)(2:2)==rnxohd%system(i)%obstype(j)(2:2).and.(rnxohd%system(i)%obstype(j)(1:1)=="C".or.&
                rnxohd%system(i)%obstype(j)(1:1)=="P").and.&
                (lack(i)%obstypelack(mp_flag(3,ix))>lack(i)%obstypelack(j).or.mp_flag(3,ix)==0))then
                    mp_flag(3,ix)=j
                endif
                if(CTRL%sysfreq(ix)(3:3)==rnxohd%system(i)%obstype(j)(2:2).and.(rnxohd%system(i)%obstype(j)(1:1)=="C".or.&
                rnxohd%system(i)%obstype(j)(1:1)=="P").and.&
                (lack(i)%obstypelack(mp_flag(4,ix))>lack(i)%obstypelack(j).or.mp_flag(4,ix)==0))then
                    mp_flag(4,ix)=j
                endif
            enddo
        enddo
		do i=1,epochsum
			if(rnxodt(i)%epochflag==0)then
				do j=1,satsum
					if(flag_all(i,j)==0.or.flag_all(i,j)==1)then
                        ix=index('GCE',rnxodt(i)%satobsdata(j)%prn(1:1))
                        call cal_freq(rnxodt(i)%satobsdata(j)%prn,obstype1(ix),f1,rnxohd%version)
                        call cal_freq(rnxodt(i)%satobsdata(j)%prn,obstype2(ix),f2,rnxohd%version)   
						alpha=f1*f1/f2/f2
						mp(i,j)%mp1=rnxodt(i)%satobsdata(j)%obs(mp_flag(3,ix))-(1+2d0/(alpha-1))*rnxodt(i)%satobsdata(j)%&
						obs(mp_flag(1,ix))*vc/f1+(2d0/(alpha-1))*rnxodt(i)%satobsdata(j)%obs(mp_flag(2,ix))*vc/f2
						mp(i,j)%mp2=rnxodt(i)%satobsdata(j)%obs(mp_flag(4,ix))-(2*alpha)/(alpha-1)*rnxodt(i)%satobsdata(j)%&
						obs(mp_flag(1,ix))*vc/f1+(2*alpha/(alpha-1)-1)*rnxodt(i)%satobsdata(j)%obs(mp_flag(2,ix))*vc/f2
						if(mp1sum(j)%num==0)then
							mp1sum(j)%num=1
							mp1sum(j)%pn=1
							mp1sum(j)%ex=mp(i,j)%mp1
							mp1sum(j)%ex2=mp(i,j)%mp1*mp(i,j)%mp1
							mp1sum(j)%dx=0
						else
							mp1sum(j)%num=mp1sum(j)%num+1
							mp1sum(j)%pn=1d0/mp1sum(j)%num
							mp1sum(j)%ex=mp1sum(j)%pn*mp(i,j)%mp1+mp1sum(j)%ex*(1-mp1sum(j)%pn)
							mp1sum(j)%ex2=mp1sum(j)%pn*mp(i,j)%mp1*mp(i,j)%mp1+mp1sum(j)%ex2*(1-mp1sum(j)%pn)
							mp1sum(j)%dx=(mp1sum(j)%ex2-mp1sum(j)%ex*mp1sum(j)%ex)**0.5
							if((flag_all(i,j)==1.or.i==epochsum).and.mp1sum(j)%num>=3.and.mp1sum(j)%dx<2)then
								mp1_num(j)=mp1_num(j)+mp1sum(j)%num
								mp1_sum(j)=mp1_sum(j)+mp1sum(j)%dx*mp1sum(j)%num
								mp1sum(j)%num=0
							endif
						endif
						if(mp2sum(j)%num==0)then
							mp2sum(j)%num=1
							mp2sum(j)%pn=1
							mp2sum(j)%ex=mp(i,j)%mp2
							mp2sum(j)%ex2=mp(i,j)%mp2*mp(i,j)%mp2
							mp2sum(j)%dx=0
						else
							mp2sum(j)%num=mp2sum(j)%num+1
							mp2sum(j)%pn=1d0/mp2sum(j)%num
							mp2sum(j)%ex=mp2sum(j)%pn*mp(i,j)%mp2+mp2sum(j)%ex*(1-mp2sum(j)%pn)
							mp2sum(j)%ex2=mp2sum(j)%pn*mp(i,j)%mp2*mp(i,j)%mp2+mp2sum(j)%ex2*(1-mp2sum(j)%pn)
							mp2sum(j)%dx=(mp2sum(j)%ex2-mp2sum(j)%ex*mp2sum(j)%ex)**0.5
							if((flag_all(i,j)==1.or.i==epochsum).and.mp2sum(j)%num>=3.and.mp2sum(j)%dx<2)then
								mp2_num(j)=mp2_num(j)+mp2sum(j)%num
								mp2_sum(j)=mp2_sum(j)+mp2sum(j)%dx*mp2sum(j)%num
								mp2sum(j)%num=0
								!if(j==4)then
								!	print*,mp2_num(j),mp2_sum(j)
								!endif
							endif
						endif
					endif
				enddo
			endif
		enddo
		!do i=1,satsum
		!	print*,i
		!	print*,mp1_sum(i)/mp1_num(i),mp1_sum(i),mp1_num(i)
		!enddo
		do i=1,satsum
			if(mp1_num(i)>=50.and.mp1_sum(i)/mp1_num(i)<0.8)then
				if(i<=40)then
					sys_sum1(1)=sys_sum1(1)+mp1_sum(i)/mp1_num(i)
					sys_num1(1)=sys_num1(1)+1
					aversum1=aversum1+mp1_sum(i)/mp1_num(i)
					avernum1=avernum1+1
				else if(i>40.and.i<=100)then
					sys_sum1(2)=sys_sum1(2)+mp1_sum(i)/mp1_num(i)
					sys_num1(2)=sys_num1(2)+1
					aversum1=aversum1+mp1_sum(i)/mp1_num(i)
					avernum1=avernum1+1
				else if(i>100)then
					sys_sum1(3)=sys_sum1(3)+mp1_sum(i)/mp1_num(i)
					sys_num1(3)=sys_num1(3)+1
					aversum1=aversum1+mp1_sum(i)/mp1_num(i)
					avernum1=avernum1+1
				endif
			endif
		enddo
		do i=1,satsum
			if(mp2_num(i)>=50.and.mp2_sum(i)/mp2_num(i)<0.8.and.i<=40)then
				sys_sum2(1)=sys_sum2(1)+mp2_sum(i)/mp2_num(i)
				sys_num2(1)=sys_num2(1)+1
				aversum2=aversum2+mp2_sum(i)/mp2_num(i)
				avernum2=avernum2+1
			else if(mp2_num(i)>=50.and.mp2_sum(i)/mp2_num(i)<0.8.and.i>40.and.i<=100)then
				sys_sum2(2)=sys_sum2(2)+mp2_sum(i)/mp2_num(i)
				sys_num2(2)=sys_num2(2)+1
				aversum2=aversum2+mp2_sum(i)/mp2_num(i)
				avernum2=avernum2+1
			else if(mp2_num(i)>=50.and.mp2_sum(i)/mp2_num(i)<0.8.and.i>100)then
				sys_sum2(3)=sys_sum2(3)+mp2_sum(i)/mp2_num(i)
				sys_num2(3)=sys_num2(3)+1
				aversum2=aversum2+mp2_sum(i)/mp2_num(i)
				avernum2=avernum2+1
			endif
		enddo

		write(*,*)"MP1"
		write(321123,*)"MP1"
		do i=1,3
			if(sys_num1(i)/=0)then
				write(*,'(1x,a1,1x,f6.3)')"GCE"(i:i),sys_sum1(i)/sys_num1(i)
				write(321123,'(1x,a1,1x,f6.3)')"GCE"(i:i),sys_sum1(i)/sys_num1(i)
			endif
		enddo
		write(*,'(1x,a1,1x,f6.3)')"A",aversum1/avernum1
		write(321123,'(1x,a1,1x,f6.3)')"A",aversum1/avernum1
		write(*,*)"MP2"
		write(321123,*)"MP2"
		do i=1,3
			if(sys_num2(i)/=0)then
				write(*,'(1x,a1,1x,f6.3)')"GCE"(i:i),sys_sum2(i)/sys_num2(i)
				!print*,sys_sum2(i),sys_num2(i)
				write(321123,'(1x,a1,1x,f6.3)')"GCE"(i:i),sys_sum2(i)/sys_num2(i)
			endif
		enddo
		write(*,'(1x,a1,1x,f6.3)')"A",aversum2/avernum2
		write(321123,'(1x,a1,1x,f6.3)')"A",aversum2/avernum2
		write(321123,'(a)')"multipath"
		write(321123,'(a)')"PRN MP1   MP2"
		do i=1,satsum
			if(mp1_num(i)/=0.or.mp2_num(i)/=0)then
				write(321123,'(a3,f6.3,f6.3)')get_satprn(i),mp1_sum(i)/mp1_num(i),mp2_sum(i)/mp2_num(i)
			endif
		enddo
		!print*,"G",G_sum1/G_num1,"C",C_sum1/C_num1,"E",E_sum1/E_num1
		!print*,"G",G_sum2/G_num2,"C",C_sum2/C_num2,"E",E_sum2/E_num2
	end subroutine
end module
