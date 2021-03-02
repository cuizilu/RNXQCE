module m_qc_slip
	use parameters
	use m_rnxqce_rinexread
	use m_rnxqce_control
	use m_qc_findlack,only:lack_result
	use m_rnxqce_frequency

	type slip_result
		real*8::LG=0
		real*8::MW=0
	endtype
	
	type slip_find
		integer::num=0
		real*8::pn=0
		real*8::ex=0
		real*8::ex2=0
		real*8::dx=0
	endtype

	contains
	subroutine qc_slip(rnxohd,rnxodt,epochsum,CTRL,flag_all,lack)
	implicit none
		type(T_rnxohd)::rnxohd
		type(T_rnxodt),allocatable::rnxodt(:)
		type(T_rnxqce_ctrl)::CTRL
		type(slip_result),allocatable::slip(:,:)
		integer*4,allocatable::flag_all(:,:)
		integer::slip_flag(4,3),epochsum,i,j,ix,memory(satsum),aversum=0,avernum=1,slipnum(3),slipsum(3)
		real*8::f1,f2,alpha,LG,MW
		character(len=3)::obstype1(3),obstype2(3)
		type(lack_result),allocatable::lack(:)
		type(slip_find)::LGfind(satsum),MWfind(satsum)
		allocate(slip(epochsum,satsum))
		slip_flag(:,:)=0
		memory(:)=0
		slipnum(:)=1
		slipsum(:)=0
        do i=1,rnxohd%systemnumber
            ix=index('GCE',rnxohd%system(i)%systemname)
            if(ix==0)then
                cycle
            endif
            do j=1,rnxohd%system(i)%obstypenumber!1,2-->phase observation 3,4-->code observation
                if(CTRL%sysfreq(ix)(2:2)==rnxohd%system(i)%obstype(j)(2:2).and.rnxohd%system(i)%obstype(j)(1:1)=="L".and.&
                (lack(i)%obstypelack(slip_flag(1,ix))>lack(i)%obstypelack(j).or.slip_flag(1,ix)==0))then
                    slip_flag(1,ix)=j
					obstype1(ix)=rnxohd%system(i)%obstype(j)
                endif
                if(CTRL%sysfreq(ix)(3:3)==rnxohd%system(i)%obstype(j)(2:2).and.rnxohd%system(i)%obstype(j)(1:1)=="L".and.&
                (lack(i)%obstypelack(slip_flag(2,ix))>lack(i)%obstypelack(j).or.slip_flag(2,ix)==0))then
                    slip_flag(2,ix)=j
					obstype2(ix)=rnxohd%system(i)%obstype(j)
                endif
                if(CTRL%sysfreq(ix)(2:2)==rnxohd%system(i)%obstype(j)(2:2).and.(rnxohd%system(i)%obstype(j)(1:1)=="C".or.&
				rnxohd%system(i)%obstype(j)(1:1)=="P").and.&
				(lack(i)%obstypelack(slip_flag(3,ix))>lack(i)%obstypelack(j).or.slip_flag(3,ix)==0))then
                    slip_flag(3,ix)=j
                endif
                if(CTRL%sysfreq(ix)(3:3)==rnxohd%system(i)%obstype(j)(2:2).and.(rnxohd%system(i)%obstype(j)(1:1)=="C".or.&
				rnxohd%system(i)%obstype(j)(1:1)=="P").and.&
				(lack(i)%obstypelack(slip_flag(4,ix))>lack(i)%obstypelack(j).or.slip_flag(4,ix)==0))then
                    slip_flag(4,ix)=j
                endif
            enddo
        enddo
		do i=1,epochsum
			if(rnxodt(i)%epochflag==0)then
				do j=1,satsum
					if(flag_all(i,j)==0)then
						ix=index('GCE',rnxodt(i)%satobsdata(j)%prn(1:1))
						call cal_freq(rnxodt(i)%satobsdata(j)%prn,obstype1(ix),f1,rnxohd%version)
						call cal_freq(rnxodt(i)%satobsdata(j)%prn,obstype2(ix),f2,rnxohd%version)
						if(rnxodt(i)%satobsdata(j)%obs(slip_flag(3,ix))==0.or.rnxodt(i)%satobsdata(j)%obs(slip_flag(4,ix))==0)then
							slip(i,j)%LG=0
							slip(i,j)%MW=0
							flag_all(i,j)=5
						else
							alpha=f1*f1/f2/f2
							LG=slip(memory(j),j)%LG
							MW=slip(memory(j),j)%MW
							slip(i,j)%LG=vc/f1*rnxodt(i)%satobsdata(j)%obs(slip_flag(1,ix))-vc/f2*rnxodt(i)%satobsdata(j)%obs(slip_flag(2,ix))
							!slip(i,j)%LG=rnxodt(i)%satobsdata(j)%obs(slip_flag(1,ix))-f1/f2*rnxodt(i)%satobsdata(j)%obs(slip_flag(2,ix))
							slip(i,j)%MW=rnxodt(i)%satobsdata(j)%obs(slip_flag(1,ix))-rnxodt(i)%satobsdata(j)%obs(slip_flag(2,ix))&
							-(f1-f2)/(f1+f2)*(rnxodt(i)%satobsdata(j)%obs(slip_flag(3,ix))*f1/vc+rnxodt(i)%satobsdata(j)%obs(slip_flag(4,ix))*f2/vc)
							if(LGfind(j)%num==0)then
								LGfind(j)%num=1
								LGfind(j)%pn=1
								LGfind(j)%ex=slip(i,j)%LG
								LGfind(j)%ex2=slip(i,j)%LG*slip(i,j)%LG
								LGfind(j)%dx=0
							else
								LGfind(j)%num=LGfind(j)%num+1
								LGfind(j)%pn=1d0/LGfind(j)%num
								LGfind(j)%ex=LGfind(j)%pn*slip(i,j)%LG+LGfind(j)%ex*(1-LGfind(j)%pn)
								LGfind(j)%ex2=LGfind(j)%pn*slip(i,j)%LG*slip(i,j)%LG+LGfind(j)%ex2*(1-LGfind(j)%pn)
								LGfind(j)%dx=(LGfind(j)%ex2-LGfind(j)%ex*LGfind(j)%ex)**0.5
								if(abs(LG-slip(i,j)%LG)>4*LGfind(j)%dx)then
									LGfind(j)%num=0
									flag_all(memory(j),j)=1
									!print*,"GF",rnxodt(memory(j))%cal,rnxodt(memory(j))%satobsdata(j)%prn
								endif
							endif
							if(MWfind(j)%num==0)then
								MWfind(j)%num=1
								MWfind(j)%pn=1
								MWfind(j)%ex=slip(i,j)%MW
								MWfind(j)%ex2=slip(i,j)%MW*slip(i,j)%MW
								MWfind(j)%dx=0
							else
								MWfind(j)%num=MWfind(j)%num+1
								MWfind(j)%pn=1d0/MWfind(j)%num
								MWfind(j)%ex=MWfind(j)%pn*slip(i,j)%MW+MWfind(j)%ex*(1-MWfind(j)%pn)
								MWfind(j)%ex2=MWfind(j)%pn*slip(i,j)%MW*slip(i,j)%MW+MWfind(j)%ex2*(1-MWfind(j)%pn)
								MWfind(j)%dx=(MWfind(j)%ex2-MWfind(j)%ex*MWfind(j)%ex)**0.5
								if(abs(MW-slip(i,j)%MW)>4*MWfind(j)%dx)then
									MWfind(j)%num=0
									flag_all(memory(j),j)=1
									!print*,"MW",rnxodt(memory(j))%cal,rnxodt(memory(j))%satobsdata(j)%prn
								endif
							endif
							memory(j)=i
						endif
					endif
				enddo
			endif
		enddo
		do i=1,epochsum
			do j=1,satsum
				if(flag_all(i,j)==1)then
					select case(j)
						case(1:40)
							slipnum(1)=slipnum(1)+1
							slipsum(1)=slipsum(1)+1
						case(41:100)
							slipnum(2)=slipnum(2)+1
							slipsum(2)=slipsum(2)+1
						case default
							slipnum(3)=slipnum(3)+1
							slipsum(3)=slipsum(3)+1
					end select
					aversum=aversum+1
					avernum=avernum+1
				endif
				if(flag_all(i,j)==0)then
					select case(j)
						case(1:40)
							slipsum(1)=slipsum(1)+1
						case(41:100)
							slipsum(2)=slipsum(2)+1
						case default
							slipsum(3)=slipsum(3)+1
					end select
					aversum=aversum+1
				endif
			enddo
		enddo
		write(*,*)"o/slps"
		write(321123,*)"o/slps    expt    slips"
		do i=1,3
			if(slipsum(i)/=0)then
				write(*,'(1x,a1,f8.2)')"GCE"(i:i),slipsum(i)*1.0/slipnum(i)
				write(321123,'(1x,a1,f7.1,2x,i5,2x,i5)')"GCE"(i:i),slipsum(i)*1.0/slipnum(i),slipsum(i),slipnum(i)
			endif
		enddo
		write(*,'(1x,a1,f8.2)')"A",aversum*1.0/avernum
		write(321123,'(1x,a1,f7.1,2x,i5,2x,i5)')"A",aversum*1.0/avernum,aversum,avernum
	end subroutine
	end module
