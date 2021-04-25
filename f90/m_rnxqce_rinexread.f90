module m_rnxqce_rinexread
use m_rnxqce_control
use coordinate
use timetrans
type obsys
    character*1::systemname!系统名称
    integer::systemsum!系统总数
    integer::obstypenumber!观测值类型数目
    character*3,allocatable::obstype(:)!观测值类型集合数组
endtype


type T_rnxohd
	real*4::version!版本号
	integer::navflag!是否有星历文件输入
	character::systemtype*1
	character::mark_name*60
	character::receivernumber*20
	character::receivertype*20
	character::receiverversion*20
	character::antennanumber*20
	character::antennatype*20
	type(cartesian)::approxpos
	type(topocentric)::antennadelta
	integer::obstypenumber!2版本观测值类型数
	character*2,allocatable::obstype(:)!2版本观测值类型
	integer,allocatable::obstypelack(:)
	integer,allocatable::obstypesum(:)
	type(obsys),allocatable::system(:)!系统集合数组
	character(len=stringlength)::systemsumstring!系统总和制成的字符串
	real*8::interval
	type(gpsTime)::starttime
	type(gpsTime)::endtime
	integer::systemnumber=0
endtype

type sat_data
	character*3::prn
	real*8,allocatable::obs(:)
	integer,allocatable::lli(:)
	integer,allocatable::signalstrength(:)
endtype

type T_rnxodt
	type(gpsTime)::epoch
	type(calendarTime)::cal
	integer::epochflag
	integer::numsat
	real*8::receiverclkbias
	integer::startline
	integer::endline
	type(sat_data)::satobsdata(satsum)! 1-40--GPS  41-100--BDS  101-150--GAL
end type

contains

subroutine read_rinex_head(path,rnxohd)
	integer::pathflag,ios,sysflag=0,flagcui1=0,jump,pp,obstype_lines=0,systemnumber=0,ii,flagcui=0
	real::version,ww
	type(calendarTime)::cal1,cal2
	character(len=stringlength)::path,path2
	character(len=tmplength)::tmp,tmp_obstype(20)
	type(T_rnxohd)rnxohd
	rnxohd%systemsumstring=""
	open(123,file=path,blank="null",status='old',iostat=ios)
	pathflag=1
	do while(path(pathflag:pathflag)/=".")
		pathflag=pathflag+1
	enddo
	path2=path(1:pathflag+2)//"s"
	open(unit=321123,file=path2,form='formatted',status='replace')
	if(ios>0)then
		print*,"undifined file, please check your input"
	else
		do while(.true.)
			read(123,'(a)',iostat=jump)tmp
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
            	    rnxohd%version=pp
            	elseif(flag==2)then
            	    rnxohd%version=ww
            	endif
            	rnxohd%systemtype=tmp(41:41)
				rnxohd%interval=30
        	endif
			if(tmp(61:71)=="MARKER NAME")then
                read(tmp,'(a60)')rnxohd%mark_name
            endif
            if(tmp(61:65)=="REC #")then
                read(tmp,'(3a20)')rnxohd%receivernumber,rnxohd%receivertype,rnxohd%receiverversion
            endif
            if(tmp(61:72)=="ANT # / TYPE")then
                read(tmp,'(2a20)')rnxohd%antennanumber,rnxohd%antennatype
            endif
            if(tmp(61:79)=="APPROX POSITION XYZ")then
                read(tmp,'(3f14.4)')rnxohd%approxpos%x,rnxohd%approxpos%y,rnxohd%approxpos%z
            endif
            if(tmp(61:80)=="ANTENNA: DELTA H/E/N")then
                read(tmp,'(3f14.4)')rnxohd%antennadelta%u,rnxohd%antennadelta%e,rnxohd%antennadelta%n
            endif
			if(tmp(61:79)=="# / TYPES OF OBSERV")then
				obstype_lines=obstype_lines+1
				tmp_obstype(obstype_lines)=tmp
			endif
			if(tmp(61:79)=="SYS / # / OBS TYPES")then
				obstype_lines=obstype_lines+1
				tmp_obstype(obstype_lines)=tmp
			endif
			if(tmp(61:68)=="INTERVAL")then
                read(tmp,'(f10.3)')rnxohd%interval
				if(rnxohd%interval==0)then
					rnxohd%interval=30
				endif
            endif
            if(tmp(61:77)=="TIME OF FIRST OBS")then
                read(tmp,'(5i6,f13.7)')cal1%year,cal1%month,cal1%day,cal1%hour,cal1%minute,cal1%second
                call caltogps(cal1,rnxohd%starttime)
            endif
            if(tmp(61:76)=="TIME OF LAST OBS")then
                read(tmp,'(5i6,f13.7)')cal2%year,cal2%month,cal2%day,cal2%hour,cal2%minute,cal2%second
                call caltogps(cal2,rnxohd%endtime)
            endif

		enddo
		if(rnxohd%version>3)then!read the line of obs types and the lines nearby
			do i=1,obstype_lines
				if(tmp_obstype(i)(61:79)=="SYS / # / OBS TYPES")then
					if(tmp_obstype(i)(1:1)/=" ")then
                        systemnumber=systemnumber+1
                    endif
                endif
            enddo      
			rnxohd%systemnumber=systemnumber
			allocate(rnxohd%system(systemnumber))
		endif
		if(rnxohd%version<3)then
			read(tmp_obstype(1),'(i6)')rnxohd%obstypenumber
			allocate(rnxohd%obstype(rnxohd%obstypenumber))
			do i=1,rnxohd%obstypenumber
				if(mod(i,9)==0)then
					flagcui=9
				else
					flagcui=0
				endif
				read(tmp_obstype(1+(i-mod(i,9)-flagcui)/9)((flagcui+mod(i,9))*6+5:(flagcui+mod(i,9))*6+6),'(a2)')rnxohd%obstype(i)
			enddo
		else
			do i=1,obstype_lines
				if(flagcui>=1)then
					flagcui=flagcui-1
					cycle
				endif
				sysflag=sysflag+1
                read(tmp_obstype(i),'(a1,3x,i2)')rnxohd%system(sysflag)%systemname,rnxohd%system(sysflag)%obstypenumber
				rnxohd%systemsumstring=rnxohd%systemsumstring(1:sysflag-1)//rnxohd%system(sysflag)%systemname
                allocate(rnxohd%system(sysflag)%obstype(rnxohd%system(sysflag)%obstypenumber))
                if(rnxohd%system(sysflag)%obstypenumber>13.and.mod(rnxohd%system(sysflag)%obstypenumber,13)/=0)then
                	flagcui=int(rnxohd%system(sysflag)%obstypenumber/13)
                elseif(rnxohd%system(sysflag)%obstypenumber>13.and.mod(rnxohd%system(sysflag)%obstypenumber,13)==0)then
                    flagcui=int(rnxohd%system(sysflag)%obstypenumber/13)-1
                endif
                do ii=1,rnxohd%system(sysflag)%obstypenumber
                    if(mod(ii,13)==0)then
                    	flagcui1=13
                    else
                        flagcui1=0
                    endif
                    read(tmp_obstype(i+(ii-mod(ii,13)-flagcui1)/13)((flagcui1+mod(ii,13))*4+4:(flagcui1+mod(ii,13))*4+7),'(a3)')&
					rnxohd%system(sysflag)%obstype(ii)
                    rnxohd%system(sysflag)%systemsum=sysflag
                enddo
			enddo
		endif
	endif
end subroutine

subroutine read_rinex_data(rnxohd,rnxodt,flag,epochsum)
	implicit none
	type(T_rnxohd)rnxohd
	type(T_rnxodt),allocatable::rnxodt(:)
	type(calendarTime)::cal
	type(jd)::jd0,jd1
	integer::i,ix,j,jump,obstypenumber(3),epochsum,epochflag,numsat,epoch,numsatflag,systemnum(3)
	integer*4,allocatable :: flag(:,:)
	integer::pos(satsum)
	real::receiverclkbias
	character(len=3)::prn
	character(len=tmplength)::tmp,tmp_near(3)
	epochsum=24*60*60/int(rnxohd%interval)
	jd0%day=-1
	allocate(rnxodt(epochsum))
	allocate(flag(epochsum,satsum))
	flag(:,:)=9
	if(rnxohd%version>3)then
		do i=1,rnxohd%systemnumber
			if(rnxohd%system(i)%systemname=="G")then
				obstypenumber(1)=rnxohd%system(i)%obstypenumber
				systemnum(1)=i
			elseif(rnxohd%system(i)%systemname=="C")then
				obstypenumber(2)=rnxohd%system(i)%obstypenumber
				systemnum(2)=i
			elseif(rnxohd%system(i)%systemname=="E")then
				obstypenumber(3)=rnxohd%system(i)%obstypenumber
				systemnum(3)=i
			endif
		enddo
	endif
	do while(.true.)
		read(123,'(a)',iostat=jump)tmp
		if(jump/=0)then
			exit
		endif
		if(rnxohd%version<3)then
			if(tmp(30:30)==" ".and.tmp(28:28)==" ".and.tmp(29:29)/=" ")then
				numsatflag=0
				read(tmp,'(1x,i2.2,4(1x,i2),f11.7,2x,i1,i3)')cal%year,cal%month,cal%day,cal%hour,cal%minute,cal%second,epochflag,numsat!get the epoch head
				if(epochflag==4)then
					cycle
				endif
				if(cal%year<100)then
                	if(cal%year<30)then
                	    cal%year=cal%year+2000
                	else
                	    cal%year=cal%year+1900
                	endif
                endif
				call caltojd(cal,jd1)
				if(jd0%day==-1)then
					jd0=jd1
				endif
				epoch=1+anint(((jd1%day-jd0%day)*24*60*60+jd1%tod%sn+anint(jd1%tod%tos)-jd0%tod%sn-anint(jd0%tod%tos))/rnxohd%interval)
				rnxodt(epoch)%cal=cal
				call caltogps(cal,rnxodt(epoch)%epoch)
				rnxodt(epoch)%epochflag=epochflag
				rnxodt(epoch)%numsat=numsat
				tmp_near(1)=tmp
				do i=1,get_lines(numsat,12)-1
					read(123,'(a)')tmp_near(i+1)
				enddo
				do i=1,numsat
					read(tmp_near(get_lines(i,12))((i-(get_lines(i,12)-1)*12)*3+30:(i-(get_lines(i,12)-1)*12)*3+32),'(a3)')prn
					if(index(rnxohd%systemsumstring(1:rnxohd%systemnumber),prn(1:1))==0)then
						rnxohd%systemsumstring=rnxohd%systemsumstring(1:rnxohd%systemnumber)//prn(1:1)
						rnxohd%systemnumber=rnxohd%systemnumber+1
					endif
					pos(i)=get_satpos(prn)
					rnxodt(epoch)%satobsdata(pos(i))%prn=prn
					if(pos(i)/=0)then
						if(.not.allocated(rnxodt(epoch)%satobsdata(pos(i))%obs))then
						allocate(rnxodt(epoch)%satobsdata(pos(i))%obs(rnxohd%obstypenumber))
						endif
						if(.not.allocated(rnxodt(epoch)%satobsdata(pos(i))%lli))then
						allocate(rnxodt(epoch)%satobsdata(pos(i))%lli(rnxohd%obstypenumber))
						endif
						if(.not.allocated(rnxodt(epoch)%satobsdata(pos(i))%signalstrength))then
						allocate(rnxodt(epoch)%satobsdata(pos(i))%signalstrength(rnxohd%obstypenumber))
						endif
					endif
				enddo
			else
				if(epochflag==4)then
					cycle
				endif
				numsatflag=numsatflag+1
				tmp_near(1)=tmp
				do i=1,get_lines(rnxohd%obstypenumber,5)-1
					read(123,'(a)')tmp_near(i+1)
				enddo
				flag(epoch,pos(numsatflag))=0
				do i=1,rnxohd%obstypenumber
					if(pos(numsatflag)/=0)then
						read(tmp_near(get_lines(i,5))((i-(get_lines(i,5)-1)*5)*16-15:(i-(get_lines(i,5)-1)*5)*16-2),'(f14.3)')&
						rnxodt(epoch)%satobsdata(pos(numsatflag))%obs(i)
						!if(rnxodt(epoch)%satobsdata(pos(numsatflag))%obs(i)==0)then
						!	flag(epoch,pos(numsatflag))=5
						!endif
						read(tmp_near(get_lines(i,5))((i-(get_lines(i,5)-1)*5)*16+15:(i-(get_lines(i,5)-1)*5)*16+15),'(i1)')&
						rnxodt(epoch)%satobsdata(pos(numsatflag))%lli(i)
						read(tmp_near(get_lines(i,5))((i-(get_lines(i,5)-1)*5)*16+16:(i-(get_lines(i,5)-1)*5)*16+16),'(i1)')&
						rnxodt(epoch)%satobsdata(pos(numsatflag))%signalstrength(i)
					endif
				enddo
			endif
		else!if rinex > 3
			if(tmp(1:1)==">")then
				if(tmp(3:3)==" ")then
					print*,"check your rinex, there may be some error after epoch",epoch
					exit
				endif
				read(tmp,'(2x,i4,4(1x,i2.2),f11.7,2x,i1,i3,6x,f15.12)')cal%year,cal%month,cal%day,cal%hour,cal%minute,&
				cal%second,epochflag,numsat,receiverclkbias
				if(epochflag==4)then
					cycle
				endif
				call caltojd(cal,jd1)
                if(jd0%day==-1)then
                    jd0=jd1
                endif
				epoch=1+anint(((jd1%day-jd0%day)*24*60*60+jd1%tod%sn+anint(jd1%tod%tos)-jd0%tod%sn-anint(jd0%tod%tos))/rnxohd%interval)
                rnxodt(epoch)%cal=cal
				call caltogps(cal,rnxodt(epoch)%epoch)
                rnxodt(epoch)%epochflag=epochflag
                rnxodt(epoch)%numsat=numsat
				rnxodt(epoch)%receiverclkbias=receiverclkbias
			else
				if(epochflag==4)then
					cycle
				endif
				pos(1)=get_satpos(tmp(1:3))
				if(pos(1)==0)then
					cycle
				endif
				rnxodt(epoch)%satobsdata(pos(1))%prn=tmp(1:3)
				ix=index('GCE',tmp(1:1))    
				if(pos(1)/=0)then
					if(.not.allocated(rnxodt(epoch)%satobsdata(pos(1))%obs))then
                    allocate(rnxodt(epoch)%satobsdata(pos(1))%obs(obstypenumber(ix)))
					endif
					if(.not.allocated(rnxodt(epoch)%satobsdata(pos(1))%lli))then
                    allocate(rnxodt(epoch)%satobsdata(pos(1))%lli(obstypenumber(ix)))
					endif
					if(.not.allocated(rnxodt(epoch)%satobsdata(pos(1))%signalstrength))then
                    allocate(rnxodt(epoch)%satobsdata(pos(1))%signalstrength(obstypenumber(ix)))
					endif
                endif
				if(ix/=0)then
					flag(epoch,pos(1))=0
					do i=1,obstypenumber(ix)
						read(tmp(4+(i-1)*16:19+(i-1)*16),'(f14.3,i1,i1)')rnxodt(epoch)%satobsdata(pos(1))%obs(i),&
						rnxodt(epoch)%satobsdata(pos(1))%lli(i),rnxodt(epoch)%satobsdata(pos(1))%signalstrength(i)
					enddo
				endif
			endif
		endif
	enddo
	close(123)
	epochsum=epoch
	if(rnxohd%version<3)then
		allocate(rnxohd%system(rnxohd%systemnumber))
		do i=1,rnxohd%systemnumber
			rnxohd%system(i)%systemname=rnxohd%systemsumstring(i:i)
			rnxohd%system(i)%obstypenumber=rnxohd%obstypenumber
			allocate(rnxohd%system(i)%obstype(rnxohd%system(i)%obstypenumber))
			do j=1,rnxohd%system(i)%obstypenumber
				rnxohd%system(i)%obstype(j)=rnxohd%obstype(j)
			enddo
		enddo
	endif
end subroutine

end module

