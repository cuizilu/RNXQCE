module m_rnxqce_rinexmodify
use m_rnxqce_rinexread
use m_rnxqce_control
use parameters

contains
subroutine modify_rinex(CTRL)
	implicit none
	character(len=tmplength)::tmp
	integer::ios,jump,i
	character(len=stringlength)::path
	character(len=20)::empty=" "
	type(T_rnxqce_ctrl) CTRL
	open(321,file=CTRL%frnxo,blank="null",status='old',iostat=ios)
	path=trim(CTRL%frnxo)//"cui"
	open(unit=3210,file=path,form='formatted',status='replace')
	do i=1,len(CTRL%anttype)
		if(CTRL%anttype(i:i)=='^')then
			CTRL%anttype(i:i)=' '
		endif
	enddo
	do i=1,len(CTRL%rectype)
		if(CTRL%rectype(i:i)=='^')then
			CTRL%rectype(i:i)=' '
		endif
	enddo
    if(ios>0)then
        print*,"undifined file, please check your input"
    else
		do while(.true.)
        	read(321,'(a)',iostat=jump)tmp
			if(jump/=0)then
    	        exit
	        endif
			if(tmp(61:72)=="ANT # / TYPE")then
				if(CTRL%anttype/='')then
					write(3210,'(3a20,a20)')tmp(1:20),adjustl(trim(CTRL%anttype))//empty,tmp(41:60),tmp(61:80)
				endif
			else if(tmp(61:65)=="REC #")then
				if(CTRL%rectype/='')then
					write(3210,'(3a20,a20)')tmp(1:20),adjustl(trim(CTRL%rectype))//empty,tmp(41:60),tmp(61:80)
				endif
			else if(tmp(61:80)=="ANTENNA: DELTA H/E/N")then
				if(CTRL%enu(1)/=1000.or.CTRL%enu(2)/=1000.or.CTRL%enu(3)/=1000)then
					write(3210,'(3f14.4,18x,a20)')CTRL%enu(1),CTRL%enu(2),CTRL%enu(3),tmp(61:80)
				endif
			else
				write(3210,'(a)')trim(tmp)
			endif
		enddo
	endif
	close(321,status="delete")
	close(3210)
	call rename(path,CTRL%frnxo,ios)
end subroutine
end module

