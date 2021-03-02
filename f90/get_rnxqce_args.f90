subroutine get_rnxqce_args(CTRL)
	  use M_rnxqce_control
	  implicit none
	  type(T_rnxqce_ctrl) CTRL
      integer*4 np,i,ix,j
	  character::cc*20

	  np = iargc()
	  if(np.le.1) goto 100

	  i=1
	  do while(i.le.np) 
         call getarg(i,cc)
		 if(cc.eq.'-rnxo') then
            call getarg(i+1,CTRL%frnxo)
		 else if(cc.eq.'-qc') then
			CTRL%lqc=.true.
		 else if(cc.eq.'-mr') then
		 	CTRL%lmr=.true.
		 else if(cc.eq.'-rnxn') then
            call getarg(i+1,CTRL%frnxn)
		 else if(cc.eq.'-freq') then
			j=1
            do while((i+j).le.np) 
			   call getarg(i+j,cc)
			   if(cc(1:1).eq.'-')then
			   		exit
			   endif
			   ix=index('GCE',cc(1:1))
			   CTRL%sysfreq(ix)=cc
			   j=j+1
			enddo
		 else if(cc.eq.'-hd.ant') then
		 	j=1
		 	do while((i+j).le.np)
				call getarg(i+j,cc)
				if(cc(1:1)=='-')then
					exit
				endif
				if(j/=1)then
					CTRL%anttype=trim(CTRL%anttype)//" "//trim(cc)
				else
					CTRL%anttype=cc
				endif
				j=j+1
			enddo
		 else if(cc.eq.'-hd.rec') then
		 	j=1
            do while((i+j).le.np)
                call getarg(i+j,cc)
                if(cc(1:1)=='-')then
                    exit
                endif
                CTRL%rectype=trim(CTRL%rectype)//" "//trim(cc)
                j=j+1
            enddo
		 else if(cc.eq.'-hd.enu') then
            call getarg(i+1,cc)
			read(cc,*) CTRL%enu(1)
            call getarg(i+2,cc)
			read(cc,*) CTRL%enu(2)
            call getarg(i+3,cc)
			read(cc,*) CTRL%enu(3)
		 else if(cc.eq.'-cutoff') then
            call getarg(i+1,cc)
			read(cc,*) CTRL%cutoff
		else

	    endif

		i=i+1    

	  enddo

    if(len_trim(CTRL%frnxo).eq.0) then
	  write(*,*) '***ERROR: RINEX O-file must be defined, check the input '
	  goto 100
	endif

	return

100 continue
    write(*,*) ''
    write(*,*) './rnxqce -rnxo bjfs1000.20o -qc -rnxn brdm100.20p -freq G12 C26 R12 E25 -cutoff 10 ' 
    
    write(*,*) './rnxqce -rnxo bjfs1000.20o -mr -hd.ant TRM59800.00 -hd.rec TRIMBLE^NETR9  -hd.enu 0.0 0.0 0.0 ' 
    stop 	
end subroutine

