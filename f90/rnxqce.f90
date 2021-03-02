! rinex quality checking and editing
! rnxqce -rnxo bjfs1000.20o -qc -rnxn brdm100.20p 
!        -freq G12 C26 R12 E25 

! rnxqce -rnxo bjfs1000.20o -hd.ant TRM59800.00
!        -hd.rec TRIMBLE NETR9  -hd.enu 0.0 0.0 0.0  

! rnxqce -help 

program rnxqce 
	
	use m_rnxqce_rinexread
	use m_rnxqce_navread
    use M_rnxqce_control	
	use vertical_angle
	use parameters
	use m_qc_findlack
	use m_qc_snr
	use m_qc_slip
	use m_qc_multipath
	use m_rnxqce_rinexmodify
    implicit none
    

	type(T_rnxqce_ctrl) CTRL 
    type(T_rnxohd) rnxohd
    type(T_rnxnhd) rnxnhd 
    type(T_rnxodt),allocatable::rnxodt(:)
    type(T_rnxndt),allocatable::rnxndt(:,:)
	integer*4,allocatable :: flag(:,:)
	real,allocatable::vz(:,:)
	integer*4::epochsum,navsum(150)
	type(lack_result),allocatable::lack(:)
! flag == 9 no data
!      == 0 good data
!      == 1 slip start
!      == 5 lack data
!      == 6 cutoff data

	call init_sysfreq(CTRL%sysfreq)
	call get_rnxqce_args(CTRL)
    
	if(CTRL%lqc) then
		call read_rinex_head(CTRL%frnxo,rnxohd)
		call read_rinex_data(rnxohd,rnxodt,flag,epochsum)
		if(CTRL%frnxn/='')then
			call read_nav_head(CTRL%frnxn,rnxnhd)
			call read_nav_data(CTRL%frnxn,rnxnhd,rnxndt,navsum)
			call calculate_vertical_angel(rnxohd,rnxnhd,rnxodt,rnxndt,epochsum,navsum,vz,flag,CTRL)
		endif
	
	    call qc_findlack(rnxohd,rnxodt,epochsum,CTRL,flag,lack) 
	    call qc_snr(rnxohd,rnxodt,epochsum,CTRL,flag,lack)
	    call qc_slip(rnxohd,rnxodt,epochsum,CTRL,flag,lack)
	    call qc_mp(rnxohd,rnxodt,epochsum,CTRL,flag,lack)
		call vzprint(vz,epochsum)
	endif
	if(CTRL%lmr) then
!	  ! modify antenna type receiver type and antenna offsets
		call modify_rinex(CTRL)
!      call write_rnxo(CTRL,rnxohd,rnxorec)
	endif

end program
