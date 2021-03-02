module m_rnxqce_control
  use parameters 
type T_rnxqce_ctrl	

	 character(len=stringlength) :: frnxo='' 
	 character(len=stringlength) :: frnxn=''
	 logical*4 :: lqc=.false.
	 logical*4 :: lmr=.false.
	 character(len=3) sysfreq(4)
	 character(len=20):: rectype=''
	 character(len=20):: anttype=''
	 real*8 :: enu(3)=1000
	 real*8 :: cutoff=10

end type
 contains
	subroutine init_sysfreq(sysfreq)
	  implicit none
	  character(len=3) sysfreq(3)
	  sysfreq(1)='G12'
	  sysfreq(2)='C26'
	  sysfreq(3)='E15'
	end subroutine


end
