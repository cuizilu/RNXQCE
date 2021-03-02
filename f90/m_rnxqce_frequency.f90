module m_rnxqce_frequency
	use parameters
	contains
	subroutine cal_freq(prn,obstype,f,version)
	implicit none
	character(len=3)::prn,obstype
	real::version
	real*8::f
	f=0
	if(prn(1:1)=="G".and.obstype(2:2)=="1")then
	    f=1575.42d6
	else if(prn(1:1)=="G".and.obstype(2:2)=="2")then
	    f=1227.6d6
	else if(prn(1:1)=="G".and.obstype(2:2)=="5")then
	    f=1176.45d6
	else if(prn(1:1)=="E".and.obstype(2:2)=="1")then
	    f=1575.42d6
	else if(prn(1:1)=="E".and.obstype(2:2)=="5")then
	    f=1176.45d6
	else if(prn(1:1)=="E".and.obstype(2:2)=="7")then
	    f=1207.14d6
	else if(prn(1:1)=="E".and.obstype(2:2)=="8")then
	    f=1191.795d6
	else if(prn(1:1)=="E".and.obstype(2:2)=="6")then
	    f=1278.75d6
	else if(prn(1:1)=="S".and.obstype(2:2)=="1")then
	    f=1575.42d6   
	else if(prn(1:1)=="S".and.obstype(2:2)=="5")then
	    f=1176.45d6
	else if(prn(1:1)=="J".and.obstype(2:2)=="1")then
	    f=1575.42d6
	else if(prn(1:1)=="J".and.obstype(2:2)=="2")then
	    f=1227.6d6
	else if(prn(1:1)=="J".and.obstype(2:2)=="5")then
	    f=1176.45d6
	else if(prn(1:1)=="J".and.obstype(2:2)=="6")then
	    f=1278.75d6
	else if(prn(1:1)=="C".and.(obstype(2:2)=="2".or.obstype(2:2)=="1").and.version<3.025.and.version>3.015)then
		f=1561.098d6
	else if(prn(1:1)=="C".and.obstype(2:2)=="2".and.version/=3.02)then
	    f=1561.098d6
	else if(prn(1:1)=="C".and.obstype(2:2)=="1".and.version/=3.02)then
	    f=1575.42d6
	else if(prn(1:1)=="C".and.obstype(2:2)=="5")then
	    f=1176.45d6
	else if(prn(1:1)=="C".and.obstype(2:2)=="7")then
	    f=1207.14d6
	else if(prn(1:1)=="C".and.obstype(2:2)=="8")then
	    f=1191.795d6
	else if(prn(1:1)=="C".and.obstype(2:2)=="6")then
	    f=1268.52d6
	else if(prn(1:1)=="I".and.obstype(2:2)=="5")then
	    f=1176.45d6
	else if(prn(1:1)=="I".and.obstype(2:2)=="9")then
	    f=2492.028d6
	endif
	end subroutine 
	end module

