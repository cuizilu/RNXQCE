module calculatesiteposition
    use m_rnxqce_navread
    use parameters
    use coordinate
    implicit none
    type siteposition
        real*8::x
        real*8::y
        real*8::z
    endtype

    contains 
    subroutine calsitepos(rec,t,deltat,sitepos)
	implicit none
    type(T_rnxndt)::rec
    type(cartesian)::sitepos
    real*8 t,deltat,tk,a,n0,mk,ek,ek1,n,vk,phik,phiuk,phirk,phiik,uk,rk,ik,xk,yk,lk,deltatr,deltatsv,start,finish
    ek=0
    tk=t-rec%toe
    if(tk>302400)then
        tk=tk-604800
    elseif(tk<-302400)then
        tk=tk+604800
    else
        tk=tk
    endif
    a=rec%sqrta*rec%sqrta
    n0=sqrt(miu/(a*a*a))
    n=n0+rec%deltan
    mk=rec%m0+n*tk
    ek1=mk
    call CPU_TIME(start)
    do while(abs(ek-ek1)>0.0000000001)
        call CPU_TIME(finish)
        if(finish-start>0.1)then
            exit
        endif
        ek=ek1
        ek1=mk+rec%e*sin(ek)
    enddo
    vk=atan2(sqrt(1-rec%e*rec%e)*sin(ek),(cos(ek)-rec%e))
    phik=vk+rec%omega1
    phiuk=rec%cus*sin(2*phik)+rec%cuc*cos(2*phik)
    phirk=rec%crs*sin(2*phik)+rec%crc*cos(2*phik)
    phiik=rec%cis*sin(2*phik)+rec%cic*cos(2*phik)
    uk=phik+phiuk
    rk=a*(1-rec%e*cos(ek))+phirk
    ik=rec%i0+phiik+rec%idot*tk
    xk=rk*cos(uk)
    yk=rk*sin(uk)
    lk=rec%omega+(rec%omegadot-omegae)*tk-omegae*rec%toe
    sitepos%x=xk*cos(lk)-yk*cos(ik)*sin(lk)
    sitepos%y=xk*sin(lk)+yk*cos(ik)*cos(lk)
    sitepos%z=yk*sin(ik)
    deltatr=-2*sqrt(a*miu)/vc/vc*rec%e*sin(ek)
    deltatsv=rec%clkbias+rec%clkdrift*(t-rec%toc%wn*7*24*60*60-rec%toc%tow%sn-rec%toc%tow%tos)+rec%clkdriftrate*(t-rec%toc%wn*7*24&
	*60*60-rec%toc%tow%sn-rec%toc%tow%tos)*(t-rec%toc%wn*7*24*60*60-rec%toc%tow%sn-rec%toc%tow%tos)+deltatr-rec%tgd
    deltat=deltatr
    endsubroutine
    endmodule
    
    
