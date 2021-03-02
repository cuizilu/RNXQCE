module timetrans
    type calendarTime
        integer year
        integer month
        integer day
        integer hour
        integer minute
        real*16 second
    endtype  
    type timeofDay
        integer*8::sn
        real*16::tos
    endtype timeofDay
    type jd
        integer*8::day
        type(timeofDay)::tod
    endtype jd 
    type timeofWeek
       integer*8::sn
       real*16::tos
    endtype timeofWeek
    type gpsTime
        integer*8::wn
        type(timeofWeek)::tow
    endtype gpsTime
    type timeofDoy
       integer*8::sn
       real*16::tos
    endtype timeofDoy
    type doytime
        integer::year
        integer::doy
        type(timeofDoy)::toy
    endtype doytime
    
    contains 
    subroutine caltojd(cal,jd1)
    type(calendarTime)::cal
    type(jd)::jd1
    integer*8::y,m
    if(cal%month<=2)then
        y=cal%year-1
        m=cal%month+12
    elseif(cal%month>2)then
        y=cal%year
        m=cal%month
    endif
    jd1%day=int(365.25*y)+int(30.6001*(m+1))+cal%day+(cal%hour+cal%minute/60.0+cal%second/3600.0)/24+1720981.5
    jd1%tod%sn=int((int(365.25*y)+int(30.6001*(m+1))+cal%day+(cal%hour+cal%minute/60.0&
	+cal%second/3600.0)/24+1720981.5-jd1%day)*24*60*60)
    jd1%tod%tos=((int(365.25*y)+int(30.6001*(m+1))+cal%day+(cal%hour+cal%minute/60.0&
	+cal%second/3600.0)/24+1720981.5-jd1%day)*24*60*60)-jd1%tod%sn
    end subroutine
    
    subroutine jdtocal(cal,jd1)
    type(calendarTime)::cal
    type(jd)::jd1
    integer*8 a,b,c,d,e,m,y
    real*16 day
    a=int(jd1%day+(jd1%tod%sn+jd1%tod%tos)/24/60/60+0.5)
    b=a+1537
    c=int((b-122.1)/365.25)
    d=int(365.25*c)
    e=int((b-d)/30.600)
    day=b-d-int(30.6001*e)+(jd1%day+(jd1%tod%sn+jd1%tod%tos)/24/60/60+0.5)-a
    m=e-1-12*int(e/14.0)
    y=c-4715-int((7+m)/10)
    n=mod(a,7)
    cal%year=y
    cal%month=m
    cal%day=int(day)
    cal%hour=int((day-int(day))*24)
    cal%minute=int((((day-int(day))*24)-cal%hour)*60)
    cal%second=int((((((day-int(day))*24)-cal%hour)*60)-cal%minute)*60)
    end subroutine
    
    subroutine jdtogps(jd1,gps)
    type(jd)::jd1
    type(gpsTime)::gps
    gps%wn=int(((jd1%day+(jd1%tod%sn+jd1%tod%tos)/24/60/60-2444244.5))/7)
    gps%tow%sn=int((jd1%day+(jd1%tod%sn+jd1%tod%tos)/24/60/60-2444244.5-7*gps%wn)*24*60*60)
    gps%tow%tos=((jd1%day+(jd1%tod%sn+jd1%tod%tos)/24/60/60-2444244.5-7*gps%wn)*24*60*60)-gps%tow%sn
    end subroutine 
    
    subroutine caltogps(cal,gps)
    type(calendarTime)::cal
    type(jd)::jd1
    type(gpsTime)::gps
    call caltojd(cal,jd1)
    call jdtogps(jd1,gps)
    end subroutine
    
    subroutine gpstojd(jd1,gps)
    type(jd)::jd1
    type(gpsTime)::gps
    jd1%day=int(gps%wn*7+(gps%tow%sn+gps%tow%tos)/86400.0+2444244.5)
    jd1%tod%sn=int(((gps%wn*7+(gps%tow%sn+gps%tow%tos)/86400.0+2444244.5)-jd1%day)*24*60*60)
    jd1%tod%tos=(((gps%wn*7+(gps%tow%sn+gps%tow%tos)/86400.0+2444244.5)-jd1%day)*24*60*60)-jd1%tod%sn
    end subroutine
    
    
    subroutine gpstocal(cal,gps)
    type(calendarTime)::cal
    type(jd)::jd1
    type(gpsTime)::gps
    call gpstojd(jd1,gps)
    call jdtocal(cal,jd1)
    end subroutine
    
    subroutine jdtodoy(jd1,doy1)
    type(calendarTime)::cal
    type(calendarTime)::cal1
    type(doytime)::doy1
    type(jd)::jd1
    type(jd)::jd2
    call jdtocal(cal,jd1)
    cal1%year=cal%year
    cal1%month=1
    cal1%day=1
    cal1%hour=0
    cal1%minute=0
    cal1%second=0
    call caltojd(cal1,jd2)
    doy1%year=cal%year
    doy1%doy=jd1%day-jd2%day+1
    doy1%toy%sn=jd1%tod%sn-jd2%tod%sn
    doy1%toy%tos=jd1%tod%tos-jd2%tod%tos
    end subroutine
    
    subroutine doytojd(jd1,doy1)
    type(doytime)::doy1
    type(calendarTime)::cal1
    type(jd)::jd2
    type(jd)::jd1
    cal1%year=doy1%year
    cal1%month=1
    cal1%day=1
    cal1%hour=0
    cal1%minute=0
    cal1%second=0
    call caltojd(cal1,jd2)
    jd1%day=doy1%doy+jd2%day-1
    jd1%tod%sn=doy1%toy%sn+jd2%tod%sn
    jd1%tod%tos=doy1%toy%tos+jd2%tod%tos
    end subroutine
    
    subroutine caltodoy(cal,doy1)
    type(calendarTime)::cal
    type(doytime)::doy1
    type(jd)::jd1
    call caltojd(cal,jd1)
    call jdtodoy(jd1,doy1)
    end subroutine
    
    subroutine gpstodoy(gps,doy1)
    type(gpsTime)::gps
    type(jd)::jd1
    type(doytime)::doy1
    call gpstojd(jd1,gps)
    call jdtodoy(jd1,doy1)
    end subroutine
    
    subroutine doytocal(cal,doy1)
    type(calendarTime)::cal
    type(doytime)::doy1
    type(jd)::jd1
    call doytojd(jd1,doy1)
    call jdtocal(cal,jd1)
    end subroutine
    
    subroutine doytogps(gps,doy1)
    type(gpsTime)::gps
    type(jd)::jd1
    type(doytime)::doy1
    call doytojd(jd1,doy1)
    call jdtogps(jd1,gps)
    end subroutine
    
    end module
