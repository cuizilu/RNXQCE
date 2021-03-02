module coordinate
    use parameters
    type cartesian!笛卡尔坐标系
        real*16::x
        real*16::y
        real*16::z
    end type
    type geodetic!大地坐标系
        real*16::b
        real*16::l
        real*16::h
    endtype
    type topocentric!站心低平坐标系
        real*16::n
        real*16::e
        real*16::u
    endtype
    type topopolar!站心极坐标系
        real*16::s
        real*16::e
        real*16::a
    endtype
    
    contains
    subroutine geotocar(car,geo)!大地坐标系转笛卡尔坐标系
    type(cartesian)::car
    type(geodetic)::geo
    real*16 n
    n=am/sqrt(1-e2*sin(geo%b)*sin(geo%b))
    car%x=(n+geo%h)*cos(geo%b)*cos(geo%l)
    car%y=(n+geo%h)*cos(geo%b)*sin(geo%l)
    car%z=(n*(1-e2)+geo%h)*sin(geo%b)
    endsubroutine
    
    subroutine cartogeo(car,geo)!笛卡尔坐标系转大地坐标系
    type(cartesian)::car
    type(geodetic)::geo
    real*16 w,n,bb
    geo%b=0
    geo%l=atan2(car%y,car%x)
    bb=atan2(car%z,sqrt(car%x*car%x+car%y*car%y))
    do while(abs(bb-geo%b)>0.00000000000001)
    geo%b=bb
    w=sqrt(1-e2*sin(geo%b)*sin(geo%b))
    n=am/w
    bb=atan2((car%z+n*e2*sin(geo%b)),sqrt(car%x*car%x+car%y*car%y))
    geo%h=car%z/sin(geo%b)-n*(1-e2)
    enddo

    end subroutine
    
    subroutine cartotop(car1,car2,top)!笛卡尔坐标系转站心地平坐标系
    type(cartesian)::car1,car2
    type(topocentric)::top
    type(geodetic)::geo
    call cartogeo(car1,geo)
    top%n=-sin(geo%b)*cos(geo%l)*(car2%x-car1%x)-sin(geo%b)*sin(geo%l)*(car2%y-car1%y)+cos(geo%b)*(car2%z-car1%z)
    top%e=-sin(geo%l)*(car2%x-car1%x)+cos(geo%l)*(car2%y-car1%y)
    top%u=cos(geo%b)*cos(geo%l)*(car2%x-car1%x)+cos(geo%b)*sin(geo%l)*(car2%y-car1%y)+sin(geo%b)*(car2%z-car1%z)
    end subroutine
    
    subroutine toptocar(car1,car2,top)!站心地平坐标系转笛卡尔坐标系
    type(cartesian)::car1,car2
    type(topocentric)::top
    type(geodetic)::geo
    call cartogeo(car1,geo)
    car2%x=-sin(geo%b)*cos(geo%l)*top%n-sin(geo%l)*top%e+cos(geo%b)*cos(geo%l)*top%u+car1%x
    car2%y=-sin(geo%b)*sin(geo%l)*top%n+cos(geo%l)*top%e+cos(geo%b)*sin(geo%l)*top%u+car1%y
    car2%z=cos(geo%b)*top%n+sin(geo%b)*top%u+car1%z
    endsubroutine
    
    subroutine toptopop(top,pop)!站心地平坐标系转站心极坐标系
    type(topocentric)::top
    type(topopolar)::pop
    pop%s=sqrt(top%n*top%n+top%e*top%e+top%u*top%u)
    pop%e=asin(top%u/pop%s)
    pop%a=atan2(top%e,top%n)
    endsubroutine
    
    subroutine poptotop(top,pop)!站心极坐标系转站心地平坐标系
    type(topocentric)::top
    type(topopolar)::pop
    top%n=pop%s*cos(pop%e)*cos(pop%a)
    top%e=pop%s*cos(pop%e)*sin(pop%a)
    top%u=pop%s*sin(pop%e)
    endsubroutine
    
    
    end module
