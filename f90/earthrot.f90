module earthrot
    use parameters
    use coordinate
    contains
    subroutine calearth(ts,tr,x1,x2)
    type(cartesian)::x1,x2
    real*8::ts,tr
    x2%x=cos(omegae*(tr-ts))*x1%x+sin(omegae*(tr-ts))*x1%y
    x2%y=-sin(omegae*(tr-ts))*x1%x+cos(omegae*(tr-ts))*x1%y
    x2%z=x1%z
    endsubroutine
    end module
