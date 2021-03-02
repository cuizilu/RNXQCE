module clockbias
    contains 
    subroutine calclock(a0,a1,a2,tsminustoc,deltas)
    real*8::a0,a1,a2,tsminustoc,deltas
    deltas=a0+a1*tsminustoc+a2*tsminustoc**2
    end subroutine
    end module