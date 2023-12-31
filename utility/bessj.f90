      FUNCTION bessj(n,x)
      implicit none
      INTEGER n,IACC
      REAL*8 bessj,x,BIGNO,BIGNI
      PARAMETER (IACC=100,BIGNO=1.d10,BIGNI=1.d-10)
!U    USES bessj0,bessj1
      INTEGER j,jsum,m
      REAL*8 ax,bj,bjm,bjp,sum,tox,bessj0,bessj1
      if(n<2) then
        if (n==1) bessj=bessj1(x)
        if (n==0) bessj=bessj0(x)
        return
      end if
      ax=abs(x)
      if(ax.eq.0.d0)then
        bessj=0.d0
      else if(ax.gt.dble(n))then
        tox=2.d0/ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
 11     continue
        bessj=bj
      else
        tox=2.d0/ax
        m=2*((n+int(dsqrt(dble(IACC*n))))/2)
        bessj=0.d0
        jsum=0
        sum=0.d0
        bjp=0.d0
        bj=1.d0
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
 12     continue
        sum=2.d0*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0.d0.and.mod(n,2).eq.1)bessj=-bessj
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
