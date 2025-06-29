      SUBROUTINE STRUC(X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/INTINIP/IINIP
      INTEGER NCOUNT
      DATA NCOUNT/0/
      SAVE NCOUNT

      IF(NCOUNT.EQ.0) THEN
      IINIP=0
      ENDIF
C      WRITE(*,1) IINIP
C    1 FORMAT(' ','IINIP =',I3)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      Q2=SCALE*SCALE
      IF(Q2.LT.qsqmin.OR.Q2.GT.qsqmax) PRINT 99
      if(X.LT.xmin.or.X.GT.xmax)       PRINT 98
  99  FORMAT('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  FORMAT('  WARNING:   X  VALUE IS OUT OF RANGE   ')


      IF(MODE.EQ.1) CALL MRST1(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.2) CALL MRST2(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.3) CALL MRST3(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.4) CALL MRST4(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)      
      IF(MODE.EQ.5) CALL MRST5(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.6) CALL MRST6(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.7) CALL MRST7(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.8) CALL MRST8(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.9) CALL MRST9(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.10) CALL GRVLO(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.11) CALL GRVHO(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.12) CALL GRV94LO(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.13) CALL GRV94HO(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.14) CALL GRV94DI(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.15) CALL GRV98(1,X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.16) CALL GRV98(2,X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF(MODE.EQ.17) CALL GRV98(3,X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     1GLU)
      IF (MODE.EQ.51) CALL CTEQ6L(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,
     1BOT,GLU)
      IF (MODE.EQ.52) CALL CTEQ6M(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,
     1BOT,GLU)
      IF (MODE.EQ.53) CALL CTEQ6D(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,
     1BOT,GLU)
      NCOUNT=NCOUNT+1
      RETURN
      END


C      subroutine mrstlo(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2001 LO parton            C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0201xxx                                  C
C                                                               C
C  There is 1 pdf set corresponding to mode = 1                 C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.220 C
C  corresponding to alpha_s(M_Z) of 0.130                       C
C  This set reads a grid whose first number is 0.02868          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
C      implicit real*8(a-h,o-z)
C      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
C      q2=q*q
C      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
C      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
C          if(mode.eq.1) then
C        call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
C      endif 
C  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
C  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
C      return
C      end

      subroutine mrst9(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='lo2002.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe5(nx,nq,xxl,qql,f1,cc1)
      call jeppe5(nx,nq,xxl,qql,f2,cc2)
      call jeppe5(nx,nq,xxl,qql,f3,cc3)
      call jeppe5(nx,nq,xxl,qql,f4,cc4)
      call jeppe5(nx,nq,xxl,qql,f6,cc6)
      call jeppe5(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe5(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe5(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe6(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe6(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe6(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe6(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe6(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe6(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe6(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe6(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 
c      subroutine jeppe5(nx,my,xx,yy,ff,cc)
c      implicit real*8(a-h,o-z)
c      dimension xx(nx),yy(my),ff(nx,my),ff1(nx,my),ff2(nx,my),
c     xff12(nx,my),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
c     xcl(16),cc(nx,my,4,4),iwt(16,16)

      subroutine jeppe5(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      PARAMETER(NNX=49,MMY=37)
      dimension xx(nx),yy(my),ff(nx,my),ff1(NNX,MMY),ff2(NNX,MMY),
     xff12(NNX,MMY),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     xcl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv2(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppe6(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      

      n=locx2(xx,nx,x)
      m=locx2(yy,my,y)

      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))

      z=0.
      do 1 l=4,1,-1
      z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     .       +cc(n,m,l,2))*u+cc(n,m,l,1)
    1 continue
      return
      end

      integer function locx2(xx,nx,x)
      implicit real*8(a-h,o-z)
      dimension xx(nx)
      if(x.le.xx(1)) then
      locx2=1
      return
      endif
      if(x.ge.xx(nx)) then 
      locx2=nx-1  
      return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
      jl=jm
      else
      ju=jm
      endif
      go to 1
    2 locx2=jl
      return
      end


      real*8 function  polderiv2(x1,x2,x3,y1,y2,y3)
      implicit real*8(a-h,o-z)
      polderiv2=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
     .(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

C      subroutine mrst2001(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2001 NLO parton           C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0110215                                  C
C                                                               C
C  There are 4 pdf sets corresponding to mode = 1, 2, 3, 4      C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.323 C
C  corresponding to alpha_s(M_Z) of 0.119                       C
C  This set reads a grid whose first number is 0.00927          C
C                                                               C
C  Mode=2 gives the set with Lambda(MSbar,nf=4) = 0.290         C
C  corresponding to alpha_s(M_Z) of 0.117                       C
C  This set reads a grid whose first number is 0.00953          C
C                                                               C
C  Mode=3 gives the set with Lambda(MSbar,nf=4) = 0.362         C
C  corresponding to alpha_s(M_Z) of 0.121                       C
C  This set reads a grid whose first number is 0.00889          C
C                                                               C
C  Mode=4 gives the set MRST2001J which gives better agreement  C
C  with the Tevatron inclusive jet data but has unattractive    C
C  gluon behaviour at large x (see discussion in paper)         C
C  This set has Lambda(MSbar,nf=4) = 0.353(alpha_s(M_Z) = 0.121 C 
C  This set reads a grid whose first number is 0.00826          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
C      implicit real*8(a-h,o-z)
C      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
C      q2=q*q
C      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
C      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
C          if(mode.eq.1) then
C        call mrst5(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
C      elseif(mode.eq.2) then
C        call mrst6(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
C      elseif(mode.eq.3) then
C        call mrst7(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
C      elseif(mode.eq.4) then
C        call mrst8(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
C      endif 
C  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
C  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
C      return
C      end

      subroutine mrst5(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='alf119.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe3(nx,nq,xxl,qql,f1,cc1)
      call jeppe3(nx,nq,xxl,qql,f2,cc2)
      call jeppe3(nx,nq,xxl,qql,f3,cc3)
      call jeppe3(nx,nq,xxl,qql,f4,cc4)
      call jeppe3(nx,nq,xxl,qql,f6,cc6)
      call jeppe3(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe3(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe3(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe4(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe4(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 
      subroutine mrst6(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='alf117.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe3(nx,nq,xxl,qql,f1,cc1)
      call jeppe3(nx,nq,xxl,qql,f2,cc2)
      call jeppe3(nx,nq,xxl,qql,f3,cc3)
      call jeppe3(nx,nq,xxl,qql,f4,cc4)
      call jeppe3(nx,nq,xxl,qql,f6,cc6)
      call jeppe3(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe3(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe3(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe4(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe4(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst7(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='alf121.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe3(nx,nq,xxl,qql,f1,cc1)
      call jeppe3(nx,nq,xxl,qql,f2,cc2)
      call jeppe3(nx,nq,xxl,qql,f3,cc3)
      call jeppe3(nx,nq,xxl,qql,f4,cc4)
      call jeppe3(nx,nq,xxl,qql,f6,cc6)
      call jeppe3(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe3(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe3(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe4(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe4(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst8(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='j121.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe3(nx,nq,xxl,qql,f1,cc1)
      call jeppe3(nx,nq,xxl,qql,f2,cc2)
      call jeppe3(nx,nq,xxl,qql,f3,cc3)
      call jeppe3(nx,nq,xxl,qql,f4,cc4)
      call jeppe3(nx,nq,xxl,qql,f6,cc6)
      call jeppe3(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe3(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe3(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe4(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe4(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe4(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine jeppe3(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),ff(nx,my),ff1(nx,my),ff2(nx,my),
     xff12(nx,my),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     xcl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv1(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv1(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv1(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppe4(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      

      n=locx1(xx,nx,x)
      m=locx1(yy,my,y)

      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))

      z=0.
      do 1 l=4,1,-1
      z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     .       +cc(n,m,l,2))*u+cc(n,m,l,1)
    1 continue
      return
      end

      integer function locx1(xx,nx,x)
      implicit real*8(a-h,o-z)
      dimension xx(nx)
      if(x.le.xx(1)) then
      locx1=1
      return
      endif
      if(x.ge.xx(nx)) then 
      locx1=nx-1  
      return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
      jl=jm
      else
      ju=jm
      endif
      go to 1
    2 locx1=jl
      return
      end


      real*8 function  polderiv1(x1,x2,x3,y1,y2,y3)
      implicit real*8(a-h,o-z)
      polderiv1=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
     .(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end


*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*    G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S     *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE :                *
*              M. GLUECK, E.REYA, A.VOGT : DO-TH 91/07            *
*               (PUBLISHED IN Z.PHYS. C53 (1992) 127)             *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE PARTON DISTRIBUTIONS   *
*   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   *
*   / HO) AND  1.E8 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   *
*   REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION IS NEG-   *
*   LIGIBLE, I.E. BELOW ABOUT 1.E-4, WERE EXCLUDED FROM THE FIT.  *
*                                                                 *
*              HEAVY QUARK THRESHOLDS  Q(H) = M(H) :              *
*         M(C)  =  1.5,  M(B)  =  4.5,  M(T)  =  100  GEV         *
*                                                                 *
*      CORRESPONDING LAMBDA(F) VALUES FOR F ACTIVE FLAVOURS :     *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,   LAMBDA(6)  =  0.082  GEV     *
*      HO :   LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131,   LAMBDA(6)  =  0.053  GEV     *
*                                                                 *
*   HO DISTRIBUTION REFER TO THE MS-BAR SCHEME OF BARDEEN ET AL.  *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS :
*
*    X   = MOMENTUM FRACTION
*    Q2  = SCALE Q**2 IN GEV**2
*
*...OUTPUT :
*
*    UDV = U(VAL) + D(VAL)
*    DV  = D(VAL)
*    GL  = GLUON
*    UDB = U(BAR) = D(BAR) = U(SEA) = D(SEA)
*    SB  = S = S(BAR)
*    CB  = C = C(BAR)
*    BB  = B = B(BAR)
*
*...LO PARAMETRIZATION :
*
      SUBROUTINE GRVLO(X,Q,UV,DV,UDB,DDB,SB,CB,BB,GL)
      IMPLICIT REAL*8(A - Z)
      MU2=0.25D0
      LAM2=0.232D0*0.232D0
      Q2=Q*Q
      S=DLOG(DLOG(Q2/LAM2)/DLOG(MU2/LAM2))
      S2=S*S
      S3=S2*S
*...X * (UV + DV) :
      NUD=0.663D0+0.191D0*S-0.041D0*S2+0.031D0*S3
      AKUD=0.326D0
      AGUD=-1.97D0+6.74D0*S-1.96D0*S2
      BUD=24.4D0-20.7D0*S+4.08D0*S2
      DUD=2.86D0+0.70D0*S-0.02D0*S2
      UDV=FV(X,NUD,AKUD,AGUD,BUD,DUD)
*...X * DV :
      ND=0.579D0+0.283D0*S+0.047D0*S2
      AKD=0.523D0-0.015D0*S
      AGD=2.22D0-0.59D0*S-0.27D0*S2
      BD=5.95D0-6.19D0*S+1.55D0*S2
      DD=3.57D0+0.94D0*S-0.16D0*S2
      DV=FV(X,ND,AKD,AGD,BD,DD)
      UV=UDV-DV
*...X * G :
      ALG=0.558D0
      BEG=1.218D0
      AKG=1.0D0-0.17D0*S
      BKG=0.0D0
      AGG=0.0D0+4.879D0*S-1.383D0*S2
      BGG=25.92D0-28.97D0*S+5.596D0*S2
      CG=-25.69D0+23.68D0*S-1.975D0*S2
      DG=2.537D0+1.718D0*S+0.353D0*S2
      EG=0.595D0+2.138D0*S
      ESG=4.066D0
      GL=FW(X,S,ALG,BEG,AKG,BKG,AGG,BGG,CG,DG,EG,ESG)
*...X * UBAR = X * DBAR :
      ALU=1.396D0
      BEU=1.331D0
      AKU=0.412D0-0.171D0*S
      BKU=0.566D0-0.496D0*S
      AGU=0.363D0
      BGU=-1.196D0
      CU=1.029D0+1.785D0*S-0.459D0*S2
      DU=4.696D0+2.109D0*S
      EU=3.838D0+1.944D0*S
      ESU=2.845D0
      UDB=FW(X,S,ALU,BEU,AKU,BKU,AGU,BGU,CU,DU,EU,ESU)
      DDB=UDB
*...X * SBAR = X * S :
      SS=0.0D0
      ALS=0.803D0
      BES=0.563D0
      AKS=2.082D0-0.577D0*S
      AGS=-3.055D0+1.024D0*S**0.67D0
      BS=27.4D0-20.0D0*S**0.154D0
      DS=6.22D0
      EST=4.33D0+1.408D0*S
      ESS=8.27D0-0.437D0*S
      SB=FWS(X,S,SS,ALS,BES,AKS,AGS,BS,DS,EST,ESS)
*...X * CBAR = X * C :
      SC=0.888D0
      ALC=1.01D0
      BEC=0.37D0
      AKC=0.0D0
      AGC=0.0D0
      BC=4.24D0-0.804D0*S
      DC=3.46D0+1.076D0*S
      EC=4.61D0+1.490D0*S
      ESC=2.555D0+1.961D0*S
      CB=FWS(X,S,SC,ALC,BEC,AKC,AGC,BC,DC,EC,ESC)
*...X * BBAR = X * B :
      SBO=1.351D0
      ALB=1.00D0
      BEB=0.51D0
      AKB=0.0D0
      AGB=0.0D0
      BBO=1.848D0
      DB=2.929D0+1.396D0*S
      EB=4.71D0+1.514D0*S
      ESB=4.02D0+1.239D0*S
      BB=FWS(X,S,SBO,ALB,BEB,AKB,AGB,BBO,DB,EB,ESB)
      RETURN
      END
*
*...HO PARAMERTRIZATION :
*
      SUBROUTINE GRVHO(X,Q,UV,DV,UDB,DDB,SB,CB,BB,GL)
      IMPLICIT REAL*8(A - Z)
      MU2=0.3D0
      LAM2=0.248D0 *0.248D0
      Q2=Q*Q
      S=DLOG(DLOG(Q2/LAM2)/DLOG(MU2/LAM2))
      DS=DSQRT(S)
      S2=S*S
      S3=S2*S
*...X * (UV + DV) :
      NUD=0.330D0+0.151D0*S-0.059D0*S2+0.027D0*S3
      AKUD=0.285D0
      AGUD=-2.28D0+15.73D0*S-4.58D0*S2
      BUD=56.7D0-53.6D0*S+11.21D0*S2
      DUD=3.17D0+1.17D0*S-0.47D0*S2+0.09D0*S3
      UDV=FV(X,NUD,AKUD,AGUD,BUD,DUD)
*...X * DV :
      ND= 0.459D0+0.315D0*DS+0.515D0*S
      AKD=0.624D0-0.031D0*S
      AGD=8.13D0-6.77D0*DS+0.46D0*S
      BD=6.59D0-12.83D0*DS+5.65D0*S
      DD=3.98D0+1.04D0*S-0.34D0*S2
      DV=FV(X,ND,AKD,AGD,BD,DD)
      UV=UDV-DV
*...X * G :
      ALG=1.128D0
      BEG=1.575D0
      AKG=0.323D0+1.653D0*S
      BKG=0.811D0+2.044D0*S
      AGG=0.0D0+1.963D0*S-0.519D0*S2
      BGG=0.078D0+6.24D0*S
      CG=30.77D0-24.19D0*S
      DG=3.188D0+0.720D0*S
      EG=-0.881D0+2.687D0*S
      ESG=2.466D0
      GL=FW(X,S,ALG,BEG,AKG,BKG,AGG,BGG,CG,DG,EG,ESG)
*...X * UBAR = X * DBAR :
      ALU=0.594D0
      BEU=0.614D0
      AKU=0.636D0-0.084D0*S
      BKU=0.0D0
      AGU=1.121D0-0.193D0*S
      BGU=0.751D0-0.785D0*S
      CU=8.57D0-1.763D0*S
      DU=10.22D0+0.668D0*S
      EU=3.784D0+1.280D0*S
      ESU=1.808D0+0.980D0*S
      UDB=FW(X,S,ALU,BEU,AKU,BKU,AGU,BGU,CU,DU,EU,ESU)
      DDB=UDB
*...X * SBAR = X * S :
      SS=0.0D0
      ALS=0.756D0
      BES=0.101D0
      AKS=2.942D0-1.016D0*S
      AGS=-4.60D0+1.167D0*S
      BS=9.31D0-1.324D0*S
      DS=11.49D0-1.198D0*S+0.053D0*S2
      EST=2.630D0+1.729D0*S
      ESS=8.12D0
      SB=FWS(X,S,SS,ALS,BES,AKS,AGS,BS,DS,EST,ESS)
*...X * CBAR = X * C :
      SC=0.820D0
      ALC=0.98D0
      BEC=0.0D0
      AKC=-0.625D0-0.523D0*S
      AGC=0.0D0
      BC=1.896D0+1.616D0*S
      DC=4.12D0+0.683D0*S
      EC=4.36D0+1.328D0*S
      ESC=0.677D0+0.679D0*S
      CB=FWS(X,S,SC,ALC,BEC,AKC,AGC,BC,DC,EC,ESC)
*...X * BBAR = X * B :
      SBO=1.297D0
      ALB=0.99D0
      BEB=0.0D0
      AKB=0.0D0-0.193D0*S
      AGB=0.0D0
      BBO=0.0D0
      DB=3.447D0+0.927D0*S
      EB=4.68D0+1.259D0*S
      ESB=1.892D0+2.199D0*S
      BB=FWS(X,S,SBO,ALB,BEB,AKB,AGB,BBO,DB,EB,ESB)
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION FV(X,N,AK,AG,B,D)
      IMPLICIT REAL*8(A - Z)
      DX=DSQRT(X)
      FV=N*X**AK*(1.0D0+AG*DX+B*X)*(1.0D0-X)**D
      RETURN
      END
*
*...FUNCTIONAL FORMS OF THE PARAMETRIZATIONS :
*
      DOUBLE PRECISION FUNCTION FW(X,S,AL,BE,AK,BK,AG,BG,
     1C,D,E,ES)
      IMPLICIT REAL*8(A - Z)
      LX=DLOG(1.0D0/X)
      FW=(X**AK*(AG+X*(BG+X*C))*LX**BK+S**AL
     1*DEXP(-E+DSQRT(ES*S**BE*LX)))*(1.0D0-X)**D
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION FWS(X,S,ST,AL,BE,AK,AG,B,
     1D,E,ES)
      IMPLICIT REAL*8(A - Z)
      DX=DSQRT(X)
      LX=DLOG(1.0D0/X)
      IF(S.LE.ST) THEN
      FWS=0.0D0
      ELSE
      FWS=(S-ST)**AL/LX**AK*(1.0D0+AG*DX+B*X)*(1.0D0-X)**D
     1*DEXP(-E+DSQRT(ES*S**BE*LX))
      ENDIF
      RETURN
      END



* file: grv94par.f
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*    G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S     *
*                                                                 *
*                         1994 UPDATE                             *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE                  *
*                   M. GLUECK, E.REYA, A.VOGT :                   *
*                   DO-TH 94/24  =  DESY 94-206                   *
*             (PUBLISHED IN Z. PHYS. C67 (1995) 433)              *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE EVOLVED PARTONS FOR    *
*        Q**2 / GEV**2  BETWEEN   0.4   AND  1.E6                 *
*             X         BETWEEN  1.E-5  AND   1.                  *
*   LARGE-X REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION   *
*   IS NEGLIGIBLY SMALL, WERE EXCLUDED FROM THE FIT.              *
*                                                                 *
*   HEAVY QUARK THRESHOLDS  Q(H) = M(H)  IN THE BETA FUNCTION :   *
*                   M(C)  =  1.5,  M(B)  =  4.5                   *
*   CORRESPONDING LAMBDA(F) VALUES IN GEV FOR  Q**2 > M(H)**2 :   *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,                                *
*      NLO :  LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131.                                *
*   THE NUMBER OF ACTIVE QUARK FLAVOURS IS  NF = 3  EVERYWHERE    *
*   EXCEPT IN THE BETA FUNCTION, I.E. THE HEAVY QUARKS C,B,...    *
*   ARE NOT PRESENT AS PARTONS IN THE Q2-EVOLUTION.               *
*   IF NEEDED, HEAVY QUARK DENSITIES CAN BE TAKEN FROM THE 1991   *
*   GRV PARAMETRIZATION.                                          *
*                                                                 *
*   NLO DISTRIBUTIONS ARE GIVEN IN MS-BAR FACTORIZATION SCHEME    *
*   (SUBROUTINE GRV94HO) AS WELL AS IN THE DIS SCHEME (GRV94DI),  *
*   THE LEADING ORDER PARAMETRIZATION IS PROVIDED BY "GRV94LO".   *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS :
*
*    X   = MOMENTUM FRACTION
*    Q2  = SCALE Q**2 IN GEV**2
*
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION) :
*
*    UV  = U(VAL) = U - U(BAR)
*    DV  = D(VAL) = D - D(BAR)
*    DEL = D(BAR) - U(BAR)
*    UDB = U(BAR) + D(BAR)
*    SB  = S = S(BAR)
*    GL  = GLUON
*
*...LO PARAMETRIZATION :
*

       SUBROUTINE GRV94LO (X,Q,UV,DV,USB,DSB,SB,CB,BB,GL)
       IMPLICIT REAL*8 (A - Z)
       MU2=0.23D0
       LAM2=0.2322D0*0.2322D0
       Q2=Q*Q
       S=DLOG(DLOG(Q2/LAM2)/DLOG(MU2/LAM2))
       DS=DSQRT(S)
       S2=S*S
       S3=S2*S
*...UV :
       NU=2.284D0+0.802D0*S+0.055D0*S2
       AKU=0.590D0-0.024D0*S
       BKU=0.131D0+0.063D0*S
       AU=-0.449D0-0.138D0*S-0.076D0*S2
       BU=0.213D0+2.669D0*S-0.728D0*S2
       CU=8.854D0-9.135D0*S+1.979D0*S2
       DU=2.997D0+0.753D0*S-0.076D0*S2
       UV=FFV(X,NU,AKU,BKU,AU,BU,CU,DU)
*...DV :
       ND=0.371D0+0.083D0*S+0.039D0*S2
       AKD=0.376D0
       BKD=0.486D0+0.062D0*S
       AD=-0.509D0+3.310D0*S-1.248D0*S2
       BD=12.41D0-10.52D0*S+2.267D0*S2
       CD=6.373D0-6.208D0*S+1.418D0*S2
       DD=3.691D0+0.799D0*S-0.071D0*S2
       DV=FFV(X,ND,AKD,BKD,AD,BD,CD,DD)
*...DEL :
       NE=0.082D0+0.014D0*S+0.008D0*S2
       AKE=0.409D0-0.005D0*S
       BKE=0.799D0+0.071D0*S
       AE=-38.07D0+36.13D0*S-0.656D0*S2
       BE=90.31D0-74.15D0*S+7.645D0*S2
       CE=0.0D0
       DE=7.486D0+1.217D0*S-0.159D0*S2
       DEL=FFV(X,NE,AKE,BKE,AE,BE,CE,DE)
*...UDB :
       ALX=1.451D0
       BEX=0.271D0
       AKX=0.410D0-0.232D0*S
       BKX=0.534D0-0.457D0*S
       AGX=0.890D0-0.140D0*S
       BGX=-0.981D0
       CX=0.320D0+0.683D0*S
       DX=4.752D0+1.164D0*S+0.286D0*S2
       EX=4.119D0+1.713D0*S
       ESX=0.682D0+2.978D0*S
       UDB=FFW(X,S,ALX,BEX,AKX,BKX,AGX,BGX,CX,DX,EX,ESX)
       USB=0.5D0*(UDB-DEL)
       DSB=0.5D0*(UDB+DEL)
       CB=0.0D0
       BB=0.0D0
*...SB :
       ALS=0.914D0
       BES=0.577D0
       AKS=1.798D0-0.596D0*S
       AS=-5.548D0+3.669D0*DS-0.616D0*S
       BS=18.92D0-16.73D0*DS+5.168D0*S
       DST=6.379D0-0.350D0*S+0.142D0*S2
       EST=3.981D0+1.638D0*S
       ESS=6.402D0
       SB=FFWS(X,S,ALS,BES,AKS,AS,BS,DST,EST,ESS)
**...GL :
       ALG=0.524D0
       BEG=1.088D0
       AKG=1.742D0-0.930D0*S
       BKG=-0.399D0*S2
       AG=7.486D0-2.185D0*S
       BG=16.69D0-22.74D0*S+5.779D0*S2
       CG=-25.59D0+29.71D0*S-7.296D0*S2
       DG=2.792D0+2.215D0*S+0.422D0*S2-0.104D0*S3
       EG=0.807D0+2.005D0*S
       ESG=3.841D0+0.316D0*S
       GL=FFW(X,S,ALG,BEG,AKG,BKG,AG,BG,CG,DG,EG,ESG)
       RETURN
       END
*
*...NLO PARAMETRIZATION (MS(BAR)) :
*
       SUBROUTINE GRV94HO(X,Q,UV,DV,USB,DSB,SB,CB,BB,GL)
       IMPLICIT REAL*8 (A - Z)
       MU2=0.34D0
       LAM2=0.248D0*0.248D0
       Q2=Q*Q
       S=DLOG(DLOG(Q2/LAM2)/DLOG(MU2/LAM2))
       DS=DSQRT(S)
       S2=S*S
       S3=S2*S
*...UV :
       NU=1.304D0+0.863D0*S
       AKU=0.558D0-0.020D0*S
       BKU=0.183D0*S
       AU=-0.113D0+0.283D0*S-0.321D0*S2
       BU=6.843D0-5.089D0*S+2.647D0*S2-0.527D0*S3
       CU=7.771D0-10.09D0*S+2.630D0*S2
       DU=3.315D0+1.145D0*S-0.583D0*S2+0.154D0*S3
       UV=FFV(X,NU,AKU,BKU,AU,BU,CU,DU)
*...DV :
       ND=0.102D0-0.017D0*S+0.005D0*S2
       AKD=0.270D0-0.019D0*S
       BKD=0.260D0
       AD=2.393D0+6.228D0*S-0.881D0*S2
       BD=46.06D0+4.673D0*S-14.98D0*S2+1.331D0*S3
       CD=17.83D0-53.47D0*S+21.24D0*S2
       DD=4.081D0+0.976D0*S-0.485D0*S2+0.152D0*S3
       DV=FFV(X,ND,AKD,BKD,AD,BD,CD,DD)
*...DEL :
       NE=0.070D0+0.042D0*S-0.011D0*S2+0.004D0*S3
       AKE=0.409D0-0.007D0*S
       BKE=0.782D0+0.082D0*S
       AE=-29.65D0+26.49D0*S+5.429D0*S2
       BE=90.20D0-74.97D0*S+4.526D0*S2
       CE=0.0D0
       DE=8.122D0+2.120D0*S-1.088D0*S2+0.231D0*S3
       DEL=FFV(X,NE,AKE,BKE,AE,BE,CE,DE)
*...UDB :
       ALX=0.877D0
       BEX=0.561D0
       AKX=0.275D0
       BKX=0.0D0
       AGX=0.997D0
       BGX=3.210D0-1.866D0*S
       CX=7.300D0
       DX=9.010D0+0.896D0*DS+0.222D0*S2
       EX=3.077D0+1.446D0*S
       ESX=3.173D0-2.445D0*DS+2.207D0*S
       UDB=FFW(X,S,ALX,BEX,AKX,BKX,AGX,BGX,CX,DX,EX,ESX)
       USB=0.5D0*(UDB-DEL)
       DSB=0.5D0*(UDB+DEL)
       CB=0.0D0
       BB=0.0D0
*...SB :
       ALS=0.756D0
       BES=0.216D0
       AKS=1.690D0+0.650D0*DS-0.922D0*S
       AS=-4.329D0+1.131D0*S
       BS=9.568D0-1.744D0*S
       DST=9.377D0+1.088D0*DS-1.320D0*S+0.130D0*S2
       EST=3.031D0+1.639D0*S
       ESS=5.837D0+0.815D0*S
       SB=FFWS(X,S,ALS,BES,AKS,AS,BS,DST,EST,ESS)
*...GL :
       ALG=1.014D0
       BEG=1.738D0
       AKG=1.724D0+0.157D0*S
       BKG=0.800D0+1.016D0*S
       AG=7.517D0-2.547D0*S
       BG=34.09D0-52.21D0*DS+17.47D0*S
       CG=4.039D0+1.491D0*S
       DG=3.404D0+0.830D0*S
       EG=-1.112D0+3.438D0*S-0.302D0*S2
       ESG=3.256D0-0.436D0*S
       GL=FFW(X,S,ALG,BEG,AKG,BKG,AG,BG,CG,DG,EG,ESG)
       RETURN
       END
*
*...NLO PARAMETRIZATION (DIS) :
*
       SUBROUTINE GRV94DI(X,Q,UV,DV,USB,DSB,SB,CB,BB,GL)
       IMPLICIT REAL*8 (A - Z)
       MU2=0.34D0
       LAM2=0.248D0*0.248D0
       Q2=Q*Q
       S=DLOG(DLOG(Q2/LAM2)/DLOG(MU2/LAM2))
       DS=DSQRT(S)
       S2=S*S
       S3=S2*S
*...UV :
       NU=2.484D0+0.116D0*S+0.093D0*S2
       AKU=0.563D0-0.025D0*S
       BKU=0.054D0+0.154D0*S
       AU=-0.326D0-0.058D0*S-0.135D0*S2
       BU=-3.322D0+8.259D0*S-3.119D0*S2+0.291D0*S3
       CU=11.52D0-12.99D0*S+3.161D0*S2
       DU=2.808D0+1.400D0*S-0.557D0*S2+0.119D0*S3
       UV=FFV(X,NU,AKU,BKU,AU,BU,CU,DU)
*...DV :
       ND=0.156D0-0.017D0*S
       AKD=0.299D0-0.022D0*S
       BKD=0.259D0-0.015D0*S
       AD=3.445D0+1.278D0*S+0.326D0*S2
       BD=-6.934D0+37.45D0*S-18.95D0*S2+1.463D0*S3
       CD=55.45D0-69.92D0*S+20.78D0*S2
       DD=3.577D0+1.441D0*S-0.683D0*S2+0.179D0*S3
       DV=FFV(X,ND,AKD,BKD,AD,BD,CD,DD)
*...DEL :
       NE=0.099D0+0.019D0*S+0.002D0*S2
       AKE=0.419D0-0.013D0*S
       BKE=1.064D0-0.038D0*S
       AE=-44.00D0+98.70D0*S-14.79D0*S2
       BE=28.59D0-40.94D0*S-13.66D0*S2+2.523D0*S3
       CE=84.57D0-108.8D0*S+31.52D0*S2
       DE=7.469D0+2.480D0*S-0.866D0*S2
       DEL=FFV(X,NE,AKE,BKE,AE,BE,CE,DE)
*...UDB :
       ALX=1.215D0
       BEX=0.466D0
       AKX=0.326D0+0.150D0*S
       BKX=0.956D0+0.405D0*S
       AGX=0.272D0
       BGX=3.794D0-2.359D0*DS
       CX=2.014D0
       DX=7.941D0+0.534D0*DS-0.940D0*S+0.410D0*S2
       EX=3.049D0+1.597D0*S
       ESX=4.396D0-4.594D0*DS+3.268D0*S
       UDB=FFW(X,S,ALX,BEX,AKX,BKX,AGX,BGX,CX,DX,EX,ESX)
       USB=0.5D0*(UDB-DEL)
       DSB=0.5D0*(UDB+DEL)
       CB=0.0D0
       BB=0.0D0
*...SB :
       ALS=0.175D0
       BES=0.344D0
       AKS=1.415D0-0.641D0*DS
       AS=0.580D0-9.763D0*DS+6.795D0*S-0.558D0*S2
       BS=5.617D0+5.709D0*DS-3.972D0*S
       DST=13.78D0-9.581D0*S+5.370D0*S2-0.996D0*S3
       EST=4.546D0+0.372D0*S2
       ESS=5.053D0-1.070D0*S+0.805D0*S2
       SB=FFWS(X,S,ALS,BES,AKS,AS,BS,DST,EST,ESS)
**...GL :
       ALG=1.258D0
       BEG=1.846D0
       AKG=2.423D0
       BKG=2.427D0+1.311D0*S-0.153D0*S2
       AG=25.09D0-7.935D0*S
       BG=-14.84D0-124.3D0*DS+72.18D0*S
       CG=590.3D0-173.8D0*S
       DG=5.196D0+1.857D0*S
       EG=-1.648D0+3.988D0*S-0.432D0*S2
       ESG=3.232D0-0.542D0*S
       GL=FFW(X,S,ALG,BEG,AKG,BKG,AG,BG,CG,DG,EG,ESG)
       RETURN
       END
*
*...FUNCTIONAL FORMS OF THE PARAMETRIZATIONS :
*
       REAL*8 FUNCTION FFV(X,N,AK,BK,A,B,C,D)
       IMPLICIT REAL*8 (A - Z)
       DX=DSQRT (X)
       FFV=N*X**AK*(1.0D0+A*X**BK+X*(B+C*DX))*(1.0D0-X)**D
       RETURN
       END
*
       REAL*8 FUNCTION FFW(X,S,AL,BE,AK,BK,A,B,C,D,E,ES)
       IMPLICIT REAL*8 (A - Z)
       LX=DLOG(1.0D0/X)
       FFW=(X**AK*(A+X*(B+X*C))*LX**BK+S**AL
     1 *DEXP(-E+DSQRT(ES*S**BE*LX)))*(1.0D0- X)**D
       RETURN
       END
*
       REAL*8 FUNCTION FFWS(X,S,AL,BE,AK,AG,B,D,E,ES)
       IMPLICIT REAL*8 (A - Z)
       DX=DSQRT(X)
       LX=DLOG(1.0D0/X)
       FFWS=S**AL/LX**AK*(1.0D0+AG*DX+B*X)*(1.0D0-X)**D
     1 *DEXP(-E+DSQRT(ES*S**BE*LX))
       RETURN
       END




* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                   *
*     G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S      *
*                                                                   *
*                          1998 UPDATE                              *
*                                                                   *
*                  For a detailed explanation see                   *
*                   M. Glueck, E. Reya, A. Vogt :                   *
*        hep-ph/9806404  =  DO-TH 98/07  =  WUE-ITP-98-019          *
*                  (To appear in Eur. Phys. J. C)                   *
*                                                                   *
*   This package contains subroutines returning the light-parton    *
*   distributions in NLO (for the MSbar and DIS schemes) and LO;    * 
*   the respective light-parton, charm, and bottom contributions    *
*   to F2(electromagnetic); and the scale dependence of alpha_s.    *
*                                                                   *
*   The parton densities and F2 values are calculated from inter-   *
*   polation grids covering the regions                             *
*         Q^2/GeV^2  between   0.8   and  1.E6 ( 1.E4 for F2 )      *
*            x       between  1.E-9  and   1.                       *
*   Any call outside these regions stops the program execution.     *
*                                                                   *
*   At Q^2 = MZ^2, alpha_s reads  0.114 (0.125) in NLO (LO); the    *
*   heavy quark thresholds, QH^2 = mh^2, in the beta function are   *
*            mc = 1.4 GeV,  mb = 4.5 GeV,  mt = 175 GeV.            *
*   Note that the NLO alpha_s running is different from GRV(94).    * 
*                                                                   *
*    Questions, comments etc to:  avogt@physik.uni-wuerzburg.de     *
*                                                                   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*
*
      SUBROUTINE GRV98 (ISET, X, Q, UV, DV, US, DS, SS, CH , BT ,GL)
*********************************************************************
*                                                                   *
*   THE PARTON ROUTINE.                                             *
*                                     __                            *
*   INPUT:   ISET =  1 (LO),  2 (NLO, MS), or  3 (NLO, DIS)         *
*            X  =  Bjorken-x        (between  1.E-9 and 1.)         *
*            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)         *
*                                                                   *
*   OUTPUT:  UV = u - u(bar),  DV = d - d(bar),  US = u(bar),       *
*            DS = d(bar),  SS = s = s(bar),  GL = gluon.            *
*            Always x times the distribution is returned.           *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINIP / IINIP , and the     *
*            integer variable  IINIP  has always to be zero when    *
*            GRV98PA is called for the first time or when  ISET     *
*            has been changed.                                      *
*                                                                   *
*   GRIDS:   1. grv98lo.grid, 2. grv98nlm.grid, 3. grv98nld.grid,   *
*            (1+1809 lines with 6 columns, 4 significant figures)   *
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NPART=6, NX=68, NQ=27, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XDEF(NX,NQ), XUDF(NX,NQ),
     1          XSF(NX,NQ), XGF(NX,NQ), PARTON (NPART,NQ,NX-1), 
     2          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      CHARACTER*80 LINE
      COMMON / INTINIP / IINIP
      SAVE XUVF, XDVF, XDEF, XUDF, XSF, XGF, NA, ARRF
      Q2=Q*Q
      CH=0.0D0
      BT=0.0D0
*
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8E0, 
     1           1.0E0, 1.3E0, 1.8E0, 2.7E0, 4.0E0, 6.4E0,
     2           1.0E1, 1.6E1, 2.5E1, 4.0E1, 6.4E1,
     3           1.0E2, 1.8E2, 3.2E2, 5.7E2,
     4           1.0E3, 1.8E3, 3.2E3, 5.7E3,
     5           1.0E4, 2.2E4, 4.6E4,
     6           1.0E5, 2.2E5, 4.6E5, 
     7           1.E6 /
       DATA XB / 1.0E-9, 1.8E-9, 3.2E-9, 5.7E-9, 
     A           1.0E-8, 1.8E-8, 3.2E-8, 5.7E-8, 
     B           1.0E-7, 1.8E-7, 3.2E-7, 5.7E-7, 
     C           1.0E-6, 1.4E-6, 2.0E-6, 3.0E-6, 4.5E-6, 6.7E-6,
     1           1.0E-5, 1.4E-5, 2.0E-5, 3.0E-5, 4.5E-5, 6.7E-5,
     2           1.0E-4, 1.4E-4, 2.0E-4, 3.0E-4, 4.5E-4, 6.7E-4,
     3           1.0E-3, 1.4E-3, 2.0E-3, 3.0E-3, 4.5E-3, 6.7E-3,
     4           1.0E-2, 1.4E-2, 2.0E-2, 3.0E-2, 4.5E-2, 0.06, 0.08,
     5           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     6           0.3, 0.325, 0.35, 0.375, 0.4,  0.45, 0.5, 0.55,
     7           0.6, 0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1. /
*
*...CHECK OF X AND Q2 VALUES : 
      IF ( (X.LT.0.99D-9) .OR. (X.GT.1.D0) ) THEN
         WRITE(6,91) 
  91     FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
         STOP
      ENDIF
      IF ( (Q2.LT.0.799) .OR. (Q2.GT.1.01E6) ) THEN
         WRITE(6,92) 
  92     FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
         STOP
      ENDIF
      IF (IINIP .NE. 0) GOTO 16
*
*...INITIALIZATION, IF REQUIRED :
*
*    SELECTION AND READING OF THE GRID : 
*    (COMMENT: FIRST NUMBER IN THE FIRST LINE OF THE GRID)
      IF (ISET .EQ. 1) THEN
        OPEN (11,FILE='grv98lo.grid',STATUS='old')   ! 7.332E-05
      ELSE IF (ISET .EQ. 2) THEN
        OPEN (11,FILE='grv98nlm.grid',STATUS='old')  ! 1.015E-04
      ELSE IF (ISET .EQ. 3) THEN
        OPEN (11,FILE='grv98nld.grid',STATUS='old')  ! 1.238E-04
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'NO OR INVALID PARTON SET CHOICE')
        STOP
      END IF
      IINIP = 1
      READ(11,89) LINE
  89  FORMAT(A80)
      DO 15 M = 1, NX-1 
      DO 15 N = 1, NQ
      READ(11,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1            PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M) 
  90  FORMAT (6(1PE10.3))
  15  CONTINUE
      CLOSE(11)
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0V = XB(IX)**0.5 
        XB0S = XB(IX)**(-0.2) 
        XB1 = 1.-XB(IX)
        XUVF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0V)
        XDVF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0V)
        XDEF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0V) 
        XUDF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0S)
        XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**7 * XB0S)
        XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0S)
  20  CONTINUE
        XUVF(NX,IQ) = 0.E0
        XDVF(NX,IQ) = 0.E0
        XDEF(NX,IQ) = 0.E0
        XUDF(NX,IQ) = 0.E0
        XSF(NX,IQ)  = 0.E0
        XGF(NX,IQ)  = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
*
*...CONTINUATION, IF INITIALIZATION WAS DONE PREVIOUSLY.
*
  16  CONTINUE
*
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      X1 = 1.- X
      XV = X**0.5
      XS = X**(-0.2)
      UV = FFINT(NARG,XT,NA,ARRF,XUVF) * X1**3 * XV
      DV = FFINT(NARG,XT,NA,ARRF,XDVF) * X1**4 * XV
      DE = FFINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
      UD = FFINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
      IF (DE. GT. UD) DE = 0.D0
      US = 0.5 * (UD - DE)
      DS = 0.5 * (UD + DE)
      SS = FFINT(NARG,XT,NA,ARRF,XSF)  * X1**7 * XS
      GL = FFINT(NARG,XT,NA,ARRF,XGF)  * X1**5 * XS 
*
 60   RETURN
      END
*
*
*
*
      FUNCTION FFINT(NARG,ARG,NENT,ENT,TABLE)
*********************************************************************
*                                                                   *
*   THE INTERPOLATION ROUTINE (CERN LIBRARY ROUTINE E104)           *
*                                                                   *
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FFINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FFINT=FFINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END


C      subroutine mrstnnlo(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the MRST 2002 NNLO parton distributionsC
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0201127                                  C
C                                                               C
C  There are 4 pdf sets corresponding to mode = 1, 2, 3, 4      C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.235 C
C  corresponding to alpha_s(M_Z) of 0.1155                      C
C  This set is the `average' of the slow and fast evolutions    C 
C  This set reads a grid whose first number is 0.00725          C
C                                                               C
C  Mode=2 gives the set with Lambda(MSbar,nf=4) = 0.235         C
C  corresponding to alpha_s(M_Z) of 0.1155                      C
C  This set is the fast evolution                               C 
C  This set reads a grid whose first number is 0.00734          C
C                                                               C
C  Mode=3 gives the set with Lambda(MSbar,nf=4) = 0.235         C
C  corresponding to alpha_s(M_Z) of 0.1155                      C
C  This set is the slow evolution                               C 
C  This set reads a grid whose first number is 0.00739          C
C                                                               C
C  Mode=4 gives the set MRSTNNLOJ which gives better agreement  C
C  with the Tevatron inclusive jet data but has unattractive    C
C  gluon behaviour at large x (see discussion in paper)         C
C  This set has Lambda(MSbar,nf=4) = 0.267(alpha_s(M_Z) =0.1180 C 
C  This set reads a grid whose first number is 0.00865          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
C      implicit real*8(a-h,o-z)
C      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
C      q2=q*q
C      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
C      if(x.lt.xmin.or.x.gt.xmax)       print 98,x
C          if(mode.eq.1) then
C        call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
C      elseif(mode.eq.2) then
C        call mrst2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
C      elseif(mode.eq.3) then
C        call mrst3(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
C      elseif(mode.eq.4) then
C        call mrst4(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
C      endif 
C  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
C  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
C      return
C      end

      subroutine mrst1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='vnvalf1155.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 
      subroutine mrst2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='vnvalf1155a.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst3(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='vnvalf1155b.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst4(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='vnvalf1180j.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine jeppe1(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      PARAMETER(NNX=49,MMY=37)
      dimension xx(nx),yy(my),ff(nx,my),ff1(NNX,MMY),ff2(NNX,MMY),
     xff12(NNX,MMY),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     xcl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppe2(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      

      n=locx(xx,nx,x)
      m=locx(yy,my,y)

      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))

      z=0.
      do 1 l=4,1,-1
      z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     .       +cc(n,m,l,2))*u+cc(n,m,l,1)
    1 continue
      return
      end

      integer function locx(xx,nx,x)
      implicit real*8(a-h,o-z)
      dimension xx(nx)
      if(x.le.xx(1)) then
      locx=1
      return
      endif
      if(x.ge.xx(nx)) then 
      locx=nx-1  
      return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
      jl=jm
      else
      ju=jm
      endif
      go to 1
    2 locx=jl
      return
      end


      real*8 function  polderiv(x1,x2,x3,y1,y2,y3)
      implicit real*8(a-h,o-z)
      polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
     .(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCJS
C    Iset is the set label; in this version, Iset = 1, 2, 3 
C                           correspond to the following CTEQ global fits:
C ---------------------------------------------------------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
C   3    CTEQ6L   Leading Order           0.118     326   226    cteq6l.tbl
C ---------------------------------------------------------------------------
      SUBROUTINE CTEQ6L(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU) 
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NCOUNT
      DATA NCOUNT/0/
      SAVE NCOUNT
c
      call SetCtq6(3)
c
      CALL CTEQ6(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
      IF(NCOUNT.EQ.0)THEN
      WRITE(6,531) X,SCALE,UPV,DNV,USEA
 531  FORMAT(5(2X,D10.5))
      WRITE(6,532) DSEA,STR,CHM,BOT,GLU
 532  FORMAT(5(2X,D10.5))
      ENDIF
      NCOUNT=NCOUNT+1
      return
      end

      SUBROUTINE CTEQ6M(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU) 
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NCOUNT
      DATA NCOUNT/0/
      SAVE NCOUNT
c
      call SetCtq6(1)
c
      CALL CTEQ6(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
      IF(NCOUNT.EQ.0)THEN
      WRITE(6,533) X,SCALE,UPV,DNV,USEA
 533  FORMAT(5(2X,D10.5))
      WRITE(6,534) DSEA,STR,CHM,BOT,GLU
 534  FORMAT(5(2X,D10.5))
      ENDIF
      NCOUNT=NCOUNT+1
      return
      end

      SUBROUTINE CTEQ6D(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU) 
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NCOUNT
      DATA NCOUNT/0/
      SAVE NCOUNT
c
      call SetCtq6(2)
c
      CALL CTEQ6(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
      IF(NCOUNT.EQ.0)THEN
      WRITE(6,535) X,SCALE,UPV,DNV,USEA
 535  FORMAT(5(2X,D10.5))
      WRITE(6,536) DSEA,STR,CHM,BOT,GLU
 536  FORMAT(5(2X,D10.5))
      ENDIF
      NCOUNT=NCOUNT+1
      return
      end
      
      
       SUBROUTINE CTEQ6(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU) 
       implicit real*8 (a-h,o-z)
C
       Q=SCALE
       xsave=X
       qsave=Q
       U =         X * Ctq6Pdf(1,X,Q)
       D =         X * Ctq6Pdf(2,X,Q)
       USEA =      X * Ctq6Pdf(-1,X,Q)
       DSEA =      X * Ctq6Pdf(-2,X,Q)
       STR =       X * Ctq6Pdf(3,X,Q)
       CHM =       X * Ctq6Pdf(4,X,Q)
       BOT =       X * Ctq6Pdf(5,X,Q)
       GLU  =      X * Ctq6Pdf(0,X,Q)
      UPV=U-USEA
      DNV=D-DSEA
      X=xsave
      Q=qsave
      return
      end

C============================================================================
C                CTEQ Parton Distribution Functions: Version 6.0
C                             January 24, 2002
C
C   Ref: "New Generation of Parton Distributions with
C         Uncertainties from Global QCD Analysis"
C   By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       hep-ph/0201195
C
C   This package contains 3 standard sets of CTEQ6 PDF's and 40 up/down sets
C   with respect to CTEQ6M PDF's. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
C     ------------------------------
C   1xx  CTEQ6M1xx  +/- w.r.t. CTEQ6M     0.118     326   226    cteq6m1xx.tbl
C    (where xx=01--40)
C ---------------------------------------------------------------------------
C   ** ALL fits are obtained by using the same coupling strength \alpha_s(Mz)=0.118;
C   and the NLO running \alpha_s formula.  For the LO fit, the evolution of the PDF
C   and the hard cross sections are calculated at LO.  More detailed discussions are
C   given in hep-ph/0201195.
C
C   The table grids are generated for 10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV).
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq6(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   Pumplin@pa.msu.edu or Tung@pa.msu.edu.
C
C===========================================================================

      Function Ctq6Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq6Pdf: ', X
        Stop
      Endif
      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq6Pdf: ', Q
        Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq6Pdf: '
     >              , Iparton
         Endif
         Ctq6Pdf = 0D0
         Return
      Endif

      Ctq6Pdf = PartonX6 (Iparton, X, Q)
      if(Ctq6Pdf.lt.0.D0)  Ctq6Pdf = 0.D0

      Return

C                             ********************
      End

      Subroutine SetCtq6 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax0=3)
      Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
      Data (Flnm(I), I=1,Isetmax0)
     > / 'cteq6m', 'cteq6d', 'cteq6l' /
      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,101,140/
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         IU= NextUn()
         If (Iset.ge.Isetmin0 .and. Iset.le.Isetmax0) Then
            Tablefile=Flnm(Iset)//'.tbl'
         Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
            write(nn,'(I3)') Iset
            Tablefile=Flnm(1)//nn//'.tbl'
         Else
            Print *, 'Invalid Iset number in SetCtq6 :', Iset
            Stop
         Endif
         Open(IU, File=Tablefile, Status='OLD', Err=100)
 21      Call RReadTbl (IU)
         Close (IU)
         Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq6!!'
      Stop
C                             ********************
      End

      Subroutine RReadTbl (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 96, MXQ = 20, MXF = 5)
      PARAMETER (MXPQX = (MXF + 3) * MXQ * MXX)
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)

      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         TV(Iq) = Log(Log (TV(Iq) /Al))
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage,
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      Parameter (MXX = 96, MXQ = 20, MXF = 5)
      Parameter (MXQX= MXQ * MXX,   MXPQX = MXQX * (MXF+3))

      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Save ientry,xvpow

c store the powers used for interpolation on first call...
      if(ientry .eq. 0) then
         ientry = 1

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf


c get the pdf function values at the lattice points...

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint (XVpow(0), Fij(1), 4, ss, Fx, Dfx)

         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)

        Fvec(it) = Fx

       Else
C                       for all interior points, use Jon's in-line function
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4

         Fvec(it) = (const5*(Upd(J1)-g1)
     &               + const6*(Upd(J1+3)-g4)
     &               + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint (TV(0), Fvec(1), 4, tt, ff, Dfq)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C                                       ********************
      End
      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        ADAPTED FROM "NUMERICAL RECIPES"
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


