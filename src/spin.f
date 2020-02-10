      subroutine spin(xxm,nclu,iangle)

c Randomly spin each collider

      implicit none
     
      integer mnat
      parameter(mnat=100) ! max number of atoms for each collider

      double precision xxm(3,mnat),ppm(3,mnat),tmp1,tmp2,tmp3
C
      logical iangle
      double precision rot(3,3),q0,q1,q2,q3
C
      integer i,j,nclu
      double precision s1,s2,s3,cs1,ss1,cs2,ss2,cs3,ss3,ranm(4)
C
      double precision pi
      parameter(pi=3.1415926536d0)

      if (.not.iangle) then
        s1=2.0d0
        s2=2.0d0
C     generate two random numbers(-1,1) with sum of square < 1
        do while(s1.ge.1.0d0)
          ranm(1) = 1.0-2.0*dble(rand())
          ranm(2) = 1.0-2.0*dble(rand())
          s1 = ranm(1)**2 + ranm(2)**2
        enddo
C     generate another two random numbers with sum of square < 1
        do while(s2.ge.1.0d0 .or. s2.eq.0.d0)
          ranm(3) = 1.0-2.0*dble(rand())
          ranm(4) = 1.0-2.0*dble(rand())
          s2 = ranm(3)**2 + ranm(4)**2
        enddo
        q0 = ranm(1)
        q1 = ranm(2)
        q2 = ranm(3)*dsqrt( (1-s1)/s2 )
        q3 = ranm(4)*dsqrt( (1-s1)/s2 )
C     Construction the transformation matrix
        rot(1,1)=q0**2 + q1**2 - q2**2 - q3**2
        rot(1,2)=2.d0*(q1*q2+q0*q3)
        rot(1,3)=2.d0*(q1*q3-q0*q2)
C
        rot(2,1)=2.d0*(q1*q2-q0*q3)
        rot(2,2)=q0**2 - q1**2 + q2**2 - q3**2
        rot(2,3)=2.d0*(q2*q3+q0*q1)
C
        rot(3,1)=2.d0*(q1*q3+q0*q2)
        rot(3,2)=2.d0*(q2*q3-q0*q1)
        rot(3,3)=q0**2 - q1**2 - q2**2 + q3**2
C     Change on angles s1(theta): 0-pi (angle to +z), s2(phi),
C     s3(posi): 0-2pi
      else
C     generate three random numbers(-1,1)
        ranm(1) = dble(rand())
        ranm(2) = dble(rand())
        ranm(3) = dble(rand())
C     Three angles: uniform on cos(theta)
        cs1 = 1.d0-2.d0*ranm(1)
        s1 = dacos(cs1)
        s2 = 2.d0*pi*(1.d0-2.d0*ranm(2))
        s3 = 2.d0*pi*(1.d0-2.d0*ranm(3))
C     cos and sin of angles s1,s2,s3
        ss1 = dsin(s1)
        cs2 = dcos(s2)
        ss2 = dsin(s2)
        cs3 = dcos(s3)
        ss3 = dsin(s3)
C
        rot(1,1)= cs2*cs3 - ss2*cs1*ss3
        rot(1,2)= ss2*cs3 + cs2*cs1*ss3
        rot(1,3)= ss1*ss3
C      
        rot(2,1)=-cs2*ss3 - ss2*cs1*cs3
        rot(2,2)=-ss2*ss3 + cs2*cs1*cs3
        rot(2,3)= ss1*cs3
C      
        rot(3,1)= ss2*ss1
        rot(3,2)=-cs2*ss1
        rot(3,3)= cs1
C      
      endif

      do i=1,nclu
        tmp1 = rot(1,1)*xxm(1,i)+rot(1,2)*xxm(2,i)+rot(1,3)*xxm(3,i)
        tmp2 = rot(2,1)*xxm(1,i)+rot(2,2)*xxm(2,i)+rot(2,3)*xxm(3,i)
        tmp3 = rot(3,1)*xxm(1,i)+rot(3,2)*xxm(2,i)+rot(3,3)*xxm(3,i)
        xxm(1,i) = tmp1
        xxm(2,i) = tmp2
        xxm(3,i) = tmp3
      enddo

      return
      end

