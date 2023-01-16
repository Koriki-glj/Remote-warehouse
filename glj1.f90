program main
implicit none
!sheng ming
!边界和内部网格
real(kind=8):: ys(2,2,200),zs(2,2,200),ya(300,300,11),za(300,300,11)
real(kind=8),allocatable :: xlast(:,:,:,:),ylast(:,:,:,:),zlast(:,:,:,:) !可分配内存数组
!分布函数
real(kind=8) :: s(100),sb(2,250)
!chi cun 
real(kind=8):: r,theta,theta1,theta2,theta3,theta4
real(kind=8),parameter:: pi=4.0*atan(1.0)   !等于180度
real(kind=8) :: zd2,hy2,sz2,temp1,temp2,temp3,xjl
!网格正交性控制和其他循环变量
real(kind=8) :: p,q
integer :: n,j,m,k,iter,i
!die dai shu
integer :: nx,ny1,nz1,ny2,nz2,ny3,nz3,ny4,nz4,ny5,nz51,nz52,ny6,nz61,nz62,ny71,ny72,ny73,nz71,nz72,nz73,ny8,nz81,nz82,ny9,nz91,nz92

!chi cun ding yi
zd2=5.0
r=5.0
theta=r/2.0
hy2=2.0*r
sz2=r+zd2
theta1=asin(theta/r)  !等于30度
theta2=pi    !等于180度
theta3=1.25*pi    !等于225度
theta4=1.75*pi   !等于315度
temp1=4.0
temp2=3.0
temp3=2.0
nx=91
ny1=31
nz1=31
ny2=31
nz2=31
ny3=31
nz3=31
ny4=31
nz4=31
ny5=31
nz51=31
nz52=31
ny6=31
nz61=31
nz62=31
ny71=31
nz71=21
ny72=31
nz72=21
ny73=31
nz73=21
ny8=31
nz81=31
nz82=31
ny9=31
nz91=31
nz92=31
xjl=50.0

!区域1半圆
m=1
n=int(ny1/2)+1
p=0.5    !从密到疏再到密
q=2.0
call strech(n,p,q,s)
do j=1,n-1
 s(n+j)=1.0-s(n-j)
end do
 do 1  j=1,ny1
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)     !sb(1，j)是分布函数
end if
1 continue
  
 ! write(*,*) (sb(1,j),j=1,ny1)
 ! pause
  
do j=1,ny1
	  theta=THETA3+(THETA4-THETA3)*sb(1,j)
	  ys(1,1,j)=   R*DCOS(THETA)             
	  zs(1,1,j)=        R*DSIN(THETA)  
	
	  ys(1,2,j)=         (-r)+2.0*r*sb(1,j)
	  zs(1,2,j)=        0.0
end do

n=int(nz1/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz1
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
do k=1,nz1
	   zs(2,1,k)=           R*DSIN(THETA3-(THETA3-THETA2)*Sb(2,k))  
	  ys(2,1,k)=      -DSQRT(R**2-zs(2,1,k)**2)		
	  zs(2,2,k)=           zs(2,1,k)
	  ys(2,2,k)=            -ys(2,1,k)
	  
end do

  !内部网格
        do  j=1,ny1
	  do  k=1,nz1   
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny1)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny1)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny1)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny1)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
		end do
		
		!区域2
		m=2
		n=int(ny2/2)+1
p=0.5
q=2.0 
call strech(n,p,q,s)
do j=1,n-1
 s(n+j)=1.0-s(n-j)
end do

 do   j=1,ny2
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)    
end if
end do


do j=1,ny2
	
		  ys(1,1,j)= -r*dcos(0.25*pi+0.5*pi*sb(1,j))
	  zs(1,1,j)= -r*dsin(0.25*pi+0.5*pi*sb(1,j))  
	  ys(1,2,j)=  -r+2*r*sb(1,j)
	  zs(1,2,j)=     -sz2
end do

n=int(nz2/2)+1
p=0.5
q=2.0 
call strech(n,p,q,s)
do k=1,n-1
 s(n+k)=1.0-s(n-k)
end do

 do   k=1,nz2
 if(k<=n)  then
 s(k)=s(k)/2.0
 sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0  
sb(2,k)=s(k)
end if
end do

do k=1,nz2
  zs(2,1,k)= -r*dsin(0.25*pi)+(r*dsin(0.25*pi)-sz2)*Sb(2,k) 
	  ys(2,1,k)= -r*dsin(0.25*pi)+(r*dsin(0.25*pi)-zd2)*Sb(2,k)
	  zs(2,2,k)=zs(2,1,k)
	  ys(2,2,k)=-ys(2,1,k)
	  
end do

  !内部网格
        do  j=1,ny2
	  do  k=1,nz2 
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny2)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny2)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny2)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny2)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
		end do
!区域3-半圆		
m=3
n=int(ny3/2)+1
p=0.5
q=2.0
call strech(n,p,q,s)


do j=1,n-1
 s(n+j)=1.0-s(n-j)
end do

 do   j=1,ny3
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)    
end if
 end do

do j=1,ny3

theta3=1.25*pi    !等于225度
theta4=1.75*pi   !等于315度
	  theta=THETA3+(THETA4-THETA3)*sb(1,j)
	  ys(1,1,j)=   R*DCOS(THETA)             
	  zs(1,1,j)=        R*DSIN(THETA)  -sz2
	
	  ys(1,2,j)=         (-r)+2.0*r*sb(1,j)
	  zs(1,2,j)=        -sz2
end do


n=int(nz3/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz3
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
do k=1,nz3
	theta2=pi    !等于180度
theta3=1.25*pi    !等于225度

	   ys(2,1,k)=           R*Dcos(THETA3-(THETA3-THETA2)*Sb(2,k))
	  zs(2,1,k)=      -DSQRT(R**2-yS(2,1,k)**2)		 -sz2 
	  zs(2,2,k)=           zs(2,1,k)
	  ys(2,2,k)=            -ys(2,1,k)
	  
end do

  !内部网格
        do  j=1,ny3
	  do  k=1,nz3 
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny3)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny3)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny3)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny3)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
		end do

!区域-4
m=4
n=int(ny4/2)+1
p=0.5
q=2.0 
call strech(n,p,q,s)
do j=1,n-1
 s(n+j)=1.0-s(n-j)
end do

 do   j=1,ny4
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)    
end if
 end do
 
do j=1,ny4
	
		  ys(1,2,j)= -r*dcos(0.25*pi+0.5*pi*sb(1,j))
	  zs(1,2,j)= -r*dsin(0.25*pi+0.5*pi*sb(1,j))  
	  ys(1,1,j)= -r*dcos(pi/4.0)+2*r*dcos(pi/4.0)*sb(1,j)
	  zs(1,1,j)=-sz2-r-temp2
end do
n=int(nz4/2)+1
p=0.5
q=2.0 
call strech(n,p,q,s)
do k=1,n-1
 s(n+k)=1.0-s(n-k)
end do

 do   k=1,nz4
 if(k<=n)  then
 s(k)=s(k)/2.0
 sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0  
sb(2,k)=s(k)
end if
end do

do k=1,nz4
  zs(2,1,k)=(-sz2-r-temp2)+(r+temp2-r*dsin(pi/4.0)) *sb(2,k)
	  ys(2,1,k)= -r*dcos(pi/4.0)
	  zs(2,2,k)=zs(2,1,k)
	  ys(2,2,k)=-ys(2,1,k)
	  
end do

  !内部网格
        do  j=1,ny4
	  do  k=1,nz4
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny4)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny4)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny4)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny4)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
		end do

		
!区域-5
m=5
n=int(ny5/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
		end do
 do   j=1,ny5
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)   
end if
 end do

	
do j=1,ny5

	  ys(1,1,j)=   r*dcos(pi/4.0)+(r-r*dcos(pi/4.0)+temp1)*sb(1,j)   
	  zs(1,1,j)=      -sz2-r-temp2 
	
	  ys(1,2,j)=   r+temp1*sb(1,j) 
	  zs(1,2,j)=   -sz2     
end do	
	
	n=int(nz51/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz51
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
		
		
	do k=1,nz51

	   ys(2,1,k)=  r*dcos(pi/4.0)       
	  zs(2,1,k)=    (-sz2-r-temp2)+(r+temp2-r*dsin(pi/4.0)) *sb(2,k)
	end do
	n=int(nz52/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz52
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
	
	
	do k=1,nz52
 ys(2,1,nz51+k-1)=	   r*dcos(pi/4.0-pi/4.0*sb(2,k))  
 zs(2,1,nz51+k-1)=-sz2-r*dsin(pi/4.0-pi/4.0*sb(2,k))  
	end do
	       
   	n=int((nz51+nz52-1)/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz51+nz52-1
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do

	  do k=1,nz51+nz52-1
	  zs(2,2,k)=    zs(1,1,ny5)+(r+temp2)*sb(2,k)
	  ys(2,2,k)=	r+temp1
	  
	end do
	
	        do  j=1,ny5
	  do  k=1,nz51+nz52-1 
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny5)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny5)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny5)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny5)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
			end do
			
!区域-6
m=6
n=int(ny6/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
		end do
 do   j=1,ny6
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)   
end if
 end do

	
do j=1,ny6

	  ys(1,1,j)=   r+temp1*sb(1,j)   
	  zs(1,1,j)=      -sz2
	  ys(1,2,j)=   r+temp1*sb(1,j) 
	  zs(1,2,j)=  0.0    
end do	
	
	n=int(nz61/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz61
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
		
		
	do k=1,nz61

	   ys(2,1,k)=  r+(r*dcos(pi/4.0)-r)     * sb(2,k)
	  zs(2,1,k)=  -sz2+(-r*dsin(pi/4.0)+sz2)*sb(2,k)
	end do
	n=int(nz62/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz62
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
	
	
	do k=1,nz62
 ys(2,1,nz61+k-1)=	   r*dcos(pi/4.0-pi/4.0*sb(2,k))  
 zs(2,1,nz61+k-1)=-r*dsin(pi/4.0-pi/4.0*sb(2,k))  
	end do
	       
   	n=int((nz61+nz62-1)/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz61+nz62-1
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do

	  do k=1,nz61+nz62-1
	  zs(2,2,k)=    zs(1,1,ny6)+(sz2)*sb(2,k)
	  ys(2,2,k)=	r+temp1
	  
	end do
	
	        do  j=1,ny6
	  do  k=1,nz61+nz62-1
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny6)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny6)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny6)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny6)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
			end do
	
		!区域-7
m=7
		n=int(ny71/2)+1
p=0.5
q=2.0 
call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
		end do
 do   j=1,ny71
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)   
end if
 end do


do j=1,ny71
	
	  ys(1,1,j)=        r+temp1*sb(1,j)
	  zs(1,1,j)=       0.0


	  ys(1,2,j)=     ys(1,1,j)       
	  zs(1,2,j)=     temp3 
end do

		n=int(nz71/2)+1
p=1.0
q=2.0 
call strech(n,p,q,s)
do k=1,n-1
 s(n+k)=1.0-s(n-k)
end do

 do   k=1,nz71
 if(k<=n)  then
 s(k)=s(k)/2.0
 sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0  
sb(2,k)=s(k)
end if
end do

do k=1,nz71
	  ys(2,1,k)=r
	zs(2,1,k)=   0.0+temp3*sb(2,k)     
	 
	      ys(2,2,k)=    r+temp1                   
	  zs(2,2,k)=         zs(2,1,k)
	
	  
end do

  !内部网格
        do  j=1,ny71
	  do  k=1,nz71 
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny71)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny71)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny71)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny71)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
		end do
		
m=10
		n=int(ny72/2)+1
p=0.5
q=2.0 
call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
		end do
 do   j=1,ny72
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)   
end if
 end do


do j=1,ny72
	
	  ys(1,1,j)=     -r+2*   r*sb(1,j)
	  zs(1,1,j)=       0.0


	  ys(1,2,j)=     ys(1,1,j)       
	  zs(1,2,j)=     temp3 
end do

		n=int(nz72/2)+1
p=1.0
q=2.0 
call strech(n,p,q,s)
do k=1,n-1
 s(n+k)=1.0-s(n-k)
end do

 do   k=1,nz72
 if(k<=n)  then
 s(k)=s(k)/2.0
 sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0  
sb(2,k)=s(k)
end if
end do

do k=1,nz72
	  ys(2,1,k)=-r
	zs(2,1,k)=   0.0+temp3*sb(2,k)     
	 
	      ys(2,2,k)=    r                  
	  zs(2,2,k)=         zs(2,1,k)
	
	  
end do

  !内部网格
        do  j=1,ny72
	  do  k=1,nz72
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny72)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny72)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny72)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny72)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
		end do
		
m=11
		n=int(ny73/2)+1
p=0.5
q=2.0 
call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
		end do
 do   j=1,ny73
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)   
end if
 end do


do j=1,ny73
	
	  ys(1,1,j)=      -  (r+temp1)+temp1*sb(1,j)
	  zs(1,1,j)=       0.0


	  ys(1,2,j)=     ys(1,1,j)       
	  zs(1,2,j)=     temp3 
end do

		n=int(nz73/2)+1
p=1.0
q=2.0 
call strech(n,p,q,s)
do k=1,n-1
 s(n+k)=1.0-s(n-k)
end do

 do   k=1,nz73
 if(k<=n)  then
 s(k)=s(k)/2.0
 sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0  
sb(2,k)=s(k)
end if
end do

do k=1,nz73
	  ys(2,1,k)=-r-temp1
	zs(2,1,k)=   0.0+temp3*sb(2,k)     
	 
	      ys(2,2,k)=    -r                   
	  zs(2,2,k)=         zs(2,1,k)
	
	  
end do

  !内部网格
        do  j=1,ny73
	  do  k=1,nz73 
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny73)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny73)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny73)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny73)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
		end do

!区域-8
m=8
n=int(ny8/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
		end do
 do   j=1,ny8
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)   
end if
 end do

	
do j=1,ny8

	  ys(1,1,j)=  -( r+temp1*sb(1,j))   
	  zs(1,1,j)=      -sz2
	  ys(1,2,j)=  -( r+temp1*sb(1,j)) 
	  zs(1,2,j)=  0.0    
end do	
	
	n=int(nz81/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz81
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
		
		
	do k=1,nz81

	   ys(2,1,k)= -( r+(r*dcos(pi/4.0)-r) * sb(2,k))
	  zs(2,1,k)=  -sz2+(-r*dsin(pi/4.0)+sz2)*sb(2,k)
	end do
	n=int(nz82/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz82
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
	
	
	do k=1,nz82
 ys(2,1,nz81+k-1)=	 -(  r*dcos(pi/4.0-pi/4.0*sb(2,k))  )
 zs(2,1,nz81+k-1)=-r*dsin(pi/4.0-pi/4.0*sb(2,k))  
	end do
	       
   	n=int((nz81+nz82-1)/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz81+nz82-1
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do

	  do k=1,nz81+nz82-1
	  zs(2,2,k)=    zs(1,1,ny8)+(sz2)*sb(2,k)
	  ys(2,2,k)=	-(r+temp1)
	  
	end do
	
	        do  j=1,ny8
	  do  k=1,nz81+nz82-1
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny8)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny8)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny8)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny8)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
			end do		
!区域-9
m=9
n=int(ny9/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
		end do
 do   j=1,ny9
 if(j<=n)  then
 s(j)=s(j)/2.0
 sb(1,j)=s(j)
else
s(j)=0.5+s(j)/2.0  
sb(1,j)=s(j)   
end if
 end do

	
do j=1,ny9

	  ys(1,1,j)=  - (r*dcos(pi/4.0)+(r-r*dcos(pi/4.0)+temp1)*sb(1,j))   
	  zs(1,1,j)=      -sz2-r-temp2 
	
	  ys(1,2,j)=  - (r+temp1*sb(1,j)) 
	  zs(1,2,j)=   -sz2     
end do	
	
	n=int(nz91/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz91
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
		
		
	do k=1,nz91

	   ys(2,1,k)=  -(r*dcos(pi/4.0))       
	  zs(2,1,k)=    (-sz2-r-temp2)+(r+temp2-r*dsin(pi/4.0)) *sb(2,k)
	end do
	
	n=int(nz92/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz92
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do
	
	
	do k=1,nz92
 ys(2,1,nz51+k-1)=	  -( r*dcos(pi/4.0-pi/4.0*sb(2,k)) ) 
 zs(2,1,nz51+k-1)=-sz2-r*dsin(pi/4.0-pi/4.0*sb(2,k))  
	end do
	       
   	n=int((nz91+nz92-1)/2)+1
        p=0.5
        q=2.0
        call strech(n,p,q,s)
do k=1,n-1
s(n+k)=1.0-s(n-k)
end do 
do k=1,nz91+nz92-1
if(k<=n)   then
s(k)=s(k)/2.0
sb(2,k)=s(k)
else
s(k)=0.5+s(k)/2.0   
sb(2,k)=s(k)
end if
end do

	  do k=1,nz91+nz92-1
	  zs(2,2,k)=    zs(1,1,ny9)+(r+temp2)*sb(2,k)
	  ys(2,2,k)=	-(r+temp1)
	  
	end do
	
	        do  j=1,ny9
	  do  k=1,nz91+nz92-1 
        ya(j,k,m)=(1-sb(2,k))*ys(1,1,j)+sb(2,k)*ys(1,2,j)+(1-sb(1,j)) &
     &       *ys(2,1,k)+sb(1,j)*ys(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *ys(1,2,ny9)+sb(1,j)*(1-sb(2,k))*ys(1,1,ny9)+sb(2,k) &
     &       *(1-sb(1,j))*ys(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*ys(1,1,1))

	za(j,k,m)=(1-sb(2,k))*zs(1,1,j)+sb(2,k)*zs(1,2,j)+(1-sb(1,j)) &
     &       *zs(2,1,k)+sb(1,j)*zs(2,2,k)-(sb(1,j)*sb(2,k) &
     &       *zs(1,2,ny9)+sb(1,j)*(1-sb(2,k))*zs(1,1,ny9)+sb(2,k) &
     &       *(1-sb(1,j))*zs(1,2,1)+(1-sb(1,j)) &
     &       *(1-sb(2,k))*zs(1,1,1)) 
	  end do
			end do
		
!open(2,file='PCHE-2d.plt')
!      WRITE(2,*) 'variables= "y","z","m"'
!
!      m=1
!
!     	WRITE(2,*) 'zone t="fluid"',' j=',nz1,'k=',ny1
!	do 21 j=1,ny1 
!	do  21 k=1,nz1  
!	write(2,*) ya(j,k,m),za(j,k,m),m
!21	continue
!	
!	      m=2
!
!     	WRITE(2,*) 'zone t="solid"',' j=',nz2,'k=',ny2
!	do 22 j=1,ny2
!	do  22 k=1,nz2  
!	write(2,*) ya(j,k,m),za(j,k,m),m
!22	continue	
!	
!	      m=3
!
!     	WRITE(2,*) 'zone t="fluid"',' j=',nz3,'k=',ny3
!	do 23 j=1,ny3
!	do  23 k=1,nz3  
!	write(2,*) ya(j,k,m),za(j,k,m),m
!23	continue
!	
!		      m=4
!
!     	WRITE(2,*) 'zone t="solid"',' j=',nz4,'k=',ny4
!	do 24 j=1,ny4
!	do  24 k=1,nz4  
!	write(2,*) ya(j,k,m),za(j,k,m),m
!24	continue
!	
!	 m=5
!
!     	WRITE(2,*) 'zone t="fluid"',' j=',nz51+nz52-1,'k=',ny5
!	do 25 j=1,ny5
!	do  25 k=1,nz51+nz52-1
!	write(2,*) ya(j,k,m),za(j,k,m),m
!25	continue
!	m=6
!
!     	WRITE(2,*) 'zone t="fluid"',' j=',nz61+nz62-1,'k=',ny6
!	do 26 j=1,ny6
!	do  26 k=1,nz61+nz62-1
!	write(2,*) ya(j,k,m),za(j,k,m),m
!26	continue
!	
!	 m=7
!
!     	WRITE(2,*) 'zone t="fluid"',' j=',nz71,'k=',ny71
!	do 271 j=1,ny71
!	do  271 k=1,nz71 
!	write(2,*) ya(j,k,m),za(j,k,m),m
!271	continue
!	
!		 m=10
!
!     	WRITE(2,*) 'zone t="fluid"',' j=',nz72,'k=',ny72
!	do 272 j=1,ny72
!	do  272 k=1,nz72 
!	write(2,*) ya(j,k,m),za(j,k,m),m
!272	continue
!	
!			 m=11
!
!     	WRITE(2,*) 'zone t="fluid"',' j=',nz73,'k=',ny73
!	do 273 j=1,ny73
!	do  273 k=1,nz73 
!	write(2,*) ya(j,k,m),za(j,k,m),m
!273	continue
!	
!	 m=8
!
!     	WRITE(2,*) 'zone t="solid"',' j=',nz81+nz82-1,'k=',ny8
!	do 28 j=1,ny8
!	do  28 k=1,nz81+nz82-1
!	write(2,*) ya(j,k,m),za(j,k,m),m
!28	continue
!
!		 m=9
!
!     	WRITE(2,*) 'zone t="solid"',' j=',nz91+nz92-1,'k=',ny9
!	do 29 j=1,ny9
!	do  29 k=1,nz91+nz92-1
!	write(2,*) ya(j,k,m),za(j,k,m),m
!29	continue
!close(2)


	iter=20
	m=1
	call  LPLS(ya,za,iter,ny1,nz1,m)
		iter=20
		m=2
	call  LPLS(ya,za,iter,ny2,nz2,m)
		iter=20
	m=3
	call  LPLS(ya,za,iter,ny3,nz3,m)
	iter=20
	m=4
	call  LPLS(ya,za,iter,ny4,nz4,m)
	iter=20
	m=5
	call  LPLS(ya,za,iter,ny5,nz51+nz52-1,m)
	iter=20
	m=6
	call  LPLS(ya,za,iter,ny6,nz61+nz62-1,m)
	iter=20
	m=7
	call  LPLS(ya,za,iter,ny71,nz71,m)
		iter=20
	m=10
	call  LPLS(ya,za,iter,ny72,nz72,m)
		iter=20
	m=11
	call  LPLS(ya,za,iter,ny73,nz73,m)
	
	iter=20
	m=8
	call  LPLS(ya,za,iter,ny8,nz81+nz82-1,m)
	iter=20
	m=9
	call  LPLS(ya,za,iter,ny9,nz91+nz92-1,m)
	
	open(2,file='PCHE-2d-lpls.plt')
      WRITE(2,*) 'variables= "y","z","m"'

      m=1

     	WRITE(2,*) 'zone t="fluid"',' j=',nz1,'k=',ny1
	do 11 j=1,ny1 
	do  11 k=1,nz1  
	write(2,*) ya(j,k,m),za(j,k,m),m
11	continue
	    m=2

     	WRITE(2,*) 'zone t="fluid"',' j=',nz2,'k=',ny2
	do 12 j=1,ny2
	do  12 k=1,nz2
	write(2,*) ya(j,k,m),za(j,k,m),m
12	continue
      m=3

     	WRITE(2,*) 'zone t="fluid"',' j=',nz3,'k=',ny3
	do 13 j=1,ny3
	do  13 k=1,nz3  
	write(2,*) ya(j,k,m),za(j,k,m),m
13	continue
	
	      m=4

     	WRITE(2,*) 'zone t="solid"',' j=',nz4,'k=',ny4
	do 14 j=1,ny4
	do  14 k=1,nz4 
	write(2,*) ya(j,k,m),za(j,k,m),m
14	continue
		 m=5

     	WRITE(2,*) 'zone t="fluid"',' j=',nz51+nz52-1,'k=',ny5
	do 15 j=1,ny5
	do  15 k=1,nz51+nz52-1
	write(2,*) ya(j,k,m),za(j,k,m),m
15	continue
	
	m=6

     	WRITE(2,*) 'zone t="fluid"',' j=',nz61+nz62-1,'k=',ny6
	do 16 j=1,ny6
	do  16 k=1,nz61+nz62-1
	write(2,*) ya(j,k,m),za(j,k,m),m
16	continue
	
	 m=7

     	WRITE(2,*) 'zone t="fluid"',' j=',nz71,'k=',ny71
	do 171 j=1,ny71
	do  171 k=1,nz71
	write(2,*) ya(j,k,m),za(j,k,m),m
171	continue
		 m=10

     	WRITE(2,*) 'zone t="fluid"',' j=',nz72,'k=',ny72
	do 172 j=1,ny72
	do  172 k=1,nz72
	write(2,*) ya(j,k,m),za(j,k,m),m
172	continue
		 m=11

     	WRITE(2,*) 'zone t="fluid"',' j=',nz73,'k=',ny73
	do 173 j=1,ny73
	do  173 k=1,nz73
	write(2,*) ya(j,k,m),za(j,k,m),m
173	continue
	
	 m=8

     	WRITE(2,*) 'zone t="solid"',' j=',nz81+nz82-1,'k=',ny8
	do 18 j=1,ny8
	do  18 k=1,nz81+nz82-1
	write(2,*) ya(j,k,m),za(j,k,m),m
18	continue

		 m=9

     	WRITE(2,*) 'zone t="solid"',' j=',nz91+nz92-1,'k=',ny9
	do 19 j=1,ny9
	do  19 k=1,nz91+nz92-1
	write(2,*) ya(j,k,m),za(j,k,m),m
19	continue
	close(2)
	pause
	
	
	
	!建立3维
n=int(nx/2)+1  
p=1.0
q=2.0 
call strech(n,p,q,s)
do i=1,n-1
s(n+i)=1.0-s(n-i)
end do 
do i=1,nx
if(i<=n)   then
s(i)=s(i)/2.0
else
s(i)=0.5+s(i)/2.0   
end if
end do

allocate(xlast(300,300,300,11))
allocate(ylast(300,300,300,11))
allocate(zlast(300,300,300,11))

!区域1--------------------------------------
      m=1
      do 211 i=1,nx
      do 211 j=1,ny1
      do 211 k=1,nz1
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
211   continue
!区域2-----------------------------------------
      m=2
     
      do 221 i=1,nx
      do 221 j=1,ny2
      do 221 k=1,nz2
       ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
221   continue
!区域3-------------------------------------
      m=3

      do 231 i=1,nx
      do 231 j=1,ny3
      do 231 k=1,nz3
         ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
231   continue
!-区域4-------------------------------------
      m=4
   
      do 241 i=1,nx
      do 241 j=1,ny4
      do 241 k=1,nz4
        ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
241   continue
!区域5----------------------------------------
      m=5
    
      do 251 i=1,nx
      do 251 j=1,ny5
      do 251 k=1,nz51+nz52-1
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
251	  continue

	 !区域6---------------------------------------
      m=6
    
      do 261 i=1,nx
      do 261 j=1,ny6
      do 261 k=1,nz61+nz62-1
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
261	  continue 
	  
		 !区域7---------------------------------------
      m=7
    
      do 2711 i=1,nx
      do 2711 j=1,ny71
      do 2711 k=1,nz71
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
2711  continue  
	  
	        m=10
    
      do 2712 i=1,nx
      do 2712 j=1,ny72
      do 2712 k=1,nz72
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
2712  continue  
	
	        m=11
    
      do 2713 i=1,nx
      do 2713 j=1,ny73
      do 2713 k=1,nz73
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
2713  continue  
	  
	  	 !区域8---------------------------------------
      m=8
    
      do 281 i=1,nx
      do 281 j=1,ny8
      do 281 k=1,nz81+nz82-1
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
281	  continue 
	  
	  
	  	  	 !区域9---------------------------------------
      m=9
    
      do 291 i=1,nx
      do 291 j=1,ny9
      do 291 k=1,nz91+nz92-1
      ylast(i,j,k,m)=ya(j,k,m)
      zlast(i,j,k,m)=za(j,k,m)
      xlast(i,j,k,m)=xjl*s(i)
291	  continue 
	  
	  
 open(2,file='PCHE-3d.plt')
      WRITE(2,*) 'variables= "x","y","z","m"'

      m=1
     
     	WRITE(2,*) 'zone t="fluid"',' i=',nz1,'j=',ny1,'k=',nx
	do 126 i=1,nx 
	do 126 j=1,ny1
      do 126 k=1,nz1
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1   !为什么有个1
126   continue
      

      m=2
      WRITE(2,*) 'zone t="fluid"',' i=',nz2,'j=',ny2,'k=',nx
	do 127 i=1,nx 
	do 127 j=1,ny2
      do 127 k=1,nz2
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
127   continue
      

      m=3
 
      WRITE(2,*) 'zone t="fluid"',' i=',nz3,'j=',ny3,'k=',nx
	do 128 i=1,nx 
	do 128 j=1,ny3
      do 128 k=1,nz3
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
128   continue
      

      m=4
 
      WRITE(2,*) 'zone t="fluid"','i=',nz4,'j=',ny4,'k=',nx
	do 129 i=1,nx 
	do 129 j=1,ny4
      do 129 k=1,nz4
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
129   continue
      

      m=5
    
      WRITE(2,*) 'zone t="fluid"',' i=',nz51+nz52-1,'j=',ny5,'k=',nx
	do 130 i=1,nx 
	do 130 j=1,ny5
      do 130 k=1,nz51+nz52-1
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
130	continue
	      m=6
    
      WRITE(2,*) 'zone t="fluid"',' i=',nz61+nz62-1,'j=',ny6,'k=',nx
	do 131 i=1,nx 
	do 131 j=1,ny6
      do 131 k=1,nz61+nz62-1
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
131	continue
	
		      m=7
    
      WRITE(2,*) 'zone t="fluid"',' i=',nz71,'j=',ny71,'k=',nx
	do 1321 i=1,nx 
	do 1321 j=1,ny71
      do 1321 k=1,nz71
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
1321 continue
	 
	 		      m=10
    
      WRITE(2,*) 'zone t="fluid"',' i=',nz72,'j=',ny72,'k=',nx
	do 1322 i=1,nx 
	do 1322 j=1,ny72
      do 1322 k=1,nz72
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
1322 continue
	 
	 		      m=11
    
      WRITE(2,*) 'zone t="fluid"',' i=',nz73,'j=',ny73,'k=',nx
	do 1323 i=1,nx 
	do 1323 j=1,ny73
      do 1323 k=1,nz73
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
1323 continue
	
	      m=8
    
      WRITE(2,*) 'zone t="fluid"',' i=',nz81+nz82-1,'j=',ny8,'k=',nx
	do 133 i=1,nx 
	do 133 j=1,ny8
      do 133 k=1,nz81+nz82-1
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
133	continue	
	
	
		      m=9
    
      WRITE(2,*) 'zone t="fluid"',' i=',nz91+nz92-1,'j=',ny9,'k=',nx
	do 134 i=1,nx 
	do 134 j=1,ny9
      do 134 k=1,nz91+nz92-1
	write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),1
134	continue
      close(2)
end program main
	
subroutine strech(n,u,v,s)
implicit none 
real(kind=8)  ::  s(100)
integer(kind=4)  :: an,al,n,l
real(kind=8)  :: deta,tqi,eta,dum,u,v
an=n-1
deta=1.0/float(an)
tqi=1.0/tanh(v)
do 1 l=1,n
   al=l-1
   eta=float(al)*deta
   dum=v*(1.0-eta)
   dum=1.0-tanh(dum)*tqi
 s(l)=u*eta+(1.0-u)*dum
1     continue
return
	end
	
	subroutine LPLS(x,y,iter,nx,ny,m)
implicit none
real(kind=8)  ::  x(300,300,m),y(300,300,m)
real(kind=8)  :: om,sumx,sumy,xxip,yxip,xetp,yetp,alfa,beta,gama,bb,xwdd,ywdd,dify,difx
integer(kind=4) :: iter,n,i,j,im,ip,jm,jp,nx,ny  ,m
om=1.0	
do 1 n=1,iter
sumx=0.0
sumy=0.0
do 2 i=2,nx-1
im=i-1
ip=i+1
do 3 j=2,ny-1
jm=j-1
jp=j+1   
xxip=(x(ip,j,m)-x(im,j,m))/2.0
yxip=(y(ip,j,m)-y(im,j,m))/2.0      
xetp=(x(i,jp,m)-x(i,jm,m))/2.0
yetp=(y(i,jp,m)-y(i,jm,m))/2.0   

alfa=xetp**2+yetp**2
beta=xxip*xetp+yxip*yetp
gama=xxip**2+yxip**2   

bb=alfa*(x(im,j,m)+x(ip,j,m))+gama*(x(i,jm,m)+x(i,jp,m))-0.5d0*beta*(x(ip,jp,m)-x(im,jp,m)-x(ip,jm,m)+x(im,jm,m))
xwdd=bb/(2.d0*(alfa+gama)+1.e-30)
difx=xwdd-x(i,j,m)
x(i,j,m)=x(i,j,m)+om*difx
sumx=sumx+difx**2   
bb=alfa*(y(im,j,m)+y(ip,j,m))+gama*(y(i,jm,m)+y(i,jp,m))-0.5d0*beta*(y(ip,jp,m)-y(im,jp,m)-y(ip,jm,m)+y(im,jm,m)) 
ywdd=bb/(2.d0*(alfa+gama)+1.e-30)  
dify=ywdd-y(i,j,m)
y(i,j,m)=y(i,j,m)+om*dify
sumy=sumy+dify**2
3    continue
2    continue
1    continue
return
end 

	
	