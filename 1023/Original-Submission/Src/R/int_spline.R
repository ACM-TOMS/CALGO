#   Restoration of function by integrals with cubic integral smoothing basis spline in R.
#   The method is based on a specially designed cubic integral basis spline with a penalty function,
# which minimizes the sum of squares of the difference between the observed integrals of the unknown
# function and the integrals of the spline being constructed, plus an additional penalty for the
# nonlinearity (roughness) of the spline.
#    Autor - Korablev Yuriy Aleksandrovich
# PhD in Economics, docent, department "System Analysis in Economics"
# Financial University under the Government of the Russian Federation
# email: yura-korablyov@yandex.ru, YUAKorablev@fa.ru

# Copyright © 2020  Korablev Yuriy Aleksandrovich
#   
#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>

#================ Spline computation ===================

int_spline = function (t,Y,s=t,m=length(s),W=rep(1,length(Y)),alpha=10^5,x=seq(t[1],t[length(t)],1),info=FALSE)
{
  # t - array of observation coordinates (beginning of integrals) 
  # Y - array of observation, where Y[i] is integral of unknow function from t[i] to t[i+1]
  # s - spline knots
  # m - number of spline knots
  # W - weight matrix
  # alpha - smoothing parameter
  # x - coordinates in which spline will be calculated
  
  n=length(t) #number of observation coordinates   
  
  # in case s is not defined
  if (m!=length(s)) #its possible when m was defined, but s wasn't
    s=seq(t[1],t[n],length=m)
  
  h=array(0,dim=m-1) #array of distance between knots
  h[1:(m-1)]=s[2:m]-s[1:(m-1)]

  #Matrix Q
  Q=matrix(0,nrow=m,ncol=m-2)
  for (i in 1:(m-2))
  {
    Q[i,i]=1/h[i];
    Q[i+1,i]=-1/h[i]-1/h[i+1];
    Q[i+2,i]=1/h[i+1]
  }

  #Matrix R
  R=matrix(0,nrow=m-2,ncol=m-2)
  for (i in 1:(m-2))
  {
    R[i,i]=1/3*(h[i]+h[i+1]);
    if (i<m-2)
    {
      R[i+1,i]=1/6*h[i+1];
      R[i,i+1]=1/6*h[i+1];
    }
  }

  #Matrix K calculation
  inv_R=solve(R)
  t_Q=t(Q)
  K=Q %*% inv_R %*% t_Q

  #Filling in V and P matrices
  V=matrix(0,nrow=n-1,ncol=m)
  P=matrix(0,nrow=n-1,ncol=m)
  k=1
  while( (s[k]<=t[1]) & (s[k+1]<t[1])) #find first k, that s[k+1]>t[1]
    k=k+1
  for (i in 1:(n-1))
  {
    #finding L, it can be 0
    for (L in 0:(m-k-1))
      if (t[i+1]<=s[k+L+1])
        break;
    l=1;
    V[i,k]=(s[k+1]-t[i])^2/h[k]/2  
    P[i,k]=h[k]^3/24-(t[i]-s[k])^2*(s[k+1]-t[i]+h[k])^2/h[k]/24    
    while (l<=L)
    {
      V[i,k+l]=(h[k+l-1]+h[k+l])/2
      P[i,k+l]=(h[k+l-1]^3+h[k+l]^3)/24
      l=l+1;
    }
    V[i,k+1]=V[i,k+1]-(t[i]-s[k])^2/h[k]/2
    P[i,k+1]=P[i,k+1]+(t[i]-s[k])^2*((t[i]-s[k])^2-2*h[k]^2)/h[k]/24
    V[i,k+L]=V[i,k+L]-(s[k+L+1]-t[i+1])^2/h[k+L]/2
    P[i,k+L]=P[i,k+L]+(s[k+L+1]-t[i+1])^2*((s[k+L+1]-t[i+1])^2-2*h[k+L]^2)/h[k+L]/24    
    V[i,k+L+1]=(t[i+1]-s[k+L])^2/h[k+L]/2
    P[i,k+L+1]=h[k+L]^3/24-(s[k+L+1]-t[i+1])^2*(t[i+1]-s[k+L]+h[k+L])^2/h[k+L]/24
    k=k+L
  }
  P=P[1:(n-1),2:(m-1)] #don't need first and last column

  #Matrix C calculation
  
  C=V-P %*% inv_R %*% t_Q

  #Calculation of g and gamma
  t_C=t(C)
  W=diag(W)                            # Weight matrix
  A=t_C %*% W %*% C + alpha * K
  g=solve(A , t_C %*% W %*% Y )
  gamma=inv_R %*% t_Q %*% g
  #After that spline is completely defined via g and gamma

  #============== Calculating and returning spline values in x coordinates  ================
  
  g2=c(0,gamma,0) #Second derivative on the edges was zero

  #x=seq(t[1],t[n],1) by default
  y=rep(0,length(x)) 

  k=1; #index of interval 
  for (j in (1:length(x)))
  {
    while (x[j]<t[n] & x[j]>s[k]+h[k])
      k=k+1;
    y[j] = ((x[j]-s[k])*g[k+1]+(s[k+1]-x[j])*g[k])/h[k] - 1/6*(x[j]-s[k])*(s[k+1]-x[j])*(g2[k+1]*(1+(x[j]-s[k])/h[k])+g2[k]*(1+(s[k+1]-x[j])/h[k])  )
  }
  if (info)
  {
    result=list(x=x,y=y,g=g,gamma=g2,s=s,h=h)
    return (result)
  }
  else 
    return(y)
}

