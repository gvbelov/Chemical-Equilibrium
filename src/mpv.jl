using LinearAlgebra
function back_subst(R,b,xrem)
 n = length(b)
 x = zeros(n)
 for i=n:-1:1
   if abs(R[i,i]) > 1.e-15 
     x[i] = (b[i] - R[i,i+1:n]'*x[i+1:n]) / R[i,i] 
   else x[i] = xrem[i]
   end
   end
 return x
end;

function CheckConvergency(p,m,k, b, Nji, F, ex, xx, xs, StandardPressure)
# check MB
  epsmb = 0
  for j in 1:m
    x=0
    for i in 1:k
      x+=ex[i]*Nji[i,j]
    end
  x=abs(b[j]-x)
  if x > epsmb epsmb = x end
  end
  SumMolarFractions=0
  SumG=0
  px=(1.e6/StandardPressure)*p
  for i in 1:k
    x=0 
    for j in 1:m
    x+=Nji[i,j]*xx[j]
    end
  SumMolarFractions+=exp(F[i]+x)/px
  SumG+=(log(px*ex[i]/xs)-F[i])*ex[i]
  end
  SumL=0
  for j in 1:m 
    SumL+=b[j]*xx[j]
  end
  epsg=abs(SumL-SumG) # sum(bj*Lamj) - sum(Gi(T,pi))
  return epsmb, epsg, SumMolarFractions, SumG
end

function MicPV(T,p,v, m,k, b, Nji, F)
      Parameters=zeros(4); EquilibriumComposition=zeros(k); bb = zeros(m+1)
      Rgas = 8.314462618
      StandardPressure = 101325   #  =1.e5, if 1 bar
      local xx, xs, eps
sum_b=0; 
for i in 1:m sum_b+=b[i] end
for i in 1:m b[i] = b[i]/sum_b end
v=v/sum_b    
    if (m == 0) | (k == 0) ier = " (m == 0) | (k == 0)"; return ier, Parameters, EquilibriumComposition,0 end 
      m1=m+1; m2=m+2; iter=0; ier="";
      x=zeros(k); ex=zeros(k); str=zeros(m2); xrem=zeros(m1)
        for i in 1:k 
          x[i] = -15.0 
          ex[i]= exp(x[i]) end
        y=0.0
while true
      eps = 0.0; S=zeros(m1,m1)
      for j in 1:m bb[j]=b[j] end
      if v <= 0.0 
	S[m1,m1]=(1.e6/StandardPressure)*p*exp(-y)
 	bb[m1]=S[m1,m1]*(1.0+y) 
      else 
        S[m1,m1]=1.0
        bb[m1]=log((Rgas/StandardPressure)*T/v) end
      for i in 1:k  
        str[m1]=-ex[i]
        str[m2]=(x[i]-F[i]-1.0)*ex[i] 
        for j in 1:m
          str[j]=Nji[i,j]*ex[i]; end
       if v<=0.0 
        for j in 1:m1  
          S[m1,j]=S[m1,j]+str[j] end;
       bb[m1]=bb[m1]+str[m2]
       end
       for jj in 1:m # do 
         for j in 1:m1 # do 
           S[jj,j] = S[jj,j]+str[j]*Nji[i,jj]; end; 
        bb[jj] = bb[jj]+str[m2]*Nji[i,jj]
        end
      end;
# find the SLE solution
      z=qr(S)   
      yy=z.Q'*bb
      xx = back_subst(z.R,yy, xrem)
# modify the results of iteraion 
      w = 0.5
      if (abs(y) > 1.0) w = 0.5*abs(y);  end
      if (abs(xx[m1]-y) < w) w=abs(xx[m1]-y);  end
      if (xx[m1] < y) w = -w;  end
      if (abs(y)>1e-10) & (abs(w/y) > eps) eps=abs(w/y); end
      y = y+w 
      xs = 0.0;
      for i in 1:k
        u=F[i]-xx[m1]              
        for j in 1:m
          u=u+Nji[i,j]*xx[j]; end
        w=3.0 
        if u>=x[i]  
          w=1.0 
          if x[i]<-7.0 w=5.0 end
        end 
        if abs(u-x[i]) < w w=abs(u-x[i]); end
        if u < x[i] w=-w; end 
        if (x[i] > -80) & (abs(x[i])>1e-10) & (abs(w/x[i])>eps) eps=abs(w/x[i]); end 
        x[i]=x[i]+w; 
        ex[i]=exp(x[i]);
        xs=xs+ex[i]
        end;
# remember the result of current iteration
      xrem = copy(xx)
      iter=iter+1;
    eps <= 1.0e-6 && break
    iter > 200 && break
    end; # while true
    if iter > 200 (ier="after 200 iterations eps="*string(eps)) end
# Prepare the results for the export
      if v <= 0.0 
        Parameters[1]=p; 
      else  
        Parameters[1]=1.e-6*Rgas*T*xs/v; 
      end;    
# check up the convergence
    epsmb, epsg, SumMolarFractions, SumG =CheckConvergency(Parameters[1],m,k, b, Nji, F, ex, xx, xs,StandardPressure)
      Parameters[2]=T; 
 println("max relative error = ", eps)
 println("material balance error = ", epsmb)
 println("Gibbs energy error = ", epsg)
 println("sum of molar fractions = ", SumMolarFractions)
      if v <= 0.0 
        Parameters[3]=(Rgas/StandardPressure)*T*sum_b/exp(y); 
      else
        Parameters[3]=v*sum_b 
      end;
      Parameters[4]=xs*sum_b;
      for i in 1:k EquilibriumComposition[i]=ex[i]*sum_b; end
      for i in 1:m b[i]=b[i]*sum_b; end
      return ier, Parameters, EquilibriumComposition, SumG*sum_b
  end 
