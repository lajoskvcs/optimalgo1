function [x,k,fval]=myBFGS2(x,e,maxit)
  
   
   if ~iscolumn(x)
    error('the initial vector has to be a column vector');
   end
   
   if nargin<3
     maxit=100;
   end
   
   
   XS=zeros(length(x),maxit+1);
   XS(:,1)=x;
   
   k=0;
   H=eye(length(x));
   
   while norm(gr(x))>e
     p=-H*gr(x);
     a=backtr(0.5,1e-4,x,p,1);
     xn=x+a*p;
     k=k+1;
     if k==maxit
       error('The number of steps exceeded maxit');
     end
     H=hessematrix(xn,x,H);
     x=xn;
     XS(:,k+1)=x;
   end
   fval=fv(x);

   
   
   
   xx=linspace(-2,2);
   yy=linspace(-2,2);
   [X,Y]=meshgrid(xx,yy);
   Z = X.^5 - 8*X + 2*Y.^3 - 3*Y - X*Y./6;
   figure('name', 'BFGS'); contour(X,Y,Z,linspace(min(min(Z)),max(max(Z)),50));
   hold on;   
   plot(XS(1,1:k+1),XS(2,1:k+1),'*-')
   
 end  
   
  % the objective function 

#  function f=fv(x)
#     f = 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
#   end
  function f=fv(x)
     f = x(1)^5 - 8*x(1) + 2*x(2)^3 - 3*x(2) - x(1)*x(2)/6;
   end
 
 % the gradient of the objective function
#   function g=gr(x)
#     g=[-400*x(1)*(x(2)-x(1)^2)-(1-x(1))*2 ;
#        200*(x(2)-x(1)^2)];
#   end
  function g=gr(x)
     g = [5*x(1)^4 - 8 - x(2)/6 ; 6*x(2)^2 - 3 - x(1)/6];
   end 

 
  % the updating formula for H 
  function H=hessematrix(xn,x,H)
    s=xn-x;
    y=gr(xn)-gr(x);
    H=(eye(size(H))-s*y'/dot(y,s))*H*(eye(size(H))-y*s'/dot(y,s))+s*s'/dot(s,y);
  end
  
  
  
 %the backtracking algorithm 
   function a=backtr(rho,c,xk,pk,a)
     while fv(xk+a*pk)>fv(xk)+a*c*dot(gr(xk),pk)  
         a=rho*a;
      end
   end
   
   

