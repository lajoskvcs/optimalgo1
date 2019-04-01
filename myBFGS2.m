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

   
   
   
   xx=linspace(-1.5,2);
   yy=linspace(-3.5,4);
   [X,Y]=meshgrid(xx,yy);
   Z=100*(Y-X.^2).^2+(1-X).^2;
   figure; contour(X,Y,Z,linspace(min(min(Z)),max(max(Z)),50));
   hold on;   
   plot(XS(1,1:k+1),XS(2,1:k+1),'*-')
   
 end  
   
  % the objective function 
   function f=fv(x)
     f=100*(x(2)-x(1)^2)^2+(1-x(1))^2;
   end
 
 % the gradient of the objective function
   function g=gr(x)
     g=[-400*x(1)*(x(2)-x(1)^2)-(1-x(1))*2 ;
        200*(x(2)-x(1)^2)];
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
   
   

