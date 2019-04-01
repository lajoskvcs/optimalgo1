function [x,k]=myNewton2(x,e,maxit)
  
  
  g=@(x) [5*x(1)^4 - 8 - x(2)/6 ; 6*x(2)^2 - 3 - x(1)/6];
  H=@(x) [20*x(1)^3, 1/6 ; 1/6, 12*x(2)];
    
     
  if ~iscolumn(x)
    error('a kezdovektornak oszlopvektornak kell lennie');
  end
  
  
  xx=linspace(-2,2);
  yy=linspace(-2,2);
  [X,Y]=meshgrid(xx,yy);1
  Z = X.^5 - 8*X + 2*Y.^3 - 3*Y - X*Y./6;
  figure('name', 'Newton'); contour(X,Y,Z,linspace(min(min(Z)),max(max(Z)),50));
  hold on;      
  
  XS=zeros(length(x),maxit+1);
  XS(:,1)=x;
  
  
  k=0;
  while norm(g(x))>e
    p=-H(x)\g(x);
    x=x+p;
    k=k+1;
    XS(:,k+1)=x;
    if k==maxit
      error('elertuk maxit-et');
    end
    
  end
  
  plot(XS(1,1:k+1),XS(2,1:k+1),'*-')
  
  
end