function [x,k]=myNewton2(x,e,maxit)
  
  
  g=@(x) [-400*x(1)*(x(2)-x(1)^2)-(1-x(1))*2 ; 200*(x(2)-x(1)^2)];
  H=@(x) [1200*x(1)^2-400*x(2)+2, -400*x(1); -400*x(1), 200];
    
     
  if ~iscolumn(x)
    error('a kezdovektornak oszlopvektornak kell lennie');
  end
  
  
  xx=linspace(-2,2);
  yy=linspace(-3.5,4);
  [X,Y]=meshgrid(xx,yy);1
  Z=100*(Y-X.^2).^2+(1-X).^2;
  figure; contour(X,Y,Z,linspace(min(min(Z)),max(max(Z)),50));
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