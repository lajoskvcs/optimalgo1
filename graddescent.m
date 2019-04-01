function [x, k] = graddescent(x, epsilon)
   c1 = 1e-4;
   alpha0 = 0.4;
   ro = 0.5;
   k = 0;
   
   xx = linspace(-2,2);
   yy = xx;

   [X, Y] = meshgrid(xx, yy);

    Z = X.^3 - 3*X + Y.^3 - 3*Y;

    figure; contour(X, Y, Z); hold on;
   
   while norm(grad_fx(x)) >= epsilon 
      pk = -1 * grad_fx(x);
      alpha = backtrack(x, c1, pk, alpha0, ro);
      x = x + alpha * pk;
      k = k + 1;
      plot(x(1), x(2), 'b*');
      pause(0.5);
   end
end

function f = fn(x)
  f = x(1)^3 + x(2)^3 - 3*x(1) - 3*x(2);
end

function g = grad_fx(x)
  g = zeros(2, 1);
  g(1) = 3 * x(1) ^ 2 - 3;
  g(2) = 3 * x(2) ^ 2 - 3;
end

function alpha = backtrack(x, c1, pk, alpha0, ro)
  alpha = alpha0;
  while fn(x + alpha * pk) > fn(x) + alpha * c1 * grad_fx(x)' * pk
    alpha = ro * alpha;
  end
end