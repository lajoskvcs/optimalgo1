# Eredeti:
# f(x) = x(1)^5 - 8*x(1) + 2*x(2)^3 - 3*x(2) - x(1)*x(2)/6

#tartomany: (2,2)^2

# deriv1: 5*x(1)^4 - 8 - x(2)/6
# deriv2: 6*x(2)^2 - 3 - x(1)/6

# deriv1x1: 20*x(1)^3
# deriv1x2: 1/6
# deriv2x1: 1/6
# deriv2x2: 12*x(2)

# hesse = [20*x(1)^3, 1/6; 1/6, 12*x(2)];

# gradv = [5*x(1)^4 - 8 - x(2)/6 ; 6*x(2)^2 - 3 - x(1)/6];

max_iterations = 30;
epsilon = 1e-5;
# x0 = [1;1] # grad: 19, newton: 5, bfgs: 6
# x0 = [-1; 0]; # grad: 21, newton: 17, bfgs: 10
x0 = [1; -0.5]; # grad: 24(1.12893, 0.72894), newton: 5(1.12039, -0.72878), bfgs: 9(1.12893, 0.72894)

xx = linspace(-2,2);
yy = xx;

[X, Y] = meshgrid(xx, yy);

Z = X.^5 - 8*X + 2*Y.^3 - 3*Y - X*Y./6;

figure('name', 'Surface'); surf(X, Y, Z); hold on;

[x_grad, k_grad] = graddescent(x0, epsilon)

[x_newton, k_newton] = myNewton2(x0, epsilon, max_iterations)

[x_bfgs, k_bfgs] = myBFGS2(x0, epsilon, max_iterations)



