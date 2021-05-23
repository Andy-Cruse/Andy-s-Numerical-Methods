format long
N = 32;
a = 0;
b = 1;
t = (b-a)/N;
LinearFunc(a, b, N, t)
%The interval [1,2] converged at a faster rate than the interval 
%[0,1] did. Even for higher values of N approximating the function 
%linearly is difficult since it's slope is approaching 0 at a fast rate. 
function LinearFunc(a, b, N, t)

for i = 1:1001
    z_i = a + ((i-1) * ( (b-a)/1000));
    domain(i) = z_i;
end

x(1) = a;
for i = 2:N+1
    x(i) = x(i-1) + t;
end

%This for loop will generate each p_i where each i is stored as 
%a row and every column is a value of x in the domain D
i = 2;
for k = 1:1001
    b = mod(k, ceil(1000/N));
    if (b == 0)
        if (k ~= 1000)
            i=i+1;
        end
    end
    p(1,k) = ( ( (domain(k)-x(i-1))*(x(i))^(1/3) ) + ( (x(i)-domain(k))*(x(i-1))^(1/3) ) ) / (x(i) - x(i-1));
    f(1,k) = domain(k)^(1/3);
end
plot(p, 'r')
hold on;
plot(f, 'b')
hold off;
dlmwrite('q_n_for_each_x_in_domain.txt', p);
error = max(abs(p-f));
disp(error)
end