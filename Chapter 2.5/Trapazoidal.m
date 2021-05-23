a = 1;
b = 3;
N = 32;
%function1(a, b, N)

a = 0;
b = 2;
N = 32;
%function2(a,b,N)

a = 0;
b = 1;
N = 32;
function3(a,b,N)

function function1(a, b, N)
    T = ( (b-a)/(2*N) )*( (a^2)*(1/exp(a)));
    for i = 2:N-1
        z_i = a + ((i-1) * ( (b-a)/N));
        T = T + ( (b-a)/(N) )*( (z_i^2)*(1/exp(z_i)));
    end
    T = T + ( (b-a)/(2*N) )*( (b^2)*(1/exp(b)));
    disp(T)
    I = (5*exp(2) - 17) / exp(3);
    disp(T - I)
end

function function2(a, b, N)
    T = ( (b-a)/(2*N) )*( (a^2)*(1/exp(a)));
    for i = 2:N-1
        z_i = a + ((i-1) * ( (b-a)/N));
        T = T + ( (b-a)/(N) )*( (z_i^2)*(1/exp(z_i)));
    end
    T = T + ( (b-a)/(2*N) )*( (b^2)*(1/exp(b)));
    disp(T)
    disp(T - (2-10*exp(-2)))
end

function function3(a, b, N)
    T = ( (b-a)/(2*N) )*( a^(1/2));
    for i = 2:N-1
        z_i = a + ((i-1) * ( (b-a)/N));
        T = T + ( (b-a)/(N) )*( z_i^(1/2) );
    end
    T = T + ( (b-a)/(2*N) )*( b^(1/2));
    disp(T)
    disp(T - (2/3))
end