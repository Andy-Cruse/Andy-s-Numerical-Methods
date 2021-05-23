x = -1;
y = 1;
n = 2;
a(1) = rightEndPoint(x,y,n);
a(2) = midpoint(x,y,n);
a(3) = trapezoid(x,y,n);
%a(4) = correctedTrapezoid(x,y,n);
%a(5) = simpsons(x,y,n);
% a(6) = gaussian(x,y,1);
% a(7) = gaussian(x,y,2);
% a(8) = gaussian(x,y,4);
%disp(a)
% a(1) = abs(-2*pi - a(1));
% a(2) = abs(-2*pi - a(2));
% a(3) = abs(-2*pi - a(3));
% a(4) = abs(-2*pi - a(4));
% a(5) = abs(-2*pi - a(5));
% disp(a(1:5))
 a(1) = abs((pi/2) - a(1));
 a(2) = abs((pi/2) - a(2));
 a(3) = abs((pi/2) - a(3));
 %a(4) = abs((pi/2) - a(4));
 %a(5) = abs((pi/2) - a(5));
% disp(a(1:5))
g(1) = gaussian(x,y,1)
g(2) = gaussian(x,y,2)
g(3) = gaussian(x,y,4)
g(4) = gaussian(x,y,8)
g(5) = gaussian(x,y,16)
disp(g)
g(1) = abs(gaussian(x,y,1) - (pi/2));
g(2) = abs(gaussian(x,y,2) - (pi/2));
g(3) = abs(gaussian(x,y,4) - (pi/2));
g(4) = abs(gaussian(x,y,8) - (pi/2));
g(5) = abs(gaussian(x,y,16) - (pi/2));
disp(g)


function fx = f(x)
   %fx = x^9; %for Q1
   %fx = x^2 * cos(x); %for Q2 & Q4
   fx = (1-x^2)^(0.5); %for Q3 & Q4
end

function Fx = rightEndPoint(a,b,n)
    h = (b-a)/n;
    for i = 0:n
        x(i+1) = a + i*h;
        fx(i+1) = f(x(i+1));
    end
    Fx = sum(fx(2:end));
    Fx = Fx*h;
end

function Fx = midpoint(a,b,n)
    h = (b-a)/n;
    for i = 1:n+1
        x(i) = a + (i-1)*h;
    end
    for i = 1:n
        fx(i) = f( (x(i) + x(i+1))/2 );
    end
    Fx = h*sum(fx(1:end));
end

function Fx = trapezoid(a,b,n)
    h = (b-a)/n;
    for i = 1:n+1
        x(i) = a + (i-1)*h;
        fx(i) = f(x(i));
    end
    Fx = (h/2)*( fx(1) + 2*sum(fx(2:n)) + fx(n+1) );
end

function Fx = correctedTrapezoid(a,b,n)
    if (n < 2)
        error('n cannot be less than 2')
    end
    
    h = (b-a)/n;
    for i = 1:n+1
        x(i) = a + (i-1)*h;
        fx(i) = f(x(i));
    end
    Fx = trapezoid(a,b,n) - (h/24)*(3*fx(n+1) - 4*fx(n) + fx(n-1) + 3*fx(1) - 4*fx(2) + fx(3));
end

function Fx = simpsons(a,b,n)
    if (n < 2 || mod(n,2) ~= 0)
        error('n cannot be odd or less than 2');
    end
    
    h = (b-a)/n;
    for i = 1:n+1
        x(i) = a + (i-1)*h;
        fx(i) = f(x(i));
    end
    Fx = fx(1) + fx(n+1);
    for i = 2:n
        if (mod(i,2) == 0)
            Fx = Fx + 4*fx(i);
        else
            Fx = Fx + 2*fx(i);
        end
    end
    Fx = (h/3)*Fx;
end

function Fx = gaussian(a,b,n)
    if (n == 1)
        z(1) = 0;
        w(1) = 2;
    elseif (n == 2)
        z(1) = 0.5773502691896257;
        w(1) = 1;
        z(2) = -z(1);
        w(2) = 1;
    elseif (n == 4)
        z(1) = 0.3399810435848562;
        w(1) = 0.6521451548625464;
        z(3) = 0.8611363115940526;
        w(3) = 0.3478548451374476;
        
        z(2) = -z(1);
        z(4) = -z(3);
        w(2) = w(1);
        w(4) = w(3);
    elseif (n == 8)
        z(1) = 0.1834346424956498;
        w(1) = 0.3626837833783620;
        z(3) = 0.5255324099163290;
        w(3) = 0.3137066458778874;
        z(5) = 0.7966664774136268;
        w(5) = 0.2223810344533745;
        z(7) = 0.9602898564975362;
        w(7) = 0.1012285362903697;
        
        z(2) = -z(1);
        z(4) = -z(3);
        z(6) = -z(5);
        z(8) = -z(7);
        w(2) = w(1);
        w(4) = w(3);
        w(6) = w(5);
        w(8) = w(7);
    elseif (n == 16)
        z(1) = 0.9501250983763744 * 10^-1;
        w(1) = 0.1894506104550685;
        z(3) = 0.2816035507792589;
        w(3) = 0.1826034150449236;
        z(5) = 0.4580167776572274;
        w(5) = 0.1691565193950024;
        z(7) = 0.6178762444026438;
        w(7) = 0.1495959888165733;
        z(9) = 0.7554044083550030;
        w(9) = 0.1246289712555339;
        z(11) = 0.865631202387831;
        w(11) = 0.9515851168249290 * 10^-1;
        z(13) = 0.9445750230732326;
        w(13) = 0.6225352393864778 * 10^-1;
        z(15) = 0.9894009349916499;
        w(15) = 0.2715245941175185 * 10^-1;
        
        z(2) = -z(1);
        z(4) = -z(3);
        z(6) = -z(5);
        z(8) = -z(7);
        z(10) = -z(9);
        z(12) = -z(11);
        z(14) = -z(13);
        z(16) = -z(15);
        
        w(2)  = w(1);
        w(4)  = w(3);
        w(6)  = w(5);
        w(8)  = w(7);
        w(10)  = w(9);
        w(12)  = w(11);
        w(14)  = w(13);
        w(16)  = w(15);
    else
        error('n must equal 1,2,4,8, or 16')
    end
    
    for i = 1:n
        Fz(i) = ((b-a)/2)*f(a + ((b-a)/2)*(z(i) + 1));
        Fz(i) = Fz(i)*w(i);
    end
    
    Fx = sum(Fz);
end
