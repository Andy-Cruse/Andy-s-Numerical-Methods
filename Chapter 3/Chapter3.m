% format long;
%root = bisection(-3, -2, 10^-6, 1, 1);
% disp(root)
%root = newton(-0.705, 10^-4, 2, 0, 1, 30);
% disp(root)
%root = secant(0, 2, 10^-6, 2, log(2), 1, 20);
% disp(root)
root = hybrid(1, -0.705, 10^-4, 2, 0)
% disp(root)
%root = chord(0, 10^-6, 1, log(2), 1, 20);
% disp(root)
%root = superHalley(0, 10^-6, 2, log(2), 1, 20);
% disp(root)

%to graph functions
dom(1) = -2;
f(1) = poly(dom(1),0);
h = (4/1000);
for i = 2:1000
     dom(i) = dom(i-1) + h;
     f(i) = poly(dom(i),0);
end
plot(dom,f)

function fc = bisection(a, b, TOL, Output_Flag, alpha)
    fa = poly(a, 0);
    fb = poly(b, 0);
    
    if (fa*fb > 0)
        error('BAD INTERVAL! NO ROOT FOUND!')
    else
        n = fix((log(b-a) - log(TOL))/log(2))+1;
        for i = 1:n
            c = a + 0.5*(b-a);
            ao(i) = a;   %for output
            bo(i) = b; 
            co(i) = c;
            fc = poly(c, 0);
            if (fa*fc < 0)
                b = c;
                fb = fc;
            elseif (fa*fc > 0)
                a = c;
                fa = fc;
            else
                alpha = c;
            end
        end
        
        ao(n+1) = a;
        bo(n+1) = b;
        co(n+1) = c;
        
        OUTPUT(:,1) = 0:n;
        OUTPUT(:,2) = ao;
        OUTPUT(:,3) = bo;
        OUTPUT(:,4) = co;
        if(Output_Flag == 1)
            OUTPUT(:,5) = 0.5*(bo - ao);
            dlmwrite('BISECTION_Q1_1.txt', OUTPUT);
        elseif(Output_Flag == 2)
            OUTPUT(:,5) = co - alpha;
            OUTPUT(:,6) = 0.5*(bo - ao);
            dlmwrite('BISECTION_OUTPUT_2.txt', OUTPUT);
        end
    end
end

function root = newton(x0, TOL, Output_Flag, alpha, Stopping_Flag, M)
    x(1) = x0;
    if (Stopping_Flag == 1)
        END = 30;
    else
        END = 1001;
    end
    
    fx(1) = abs(poly(x(1),0));
    error(1) = abs(x(1) - alpha);
    
    i = 2;
    while (i <= END)
        
        
        x(i) = x(i-1) - (poly(x(i-1), 0) / poly(x(i-1), 1));
%         if (abs( poly(x(i), 1) ) < TOL)
%             error('DERIVATIVE TOO CLOSE TO 0!')
%         end
        
        %OUTPUT INFO
        fx(i) = abs(poly(x(i-1),0));
        xdiff(i) = abs(x(i) - x(i-1));
        error(i) = abs(x(i) - alpha);
        
        if (i >= 4)
            rate(i) = log( abs( (x(i) - x(i-1))/(x(i-1) - x(i-2)) ) ) / log( abs( (x(i-1) - x(i-2))/(x(i-2) - x(i-3))));
        end
        
        if ((abs(x(i) - x(i-1)) + abs(poly(x(i), 0))) <= (TOL/5))
            root = x(i);
            END = i;
        else
            root = x(i);
            i = i + 1;
            if (i > 1000)
                error('>1000 ITERATIONS! CHOOSE BETTER X0')
            end
        end
    end
    i = i-1;
    x(i+1) = x(i) - (poly(x(i), 0) / poly(x(i), 1));
    rate(i+1) = log( abs( (x(i+1) - x(i))/(x(i) - x(i-1)) ) ) / log( abs( (x(i) - x(i-1))/(x(i-1) - x(i-2))));
    error(i+1) = abs(x(i+1) - alpha);
    
    OUTPUT(:,1) = 0:i ;
    OUTPUT(2:end,2) = xdiff;
    OUTPUT(:,3) = x;
    if (Output_Flag == 1)
        OUTPUT(:,4) = rate;
        dlmwrite('NEWTON_Q1_3.txt',OUTPUT);
    elseif (Output_Flag == 2)
        OUTPUT(:,4) = error;
        OUTPUT(:,5) = rate;
        dlmwrite('NEWTON_Q3_3.txt',OUTPUT);
    end
    
        
end

function root = secant(x0, x1, TOL, Output_Flag, alpha, Stopping_Flag, M)
    x(1) = x0;
    x(2) = x1;
    fx0 = poly(x0, 0);
    fx1 = poly(x1, 0);
    
    if (Stopping_Flag == 1)
        END = M;
    else
        END = 1000;
    end
    
    fx(1) = abs(poly(x(1),0));
    error(1) = abs(x(1) - alpha);
    
    i = 3;
    while (i < END)
        if (abs(poly(x(i-1), 0) - poly(x(i-2), 0)) < TOL*abs(x(i-1) - x(i-2)))
            error('DERIVATIVE TOO CLOSE TO 0!')
        end
        x(i) = x(i-1) - poly(x(i-1),0)*( (x(i-1) - x(i-2)) / (poly(x(i-1),0) - poly(x(i-2),0)));
        
        %OUTPUT INFO
        fx(i) = abs(poly(x(i-1),0));
        xdiff(i) = abs(x(i) - x(i-1));
        error(i) = abs(x(i) - alpha);
        if (i >= 4)
            rate(i) = log( abs( (x(i) - x(i-1))/(x(i-1) - x(i-2)) ) ) / log( abs( (x(i-1) - x(i-2))/(x(i-2) - x(i-3))));
        end
        
        if ((abs(x(i) - x(i-1)) + abs(poly(x(i), 0))) <= (TOL/5))
            root = x(i);
            END = i;
        else
            root = x(i);
            i = i + 1;
            if (i > 1000)
                error('>1000 ITERATIONS! CHOOSE BETTER X0')
            end
        end
    end
    
    x(i+1) = x(i) - (poly(x(i), 0) / poly(x(i), 1));
    rate(i+1) = log( abs( (x(i+1) - x(i))/(x(i) - x(i-1)) ) ) / log( abs( (x(i) - x(i-1))/(x(i-1) - x(i-2))));
    error(i+1) = abs(x(i+1) - alpha);
    
    OUTPUT(:,1) = 0:i ;
    OUTPUT(2:end,2) = xdiff;
    OUTPUT(:,3) = x;
    if (Output_Flag == 1)
        OUTPUT(:,4) = rate;
        dlmwrite('SECANT_Q1_4.txt',OUTPUT);
    elseif (Output_Flag == 2)
        OUTPUT(:,4) = error;
        OUTPUT(:,5) = rate;
        dlmwrite('SECANT_Q2_3.txt',OUTPUT);
    end
end


function root = hybrid(a, b, TOL, Output_Flag, alpha)
    x(1) = a;
    x(2) = b;
    aa(1) = a;
    bb(1) = b;
    k = 1;
    
    if (poly(a, 0)*poly(b, 0) > 0)
        error('BAD INTERVAL')
    end
    while ( (abs(x(k+1) - x(k)) + abs(poly(x(k+1), 0))) > (TOL/5) && k < 10)
    
        c = x(k+1) - ((poly(x(k+1), 0)*(x(k+1) - x(k))) / (poly(x(k+1), 0) - poly(x(k), 0)));
        
        if (c < aa(k) || c > bb(k))
            c = aa(k) + 0.5*(bb(k)-a);
            if (poly(c, 0)*poly(aa(k), 0) < 0)
                aa(k+1) = aa(k);
                bb(k+1) = c;
                x(k+2) = c;
                x(k+1) = aa(k);
            elseif (poly(c, 0)*poly(bb(k), 0) < 0)
                aa(k+1) = c;
                bb(k+1) = bb(k);
                x(k+2) = c;
                x(k+1) = bb(k);
            end
            
        else
            if(poly(aa(k),0)*poly(c,0) < 0)
                aa(k+1) = aa(k);
                bb(k+1) = c;
            end
            if(poly(c,0)*poly(bb(k),0) < 0)
                aa(k+1) = c;
                bb(k+1) = bb(k);
            end
            x(k+2) = c;
        end
        k = k+1;
    end
    root = x(k);
    
    OUTPUT(:,1) = 0:k-1 ;
    OUTPUT(:,2) = aa;
    OUTPUT(:,3) = bb;
    OUTPUT(:,4) = x(1:end-1);
    if (Output_Flag == 1)
        OUTPUT(:,5) = (1/2)*(bb-aa);
        dlmwrite('HYBRID_Q1_5.txt', OUTPUT);
    else
        OUTPUT(:,5) = abs(x(1:end-1) - alpha);
        OUTPUT(:,6) = (1/2)*(bb-aa);
        dlmwrite('HYBRID_Q3_4.txt', OUTPUT);
    end
end

function root = chord(x0, TOL, Output_Flag, alpha, Stopping_Flag, M)
    x(1) = x0;
    if (Stopping_Flag == 1)
        END = M;
    else
        END = 1001;
    end
    
    i = 2;
    while (i < END)
        
        x(i) = x(i-1) - (poly(x(i-1), 0) / poly(x(1), 1));
        if (abs(poly(x(i),1)) < TOL)
            error('DERIVATIVE TOO CLOSE TO 0!')
        end
        if ((abs(x(i) - x(i-1)) + abs(poly(x(i), 0))) <= (TOL/5))
            root = x(i);
            i = END;
        else
            root = x(i);
            i = i + 1;
            if (i > 1000)
                error('>1000 ITERATIONS! CHOOSE BETTER X0');
            end
        end
    end
    n = size(x);
    n = n(1,2);
    for i = 1:n
        fun(i) = abs(poly(x(i),0));
    end
    for i = 2:n
        xx(i) = abs(x(i) - x(i-1));
    end
    OUTPUT(:,1) = 0:n-1;
    OUTPUT(:,2) = fun;
    OUTPUT(:,3) = xx;
    OUTPUT(:,4) = x;
    OUTPUT(:,5) = abs(x-alpha);
    dlmwrite('CHORD_Q2_2.txt',OUTPUT);
end

function root = superHalley(x0, TOL, Output_Flag, alpha, Stopping_Flag, M)
    x(1) = x0;
    w(1) = (poly(x0, 0)*poly(x0, 2)) / (poly(x0, 1)^2);

    if (Stopping_Flag == 1)
        END = M;
    else
        END = 1001;
    end
    
    i = 2;
    while (i < END)
        x(i) = x(i-1) - ( 1 + (1/2)*(w(i-1)/(1-w(i-1))))*(poly(x(i-1),0)/poly(x(i-1),1));
        disp(abs(poly(x(i),1)))
        if(abs(poly(x(i), 1)) < TOL)
            disp(i);
            error('DERIVATIVE TOO CLOSE TO 0!')
        end
        w(i) = (poly(x(i), 0)*poly(x(i), 2)) / (poly(x(i), 1)^2);
        if( w(i) - 1 == 0)
            error('>1000 ITERATIONS! CHOOSE BETTER X0!')
        end
        if ((abs(x(i) - x(i-1)) + abs(poly(x(i), 0))) <= (TOL/5))
            root = x(i);
            i = END;
        else
            root = x(i);
            i = i + 1;
        end
    end
    n = size(x);
    n = n(1,2);
    for i = 1:n
        fun(i) = abs(poly(x(i),0));
    end
    for i = 2:n
        xx(i) = abs(x(i) - x(i-1));
    end
    for i = 4:n
        rate(i) = log( abs( (x(i) - x(i-1))/(x(i-1) - x(i-2)) ) ) / log( abs( (x(i-1) - x(i-2))/(x(i-2) - x(i-3))));
    end
    OUTPUT(:,1) = 0:n-1;
    OUTPUT(:,2) = fun;
    OUTPUT(:,3) = xx;
    OUTPUT(:,4) = x;
    
    if(Output_Flag == 1)
        OUTPUT(:,5) = rate;
        dlmwrite('HALLEY_Q1_6.txt',OUTPUT);
    else
        OUTPUT(:,5) = abs(x - alpha);
        dlmwrite('HALLEY_Q2_5.txt', OUTPUT);
    end
end

function fx = poly(x, p)
    if (p == 0)                       %%%% ORIGINAL FUNCTION
       % fx = 2-exp(x);
       fx = 10*x*(1/exp(x^2));
    elseif (p == 1)                   %%%% FIRST DERIVATIVE 
       % fx = (-1)*exp(x);
       fx = 10*(exp(-x^2) - 2*exp(-x^2)*x^2);
    elseif (p == 2)                   %%%% SECOND DERIVATIVE
        %fx = (-1)*exp(x);
       fx = 10*(4*exp(-x^2)*x^3 - 6*exp(-x^2)*x);
    end  
end

function func = func(x, p)
    X = [-525, 270, 61, -44 , 4];
    X = (1/100)*X;
    n = size(X);
    n = n(1,2);
    f(1) = X(n);
    for i = 2:n
        f(i) = x*f(i-1) + X(n-i+1);
    end
    if (p >= 1)
        df(1) = f(1);
        for i = 2:n-1
            df(i) = x*df(i-1) + f(i);
        end
        if (p >= 2)
            d2f = (n-1)*(n-2)*X(n);
            for k = n-1:-1:3
                d2f = (k-1)*(k-2)*X(k) + d2f*x;
            end
        end
    end 
    
    if (p == 0)
        func = f(n);
    elseif (p == 1)
        func = df(n-1);
    elseif (p == 2)
        func = d2f;
    end
end