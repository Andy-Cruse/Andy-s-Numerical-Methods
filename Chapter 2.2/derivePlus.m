format long
x = 1;
h = [1/2 1/4 1/8 1/16 1/32 1/64];
fx1 = (x+1)^(1/2);
fx2 = exp(x);
deriveOne(x, h, fx1, fx2);
deriveTwo(x, h, fx1, fx2);
deriveThree(x, h, fx1, fx2);
deriveFour(x, h, fx1, fx2);

% Function to compute the D+f(x) derivative approximation
% and the error rates of it compared to the original 
% derivative
function derivePlus = deriveOne(x, h, fx1, fx2)
disp("D^+(x+1)^(1/2): ") 
fprime1 = ( (2)^(1/2) / 4 ); %actual derivative value
for i = 1:6
    disp(h(1))
    d_fx1 = ( (x+h(i)+1)^(1/2) - fx1)/h(i); %D+f(x) approximation
    disp(d_fx1)
    error =  fprime1 - (d_fx1); %substracts D^+ from f'(1)
    disp(error)
    d_2h_fx1 = ( (x+(2*h(i))+1)^(1/2) - fx1)/(2*h(i));
    errorRate = log(abs( (fprime1 - d_2h_fx1) / error)) / log(2);
    disp(errorRate)
end

disp("D^+e^x: ")
fprime2 = exp(1);
for i = 1:6
    d_fx2 = ( exp(x+h(i)) - fx2 ) / h(i);
    disp(d_fx2)
    error = fprime2 - d_fx2;
    disp(error)
    d_2h_fx2 = ( exp(x+(2*h(i))) - fx2 ) / (2*h(i));
    errorRate = log(abs( (fprime2 - d_2h_fx2) / error)) / log(2);
    disp(errorRate)
end
end

% Function to compute the D-f(x) derivative approximation
% and the error rates of it compared to the original 
% derivative
function deriveMinus = deriveTwo(x, h, fx1, fx2)
disp("D^-(x+1)^(1/2): ")
fprime1 = ( (2)^(1/2) / 4 );
for i = 1:6
    disp(h(i))
    d_fx1 = ( fx1 - (x-h(i)+1)^(1/2))/h(i);
    disp(d_fx1)
    error =  fprime1 - (d_fx1); 
    disp(error)
    d_2h_fx1 = (fx1 - (x-(2*h(i))+1)^(1/2))/(2*h(i));
    errorRate = log(abs( (fprime1 - d_2h_fx1) / error)) / log(2);
    disp(errorRate)
end

disp("D^-e^x: ")
fprime2 = exp(1);
for i = 1:6
    d_fx2 = (fx2 - exp(x-h(i)) ) / h(i);
    disp(d_fx2)
    error = fprime2 - d_fx2;
    disp(error)
    d_2h_fx2 = (fx2 - exp(x-(2*h(i))) ) / (2*h(i));
    errorRate = log(abs( (fprime2 - d_2h_fx2) / error)) / log(2);
    disp(errorRate)
end
end

% Function to compute the Df(x) derivative approximation
% and the error rates of it compared to the original 
% derivative
function deriveMid = deriveThree(x, h, fx1, fx2)
disp("D(x+1)^(1/2): ")
fprime1 = ( (2)^(1/2) / 4 );
for i = 1:6
    disp(h(i))
    d_fx1 = ( (x+h(i)+1)^(1/2) - (x-h(i)+1)^(1/2))/(2*h(i));
    disp(d_fx1)
    error =  fprime1 - (d_fx1); 
    disp(error)
    d_2h_fx1 = (( (x+(2*h(i))+1)^(1/2) - (x-(2*h(i))+1)^(1/2))/(4*h(i)));
    errorRate = log(abs( (fprime1 - d_2h_fx1) / error)) / log(2);
    disp(errorRate)
end

disp("De^x: ")
fprime2 = exp(1);
for i = 1:6
    d_fx2 = (exp(x+h(i)) - exp(x-h(i)) ) / (2*h(i));
    disp(d_fx2)
    error = fprime2 - d_fx2;
    disp(error)
    d_2h_fx2 = (exp(x+(2*h(i))) - exp(x-(2*h(i))) ) / (4*h(i));
    errorRate = log(abs( (fprime2 - d_2h_fx2) / error)) / log(2);
    disp(errorRate)
end
end

% Function to compute the D~+f(x) derivative approximation
% and the error rates of it compared to the original 
% derivative
function derive2h = deriveFour(x, h, fx1, fx2)
disp("D(x+1)^(1/2): ")
fprime1 = ( (2)^(1/2) / 4 );
for i = 2:6
    d2_fx1 = 2*( (x+h(i)+1)^(1/2) - fx1)/h(i);
    d_2h_fx1 = ( (x+(2*h(i))+1)^(1/2) - fx1)/(2*h(i));
    d_tilde = d2_fx1 - d_2h_fx1;
    disp(d_tilde)
    error = fprime1 - d_tilde;
    disp(error)
    errorRate = log(abs( (fprime1 - (2*( (x+(2*h(i))+1)^(1/2) - fx1)/(2*h(i)) - ( (x+(4*h(i))+1)^(1/2) - fx1)/(4*h(i)) )) / error)) / log(2);
    disp(errorRate)
end

disp("f(x) = e^x")
fprime2 = exp(1); %actual derivative value at x=1
for i = 2:6 %h = 1/2 not needed
    d2_fx2 = (( exp(x+h(i)) - fx2 ) / h(i))*2; %2D+f(x)
    d_2h_fx2 = ( exp(x+(2*h(i))) - fx2 ) / (2*h(i)); %D+2hf(x)
    d_fx2 = d2_fx2 - d_2h_fx2; %D~+f(x)
    disp(d_fx2)
    error = fprime2 - d_fx2;
    disp(error)
    errorRate = log(abs( (fprime2 - ( ((( exp(x+(2*h(i))) - fx2 ) / (2*h(i)))*2) - (( exp(x+(4*h(i))) - fx2 ) / (4*h(i))) )  ) / error)) / log(2);
    disp(errorRate)
end
end
