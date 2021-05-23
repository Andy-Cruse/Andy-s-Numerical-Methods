format long
%h = 1/2
%h = 1/4
%h = 1/8
%h = 1/16
h = 1/32
N = 1/h;
y = 2; %y_0
t = 1; %t_0
y_values(1) = y;
t_values(1) = t;
f = func(t, y);
y_t = solution(t);
y_t_values(1) = y_t;
disp(y_t)

eulerMethod(t, y, h, N);
RK4Method(t, y, h, N);

function g = eulerMethod(t, y, h, N)
    for k = 1:N
        yn = y + h*func(t, y);
        y = yn;
        t = t+h;
        y_values(k+1) = y; %y added to an array 
        t_values(k+1) = t; %t added to an array
        disp(y)
        disp(t)
        y_t = solution(t);
        y_t_values(k+1) = y_t; %actual solution added to an array
        disp(y_t)
        error = y_t - y;
        disp(error)
    end
    
    DATA(:,1) = t_values;
    DATA(:,2) = y_values;
    DATA(:,3) = y_t_values;
    DATA(:,4) = y_t_values - y_values;
    dlmwrite("DATAFILE.txt", DATA);
    plot(t_values, y_t_values, 'r') %Plots the real function in red 
    hold on 
    plot(t_values, y_values, 'b') %Plots the approximation in blue
    hold off
end

%RK4 Method written from pseudocode
function g_2 = RK4Method(t, y, h, N)
    for k = 1:N
        k1 = func(t, y);
        k2 = func(t+h/2, y+h*k1/2);
        k3 = func(t+h/2, y+h*k2/2);
        k4 = func(t+h, y+h*k3);
        yn = y + h*(k1+2*k2+2*k3+k4)/6;
        y = yn;
        t = t + h;
        disp(y)
        disp(t)
        y_t = solution(t);
        error = y_t - y;
        disp(error)
    end
end

%Actually solution too the problem
function y_t = solution(t)
    y_t = (1/2)*t^(-2)*(4+cos(2)-cos(2*t));
end

%Used in both Euler's Method and the RK4 Method to approximate solution
function f = func(t, y)
    f = t^(-2)*(sin(2*t) - 2*t*y);
end