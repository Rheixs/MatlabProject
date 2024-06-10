% Define the function
str_f = input('Enter the function (e.g., "exp(-x) - x"): ', 's');
f = str2func(['@(x)', str_f]);

xl = input('Enter the lower bound (e.g., -2): ');
if ~isnumeric(xl)
    error('The lower bound must be a number.');
end

delta_x = input('Enter the delta x (e.g., 0.01): ');
if ~isnumeric(delta_x) || delta_x <= 0
    error('Delta x must be a positive number.');
end

x = xl:delta_x:2;


try 
    disp('Options:');
    disp('1. Bisection Method');
    disp('2. Incremental Method');
    disp('3. Graphical Method');
    disp('4. False Position Method');
    disp('5. Simple Fixed Point Method');
    disp('6. Newton-Raphson Method');
    disp('7. Secant Method');
    
    method = input('Enter the number of the method to use: ');
    if ~ismember(method, 1:7)
        error('Method must be a number from 1 to 7.');
    end

root = NaN; 

if method == 1
    % Bisection Method
    a = xl; 
    b = 2; 
    tol = 1e-6; 
    i = 0; 
    
    A = [];
    B = [];
    DeltaX = [];
    FXL = [];
    FXMU = [];
    FXLFXMU = [];

    while abs(b-a) > tol
        c = (a + b) / 2;
        deltaX = abs(b - a);
        fxL = f(a);
        fxMU = f(c);
        fxlfxmu = fxL * fxMU;
        A = [A; a];
        B = [B; b];
        DeltaX = [DeltaX; deltaX];
        FXL = [FXL; fxL];
        FXMU = [FXMU; fxMU];
        FXLFXMU = [FXLFXMU; fxlfxmu];
        if f(c) == 0
            break;
        elseif f(a)*f(c) < 0
            b = c;
        else
            a = c;
        end
        i = i + 1;
    end

    T = table(A, B, DeltaX, FXL, FXMU, FXLFXMU);
    disp(T);

    root = (a + b) / 2;
    disp(['Root found with Bisection Method: ', num2str(root)]);

elseif method == 2
    % Incremental Method
    if delta_x > 0.01
        disp('Warning: delta x is too large for Incremental Method. Setting delta x to 0.01.');
        delta_x = 0.01;
    end
    [min_val, min_index] = min(abs(f(xl:delta_x:2)));
    root = xl + delta_x * (min_index - 1);
    if isnan(root)
        disp('No root found with Incremental Method');
    else
        disp(['Root found with Incremental Method: ', num2str(root)]);
    end

elseif method == 3
    % Graphical Method
    [min_val, min_index] = min(abs(f(x)));
    root = x(min_index);
    disp(['Estimated root from Graphical Method: ', num2str(root)]);

elseif method == 4
    % False Position Method
    a = xl;
    b = 2;
    tol = 1e-6;
    while abs(b-a) > tol
        c = a - ((b-a)/(f(b)-f(a)))*f(a);
        if f(c) == 0
            break;
        elseif f(a)*f(c) < 0
            b = c;
        else
            a = c;
        end
    end
    root = a - ((b-a)/(f(b)-f(a)))*f(a);
    disp(['Root found with False Position Method: ', num2str(root)]);

elseif method == 5
    % Simple Fixed Point Method
    g = @(x) x - f(x);
    x0 = xl;
    tol = 1e-6;
    while abs(g(x0) - x0) > tol
        x0 = g(x0);
    end
    root = g(x0);
    disp(['Root found with Simple Fixed Point Method: ', num2str(root)]);

elseif method == 6
    % Newton-Raphson Method
    df = @(x) (f(x+tol) - f(x)) / tol; % derivative of f
    x0 = xl;
    while abs(f(x0)) > tol
        x0 = x0 - f(x0)/df(x0);
    end
    root = x0;
    disp(['Root found with Newton-Raphson Method: ', num2str(root)]);

elseif method == 7
    % Secant Method
    x0 = xl;
    x1 = 2;
    tol = 1e-6;
    while abs(x1 - x0) > tol
        x_temp = x1 - f(x1)*((x1 - x0)/(f(x1) - f(x0)));
        x0 = x1;
        x1 = x_temp;
    end
    root = x1;
    disp(['Root found with Secant Method: ', num2str(root)]);
end
catch 
    fprintf('An error occurred while finding the root: %s\n', ME.message);
end

if ~isnan(root)
    figure;
    plot(x, f(x), 'b', root, f(root), 'ro');
    title('Root Finding Method');
    xlabel('x');
    ylabel('f(x)');
    grid on;
end