%%% sampling two-dimensional problem
WORKSPACE.VARIABLES = {    'x'     ;    'y';    };   % variables
WORKSPACE.GRIDPARAM = { [-1,1,201] ; [ 0,1,101] };   % [min, max, length]

WORKSPACE.PARAMETERS  = { 'a' , 'b' , 'c' , 'd' };   % parameters
WORKSPACE.PARSAMPVAL  = { -1  ,  0  , 1.3 , 1.3 };   % sampling values 

syms x y real       
syms a b c d real 

    
f1 = a*x + b*y;

f2 = a*x^2 + b*x*y + c*y^2;

f3 = sin(f1)*cos(d*f2);

f4 = atan(f1/f2); %in numerics, function 'atan' is ever replaced by 'atan2'

f5 = tanh(f1*f2);

f6 = log(1+(a*x+b*y)^2);

f7 = exp(-c*x*y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   T H E   E N D   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%