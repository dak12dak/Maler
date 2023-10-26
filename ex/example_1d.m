%%% sampling one-dimensional problem 
WORKSPACE.VARIABLES = {      'x'       };   % variables
WORKSPACE.GRIDPARAM = { [-pi, pi, 101] };   % domain and the vector length

WORKSPACE.PARAMETERS = {  'a' ,  'b' };     % parameters
WORKSPACE.PARSAMPVAL = { 1.32 , 1.06 };     % sampling values of parameters

syms x real       
syms a b real
     
f1 = sin(a*x+b);

f2 = diff(f1,x);

f3 = diff(f1,a);

f4 = diff(f1,b);

f5 = sqrt(f1^2 + f2^2);

% to use matrix operations specify each matrix element is a scalar function
f611 = x;
f612 = 1 + a*x;
f711 = x;
f721 = 1 - b*x;
f8   = simplify( [f611 , f612] * [f711 ; f721] ) ;

%%% functions definitions in a text form are also allowed
g = 'a*x^2+b';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   T H E   E N D   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%