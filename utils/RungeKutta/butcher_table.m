function T = butcher_table(ID, num_flag)
%BUTCHER_TABLE Get RK Butcher table
%
%INPUT:
%   id          :   String: Identification of the RK method
%   num_flag    :   INT. OPTIONAL. 
%                       Numerically evaluate tables:
%                           1: Numerically evaluate Butcher table
%                           0: Symbolically evaluate Butcher table
%
%OUTPUT:
%   T   :   Butcher table in the form of a struct with fields A, b, c.
%
%NOTES:
% I only have the symbolic form implemented for some schemes. Ones that
% have only rational components can be converted into symbolic form after
% being returned numerically, i.e. by using sym(T.A), for example. This is
% not true for those with irrational components, like SDIRK2, for example.

% By default, numerically evaluate Butcher tables.
if nargin == 1
    num_flag = 1;
else
    if num_flag ~= 0
        num_flag = 1;
    end
end

%% --- Explicit methods --- %%

% Forward/explicit Euler
if strcmp(ID, 'FE')
    T.c = [0];
    T.b = [1];
    T.A = [0];
    
% This is a 2-stage method. It's always at least first-order accurate, and
% it's second-order accurate i.f.f b = (1 - 1/(2*a)), for any value of a.
elseif strcmp(ID, 'ERK1(2)')
    a = 1;
    b = 0.75*(1 - 1/(2*a));
    T.c = [0; a];
    T.b = [b; 1-b];
    T.A = [0, 0;
           a, 0];    

% Optimal 2nd-order strong stability preserving scheme.     
elseif strcmp(ID, 'SSPERK2')
    T.c = [0; 1];
    T.b = [1/2; 1/2];
    T.A = [0,  0;
           1,  0];    
    
% 2nd-order variant of Ralston's method    
elseif strcmp(ID, 'Ralston2')
    T.c = [0; 2/3];
    T.b = [1/4; 3/4];
    T.A = [0,  0;
          2/3, 0];
    
% 2nd-order variant of Heun's method      
elseif strcmp(ID, 'Heun2')
    T.c = [0; 1];
    T.b = [1/2; 1/2];
    T.A = [0, 0;
           1, 0]; 
       
elseif strcmp(ID, 'Kutta3')
    T.c = [0; 1/2; 1];
    T.b = [1/6; 2/3; 1/6];
    T.A = [ 0,   0,  0;
           1/2,  0,  0;
           -1,   2,  0]; 
       
% 3rd-order variant of Huen's method.       
elseif strcmp(ID, 'Heun3')
    T.c = [0; 1/3; 2/3];
    T.b = [1/4; 0; 3/4];
    T.A = [ 0,    0,   0;
           1/3,   0,   0;
            0,   2/3,  0]; 
       
% 3rd-order variant of Ralston's method.        
elseif strcmp(ID, 'Ralston3')
    T.c = [ 0;   1/2; 3/4];
    T.b = [2/9;  1/3; 4/9];
    T.A = [ 0,    0,   0;
           1/2,   0,   0;
            0,   3/4,  0]; 
       
% Optimal 3rd-order strong stability preserving scheme.         
elseif strcmp(ID, 'SSPERK3')
    T.c = [ 0;   1;  1/2];
    T.b = [1/6; 1/6; 2/3];
    T.A = [ 0,    0,   0;
            1,    0,   0;
           1/4,  1/4,  0]; 
       
% Optimal 3rd-order strong stability preserving scheme with 4 stages.         
% See (23) of https://arxiv.org/pdf/1809.04807.pdf
elseif strcmp(ID, 'SSPERK3(4)')
    T.c = [ 0;   1/2;  1; 1/2];
    T.b = [1/6; 1/6; 1/6; 1/2];
    T.A = [ 0,   0,    0,   0;
           1/2,  0,    0,   0;
           1/2,  1/2,  0,   0
           1/6,  1/6,  1/6, 0];        
    
% Runge's 3rd-order, 4-stage method. See: https://math.stackexchange.com/questions/3675436/intuitive-understanding-of-third-order-runge-kutta-method      
elseif strcmp(ID, 'Runge3(4)')
    T.c = [ 0;   1/2;  1; 1];
    T.b = [1/6; 4/6; 0; 1/6];
    T.A = [ 0,      0,   0,   0;
            1/2,    0,   0,   0;
            1,      0,   0,   0
            0,      0,   1,   0];        
       
% Classic RK4 method.   
elseif strcmp(ID, 'RK4classic')
    T.c = [0;   1/2; 1/2;   1];
    T.b = [1/6; 1/3; 1/3; 1/6];
    T.A = [0,   0,  0, 0;
          1/2,  0,  0, 0;
           0,  1/2, 0, 0;
           0,   0,  1, 0];
       
%6-stage, 5th-order method. (See Butcher (2008), p.99, eqn. (236a)). For order 5 schemes, at least 6 stages are required.       
elseif strcmp(ID, 'RK5')       
    T.c = [0;  1/4;   1/4;  1/2;    3/4;    1];
    T.b = [7/90; 0; 32/90; 12/90; 32/90; 7/90];
    T.A = [ 0,     0,   0,    0,    0,  0;
           1/4,    0,   0,    0,    0,  0;
           1/8,   1/8,  0,    0,    0,  0;
            0,     0,  1/2,   0,    0,  0;
           3/16, -3/8, 3/8,  9/16,  0,  0;
           -3/7,  8/7, 6/7, -12/7, 8/7, 0];
       
% %7-stage, 6th-order method. (See Butcher (2008), p.194)
% elseif strcmp(ID, 'RK6')       
%     T.c = [0;         1/3;   2/3;     1/3;     5/6;   1/6;   1];
%     T.b = [13/200;     0;   11/40;   11/40;   4/25;  4/25; 13/200];
%     T.A = [  0,        0,      0,      0,       0,     0,    0;
%             1/3,       0,      0,      0,       0,     0,    0;
%              0,       2/3,     0,      0,       0,     0,    0;
%             1/12,     1/3,   -1/12,    0,       0,     0,    0;
%            25/48,   -55/24,  35/48,   15/8,     0,     0,    0;
%             3/20,   -11/24,  -1/8,     1/2,    1/10,   0,    0;       
%           -261/260,  33/13,  43/156, -118/39, 32/195, 80/39, 0];

%7-stage, 6th-order method. (See Butcher (2008), p.195)
elseif strcmp(ID, 'RK6')       
    T.c = [  0;          2/5;     4/5;         2/9;     8/15;   0;   1];
    T.b = [  0;           0;    1375/4992;   6561/20384;  3375/12544;  53/768; 19/294];
    T.A = [  0,           0,        0,          0,            0,          0,      0;
            2/5,          0,        0,          0,            0,          0,      0;
             0,          4/5,       0,          0,            0,          0,      0;
           169/1458,   110/729,  -65/1458,      0,            0,          0,      0;
           -44/675,    -88/135,   76/351,     336/325,        0,          0,      0;
            21/106,       0,     -105/689,    -324/689,     45/106,       0,      0;       
          -2517/4864,  -55/38,  10615/31616,  567/7904,    7245/4864,  2597/2432, 0];
       
%9-stage, 7th-order method. (See Butcher (2008), p.196)
elseif strcmp(ID, 'RK7')        
    
    T.c = [0; 1/6; 1/3; 1/2;    2/11;            2/3;      6/7;         0;       1     ];
    T.b = [0; 0;   0;   32/105; 1771561/6289920; 243/2560; 16807/74880; 77/1440; 11/270];
    
    T.A = [  0,         0,    0,         0,            0,            0,         0,         0,     0; ...
             1/6,       0,    0,         0,            0,            0,         0,         0,     0; ...
             0,         1/3,  0,         0,            0,            0,         0,         0,     0; ...
             1/8,       0,    3/8,       0,            0,            0,         0,         0,     0; ...
             148/1331,  0,    150/1331, -56/1331,      0,            0,         0,         0,     0; ...
            -404/243,   0,   -170/27,    4024/1701,    10648/1701,   0,         0,         0,     0; ...
             2466/2401, 0,    1242/343, -19176/16807, -51909/16807,  1053/2401, 0,         0,     0; ...
             5/154,     0,    0,         96/539,      -1815/20384,  -405/2464,  49/1144,   0,     0; ...
            -113/32,    0,   -195/22,    32/7,         29403/3584,  -729/512,   1029/1408, 21/16, 0];
    
    
      
%% --- General implicit methods --- %%  

% Backward/implicit Euler
elseif strcmp(ID, 'BE') || strcmp(ID, 'SDIRK1')
    T.c = [1];
    T.b = [1];
    T.A = [1];
    
% Implicit midpoint rule    
elseif strcmp(ID, 'midpoint') 
    T.c = [1/2];
    T.b = [1];
    T.A = [1/2];
    
% Crank--Nicolson AKA trapezoidal rule
elseif strcmp(ID, 'CN')
    T.c = [0; 1];
    T.b = [1/2; 1/2];
    T.A = [ 0,   0;
           1/2, 1/2];
    
%% --- DIRK methods --- %%       
       
% % SDIRK2 from Dobrev et al. (2017)
% elseif strcmp(ID, 'SDIRK2')
%     if num_flag
%         alpha = 1/sqrt(2);
%     else
%         alpha = 1/sqrt(sym(2));
%     end
%     T.c = [1-alpha; alpha];
%     T.b = [  1/2;     1/2];
%     T.A = [ 1-alpha,     0;
%            2*alpha-1, 1-alpha];
% SDIRK2 from Dobrev et al. (2017)
% See Butcher (2008), p. 261
elseif strcmp(ID, 'SDIRK2')
    if num_flag
        alpha = sqrt(2)/2;
    else
        alpha = 1/sqrt(sym(2));
    end
    T.c = [1-alpha; 1];
    T.b = [  alpha; 1-alpha];
    T.A = [1-alpha,     0;
             alpha, 1-alpha];

% 3rd-order, 2-stage SDIRK. Is only A-stable. From Kennedy and Carpenter
% (223)
elseif strcmp(ID, 'SDIRK3(2)')
    gamma = (3 + sqrt(3))/6;
    T.c = [gamma;     1-gamma];
    T.b = [0.5;       0.5];
    T.A = [gamma,     0,  ;
           1-2*gamma, gamma];         
         
% SDIRK3 from Dobrev et al. (2017). Updated numbers from WIKI page.
% Actually, see Butcher (2008), p. 262
elseif strcmp(ID, 'SDIRK3')
    x = 0.43586652150845899942;
    
    p = 0.5*(1+x);
    q = 0.5*(1-x);
    y = -3/2*x^2 + 4*x - 1/4;
    z =  3/2*x^2 - 5*x + 5/4;
    
    T.c = [x;  p;  1];
    T.b = [y;  z;  x];
    T.A = [x,  0,  0;
           q,  x,  0;
           y,  z,  x];
          
% SDIRK4: L-stable order 4. See table 6.5, pg. 100 of Hairer and Wanner.
elseif strcmp(ID, 'SDIRK4')

    T.c = [  1/4;    3/4;  11/20;    1/2;   1];
    T.b = [25/24; -49/48; 125/16; -85/12; 1/4];
    T.A = [   1/4,        0,      0,       0,    0;
              1/2,       1/4,     0,       0,    0;
             17/50,     -1/25,   1/4,      0,    0;
           371/1360, -137/2720, 15/544,   1/4,   0;
             25/24,    -49/48,  125/16, -85/12, 1/4];       
 
% SDIRK5: See Kennedy and Carpenter, Table 24, page 98. The original
% authors claim this is L-stable, but Kennedy and Carpenter seem to
% disagree with this.
elseif strcmp(ID, 'SDIRK5')     
    
    d = 4024571134387/14474071345096;
    
    T.c = [ 4024571134387/14474071345096; 5555633399575/5431021154178;  5255299487392/12852514622453; 3/20;                         10449500210709/14474071345096];
    T.b = [-2522702558582/12162329469185; 1018267903655/12907234417901; 4542392826351/13702606430957; 5001116467727/12224457745473; 1509636094297/3891594770934];
    T.A = [ d,                             0,                            0,                            0,                           0; ...
            9365021263232/12572342979331,  d,                            0,                            0,                           0;
            2144716224527/9320917548702,  -397905335951/4008788611757,   d,                            0,                           0; ...
           -291541413000/6267936762551,    226761949132/4473940808273,  -1282248297070/9697416712681,  d,                           0; ...
           -2481679516057/4626464057815,  -197112422687/6604378783090,   3952887910906/9713059315593,  4906835613583/8134926921134, d]; 
         
% I don't quite think this is right... Also looking at the paper this came 
% from, I don't understand what's happening. But this is not a 5th-order method...         
% % SDIRK5: A-stable order 5. See exercise 4, pg. 101 of Hairer and Wanner.
% elseif strcmp(ID, 'SDIRK5')
% 
%     b = sqrt(6);
%     a = (6-b)/10;
%     
%     T.c = [a; (6+9*b)/35; 1; (4-b)/10; (4+b)/10];
%     T.b = [0; 0; 1/9; (16-b)/36; (16+b)/36];
%     T.A = [  a,                    0,                     0,                 0,               0;
%              (-6+5*b)/14,          a,                     0,                 0,               0;
%              (888+607*b)/2850,     (126-161*b)/1425,      a,                 0,               0;
%              (3153-3082*b)/14250,  (3213+1148*b)/28500,   (-267+88*b)/500    a,               0;
%              (-32583+14638)/71250, (-17199+364*b)/142500, (1329-544*b)/2500, (-96+131*b)/625, a];                
         
else
    error('RK butcher table id %s not recognised', ID)

end