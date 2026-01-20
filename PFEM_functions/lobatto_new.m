%% matlab symbolic calc
% F = 1/sqrt(2/(2*k-1)) * int(sin(2*x), -1, pi/2) 



% get L_k function
clear L l; 
syms x 
L_0 = 1;
L{1} = x; 
l0 = (1-x)/2;
l{1} = (1+x)/2;
for k = 2:40
    Lk_1 = L{k-1} ;
    if k > 2
        Lk_2 = L{k-2} ;
    else
        Lk_2 = L_0 ;
    end
    
    L{k} = simplify((2*k-1)/k*x*Lk_1 - (k-1)/k*Lk_2);
    
    l{k} = simplify(1/sqrt(2/(2*k-1)) * int(L{k-1}, -1, x)) ;
end

%%
clc
for ii = 22 : 40
    
%     simplify(l{ii})
    simplify(diff(l{ii},x))
%     simplify(dlobatto ( ii , x ))
    
%    simplify(lobatto ( ii , x )-l{ii}) 
%    simplify(dlobatto ( ii , x )-diff(l{ii},x)) 
%    pause
end

