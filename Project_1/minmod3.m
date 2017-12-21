function psi = minmod3(a, b, c)
% function psi = minmod(a, b)
% Purpose: Implement the minmod function, a,b are vectors containing 
% the values of q, i.e. we expect them to have length/N 2 for project 2
% Set up parameters and arrays
N = size(a,1); m = 3;  psi = zeros(N,1);
% Determine which entries have the same sign
s = (sign(a)+sign(b) + sign(c))/m; 
ids = find(abs(s)==1);
% If there are entries that have same signs
if(~isempty(ids))
  psi(ids) = s(ids).*min([abs(a(ids)), abs(b(ids)), abs(c(ids))],[], 2); 
end
return;