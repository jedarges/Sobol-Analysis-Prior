function y = leg_eval(n,x)
%% 1-D Legendre polynomials
% Implemented by John Darges
% n - order of Legendre polynomial
% x - Value to evaluate Legendre polynomial at
switch n
    case 0
        y = 1;
    case 1
        y = x;
    case 2
        y = (1/2)*(3*(x^2) - 1);
    case 3
        y = (1/2) * (5*(x^3) - 3*x);
    case 4
        y = (1/8) * (35*(x^4) - 30*(x^2) + 3);
    case 5
        y = (1/8) * (63*(x^5) - 70*(x^3) + 15*x);
    case 6
        y = (1/16) * (231*(x^6) - 315*(x^4) + 105*(x^2) - 5);
    case 7
        y = (1/16) * (429*(x^7) - 693*(x^5) + 315*(x^3) - 35*x);
    case 8
        y = (1/128) * (6435*(x^8) - 12012*(x^6) + 6930*(x^4) - 1260*(x^2) + 35);
    case 9
        y = (1/128) * (12155*(x^9) - 25740*(x^7) + 18018*(x^5) - 4620*(x^3) + 315*x);
    case 10
        y = (1/256) * (46189*(x^10) - 109395*(x^8) + 90090*(x^6) - 30030*(x^4) + 3465*(x^2) - 63);
end