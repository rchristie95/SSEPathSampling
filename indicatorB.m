function y = indicatorB(x,U,L)
  y = (x > L) & (x < U);  % 1 if x is between 0 and 10, 0 otherwise
end