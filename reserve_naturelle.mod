param n:= 10;
param p:= 6;
param c {k in 1..n, j in 1..n};
param proba {k in 1..4, j in 1..}

var x {i in 1..n} >= 0, <= 15;
var y >= 0;

minimize obj: y + 0.5 * (0.4 * (x[1] + x[2]) + 0.2 * x[3]);
subject to const_1: y + x[1] >= 10;
subject to const_2: y + x[2] >= 20;
subject to const_3: y + x[3] >= 30;



