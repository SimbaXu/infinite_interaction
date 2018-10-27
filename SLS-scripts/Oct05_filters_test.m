[b4, a4] = butter(3, 80, 'high', 's');
[b5, a5] = cheby1(3, 3, 80, 'high', 's')

bode(tf(b4, a4), tf(b5, a5))
legend('butter', 'cheby1')
grid