P = Oct02_SISO_model();
K = tf([1], [1, 6, 5]);
P = ss(P);
K = ss(K);

pk = lft(P, K);

A = P.A;
B = P.B;
C = P.C;
D = P.D;

Ak = K.A;
Bk = K.B;
Ck = K.C;
Dk = K.D;

Alft = pk.A;
Blft = pk.B;
Clft = pk.C;
Dlft = pk.D;

save('Oct24_testlft_data.mat', 'A', 'B', 'C', 'D', 'Ak', 'Bk', 'Ck', 'Dk', 'Alft', 'Blft', 'Clft', 'Dlft')

[mag, phase] = bode(pk, [1])