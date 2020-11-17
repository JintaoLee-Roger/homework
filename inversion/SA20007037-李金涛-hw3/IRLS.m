function m = IRLS(t, G)
% copyright
% author: Jintao Li


condition = true;
threshold = 0.0001;
epsilon = 0.0001;

m0 = inv(G' * G) * G' * t;

r = abs(t - G * m0);
r(r < epsilon) = epsilon;
R = diag(r.^-1);

i = 0;
while condition
    i = i+1;
    m1 = (G' * R * G) \ (G' * R * t);
    condition = (norm(m1 - m0) ./ (1 + norm(m1))) > threshold;
    m0 = m1;
    r = abs(t - G * m0);
    r(r < epsilon) = epsilon;
    R = diag(r.^-1);
end

m = m0;