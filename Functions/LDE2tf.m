function [G11 G12 G21 G22] = LDE2tf(A,B,C)
    s = tf('s');
    I = eye(4,4);
    G_s = C*(s*I - A)^-1*B;
    % G(s): From output i to input j
    G11 = tf(G_s.Numerator{1,1}(:,:), G_s.Denominator{1,1}(:,:));
    G12 = tf(G_s.Numerator{1,2}(:,:), G_s.Denominator{1,2}(:,:));
    G21 = tf(G_s.Numerator{2,1}(:,:), G_s.Denominator{2,1}(:,:));
    G22 = tf(G_s.Numerator{2,2}(:,:), G_s.Denominator{2,2}(:,:));
end