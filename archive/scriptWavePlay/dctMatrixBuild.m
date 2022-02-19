d = 8;

A = zeros(d^2);
cnt = 0;
for i = 1:d
    for j = 1:d
        cnt = cnt+1;
        x = zeros(d);
        x(j,i) = 1;
        x = dct(x,[],1);
        x = dct(x,[],2);
        A(:,cnt) = x(:);
    end
end

B = zeros(d);
for i = 1:d
    x = zeros(d,1);
    x(i) = 1;
    x = dct(x);
    B(:,i) = x;
end

figure(901);
subplot(1,2,1);imagesc(A);
subplot(1,2,2);imagesc(B);        