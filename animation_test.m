

for i = 1:100
    
    y(i) = rand;
    
end

x_vec = 1:100;
figure
xlim([0 100])
ylim([0 1])
for i = 1:100
    plot(x_vec(1:i),y(1:i))
    xlim([0 100])
    ylim([0 1])
    pause(.1)
end
