
%%
subplot(1,2,1)
x = linspace(min(distances),max(distances),100);

k = 8:2:16;
tar_size = 20;
hold on
for i = 1:length(k)
    y = log2(2 .* x ./ tar_size) ./ k(i);
    plot(x,y,'--','LineWidth',2)
end

scatter(distances, durations,'o','MarkerEdgeColor','none','MarkerFaceColor', 'k', 'MarkerFaceAlpha',0.3)

hold off

legend('k = 8','k = 10','k = 12','k = 14','k = 16','Trials','Location', 'northeast');

subplot(1,2,2)

x = linspace(min(distances),max(distances),100);

k = 10;
tar_size = 12:4:32;
hold on
for i = 1:length(tar_size)
    y = log2(2 .* x ./ tar_size(i)) ./ k;
    plot(x,y,'--','LineWidth',2)
end

scatter(distances, durations,'o','MarkerEdgeColor','none','MarkerFaceColor', 'k', 'MarkerFaceAlpha',0.3)

hold off

legend('size = 6','size = 10','size = 14','size = 18','size = 22','Trials','Location', 'northeast');




%%

