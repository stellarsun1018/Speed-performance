for cond = 1:2
    if cond == 1
        part = copy(1:60,:);
    elseif cond == 2
        part = copy(61:120,:);
    end
    spdlandscape = NaN(5,3);
    sizes = unique(part(:,15));
    for i = 1:length(sizes)
        for j = 1:3
%             spdlandscape(i,j) = mean(part(part(:,15)==sizes(i)&round(part(:,10),-2)==j*100,17));
        	spdlandscape(i,j) = sum(0 == (part(part(:,15)==sizes(i)&round(part(:,10),-2)==j*100,26)));

        end
    end
    figure
    heatmap(spdlandscape);
    ylabel('Target radius mm')
    xlabel('Distance mm')
    
    title('Speed mm/s')
end
%%
for i = 1:120
    plot(validTraX(i,:),validTraY(i,:),'o-')
    hold on
end