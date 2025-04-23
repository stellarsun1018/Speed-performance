clear all
temp_res = 0.001;
range = 20;
x = -range:temp_res:range;
big_x = [x x(end)+temp_res:temp_res:x(end) + range*3];
hipspeedfit = -2;
viewscale = 1.5; % unused  becasue we are using percent distance
handscale = 3;
hip = 1./(1 + exp(hipspeedfit.*x));
slience_time = 1;
slience_duration = round(slience_time ./temp_res);
start_point = 1;
target_radius = 200; % doesn't change whether scope or not becasue passives are always in hip scale and motor in hand scale
hip_passive = 25;
scope_passive = 15;
sigma_motor_slope = 1;
switchtimes = -2:0.2:10;
non_reward_period = 10;
switchgain = NaN(1,length(switchtimes));
hip_gains = NaN(length(switchtimes),length(hip));
switch_gains = NaN(length(switchtimes),length(hip));

for i = 1:length(switchtimes)
    tic
    hswitch_time = switchtimes(i); 
    hswitch_point = round((range + hswitch_time)./temp_res +1);
    beforeslope = (hip(hswitch_point) - hip(hswitch_point-1)) ./ temp_res;
    viewaftery = hip(hswitch_point);
    viewafterslope = beforeslope./handscale;
    scopespeedfit = -viewafterslope./(viewaftery .*(1-viewaftery));
    sswitch_time = log(viewaftery/(1-viewaftery)) ./ -scopespeedfit;
    sswitch_point = round((range + sswitch_time)./temp_res +1);
    viewswitchpart = 1./(1 + exp(scopespeedfit.*big_x(sswitch_point:end) ));
    viewcomb = [hip(1:hswitch_point-1), viewswitchpart];
    imagination = 1./(1 + exp(scopespeedfit.*(x(1:hswitch_point)-(hswitch_time - sswitch_time))));
    handcomb = [hip(1:hswitch_point-1), viewswitchpart.*handscale - viewaftery.*(handscale-1)];
    handimagn = handscale./(1 + exp(scopespeedfit.*(x(1:hswitch_point)-(hswitch_time - sswitch_time))))- viewaftery.*(handscale-1);
%     handswitchpart = handscale./(1 + exp(scopespeedfit.*big_x(sswitch_point:end)));
%     afterslope = (handswitchpart(2) - handswitchpart(1)) ./ temp_res; 
%     comments above are error checking to ensure the math works correct up to this part
    
    hip_max_speed = abs(hipspeedfit)*0.25; % 0.25 is the max speed for standard sigmoid, an the derivative of a sigmoid is 0.25 multiplied by this parameter at its y midpoint
    imagn_max_speed = abs(scopespeedfit)*0.25*handscale;
    
    phit_hip = compute_phit(target_radius, hip_passive + hip_max_speed .*sigma_motor_slope, hip);
    phit_scope = compute_phit(target_radius, scope_passive + imagn_max_speed .*sigma_motor_slope, viewcomb);
    phit_switch = [phit_hip(1:hswitch_point) zeros(1,slience_duration-1) phit_scope(sswitch_point+slience_duration:end)];
    
    temp_gain = [linspace(1,0,(range*2-non_reward_period)/temp_res) zeros(1,range/temp_res*8)];
    e_gain_hip = phit_hip.*temp_gain(1:length(phit_hip));
    e_gain_switch = phit_switch.*temp_gain(1:length(phit_switch));
    switchgain(i) = max(e_gain_switch) - max(e_gain_hip);
    hip_gains(i,:) = e_gain_hip;
    switch_gains(i,:) = e_gain_switch(1:length(e_gain_hip));
%     plot(hip)
%     hold on
%     plot(viewcomb)
%     plot(imagination,'--')
%     hold off
%     xlabel('time')
%     ylabel('% screen distance travelled')
    i
    toc

end

%%
figure(1)
plot(hip,'LineWidth',2)
hold on
plot(viewcomb,'LineWidth',2)
plot(imagination,'--','LineWidth',2)
hold off
xlim([10000 40000])
    xticks(0:5000:40000);
    xticklabels(0:100:800);
xlabel('time')
ylabel('% screen distance travelled')
legend('R_{Hip}','R_{Switch}','R_{iADS}')

%%
figure(2)
plot(hip,'LineWidth',2)
hold on
plot(handcomb,'LineWidth',2)
plot(handimagn,'--','LineWidth',2)
hold off
xlim([10000 40000])
    xticks(0:5000:40000);
    xticklabels(0:100:800);
xlabel('time')
ylabel('hand distance travelled')
title('True Distance Travelled by Hand')
legend('R_{Hip}','R_{Switch}','R_{iADS}')
%%
figure(3)
plot(phit_hip,'o')
hold on
plot(phit_switch,'o')
hold off
xlim([20000 30000])
    xticks(20000:5000:30000);
    xticklabels(400:50:600);
xlabel('Time (ms)')
ylabel('Probability of Hit')
% title('P(Hit)/Time')
legend('P_{Hip}','P_{Switch}')
%%
figure(5)
pause(2)
for i = 1:size(hip_gains,1)
    plot(hip_gains(i,:),'o')
    hold on
    plot(switch_gains(i,:),'o')
    hold off
    xlabel('Time')
    xticks(0:5000:40000);
    xticklabels(0:100:800);
    ylabel('Expected Gain')
    ylim([0,0.5]);
    xline((switchtimes(i)+range)/temp_res,'--');
    title(['switch time = ' num2str((switchtimes(i)+range)./0.05)])
    legend('Hipshot','Scope')
    pause(0.1)
end
%%
figure(6)
plot(switchgain,'LineWidth',2)
yline(0,'--')
xlabel('Time of Switching')
ylabel('Switch Advantage')
title('Switch Advantage = max(E(scope)) - max(E(hipshot))')
ylim([-0.3,0.12])
xlim([-9,71])
xticks(1:10:81)
xticklabels(360:40:720)

