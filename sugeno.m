%Programmer: Yousef Sharafi

clc;
clear all;
close all;

data1=xlsread('mackey-glass.xlsx');

max_data=max(data1(1,:));
min_data=min(data1(1,:));
% data = (data - min(data))/(max(data)-min(data));
data=[data1(:,1) data1(:,2) data1(:,3) data1(:,4) data1(:,5)] ;
target=[data1(:,6)];
epoch=7;
n=size(data,1);
number_train=round(0.75*n);
number_test=n-number_train;
number_mf=3;
Mean_mf_x=linspace(0.2,0.5,number_mf);
Mean_mf_y=linspace(0.2,0.5,number_mf);
Mean_mf_z=linspace(0.2,0.5,number_mf);
Mean_mf_w=linspace(0.2,0.5,number_mf);
Mean_mf_r=linspace(0.2,0.5,number_mf);
Sigma_mf_x=unifrnd(0.1,0.2,[1,number_mf]);
Sigma_mf_y=unifrnd(0.1,0.2,[1,number_mf]);
Sigma_mf_z=unifrnd(0.1,0.2,[1,number_mf]);
Sigma_mf_w=unifrnd(0.1,0.2,[1,number_mf]);
Sigma_mf_r=unifrnd(0.1,0.2,[1,number_mf]);
rule=zeros(1,number_mf*number_mf*number_mf*number_mf*number_mf);
eta=0.16;
eta_m=0.01;
% number_parameter=3*number_mf;
f1=zeros(1,number_mf*number_mf*number_mf*number_mf*number_mf);
parameter_c=unifrnd(-0.01,0.01,[number_mf*number_mf*number_mf*number_mf*number_mf,6]);
parameter_c_jacobe=zeros(number_train,number_mf*number_mf*number_mf*number_mf*number_mf*6);
I=eye(number_mf*number_mf*number_mf*number_mf*number_mf*6);
error_jacobe=zeros(number_train,1);
for iter=1:epoch
    parameter_c_temp=zeros(number_mf*number_mf*number_mf*number_mf*number_mf,6);
   for i=1:number_train     
       x=data(i,1);
       y=data(i,2);
       z=data(i,3);
       w=data(i,4);
       rr=data(i,5);
       x_input=unifrnd(x,x,[1,number_mf]);
       y_input=unifrnd(y,y,[1,number_mf]);
       z_input=unifrnd(z,z,[1,number_mf]);
       w_input=unifrnd(w,w,[1,number_mf]);
       r_input=unifrnd(rr,rr,[1,number_mf]);
       sub_rule_x=exp(-0.5*((x_input-Mean_mf_x)./(Sigma_mf_x)).^2); 
       sub_rule_y=exp(-0.5*((y_input-Mean_mf_y)./(Sigma_mf_y)).^2);
       sub_rule_z=exp(-0.5*((z_input-Mean_mf_z)./(Sigma_mf_z)).^2);
       sub_rule_w=exp(-0.5*((w_input-Mean_mf_w)./(Sigma_mf_w)).^2);
       sub_rule_r=exp(-0.5*((r_input-Mean_mf_r)./(Sigma_mf_r)).^2);
       c=1;
       for p=1:number_mf
         for r=1:number_mf 
             for i3=1:number_mf 
                 for i4=1:number_mf 
                     for i5=1:number_mf 
                        rule(c)= sub_rule_x(p)*sub_rule_y(r)*sub_rule_z(i3)*sub_rule_w(i4)*sub_rule_r(i5);
                        c=c+1;
                     end
                 end
             end
         end
       end
       sum_rule=sum(rule);
       rule=rule/sum_rule;
       output=0;
       ss=1;
       for tr=1:number_mf*number_mf*number_mf*number_mf*number_mf
           f1(ss)=parameter_c(ss,1)+parameter_c(ss,2)*x+parameter_c(ss,3)*y+parameter_c(ss,4)*z+parameter_c(ss,5)*w+parameter_c(ss,6)*rr;
           output=output+f1(ss)*rule(ss);
           ss=ss+1;
       end
       
       output_final=output;
       %****************************
       error=target(i)-output_final;
       ss=1;
       for tr=1:number_mf*number_mf*number_mf*number_mf*number_mf
          parameter_c_temp(ss,1)=-1*rule(ss);
          parameter_c_temp(ss,2)=-1*rule(ss)*x;
          parameter_c_temp(ss,3)=-1*rule(ss)*y;
          parameter_c_temp(ss,4)=-1*rule(ss)*z;
          parameter_c_temp(ss,5)=-1*rule(ss)*w;  
          parameter_c_temp(ss,6)=-1*rule(ss)*rr; 
          ss=ss+1;
       end
        parameter_c_reshape_asli=reshape(parameter_c,numel(parameter_c),1)';
        parameter_c_reshape= reshape(parameter_c_temp,numel(parameter_c_temp),1)';
        %        temp_error=temp_error+eta*error*rule;
        parameter_c_jacobe(i,:)=parameter_c_reshape;
        error_jacobe(i)=error;
        %        w_end_jacobe
    end
    miu = 0.01 * ( error_jacobe' * error_jacobe);
    parameter_c_reshape_asli   = ( parameter_c_reshape_asli' - inv( parameter_c_jacobe' * parameter_c_jacobe + miu * I) * parameter_c_jacobe' * error_jacobe)';
    parameter_c=reshape(parameter_c_reshape_asli,number_mf*number_mf*number_mf*number_mf*number_mf,6);%        for t=1:number_mf
     
   
   
 for i=1:number_train     
       x=data(i,1);
       y=data(i,2);
       z=data(i,3);
       w=data(i,4);
       rr=data(i,5);
       x_input=unifrnd(x,x,[1,number_mf]);
       y_input=unifrnd(y,y,[1,number_mf]);
       z_input=unifrnd(z,z,[1,number_mf]);
       w_input=unifrnd(w,w,[1,number_mf]);
       r_input=unifrnd(rr,rr,[1,number_mf]);
       sub_rule_x=exp(-0.5*((x_input-Mean_mf_x)./(Sigma_mf_x)).^2); 
       sub_rule_y=exp(-0.5*((y_input-Mean_mf_y)./(Sigma_mf_y)).^2);
       sub_rule_z=exp(-0.5*((z_input-Mean_mf_z)./(Sigma_mf_z)).^2);
       sub_rule_w=exp(-0.5*((w_input-Mean_mf_w)./(Sigma_mf_w)).^2);
       sub_rule_r=exp(-0.5*((r_input-Mean_mf_r)./(Sigma_mf_r)).^2);
       c=1;
       for p=1:number_mf
         for r=1:number_mf 
             for i3=1:number_mf 
                 for i4=1:number_mf 
                     for i5=1:number_mf 
                        rule(c)= sub_rule_x(p)*sub_rule_y(r)*sub_rule_z(i3)*sub_rule_w(i4)*sub_rule_r(i5);
                        c=c+1;
                     end
                 end
             end
         end
       end
       sum_rule=sum(rule);
       rule=rule/sum_rule;
       output1=0;
       ss=1;
       for tr=1:number_mf*number_mf*number_mf*number_mf*number_mf
           f1(ss)=parameter_c(ss,1)+parameter_c(ss,2)*x+parameter_c(ss,3)*y+parameter_c(ss,4)*z+parameter_c(ss,5)*w+parameter_c(ss,6)*rr;
           output1=output1+f1(ss)*rule(ss);
           ss=ss+1;
       end
       output(i)=output1;
 end


figure(1);
subplot(1,2,1),plot(output(1:number_train),'-b');
hold on;
subplot(1,2,1),plot(target(1:number_train),'-r');
hold off;
mse1=mse(output(1:number_train)-target(1:number_train)');
title(sprintf('Sugeno Train - Mackey Glass\nEpoch = %d    MSE = %.10f ',iter,mse1),'fontsize',10,'fontweight','b');
legend('Sugeno Train','Target');

 for i=1:number_test     
       x=data(number_train+i,1);
       y=data(number_train+i,2);
       z=data(number_train+i,3);
       w=data(number_train+i,4);
       rr=data(number_train+i,5);
       x_input=unifrnd(x,x,[1,number_mf]);
       y_input=unifrnd(y,y,[1,number_mf]);
       z_input=unifrnd(z,z,[1,number_mf]);
       w_input=unifrnd(w,w,[1,number_mf]);
       r_input=unifrnd(rr,rr,[1,number_mf]);
       sub_rule_x=exp(-0.5*((x_input-Mean_mf_x)./(Sigma_mf_x)).^2); 
       sub_rule_y=exp(-0.5*((y_input-Mean_mf_y)./(Sigma_mf_y)).^2);
       sub_rule_z=exp(-0.5*((z_input-Mean_mf_z)./(Sigma_mf_z)).^2);
       sub_rule_w=exp(-0.5*((w_input-Mean_mf_w)./(Sigma_mf_w)).^2);
       sub_rule_r=exp(-0.5*((r_input-Mean_mf_r)./(Sigma_mf_r)).^2);
       c=1;
       for p=1:number_mf
         for r=1:number_mf 
             for i3=1:number_mf 
                 for i4=1:number_mf 
                     for i5=1:number_mf 
                        rule(c)= sub_rule_x(p)*sub_rule_y(r)*sub_rule_z(i3)*sub_rule_w(i4)*sub_rule_r(i5);
                        c=c+1;
                     end
                 end
             end
         end
       end
       sum_rule=sum(rule);
       rule=rule/sum_rule;
       output1=0;
       ss=1;
       for tr=1:number_mf*number_mf*number_mf*number_mf*number_mf
           f1(ss)=parameter_c(ss,1)+parameter_c(ss,2)*x+parameter_c(ss,3)*y+parameter_c(ss,4)*z+parameter_c(ss,5)*w+parameter_c(ss,6)*rr;
           output1=output1+f1(ss)*rule(ss);
           ss=ss+1;
       end
       output(number_train+i)=output1;
 end

subplot(1,2,2),plot(output(number_train+1:n),'-b');
hold on;
subplot(1,2,2),plot(target(number_train+1:n),'-r');
hold off;
mse1=mse(output(number_train+1:n)-target(number_train+1:n)');
title(sprintf('Sugeno Test - Mackey Glass\n    MSE = %.10f ',mse1),'fontsize',10,'fontweight','b');
legend('Sugeno Test','Target');
end


