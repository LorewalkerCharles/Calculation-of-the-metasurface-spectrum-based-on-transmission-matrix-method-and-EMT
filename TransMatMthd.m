function [r,t] = TransMatMthd(Freq,n2,n3,d2,d3)
%% ����4�㱡Ĥ�ķ����͸�䣬���ڴ�����󷨵ķֲ�ṹ�������Ч����ģ��
% Author: cuizijian_harbin@163.com
% Data: 20240522����
% ����ʹ�ã��������������: [10.1364/PRJ.450017]
% ��4��ϵͳ����������Ƶ����ص�������ʽ���뺯������������ķ����͸��
% ���⻹�����2��3��ĺ��
% ���䷽��: Incident->n1->n2->n3->n4 ����ͨ��Ϊ���ɿռ����䣬����n1��n4һ��Ϊ������ֵΪ1
% Since it is usually incident in free space, n1 and n4 are generally air and have values of 1

% ����                         % Input       
% Freq: Ƶ������               % Frequency
% n2: ������ṹ���Ч������   % The effective refractive index of the structured-layer
% n3: �ĵײ�������             % The refractive index of the substrate
% d2: �ṹ����               % The thickness of the structured-layer
% d3: �ĵײ���               % The thickness of the substrate

% ���            % Return
% r: ����         % Reflection
% t: ͸��         % Transmission

%% ������󷨼���
% Light speed
c_light = 3e8; 
k = 2*pi*Freq/c_light;
n1 = 1; n4 = 1;
% Fresnel Equations
r12 = (n1-n2)./(n1+n2); t12 = (2*n1)./(n1+n2); 
r23 = (n2-n3)./(n2+n3); t23 = (2*n2)./(n2+n3); delta2 = n2.*k.*d2;
r34 = (n3-n4)./(n3+n4); t34 = (2*n3)./(n3+n4); delta3 = n3.*k.*d3;
r = (exp(-1i*delta3).*(r12.*exp(-1i*delta2)+r23.*exp(1i*delta2))+...
    r34.*exp(1i*delta3).*(r12.*r23.*exp(-1i*delta2)+exp(1i*delta2)))./...
    (exp(-1i*delta3).*(exp(-1i*delta2)+r12.*r23.*exp(1i*delta2))+...
    r34.*exp(1i*delta3).*(r23.*exp(-1i*delta2)+r12.*exp(1i*delta2)));
t = (t12.*t23.*t34)./...
    (exp(-1i*delta3).*(exp(-1i*delta2)+r12.*r23.*exp(1i*delta2))+...
    r34.*exp(1i*delta3).*(r23.*exp(-1i*delta2)+r12.*exp(1i*delta2)));
end

