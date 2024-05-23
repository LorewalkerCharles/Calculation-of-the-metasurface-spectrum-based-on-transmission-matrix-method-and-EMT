function [r,t] = TransMatMthd(Freq,n2,n3,d2,d3)
%% 计算4层薄膜的反射和透射，基于传输矩阵法的分层结构超表面等效介质模型
% Author: cuizijian_harbin@163.com
% Data: 20240522整理
% 如需使用，请引用相关文献: [10.1364/PRJ.450017]
% 将4层系统的折射率以频率相关的向量形式传入函数，计算整体的反射和透射
% 此外还需给出2、3层的厚度
% 入射方向: Incident->n1->n2->n3->n4 由于通常为自由空间入射，所以n1和n4一般为空气，值为1
% Since it is usually incident in free space, n1 and n4 are generally air and have values of 1

% 输入                         % Input       
% Freq: 频率序列               % Frequency
% n2: 超表面结构层等效折射率   % The effective refractive index of the structured-layer
% n3: 衬底层折射率             % The refractive index of the substrate
% d2: 结构层厚度               % The thickness of the structured-layer
% d3: 衬底层厚度               % The thickness of the substrate

% 输出            % Return
% r: 反射         % Reflection
% t: 透射         % Transmission

%% 传输矩阵法计算
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

