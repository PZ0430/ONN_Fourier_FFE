% function aligned_output = align_signals_by_max_cross_correlation(input_signal, output_signal)
%     % 最大互相关对齐函数
%     % 输入：
%     %   input_signal - 原始输入信号
%     %   output_signal - 均衡后的输出信号
%     % 输出：
%     %   aligned_output - 与输入信号对齐后的输出信号
% 
%     % 计算互相关
%     % output_signal = yd;
%     % input_signal = dataTXX;
%     [cross_corr, lag] = xcorr(output_signal, input_signal);
% 
%     % 找到互相关的最大值位置
%     [~, max_index] = max(abs(cross_corr));
%     optimal_lag = lag(max_index);
%     % disp(optimal_lag);
%     % 根据最优延迟对齐信号
%     if optimal_lag > 0
%         % 如果输出信号滞后于输入信号，补偿延迟
%         aligned_output = output_signal(optimal_lag + 1 : end);
%     elseif optimal_lag < 0
%         % 如果输出信号超前于输入信号，截掉超前部分
% 
%         % aligned_output = [zeros(abs(optimal_lag), 1) output_signal];
%         aligned_output = [zeros(1, abs(optimal_lag)) output_signal];
%         aligned_output = aligned_output(1 : length(input_signal));
%     else
%         % 如果没有延迟，则直接对齐
%         aligned_output = output_signal;
%     end
% end

function [aligned_input, aligned_output, max_CC] = align_signals_by_max_cross_correlation(input_signal, output_signal)
% 对齐信号使最大互相关
% 输入：
%   input_signal  - 原始输入信号
%   output_signal - 输出信号
% 输出：
%   aligned_input  - 对齐后的输入信号
%   aligned_output - 对齐后的输出信号

    % 确保信号是列向量
    input_signal = input_signal(:);
    output_signal = output_signal(:);

    % 计算互相关
    [cross_corr, lag] = xcorr(output_signal, input_signal);

    % 找到最大互相关对应的延迟
    [max_CC, max_index] = max(abs(cross_corr));
    optimal_lag = lag(max_index);

    % 对齐两个信号
    if optimal_lag > 0
        % 输出信号滞后，裁剪输入信号前面部分
        aligned_input = input_signal(1:end - optimal_lag);
        aligned_output = output_signal(optimal_lag + 1:end);
    elseif optimal_lag < 0
        % 输出信号超前，裁剪输出信号前面部分
        aligned_input = input_signal(abs(optimal_lag) + 1:end);
        aligned_output = output_signal(1:end + optimal_lag);
    else
        % 无延迟，长度对齐
        min_len = min(length(input_signal), length(output_signal));
        aligned_input = input_signal(1:min_len);
        aligned_output = output_signal(1:min_len);
    end

    % 最终截取为相同长度
    min_len = min(length(aligned_input), length(aligned_output));
    aligned_input = aligned_input(1:min_len);
    aligned_output = aligned_output(1:min_len);
end