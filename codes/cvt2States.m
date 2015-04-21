function [out, uni, Nuni] = cvt2States(Data)

out = zeros(size(Data,1), size(Data,2));
% out = (Data < 0) * 0;
% out_tmp1 = Data > 0;
% out_tmp2 = Data < 1;
% out_tmp3 = out_tmp1 .* out_tmp2;
% out = out + out_tmp3 .* 1;
% out_tmp3 = Data > 1;
% out = out + out_tmp3 .* 2;

for i = 1: size(Data,2)
    uni{i} = unique(Data(:,i));
       
    if ~(uni{i}(1) < 0)
        disp(['needs check! - val(state0) > 0...?' num2str(i)]);
    end
    
    if size(uni{i},1) == 2
%         out(:,i) = out(:,i) + (Data(:,i) == uni{i}(2)) .* 1;
        if (uni{i}(1) < 0) && (uni{i}(2) - uni{i}(1) < 1.1)
            out(:,i) = out(:,i) + (Data(:,i) == uni{i}(2)) .* 1;
        elseif (uni{i}(1) < 0) && (uni{i}(2) - uni{i}(1) < 2.1)
            out(:,i) = out(:,i) + (Data(:,i) == uni{i}(2)) .* 2;
        else
            disp(['needs check! - 2 states...?' num2str(i)]);
        end
        
    elseif size(uni{i},1) == 3
        out(:,i) = out(:,i) + (Data(:,i) == uni{i}(2)) .* 1;
        out(:,i) = out(:,i) + (Data(:,i) == uni{i}(3)) .* 2;
    elseif size(uni{i},1) > 3
        disp('error! - # of unique val > 3');
    end
    
    Nuni(i) = size(uni{i},1);
end

