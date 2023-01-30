% https://ch.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval

function CI = confidenceInterval(x)

SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals