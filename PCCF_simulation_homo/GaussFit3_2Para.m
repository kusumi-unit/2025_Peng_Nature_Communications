function [Sum] = GaussFit3_2Para(A, B, IniC, IniW)

% curve fitting toolbox required

%CREATEFIT(A,B)
%  近似を作成します。
%
%  '新規近似 1' に対するデータを近似:
%      X 入力:  A
%      Y 出力:  B
%  出力:
%      fitresult: 近似を表す fit オブジェクト。
%      gof: 適合性情報をもつ構造体。
%
%  参考 FIT, CFIT, SFIT.

%  MATLAB からの自動生成日: 24-Apr-2021 04:27:10


%% 近似: '新規近似 1'。
try
warning('off', 'curvefit:prepareFittingData:sizeMismatch'); % ignore unknown alert 

[xData, yData] = prepareCurveData( A, B );

% 近似タイプとオプションを設定します。
ft = fittype( 'C * exp(-(x^2)/(2*W^2)) + 1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [IniC, IniW];

% モデルをデータに近似します。
[fitresult, gof] = fit( xData, yData, ft, opts );

ResultC = fitresult.C;
ResultW = fitresult.W;
ci = confint(fitresult,0.95);
ConfC = (ci(2,1) - ci(1,1))/2;
ConfW = (ci(2,2) - ci(1,2))/2;

catch
    % if the fitting failed, return Nan
    ResultC = nan;
    ResultW = nan;
    ConfC = nan;
    ConfW = nan;
end

Sum(1,1) = ResultC;
Sum(2,1) = ResultW;
Sum(1,2) = ConfC;
Sum(2,2) = ConfW;


% % データの近似をプロットします。
% figure( 'Name', '新規近似 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'B vs. A', '新規近似 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % ラベル Axes
% xlabel( 'A', 'Interpreter', 'none' );
% ylabel( 'B', 'Interpreter', 'none' );
% grid on


end

