function [Sum] = GaussFit3_2Para(A, B, IniC, IniW)

% curve fitting toolbox required

%CREATEFIT(A,B)
%  �ߎ����쐬���܂��B
%
%  '�V�K�ߎ� 1' �ɑ΂���f�[�^���ߎ�:
%      X ����:  A
%      Y �o��:  B
%  �o��:
%      fitresult: �ߎ���\�� fit �I�u�W�F�N�g�B
%      gof: �K�����������\���́B
%
%  �Q�l FIT, CFIT, SFIT.

%  MATLAB ����̎���������: 24-Apr-2021 04:27:10


%% �ߎ�: '�V�K�ߎ� 1'�B
try
warning('off', 'curvefit:prepareFittingData:sizeMismatch'); % ignore unknown alert 

[xData, yData] = prepareCurveData( A, B );

% �ߎ��^�C�v�ƃI�v�V������ݒ肵�܂��B
ft = fittype( 'C * exp(-(x^2)/(2*W^2)) + 1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [IniC, IniW];

% ���f�����f�[�^�ɋߎ����܂��B
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


% % �f�[�^�̋ߎ����v���b�g���܂��B
% figure( 'Name', '�V�K�ߎ� 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'B vs. A', '�V�K�ߎ� 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % ���x�� Axes
% xlabel( 'A', 'Interpreter', 'none' );
% ylabel( 'B', 'Interpreter', 'none' );
% grid on


end

