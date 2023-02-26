function myeps(fname,s)
% save eps figure

    disp('MYEPS: Saving figure to EPS has been disabled!')
    return

    set(gcf,'PaperPositionMode','auto');
    if nargin < 2,
        s = 2; % scaling (2 on windows, 1 on Linux)
    end
    f = 1; % height/width
    set(gcf,'PaperSize',s*[9,6.2*f]);
    set(gcf,'PaperPosition',[-.5 -.1 s*10 s*6.4*f]);
    print(gcf,'-depsc2', fname);
end

