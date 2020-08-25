function r=ShowStatus(t,y,s)
%global TT YY 
global g
r=0;
if isempty(s)
%      TT(end+1)=t;
%      YY(1:length(y),end+1)=y;
%    
    p = round( (t(end)/g.t_end)^0.1 * 100 );
    
     if g.progress~=p
         g.progress=p;
         disp( [ 8 8 8 8 8 8 8 8  num2str(g.progress) ' %  ]'] );
     end;
    %[~, idx] = sort(y,'descend');
%      figure(6)
%      cla
%      loglog(TT,YY); %(idx(1:100),:));
%      hold on
%      loglog(TT,sum(YY),'k','LineWidth',2);%sum(YY(idx(1:100),:)),'k','LineWidth',2);
%      ylim([1e-30 max( sum(YY) )])
%      drawnow;
%      refresh;
end
