function [aa,SVset]=ExportMSAResults (fn,titl,alpha,roi_flag,calc_bs,varargin)

dat = load (fn);
if isempty(varargin)
    ROIs = dat.RegionNames;
else
    ROIs = varargin{1};
end

LW=1;%alpha = 0.05;
ParNm = {'rank','SV','sd','FDRp'};

AsteriskFact = 1.1;
nxt = 1;    
%maxScale = [30,24,66];
Nreg = size(dat.SV,2);
BigArr = zeros(Nreg,4);
%VarNm = cell(1,30);
%f=figure;      

xvec=1:Nreg;

if calc_bs == 1    
    CBS = dat.Bset{end}.CalibBootstraps;
    f = figure;
    Nbins = dat.nBS/100; bin_lim = [-50,50];
    sv = dat.Bset{end}.CIcalibmix(:,2);
    svci1 = dat.Bset{end}.CIcalibmix(:,3);
    svci0 = dat.Bset{end}.CIcalibmix(:,1);
    % histograms
    nx = 1;
    for roi=1:size(CBS,1)
        %figure;
        ah(roi)=subplot(5,6,nx); nx = nx+1; %change according to number of regions
        %ah(roi) = axes(); 
        histogram(CBS(roi,:),Nbins,'BinLimits',bin_lim); 
        frq=histcounts(CBS(roi,:),Nbins,'BinLimits',bin_lim); 
        mx = max(frq)*1.1;
        hold on       
        line ([sv(roi) sv(roi)],[0 mx],'Color','r','LineStyle','-','LineWidth',1);
        line ([svci0(roi) svci0(roi)],[0 mx],'Color','m','LineStyle','--','LineWidth',0.5);
        line ([svci1(roi) svci1(roi)],[0 mx],'Color','m','LineStyle','--','LineWidth',0.5);
        ah(roi).YLim = [0 mx];
        title(ROIs{roi});
    end
end
figure;
SVset.BS = CBS;
SV = dat.Bset{end}.CIcalibmix(:,2);
SVci = dat.Bset{end}.CIcalibmix(:,3);
SVsd = SVci - SV;
hiLim = SVci;
FDRpval = mafdr(dat.Bset{end}.pvalestmix,'BHFDR',true);
[~,SVrank] = sort(SV,'descend');
SVset.SV = SV;
SVset.pval = FDRpval;
SVset.Z = dat.Bset{end}.Zscoreestmix;
hold on
ast = ones(1,Nreg).*NaN;  
bigsv = ones(1,Nreg).*NaN;
ast(FDRpval<alpha)=(hiLim(FDRpval<alpha))*AsteriskFact;
bigsv(dat.Bset{end}.Zscoreestmix>0) = 1;
ast = ast.*bigsv;
h=line([0 Nreg+1],[0 0],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',LW);
aa = h.Parent;
bar (SV,'FaceColor',[68/256, 114/256, 196/256],'LineWidth',LW);
errorbar (xvec,SV,SVsd,'k','LineWidth',LW,'LineStyle','none');               
a = h.Parent;        
plot (xvec,ast,'LineStyle','none','Marker','*','MarkerSize',12,'MarkerEdgeColor','k');  
tith = title(titl);
tith.FontSize = 14;
tith.FontWeight = 'normal';
%title ([fn ' ' calib_str]);
a.LineWidth = LW * 1.25;
a.XLim= [0 Nreg+1];  
%a.YLim = [-4 14];        
a.XTick = xvec;
a.TickDir = 'out';
EmptyLabels = cell(1,Nreg);
for j=1:Nreg,EmptyLabels{j}='';end
if roi_flag
    a.XTickLabel = ROIs;
    a.XTickLabelRotation = 90;
else
    a.XTickLabel = EmptyLabels;
end        
%Add Panel label:
%th = text(-2,16,panels{nx-1});
th.FontSize = 18;
%th.FontWeight = 'bold';
fprintf('\n\n');        
%         disp (['---' side{sd}  '-' behav{bh} '---']);
denom = sqrt((Nreg-1)/Nreg^2);
Lsv = std(SV)/denom;
%Lci1 = std(CI(1,:))/denom;
%Lci2 = std(CI(2,:))/denom;
disp (titl);
disp (['L [Aharonov] = ',num2str(Lsv,3)]);
for j=1:Nreg
    if FDRpval(j)<alpha && dat.Bset{end}.Zscoreestmix(j)>0 > 0
    if SV(j) > 0
        Prefix = '***';
    else
        Prefix = '';
    end
        disp ([Prefix ROIs{j} ': ' num2str(SV(j),2)]);
        %disp (Prefix);
    end
end
%VarTitle = [side{sd} '_' behav{bh}];        
%         for k=1:length(ParNm)
%             VarNm{nxt+k-1} = [VarTitle '_' ParNm{k}];
%         end
BigArr(:,nxt) = SVrank';
BigArr(:,nxt+1) = SV';
BigArr(:,nxt+2) = SVsd;
BigArr(:,nxt+3) = FDRpval';
%nxt = nxt +4;        
 

%tabx = array2table(BigArr,'VariableNames',VarNm,'RowNames',ROIs);
% tabx = array2table(BigArr,'VariableNames',ParNm,'RowNames',ROIs);
% writetable(tabx,'MSA_Summary_Revised7-12-12xlsx','WriteRowNames',true);

end

