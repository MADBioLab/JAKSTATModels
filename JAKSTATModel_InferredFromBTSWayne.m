function JAKSTATModel_InferredFromBTSWayne(parameterValues)
    parameterNames = {'P1','K2','K3','E3','mEC50','K4','K6','K7',...
                      'm8','m9','K10','m11','K12','K13','m14','m15',...
                      'P16','K17','m18','K19','K22',...
                      'K23','K24','K25','K26','IL4'};

    %% load experimental data
    FittingData=packageData('data_CtrlAndIL4Replicates.csv',14,1);
    
    %% state variables 
    variableNames={'JAK','JAKp', ...
                   'STATCyt','STATpCyt','STATpDimerCyt', 'STATdimerNuc',...
                   'RE', 'luciferase','activated RE', 'luminescence',...
                   'STATNuc', 'SOCS'};
    Nvariables=length(variableNames);

    %% initial conditions
    Y0 = zeros(1,Nvariables);  
    Y0(7) = 100; % Response element

    %% load parameter values into struct
    for i=1:length(parameterNames)
        P.(parameterNames{i}) = parameterValues(i);
    end

    %% evaluate params that are scaled relative to each other
    %SOCS mediated JAK turnover, expressed as a multiple of JAK turnover
    P.K18=P.m18 * P.K2;
    %Response element deactivation, expressed as a multiplier of response element
    %activation
    P.K11=P.m11 * P.K10;
    %nuclear STAT export, expressed as a multiplier of nuclear STAT degradation
    P.K14=P.m14 * P.K13;
    %nuclear translocation of cytosolic STAT, expressed as a 
    %multiplier of nuclear STAT degradation
    P.K15=P.m15 * P.K13;
    %STAT dimer dissociation, expressed as a multiplier of STATp dimerization
    P.K8=P.m8 * P.K7;
    %nuclear translocation of STAT dimer, expressed as a multiplier of STATp dimerization
    P.K9=P.m9*P.K7;
    %Half saturation constant for activation of JAK by IL4, expressed as 
    %multiplier of IL4 concentration    
    P.EC50=P.mEC50*P.IL4;
    
    %conversion factor between dimensionless simulation time and expt time
    %(in timepoints per day)
    P.Tscale=24/7.4; 
    save('P','P');
    
    %% solve for steady state in absence of IL4
    tmax = 500; 
    P.IL4Level=0; 
    P.ATPFlag=1;
    P.dilutionPeriod=P.Tscale; %media replacement, once per day 
    options = odeset('NonNegative',1:Nvariables);
    [~,y] = ode15s(@(t,Y,P)JAKSTAT_ODEs(t,Y,P),[0 tmax],Y0',options,P);
    Y0SS = y(end,:);
    
    %% run simulations starting from ctrl steady state
    Selected=10; %index of the state variable of interest; 10 is luminescence
    %% timespan over which to evaluate outcomes
    dataTspan=linspace(0,14,1000)*P.Tscale; %14 days, converted to dimensionless simulation time
    
    %% ctrl
    [tctrl,yctrl] = ode15s(@(t,Y,P)JAKSTAT_ODEs(t,Y,P),dataTspan,Y0SS,options,P);

    %% IL4 1
    P.IL4Level=P.IL4;
    P.ATPFlag=1;    
    [~,y1] = ode15s(@(t,Y,P)JAKSTAT_ODEs(t,Y,P),dataTspan,Y0SS,options,P);
    y1=y1./yctrl;

    %% IL4 20
    P.IL4Level=20*P.IL4;
    P.ATPFlag=1;
    [~,y20] = ode15s(@(t,Y,P)JAKSTAT_ODEs(t,Y,P),dataTspan,Y0SS,options,P);
    y20=y20./yctrl;

    %% IL4 40
    P.IL4Level=40*P.IL4;
    P.ATPFlag=1;
    [~,y40] = ode15s(@(t,Y,P)JAKSTAT_ODEs(t,Y,P),dataTspan,Y0SS,options,P);  
    y40=y40./yctrl;

    %% IL4 80
    P.IL4Level=80*P.IL4;
    P.ATPFlag=1;
    [~,y80] = ode15s(@(t,Y,P)JAKSTAT_ODEs(t,Y,P),dataTspan,Y0SS,options,P);
    y80=y80./yctrl;

    %% IL4 100
    P.IL4Level=100*P.IL4;
    P.ATPFlag=1;
    [~,y100] = ode15s(@(t,Y,P)JAKSTAT_ODEs(t,Y,P),dataTspan,Y0SS,options,P);
    y100=y100./yctrl;
   
    %% output variables for plotting
    solnsT=tctrl./P.Tscale; %convert dimensionless time to days
    solnsY=[y1(:,Selected),...
            y20(:,Selected),...
            y40(:,Selected),...
            y80(:,Selected),...
            y100(:,Selected)];

    makeVariablePlot(solnsT,solnsY,FittingData,2)        
    
end

function odes = JAKSTAT_ODEs(t,Y,P )
    odes = zeros(size(Y));
    JAK = Y(1);
    JAKp = Y(2); 
    STATCyt = Y(3);
    STATpCyt =Y(4);
    STATpDimerCyt = Y(5); 
    STATdimerNuc = Y(6);
    RE = Y(7);
    luciferase = Y(8);
    RE_STAT = Y(9);
    luminescence = Y(10);
    STATNuc = Y(11);
    SOCS = Y(12);
    
    %https://bionumbers.hms.harvard.edu/search.aspx?trm=membrane+volume+
    %membrane 450 um2 surface ares; thickness 10 nm (monocyte,eukaryote);
    %68.1% of volume is cytoplasm; 25.9% nucleus; 
    %diameter is ~9um -> 380 um3 total; 260 um3 cytoplasm; 99 um3 nucleus
    %4.5 um3 membrane; cyt2mem ratio 58; cyt2nuc ratio 2.6
    VmemRatio = 58;
    VnucRatio = 2.6;
    
    R1 = VmemRatio*P.P1; %JAK synthesis
    R2 = P.K2*VmemRatio*JAK; %JAK turnover
    %% IL4 signaling
    interval=floor(t/P.dilutionPeriod);
    IL4=P.IL4Level*(14/15)^(interval);          %IL4 dilution by media replacement
    E = P.E3*IL4/(P.EC50+IL4);                  %Effect of IL4 signaling on JAK activation
    R3 = P.K3*VmemRatio*(1+E)*JAK*P.ATPFlag;    %JAK phosphorylation
    
    %%
    R4 = P.K4*VmemRatio*JAKp;                   %JAKp turnover
    R5 = VmemRatio*JAKp*STATCyt;                %STAT phosphorylation
    R6 = P.K6*STATpCyt;                         %STATpCyt turnover
    R7 = P.K7*STATpCyt^2;                       %STATp dimerization
    R8 = P.K8*STATpDimerCyt;                    %Stat dimer dissociation
    R9 = P.K9*STATpDimerCyt;                    %Stat dimer nuclear translocation
    R10 = P.K10*VnucRatio.^2*RE*STATdimerNuc;   %Response element activation
    R11 = P.K11*VnucRatio*RE_STAT;              %Response element deactivation
    R12 = P.K12*VnucRatio*STATdimerNuc;         %nuclear dimer dissocation and dephos'n
    R13 = P.K13*VnucRatio*STATNuc;              %nuclear STAT degradation
    R14 = P.K14*VnucRatio*STATNuc;              %nuclear STAT export
    R15 = P.K15*STATCyt;                        %STAT nuclear translocation
    R16 = P.P16;                                %STAT synthesis
    R17 = P.K17*STATCyt;                        %STATcyt degradation
    R18 = P.K18*VmemRatio*SOCS*JAK;             %SOCS mediated JAK degradation
    R19 = P.K19*VnucRatio*RE_STAT;              %SOCS synthesis
    R20 = VnucRatio*RE_STAT;                    %luciferase synthesis
    R21 = luciferase*P.ATPFlag;                 %luminescence production
    R22 = P.K22*luciferase;                     %luciferase degradation
    R23 = P.K23*luminescence;                   %luminescence decay
    R24 = P.K24*SOCS;                           %SOCS degradation
    R25 = P.K25*STATpCyt;                       %pSTAT dephosphorylation
    R26 = P.K26*VmemRatio*JAKp;                 %JAKp internalization
    
    % JAK = Y(1);
    odes(1)=R1+R4+R5-(R2+R3+R18);
    % JAKp = Y(2); 
    odes(2)=R3-(R4+R5+R26);
    % STATCyt = Y(3);
    odes(3)=R16+R14+R25-(R5+R15+R17);
    % STATpCyt = Y(4);
    odes(4)=R5+2*R8-(R6+2*R7+R25);
    % STATpDimerCyt = Y(5);
    odes(5)=R7-(R8+R9);
    % STATdimerNuc = Y(6);
    odes(6)=R9+R11-(R10+R12);
    % RE = Y(7);
    odes(7)=R11-R10;
    % luciferase = Y(8);
    odes(8)=R20-R22;
    % RE_STAT = Y(9);
    odes(9)=R10-R11;
    % luminescence = Y(10);
    odes(10)=R21-R23;
    % STATNuc = Y(11);
    odes(11)=2*R12+R15-(R13+R14);
    % SOCS = Y(12);
    odes(12)=R19-R24;

end

function makeVariablePlot(T,Y,FittingData,n)
    %% concentration time plots for state variables
    h=figure(n); 
    colorSelection={"#17202A","#566573","#99A3A4","#17202A","#AEB6BF","#566573","#99A3A4"};
    LineSelection={'-','--',':','-.','-','--',':','-.'};
    MarkerSelection={'o','^','x','*','+'};
    IL4concentrationLabels={'20 ng/ml','40 ng/ml','80 ng/ml','100 ng/ml'};
    for i=1:4
        p=plot(T,Y(:,i+1),'linewidth',2); p.Color=colorSelection{i};
        p.LineStyle=LineSelection(i);
        hold on;
    end
    for i=1:4
        p=errorbar(FittingData.T,FittingData.data(:,i),...
                   FittingData.std(:,i),MarkerSelection{i},'linewidth',1);
        p.Color=colorSelection{i};
    end
    hold off;
    legend([IL4concentrationLabels,IL4concentrationLabels]);
    xlabel('time [days]'); ylabel('RLU_{IL4}/RLU_{control} [-]')
    set(gca,'fontsize',12,'FontWeight','bold')
    axis square

    saveas(h,'ModelPredictions.fig')
    saveas(h,'ModelPredictions.png')
end

function Package=packageData(filename,tmax,n)
    %filename - path to file with exptl data
    %tmax - upper limit of time in days to read from file
    rawData=importdata(filename);

    %select t<tmax days
    Package.T=rawData.data(:,1); %time
    selected=(Package.T<=tmax);
    Package.T=Package.T(selected);
    NtimePoints = length(Package.T);
    data=rawData.data(selected,2:end); %RLU data

    Nreplicates=3; %number of experimental replicates

    %columns should be time followed by a column of data for each
    %replicate. Repeated for each concentration, giving
    %replicates*concentrations+1 columns
    Package.Nconcentrations=size(data,2)/Nreplicates-1; 
    %control replicates
    ctrlData=mean(data(:,1:Nreplicates),2);
    %IL4-treated replicates
    data=data(:,Nreplicates+1:end)./ctrlData;

    %evaluate mean and standard deviation across replicates for each IL4
    %concentration
    avgData=zeros(NtimePoints,Package.Nconcentrations);
    sigma=zeros(NtimePoints,Package.Nconcentrations);
    for i=1:Package.Nconcentrations
        a=(i-1)*Nreplicates+1;
        b=i*Nreplicates;
        avgData(:,i)=mean(data(:,a:b),2);
        sigma(:,i)=std(data(:,a:b),[],2);
    end
    Package.initial=avgData(1,:)-1;
    Package.data=avgData; 
    Package.std=sigma;

    %plot timecourses for means of experimental data
    figure(n);errorbar(Package.T,Package.data,Package.std,'linewidth',2);
    xlabel('time [days]'); ylabel('RLU [-]')
    set(gca,'fontsize',14)
end