function [] = startFilterMAFFile(mafFilename, outputMAFFilename, outputDir, isCheckForMatFile, isGeneratingPlots, globalPoxoG, artifactThresholdRate, lod0Thresh, oxoQP1,oxoQP2,artifactThresholdCount)
% [] = startFilterMAFFile(mafFilename, outputMAFFilename, outputDir, isCheckForMatFile, isGeneratingPlots, globalPoxoG, artifactThresholdRate, lod0Thresh, oxoQP1,oxoQP2,artifactThresholdCount)
%
%   Entry point for running D-ToxoG.  
%       The third incarnation of the C>A/G>T artifact filter.
%
%   This script will perform the filtering and produce the following:
%   	MAF Files:  Filtered and Unfiltered
%           These will be named as <outputMAFFilename> and
%           <outputMAFFilename>.reject.maf, respectively.
%       Count Files:  Only contain the values of the total number of pass
%           and rejected mutations for the input maf file.
%           <outputMAFFilename>.pass_count.txt
%           <outputMAFFilename>.reject_count.txt
%           
%       Figures
%       Tables
%
%       The latter two can be used for further reporting downstream.
%
%   Note:  This script will preserve the order of mutations with the
%       following exception:
%
%       If multiple cases are included in the input maf file and the
%       mutations of cases are interleaved.  In this case, this script will
%       produce maf files with the mutations appearing clustered by case.
%
%       Simple example: 
%           Case input (3 mutations):
%        1 ...   i1-Tumor  i1-Normal ...
%        2 ...   i2-Tumor  i2-Normal ...
%        3 ...   i1-Tumor  i1-Normal ...
%
%           Case output (3 mutations):
%        1 ...   i1-Tumor  i1-Normal ...
%        3 ...   i1-Tumor  i1-Normal ...
%        2 ...   i2-Tumor  i2-Normal ...
% 
% mafFilename -- input maf filename to be filtered.  This file will not be
%   modified.
%   Required Headers:
%     Chromosome -- Contig number without a prefix.  E.g. "3" or "X" (without quotes)
%     Start_position -- Position of the SNV
%     End_position -- Should be the same as Start_position for SNV
%     Reference_Allele -- The allele found in the reference genome at this position.  Single character representing the base, e.g. "G".
%     Tumor_Seq_Allele2 -- alternate allele.  Single character representing the base, e.g. "G".
%     Tumor_Sample_Barcode -- name of the tumor (or "case") sample.  This is used to generate file names and plots.
%     Matched_Norm_Sample_Barcode -- name of the normal (or "control") sample.  This is used to generate file names and plots.
%     ref_context -- Small window into the reference at the SNV.  The center position should be the same as Reference_Allele.  The total string should be of odd length and have a minimum length of 3. 
%       For example: Reference Allele is G, Chromosome is 1, Start_position and End_position are 120906037:  ref_context is CTTTTTTCGCGCAAAAATGCC  (string size is 21, in this case)
%     i_t_ALT_F1R2 -- the number of reads with pair orientation of F1R2 and with the alternate allele (Tumor_Seq_Allele2).
%     i_t_ALT_F2R1 -- the number of reads with pair orientation of F2R1 and with the alternate allele (Tumor_Seq_Allele2).
%     i_t_REF_F1R2 -- the number of reads with pair orientation of F1R2 and with the reference allele (Reference_Allele).
%     i_t_REF_F2R1 -- the number of reads with pair orientation of F2R1 and with the reference allele (Reference_Allele).
%     i_t_Foxog -- Foxog, as described in the methods.  Depending on the nature of the reference and alternate alleles, either i_t_ALT_F1R2/(i_t_ALT_F1R2 + i_t_ALT_F2R1)  or i_t_ALT_F2R1/(i_t_ALT_F1R2 + i_t_ALT_F2R1).
%         C>anything:  numerator is i_t_ALT_F2R1
%         A>anything:  numerator is i_t_ALT_F2R1
%         G>anything:  numerator is i_t_ALT_F1R2
%         T>anything:  numerator is i_t_ALT_F1R2
%     Variant_Type -- "SNP" (without quotes)
%
% outputMAFFilename -- filtered maf filename to be generated.  Without path
% information.
%
% outputDir (default: './') -- outputDir for filter results that are 
%   used for the report.
%   Note:  outputMAFFilename, figures, and tables will be located 
%       in outputDir.  The mat files will not be located in the outputDir.
% 
% globalPoxoG (default: .96) -- poxoG value to use in all estimations of
%   the given mafFile.  
%
% artifactThresholdRate (default .01) -- Expected proportion of mutations
%   that are artifacts.  e.g. .01 --> 1%  
%
% isCheckForMatFile (default: 0) -- looks for a mat file that contains the
%   mafTable already.  This is just to save on load time.  If value = 1, 
%   check for a mat file with name  mat/<mafFilename>.mat to load mafTable.
%   This option can be used to speed the running of this filter.
% 
% isGeneratingPlots (default: 0) -- true/false for whether this script
%   should generate the orchestra and lego png files.  These will be put in
%   the [outputDir]/figures/ .
%
% lod0Thresh (default: -1 --> No filtering) -- binomial mixture negative
%   log likelihood threshold.
%       The number being thresholded is the difference in the binomial
%       mixture model between the estimated weight of the oxoG component
%       and zero (the subtraction of the logs of each).  In other words, 
%       the log odds ratio of no oxoG contamination versus the oxoG
%       contamination (lod0).
%       If this difference is small, then this implies that there is not 
%       much OxoG contamination.  A large number implies lots of OxoG 
%       effect.  The difference can never be negative, since the oxoG 
%       estimate is the minimum.
%       Any cases where the difference is below this threshold are not
%           filtered.
%       The functionality around this parameter is not well tested.  
%       Use with caution.
%
%
% OxoQP1 (default 1e8  - filter independent of OxoQ) 
% OxoQP2 (default 0   - filter independent of OxoQ) 
%   Two fitting paramters to estimate the number of OxoG mutations from
%   oxoQ value. OxoQP1=1e8 and OxoQP2=0.5 result in a filter that remains 
%   like the old filter below oxoQ=35, becomes less
%   agressive above oxoQ~35 and does no filtering above oxoQ~45. 
%
%   OxoQP2=0 forces a sharp cutoff for filtering only OxoG<oxoQP1 samples.   
%
% artifactThresholdCount (default 100) -- Max count of allowed artifact mutations. 
%   used to lower the artifactThresholdRate in cases with lots of mutations and at least 100 estimated oxoG artifacts. 
%
% See also load_table

[~, filename, ext] = fileparts(mafFilename);
if ~exist('mat', 'dir')
   mkdir('.', 'mat');
end
matFilename = ['mat/' filename ext '.mat'];

%% Parse arguments and determine which should use default values.

if (nargin < 11 || ~exist('artifactThresholdCount') || isempty(artifactThresholdCount)), artifactThresholdCount = 100; end
if ischar(artifactThresholdCount),artifactThresholdCount = str2num(artifactThresholdCount); end
disp(['artifactThresholdCount: ' num2str(artifactThresholdCount) ])


if (nargin < 10 || ~exist('oxoQP2') || isempty(oxoQP2)), oxoQP2 = 0; end
if ischar(oxoQP2),oxoQP2 = str2num(oxoQP2); end
disp(['oxoQ P2: ' num2str(oxoQP2) ])

if (nargin < 9 || ~exist('oxoQP1') || isempty(oxoQP1)), oxoQP1 = 1e8; end
if ischar(oxoQP1),oxoQP1 = str2num(oxoQP1); end
disp(['oxoQ P1: ' num2str(oxoQP1) ])

if nargin < 8 || ~exist('lod0Thresh') || isempty(lod0Thresh)
    lod0Thresh = -1;
end
if ischar(lod0Thresh)
    lod0Thresh = str2num(lod0Thresh);
end
disp(['lod0 Threshold: ' num2str(lod0Thresh) ])


if nargin < 7 || ~exist('artifactThresholdRate') || isempty(artifactThresholdRate)
    artifactThresholdRate = .01;
end
if ischar(artifactThresholdRate)
    artifactThresholdRate = str2num(artifactThresholdRate);
end

if nargin < 6 || ~exist('globalPoxoG') || isempty(globalPoxoG)
    globalPoxoG = .96;
end
if ischar(globalPoxoG)
    globalPoxoG = str2num(globalPoxoG);
end
   

isGeneratingPlots = 0;
%if nargin < 5 || ~exist('isGeneratingPlots') || isempty(isGeneratingPlots)
%    isGeneratingPlots = 0;
%end

if nargin < 4 || ~exist('isCheckForMatFile') || isempty(isCheckForMatFile)
    isCheckForMatFile = 0;
end

if nargin < 3 || ~exist('outputDir') || isempty(outputDir)
    outputDir = './';
end

if ~exist(mafFilename, 'file') && (isCheckForMatFile && ~exist(matFilename, 'file'))
   error(['Given maf file does not exist: '  mafFilename])
end

%% Load the maf file
if ~isCheckForMatFile
    disp(['Loading ' mafFilename ' ...'])
    [mafTable] = loadMAFTable(mafFilename);
else
    if ~exist(matFilename, 'file')
        disp(['Loading ' mafFilename ' ...'])
        [mafTable] = loadMAFTable(mafFilename);
        disp(['Saving ' matFilename ' ... '])
        save(matFilename, 'mafTable');
    end
    disp(['Loading ' matFilename ' ...'])
    load(matFilename);
end

%% Retrieve all cases (one entry for each unique case, not one entry per
%   mutation)
disp(['Retrieving cases for ' mafFilename ' ...'])
[pairs] = retrieveUniqueCaseControlBarcodes(mafTable);
disp(['Retrieved ' num2str(length(pairs)) ' cases'])

% Target rate for artifacts that escape filtering
fdrThresh = artifactThresholdRate;  
fdrThreshCount = artifactThresholdCount;  

% Measured p used in the B_A(ac,p) distibution in the binomial mixture
% model.
PoxoG = globalPoxoG;

% Alt Allele Counts to use
acs = [3:50];

% Null distribution (real mutation) binomial parameter, p
p0 = .5;

%% Create the output directories
if ~strcmp(outputDir(end), '/') 
    outputDir = [outputDir '/'];
end

if ~exist(outputDir, 'dir')
   mkdir(outputDir); 
end

outputFigureDir =[ outputDir '/figures/'];
if ~exist(outputFigureDir, 'dir')
   mkdir('.', outputFigureDir); 
end

% Write the figure location to a file.
tmpFid = fopen([outputDir '/figureLoc.txt'], 'w');
fprintf(tmpFid, '%s', GetFullPath(outputFigureDir));
fclose(tmpFid);

outputTableDir = outputDir;
if ~exist(outputTableDir, 'dir')
   mkdir('.', outputTableDir); 
end

%% Create a file for storing the lego information
legoFidBefore=fopen([outputDir '/legoCasesBefore.txt'], 'w');
legoFid=fopen([outputDir '/legoCases.txt'], 'w');
if ~isGeneratingPlots
   disp('Lego raw data files will be empty when not generating plots.') 
end

%% Initialize and write a file for the parameters as specified
paramsFilename = [outputTableDir 'params.tsv'];
pfid = fopen( paramsFilename, 'w');
pHeaders = [{'PoxoG'};{'artifactThreshold'};{'lod0Threshold'};{'artifactThresholdCount'}];
fprintf(pfid, [strjoin('\t', pHeaders) '\n'])
fprintf(pfid,'%0.3f\t%0.2f\t%0.6f\t%d\n', PoxoG, fdrThresh, lod0Thresh,fdrThreshCount);
fclose(pfid)

%% Initialize a file to list the cases in the input maf
caseFid = fopen([outputDir '/allCases.txt'], 'w');

%% Initialize a table for summary data where each case is a row in a table.
tableFilename = [outputTableDir 'caseTableData.tsv'];
tfid = fopen( tableFilename, 'w');
headers = [{'case'}; {'control'}; {'N'}; {'N_A'}; {'N_nA'};{'N_filtered'};{'N_oxo'}; {'N_oxo_CI_low'}; {'N_oxo_CI_high'};{'N_mut'}; {'RR'}; {'RR_CI_low'}; {'RR_CI_hi'}; {'PA'}; {'PA_CI_low'}; {'PA_CI_hi'};{'PoxoG'};{'PoxoG_Est_CI_low'};{'PoxoG_Est_CI_hi'};{'Proportion_N_filt'};{'lod0'};{'oxoQ'}];
fprintf(tfid, [strjoin('\t', headers) '\n'])

totalRejectCount = 0;
totalPassCount = 0;

%% Initialize a map for all of the mutations
% This map will be used to add new columns to the maf files.
mutationMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

%% Initialize an array for storage of Passing Artifact & Rejected Real counts (and hi/low CIs) in a summary view
allPAs = zeros(length(pairs),3);
allRRs = zeros(length(pairs),3);
Ns = zeros(length(pairs),1);

%% For each unique case, run the filter
for i = 1:length(pairs)

    pair = pairs(i);    
    unameIndices = retrieveMutationsForGivenCase(pair.case, pair.control, mafTable);
    unameMafTable = trimStruct(mafTable, unameIndices);
    
    caseName = [pair.case ' v. ' pair.control];
    fprintf(caseFid, [caseName '\n']);
    
    disp(['Case ' num2str(i) ': ' caseName])
    
    % Create a list of the mutations segregated by artifact mode
    %  M_A are the artifact mode mutations
    %  M_nA are the non-artifact mode mutations
    [M, M_A, M_nA] = generateMutations(mafTable, unameIndices);
    
    N1=length(M.foxog)
    
    fdrThresh1 = min([fdrThresh fdrThreshCount/N1])
    
    % fetch oxoQ value from first mutation in MAF for this sample
    oxoQ=M.i_picard_oxoQ(1);    
    % estimate expected number of oxoG mutations from oxoQ fit 
    %eQ=oxoQP1*exp(-1*oxoQP2*oxoQ);
    % cut suppression factor based on oxoQ
    fQ=1./(1+exp(oxoQP2*(oxoQ-oxoQP1)));
    % if oxoQP2 is not set, then don't correct fdr in cutLine
    if (oxoQP2==0)
        fQ=1;
    end
    % Create a cut using False Discovery Rate 
    %fdrStructs(i) = cutLineFDR(M, fdrThresh1, PoxoG, fQ);
    fdrStructs(i) = cutLineFDR(M, fdrThresh, PoxoG, fQ);  % use original fdrThresh band-aid
    
    %% Estimate num oxo artifacts
    [Noxo,Nmut,alf,alfci,NoxoCI, lod0] = oxogBinomialMixture(M_A.alt_read_count, round(M_A.foxog .* M_A.alt_read_count), acs, PoxoG, p0);

    %% If we are 95% certain that this case has less than 1% artifacts, do not filter.  Do this by setting 
    isFilterThisCase = 1;
    if lod0 < lod0Thresh
        disp('The logarithmic odds ratio is not large enough to warrant filtering this case.  Doing no filtering.')
        isFilterThisCase = 0;
    end
    
    if (oxoQP2==0)&(oxoQ>oxoQP1)        
        disp(sprintf('OxoQ exceeds threshold: %f>%f   Doing no filtering.',oxoQ>oxoP1))
        isFilterThisCase = 0;
    end
    
    if NoxoCI(2) < (fdrThresh1*length(M.alt_read_count))
       disp('We are 95% certain that this case has less than 1% artifacts.  Doing no filtering.')
       isFilterThisCase = 0;
    end
    if ~isFilterThisCase
       fdrStructs(i).cut = zeros(size(fdrStructs(i).cut));
    end
    
    %% Update pass and reject counts after finalizing fdrStructs(i).cut
    totalPassCount = totalPassCount + sum(~fdrStructs(i).cut);
    totalRejectCount = totalRejectCount + sum(fdrStructs(i).cut);

    
    %% Estimate the a cutline for FoxoG x alt count
    % acVal ./ acs is FoxoG.  The line is (FoxoG(i), acs(i))
    [acVal] = estimateCutLineFromFdrOutput(fdrStructs(i), acs, PoxoG);
    cutLineFoxoGs = (acVal' ./ acs);

    %% Plotting figures
    if isGeneratingPlots
        
        %% Plot the artifact mode mutations
        disp('Plotting artifact mode (with patch)...')
        plotArtifactMode(outputFigureDir, caseName, M_A, fdrStructs(i).cut, cutLineFoxoGs, acs)

        %% Plot the non-artifact mode mutations
        disp('Plotting non-artifact mode ...')
        plotNonArtifactMode(outputFigureDir, caseName, M_nA, acs)

        %% Generate Lego plots (before and after)
        disp('Making lego plots...')
        Pbefore.label = caseName;
        figure('visible','off') 
        legoStructsBefore = plotMutationSpectrumCategLego1(unameMafTable,'unit',[outputFigureDir caseName '_lego_before'], Pbefore);        
        close(gcf);
        saveLegoStructToFile(caseName, legoStructsBefore, legoFidBefore, i==1);
        
        Pafter.label = [caseName ' Filtered'];
        figure('visible','off') 
        legoStructsAfter = plotMutationSpectrumCategLego1(trimStruct(unameMafTable, ~fdrStructs(i).cut),'unit',[outputFigureDir caseName '_lego_after'], Pafter);
        close(gcf);
        saveLegoStructToFile(caseName, legoStructsAfter, legoFid, i==1);
    else
        disp('Figures are not being generated.')
    end
    
    %% Estimate false positive (i.e. rejected real)
    [FP,FPci,FN,FNci] = estimateFalseRatesForOxoGfdr(M.alt_read_count, M.isArtifactMode, fdrStructs(i), Noxo, NoxoCI, PoxoG);
    allPAs(i,:) = [FN FNci(1) FNci(2)]; 
    allRRs(i,:) = [FP FPci(1) FPci(2)]; 
    Ns(i) = length(M.alt_read_count);
    
    %% Estimate PoxoG error.  
    % Note that this estimate is not actually used.  This estimate is
    %  recorded for later analysis.
    [PoxoGEst,PoxoGEstci] = PoxoGestimate(M.alt_read_count, round(M.foxog .* M.alt_read_count),fdrStructs(i).iART, [], Noxo, NoxoCI);
    
    %% Write stdout
    tbl_entry = [strjoin('\t', [{pair.case}; {pair.control};{num2str(length(M.foxog))}; {num2str(length(M_A.foxog))}; ...
        {num2str(length(M_nA.foxog))}; {num2str( sum(fdrStructs(i).cut))}; {num2str(Noxo)}; {num2str(NoxoCI(1))}; {num2str(NoxoCI(2))}; {num2str(Nmut)}; ...
        {num2str(FP)}; {num2str(FPci(1))}; {num2str(FPci(2))};{num2str(FN)}; {num2str(FNci(1))}; {num2str(FNci(2))}; ...
        {num2str(PoxoGEst)}; {num2str(PoxoGEstci(1))}; {num2str(PoxoGEstci(2))} ;{num2str( sum(fdrStructs(i).cut)/length(M.foxog) )};{num2str(lod0)};...
        {num2str(oxoQ)}...
        ]) '\n'];
    fprintf(1,  tbl_entry)
    
    %% Write table entry to table file
    fprintf(tfid, tbl_entry )
    
    %% Generate temporary table for later use in appending (column-wise) to output maf files
    caseMutMap = addMutationEntries(pair, unameMafTable, M, fdrStructs(i), acVal, acs, PoxoG);
    mutationMap = [mutationMap; caseMutMap];
end
fclose(legoFidBefore);
fclose(legoFid);
fclose(tfid);
fclose(caseFid);

%% Write pass and reject counts to a file.
writeFileWithSingleNumber([outputDir '/' outputMAFFilename '.pass_count.txt'], totalPassCount);
writeFileWithSingleNumber([outputDir '/' outputMAFFilename '.reject_count.txt'], totalRejectCount);

%% Write out new maf files (one for pass and one for reject)
%   Additionally, writeFilterMAFFiles appends new columns to the mafFile,
%       so overwrite the mafTable variable to include these.
disp('Writing MAF files and populating table structure with new columns...')
mafTable = writeFilterMAFFiles(mafTable, [outputDir '/' outputMAFFilename], mutationMap);

%% Plot summary data across all cases
if (isGeneratingPlots) && (length(pairs) > 0)
    disp(' Generating summary plots for all cases...')
    [dummy, basename, ext] = fileparts(outputMAFFilename);
    plotFalseRates([outputFigureDir '/' outputMAFFilename '.RR'], [basename ' Rejected Reals'], 'Rejected Real Count', Ns, allRRs);
    plotFalseRates([outputFigureDir '/' outputMAFFilename '.PA'], [basename ' Passing Artifacts'], 'Passing Artifact Count', Ns, allPAs);
    disp('False rate plots completed.')
    
    % Plot summary orchestra plot of all artifact mode mutations
    allFdrStructs = fdrStructs(1);
    for k = 2:length(fdrStructs)
        allFdrStructs = mergeStruct(allFdrStructs, fdrStructs(k));
    end
    [M_summary, M_A_summary, M_nA_summary] = generateMutations(mafTable, [1:length(mafTable.Start_position)]);
    plotArtifactMode(outputFigureDir, [outputMAFFilename ' All Artifact Mode'], M_A_summary, allFdrStructs.cut , [], acs);
    plotNonArtifactMode(outputFigureDir, [outputMAFFilename ' All Non-Artifact Mode'], M_nA_summary, acs)
    disp('Summary Orchestra plots completed.')
    
    Pbefore.label = basename;
    figure('visible','off') 
    [legoStructBeforeSet] = plotMutationSpectrumCategLego1(mafTable,'unit',[outputFigureDir outputMAFFilename '_lego_before'], Pbefore);
    close(gcf);
    Pafter.label = [basename ' Filtered'];
	figure('visible','off') 
    [legoStructAfterSet] = plotMutationSpectrumCategLego1(trimStruct(mafTable, find(~mafTable.oxoGCut)),'unit',[outputFigureDir outputMAFFilename '_lego_after'], Pafter)
    close(gcf);
    disp('Summary Lego plots completed.')
    
    summaryTable = load_table(GetFullPath(tableFilename), char(9));
    figure('visible','off') 
    hist(summaryTable.Proportion_N_filt, [.005:.01:.995])
    xlabel('N_Filtered/N', 'FontSize', 14, 'Interpreter', 'None')
    ylabel('Count', 'FontSize', 14, 'Interpreter', 'None')
    title(['Proportion of Filtered Mutations'], 'FontSize', 16, 'Interpreter', 'None')
    print('-dpng', [outputFigureDir outputMAFFilename '_filteredProportion_hist.png']);
    print(gcf,'-depsc','-loose',[outputFigureDir outputMAFFilename '_filteredProportion_hist.eps'])

    close(gcf)
    disp('Proportion of Filtered Mutations plot completed.')
    
    % oxoQ summary
    %kpass0=find(mafTable.qox<fdrThresh1);
    kpass0=find(mafTable.qox<fdrThresh); % original fdrThreshold band-aid
    kpass=find(~mafTable.oxoGCut);
    tmp=tabulate(sort(mafTable.Tumor_Sample_Barcode));
    t=[];
    t.sample=tmp(:,1);
    t.ntot=cell2mat(tmp(:,2));
    [i m]=ismember(t.sample,mafTable.Tumor_Sample_Barcode);
    t.oxoQ=NaN*t.ntot;
    t.oxoQ(i,1)=mafTable.i_picard_oxoQ(m(m>0));
    i=t.oxoQ>50; % spurious oxoQ when error rate is negative
    t.oxoQ(i)=rand(sum(i),1)*2+40;
    tmp=tabulate(sort(mafTable.Tumor_Sample_Barcode(kpass0)));
    [i m]=ismember(t.sample,tmp(:,1));
    t.npass0=0*t.ntot;
    t.npass0(i)=cell2mat(tmp(m(m>0),2));
    tmp=tabulate(sort(mafTable.Tumor_Sample_Barcode(kpass)));
    [i m]=ismember(t.sample,tmp(:,1));
    t.npass=0*t.ntot;
    t.npass(i)=cell2mat(tmp(m(m>0),2));
	figure('visible','off') 
    subplot('position',[0.1 0.53 0.85 0.39])
    h=semilogy(t.oxoQ,t.ntot*0,'go',t.oxoQ,t.ntot-t.npass0+1,'r+',t.oxoQ,t.ntot-t.npass+1,'bx','markersize',8,'linewidth',2);
    set(h(1),'color',0.5*[0.1 1 0.1])
    set(gca,'xticklabel',[])
    ylabel('Mutations cut +1', 'FontSize', 14, 'Interpreter', 'None')
    legend('No oxoG filter','filter on oxoG q-value only','final filter using oxoQ','location','NE')
    title(['Filter dependence on picard oxoQ'], 'FontSize', 16, 'Interpreter', 'None')
    grid on
    subplot('position',[0.1 0.12 0.85 0.39])
    h=semilogy(t.oxoQ,t.ntot+1,'go',t.oxoQ,t.npass0+1,'r+',t.oxoQ,t.npass+1,'bx','markersize',8,'linewidth',2);
    set(h(1),'color',0.5*[0.1 1 0.1])
    xlabel('picard oxoQ', 'FontSize', 14, 'Interpreter', 'None')    
    ylabel('Mutations passed +1', 'FontSize', 14, 'Interpreter', 'None')
    xL=get(gca,'xlim');    xT=get(gca,'xtick');    xTL=cellstr(get(gca,'xticklabel')); 
    k=find(xT>40); 
    if ~isempty(k)
        xTL{k(end)}='~ 80';
        set(gca,'xticklabel',xTL);
    end
    ylim([min([10; t.npass/2]) max([2e3; t.npass*2])])
    yL=get(gca,'ylim');
    line(40*[1 1],yL,'color',0.5*[.5 0.5 1],'linestyle',':')
    line(xL,xL*0+median(t.npass),'color',0.5*[.5 0.5 1],'linestyle','--')
    text(0.99,0.95,sprintf(' median mutations passed = %d ',round(median(t.npass))),...
        'units','normalized', 'FontSize', 14, 'HorizontalAlignment','right',... 
        'BackgroundColor','w','edgecolor','k')
    
    grid on
    print('-dpng', [outputFigureDir outputMAFFilename '_oxoQ.png']);    
    print(gcf,'-depsc','-loose',[outputFigureDir outputMAFFilename '_oxoQ.eps'])
    close(gcf)
    disp('Picard OxoQ plot completed.')
    
    % allele frequency
    plotAllelicFractionSummaries(outputFigureDir, basename, [outputMAFFilename '_AF'], mafTable)
    disp('Summary AF plots completed.')
    

    
end

%% Clean up
close all    % close all figures so xvnc can terminate
disp('Done.')

    
