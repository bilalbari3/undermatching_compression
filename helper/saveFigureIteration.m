function saveFigureIteration(fHandle, saveFigLoc, figureName, varargin)

p = inputParser;
p.addParameter('figureSubName', '')
p.addParameter('FigureSize','')
p.addParameter('AutoSave_Flag', false)
p.parse(varargin{:});

if p.Results.AutoSave_Flag == false
    answer = questdlg('Save figure?');
else
    answer = 'Yes';
end
switch answer
    case 'Yes'

        if strcmp(p.Results.FigureSize, 'max') | strcmp(p.Results.FigureSize, 'full')
        %     set(fHandle, 'Position', get(0,'Screensize'))
            set(fHandle, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        elseif strcmp(p.Results.FigureSize, 'half')
            set(fHandle, 'units','normalized','outerposition',[0 0 0.5 1]);
        end
        fHandle.Renderer = 'Painters';

        tmp = dir(saveFigLoc);
        nextVersion = sum(contains({tmp.name}, '.pdf') & contains({tmp.name}, [figureName '_' p.Results.figureSubName]));
        if nextVersion >= 10
            version = int2str(nextVersion);
        else
            version = ['0' int2str(nextVersion)];
        end
        if isempty(p.Results.figureSubName)
            saveName = fullfile(saveFigLoc, [figureName '_v' version]);
        else
            saveName = fullfile(saveFigLoc, [figureName '_' p.Results.figureSubName '_v' version]);
        end

        saveFigurePDF(fHandle, [saveName '.pdf']);
        savefig(fHandle, saveName);

        fprintf('\nSaved as %s\n', saveName);
    otherwise
        fprintf('Figure not saved\n')
end