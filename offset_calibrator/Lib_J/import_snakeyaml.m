[pth,~,~] = fileparts(mfilename('fullpath'));       
try
    import('org.yaml.snakeyaml.*');
    javaObject('Yaml');
catch
    dp = [pth filesep 'external' filesep 'snakeyaml-1.9.jar'];
    if not(ismember(dp, javaclasspath ('-dynamic')))
        javaaddpath(dp); % javaaddpath clears global variables...!?
    end
    import('org.yaml.snakeyaml.*');
end;