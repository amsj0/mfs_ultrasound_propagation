function run_from_inputfile(path_to_input)

eval(fileread([path_to_input,'inputfile.txt']));

GRIDX_VEC = str2num(MODEL_SCL)*linspace(str2num(GRIDX_INI),str2num(GRIDX_FIN),str2num(GRIDX_DEL));
GRIDY_VEC = str2num(MODEL_SCL)*linspace(str2num(GRIDY_INI),str2num(GRIDY_FIN),str2num(GRIDY_DEL));

CENTX_VEC = linspace(str2num(CENTX_INI),str2num(CENTX_FIN),str2num(CENTX_DEL));
CENTY_VEC = linspace(str2num(CENTY_INI),str2num(CENTY_FIN),str2num(CENTY_DEL));

RATSP = [str2num(RATSP_INI),str2num(RATSP_DEL),str2num(RATSP_FIN)];
RATDS = [str2num(RATDS_INI),str2num(RATDS_DEL),str2num(RATDS_FIN)];
RATPS = [str2num(RATPS_INI),str2num(RATPS_DEL),str2num(RATPS_FIN)];
    
preamble

g.path_to_input = path_to_input;
% parpool

% parfor ii=1:length(GRIDX_VEC)
for ii=1:length(GRIDX_VEC)
    ir = GRIDX_VEC(ii);
    assignin('base','g',g);
	for jj=GRIDY_VEC
        str = strcat(RATSP_INI,'_',RATSP_DEL,'_',RATSP_FIN,'_',...
        RATDS_INI,'_',RATDS_DEL,'_',RATDS_FIN,'_',...
        RATPS_INI,'_',RATPS_DEL,'_',RATPS_FIN,'_',...
        int2str(jj),'_',int2str(ir),'_',int2str(str2num(RADIA_SIZ)*100),'_',ELEME_WAV,'.mat');
        disp(str)
            try
            %if(~exist(str,'file'))
                multisource_shrink(RATSP,RATDS,RATPS,jj,ir,str2num(RADIA_SIZ),str2num(ELEME_WAV));
            %end
            catch ME
            rethrow(ME)
            end
	end
end

%delete(loopobj)

%exit;

end