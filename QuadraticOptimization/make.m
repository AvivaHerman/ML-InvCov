function make
%MAKE Build the CD package.
%  MAKE compiles all CD MEX functions, using Matlab's default MEX
%  compiler. If the MEX compiler has not been set-up before, please run
%
%    mex -setup
%
%  before using this MAKE file.

clc

% detect platform 
compstr = computer;
is64bit = strcmp(compstr(end-1:end),'64');


% compilation parameters
compile_params = cell(0);
if (is64bit)
  compile_params{1} = '-largeArrayDims';
end


% Compile files %
cdsources = {'utility.c', 'sort.c', 'coordinate_descent.c'};

disp('Compiling solveL2L1_CD_mex...');
mex('solveL2L1_CD_mex.c', cdsources{:}, compile_params{:});

disp('Compiling solveL2L1_MLCD_mex...');
mex('solveL2L1_MLCD_mex.c', cdsources{:}, compile_params{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdsources = {'utility.c', 'sort.c', 'm_coordinate_descent.c'};

disp('Compiling solve_M_CD_mex...');
mex('solve_M_CD_mex.c', cdsources{:}, compile_params{:});

disp('Compiling solve_M_MLCD_mex...');
mex('solve_M_MLCD_mex.c', cdsources{:}, compile_params{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lssources = {'utility.c', 'sort.c'};

disp('Compiling ExactLinesearchFull...');
mex('ExactLinesearchFull.c', lssources{:}, compile_params{:});
