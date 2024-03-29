function klu_make
%KLU_MAKE compiles the KLU mexFunctions
%
% Example:
%   klu_make
%
% KLU relies on AMD, COLAMD, and BTF for its ordering options, and can
% optionally use CHOLMOD, CCOLAMD, CAMD, and METIS as well.  METIS must
% be placed in ../../metis-4.0 alongside KLU, or it will not be used.
% See http://www-users.cs.umn.edu/~karypis/metis for a copy of METIS 4.0.1.
%
% You must type the klu_make command while in the KLU/MATLAB directory.
%
% See also klu.

% Copyright 2004-2009, Univ. of Florida

metis_path = '../../metis-4.0' ;
with_cholmod = exist ([metis_path '/Lib'], 'dir') ;

details = 0 ;       % if 1, print details of each command

d = '' ;
if (~isempty (strfind (computer, '64')))
    % 64-bit MATLAB
    d = '-largeArrayDims' ;
end

fprintf ('Compiling KLU ') ;
kk = 0 ; 

include = '-I. -I../../AMD/Include -I../../COLAMD/Include -I../Include -I../../SuiteSparse_config -I../../BTF/Include' ;

if (with_cholmod)
    include = [include ' -I../../CCOLAMD/Include -I../../CAMD/Include -I../../CHOLMOD/Include -I../../SuiteSparse_config -I' metis_path '/Lib -I../User'] ;
end

% do not attempt to compile CHOLMOD with large file support (not needed)
include = [include ' -DNLARGEFILE'] ;

% fix the METIS 4.0.1 rename.h file
if (with_cholmod)
    fprintf ('with CHOLMOD, CAMD, CCOLAMD, and METIS\n') ;
    f = fopen ('rename.h', 'w') ;
    if (f == -1)
        error ('unable to create rename.h in current directory') ;
    end
    fprintf (f, '/* do not edit this file; generated by klu_make */\n') ;
    fprintf (f, '#undef log2\n') ;
    fprintf (f, '#include "%s/Lib/rename.h"\n', metis_path) ;
    fprintf (f, '#undef log2\n') ;
    fprintf (f, '#define log2 METIS__log2\n') ;
    fprintf (f, '#include "mex.h"\n') ;
    fprintf (f, '#define malloc mxMalloc\n') ;
    fprintf (f, '#define free mxFree\n') ;
    fprintf (f, '#define calloc mxCalloc\n') ;
    fprintf (f, '#define realloc mxRealloc\n') ;
    fclose (f) ;
    include = ['-DNSUPERNODAL -DNMODIFY -DNMATRIXOPS -DNCHECK ' include] ;
else
    fprintf ('without CHOLMOD, CAMD, CCOLAMD, and METIS\n') ;
    include = ['-DNCHOLMOD ' include] ;
end

include = strrep (include, '/', filesep) ;

amd_src = { ...
    '../../AMD/Source/amd_1', ...
    '../../AMD/Source/amd_2', ...
    '../../AMD/Source/amd_aat', ...
    '../../AMD/Source/amd_control', ...
    '../../AMD/Source/amd_defaults', ...
    '../../AMD/Source/amd_dump', ...
    '../../AMD/Source/amd_global', ...
    '../../AMD/Source/amd_info', ...
    '../../AMD/Source/amd_order', ...
    '../../AMD/Source/amd_postorder', ...
    '../../AMD/Source/amd_post_tree', ...
    '../../AMD/Source/amd_preprocess', ...
    '../../AMD/Source/amd_valid' } ;

colamd_src = {
    '../../COLAMD/Source/colamd', ...
    '../../COLAMD/Source/colamd_global' } ;

if (with_cholmod)

    camd_src = { ...
        '../../CAMD/Source/camd_1', ...
        '../../CAMD/Source/camd_2', ...
        '../../CAMD/Source/camd_aat', ...
        '../../CAMD/Source/camd_control', ...
        '../../CAMD/Source/camd_defaults', ...
        '../../CAMD/Source/camd_dump', ...
        '../../CAMD/Source/camd_global', ...
        '../../CAMD/Source/camd_info', ...
        '../../CAMD/Source/camd_order', ...
        '../../CAMD/Source/camd_postorder', ...
        '../../CAMD/Source/camd_preprocess', ...
        '../../CAMD/Source/camd_valid' } ;

    ccolamd_src = {
        '../../CCOLAMD/Source/ccolamd', ...
        '../../CCOLAMD/Source/ccolamd_global' } ;

    metis_src = {
        'Lib/balance', ...
        'Lib/bucketsort', ...
        'Lib/ccgraph', ...
        'Lib/coarsen', ...
        'Lib/compress', ...
        'Lib/debug', ...
        'Lib/estmem', ...
        'Lib/fm', ...
        'Lib/fortran', ...
        'Lib/frename', ...
        'Lib/graph', ...
        'Lib/initpart', ...
        'Lib/kmetis', ...
        'Lib/kvmetis', ...
        'Lib/kwayfm', ...
        'Lib/kwayrefine', ...
        'Lib/kwayvolfm', ...
        'Lib/kwayvolrefine', ...
        'Lib/match', ...
        'Lib/mbalance2', ...
        'Lib/mbalance', ...
        'Lib/mcoarsen', ...
        'Lib/memory', ...
        'Lib/mesh', ...
        'Lib/meshpart', ...
        'Lib/mfm2', ...
        'Lib/mfm', ...
        'Lib/mincover', ...
        'Lib/minitpart2', ...
        'Lib/minitpart', ...
        'Lib/mkmetis', ...
        'Lib/mkwayfmh', ...
        'Lib/mkwayrefine', ...
        'Lib/mmatch', ...
        'Lib/mmd', ...
        'Lib/mpmetis', ...
        'Lib/mrefine2', ...
        'Lib/mrefine', ...
        'Lib/mutil', ...
        'Lib/myqsort', ...
        'Lib/ometis', ...
        'Lib/parmetis', ...
        'Lib/pmetis', ...
        'Lib/pqueue', ...
        'Lib/refine', ...
        'Lib/separator', ...
        'Lib/sfm', ...
        'Lib/srefine', ...
        'Lib/stat', ...
        'Lib/subdomains', ...
        'Lib/timing', ...
        'Lib/util' } ;

    for i = 1:length (metis_src)
        metis_src {i} = [metis_path '/' metis_src{i}] ;
    end

    cholmod_src = {
        '../../CHOLMOD/Core/cholmod_aat', ...
        '../../CHOLMOD/Core/cholmod_add', ...
        '../../CHOLMOD/Core/cholmod_band', ...
        '../../CHOLMOD/Core/cholmod_change_factor', ...
        '../../CHOLMOD/Core/cholmod_common', ...
        '../../CHOLMOD/Core/cholmod_complex', ...
        '../../CHOLMOD/Core/cholmod_copy', ...
        '../../CHOLMOD/Core/cholmod_dense', ...
        '../../CHOLMOD/Core/cholmod_error', ...
        '../../CHOLMOD/Core/cholmod_factor', ...
        '../../CHOLMOD/Core/cholmod_memory', ...
        '../../CHOLMOD/Core/cholmod_sparse', ...
        '../../CHOLMOD/Core/cholmod_transpose', ...
        '../../CHOLMOD/Core/cholmod_triplet', ...
        '../../CHOLMOD/Cholesky/cholmod_amd', ...
        '../../CHOLMOD/Cholesky/cholmod_analyze', ...
        '../../CHOLMOD/Cholesky/cholmod_colamd', ...
        '../../CHOLMOD/Cholesky/cholmod_etree', ...
        '../../CHOLMOD/Cholesky/cholmod_postorder', ...
        '../../CHOLMOD/Cholesky/cholmod_rowcolcounts', ...
        '../../CHOLMOD/Partition/cholmod_ccolamd', ...
        '../../CHOLMOD/Partition/cholmod_csymamd', ...
        '../../CHOLMOD/Partition/cholmod_camd', ...
        '../../CHOLMOD/Partition/cholmod_metis', ...
        '../../CHOLMOD/Partition/cholmod_nesdis' } ;
end

btf_src = {
    '../../BTF/Source/btf_maxtrans', ...
    '../../BTF/Source/btf_order', ...
    '../../BTF/Source/btf_strongcomp' } ;

klu_src = {
    '../Source/klu_free_symbolic', ...
    '../Source/klu_defaults', ...
    '../Source/klu_analyze_given', ...
    '../Source/klu_analyze', ...
    '../Source/klu_memory' } ;

if (with_cholmod)
    klu_src = [klu_src { '../User/klu_l_cholmod' }] ;                       %#ok
end

klu_zlsrc = {
    '../Source/klu', ...
    '../Source/klu_kernel', ...
    '../Source/klu_dump', ...
    '../Source/klu_factor', ...
    '../Source/klu_free_numeric', ...
    '../Source/klu_solve', ...
    '../Source/klu_scale', ...
    '../Source/klu_refactor', ...
    '../Source/klu_tsolve', ...
    '../Source/klu_diagnostics', ...
    '../Source/klu_sort', ...
    '../Source/klu_extract', ...
    } ;

klu_lobj = {
    'klu_l', ...
    'klu_l_kernel', ...
    'klu_l_dump', ...
    'klu_l_factor', ...
    'klu_l_free_numeric', ...
    'klu_l_solve', ...
    'klu_l_scale', ...
    'klu_l_refactor', ...
    'klu_l_tsolve', ...
    'klu_l_diagnostics', ...
    'klu_l_sort', ...
    'klu_l_extract', ...
    } ;

klu_zlobj = {
    'klu_zl', ...
    'klu_zl_kernel', ...
    'klu_zl_dump', ...
    'klu_zl_factor', ...
    'klu_zl_free_numeric', ...
    'klu_zl_solve', ...
    'klu_zl_scale', ...
    'klu_zl_refactor', ...
    'klu_zl_tsolve', ...
    'klu_zl_diagnostics', ...
    'klu_zl_sort', ...
    'klu_zl_extract', ...
    } ;

try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end

if (pc)
    % Windows does not have drand48 and srand48, required by METIS.  Use
    % drand48 and srand48 in CHOLMOD/MATLAB/Windows/rand48.c instead.
    obj_extension = '.obj' ;
    cholmod_src = [cholmod_src {'../../CHOLMOD/MATLAB/Windows/rand48'}] ;
    include = [include ' -I../../CHOLMOD/MATLAB/Windows'] ;
else
    obj_extension = '.o' ;
end

% compile each library source file
obj = ' ' ;

source = [amd_src btf_src klu_src colamd_src] ;
if (with_cholmod)
    source = [metis_src ccolamd_src camd_src cholmod_src source] ;
end

for f = source
    fs = strrep (f {1}, '/', filesep) ;
    slash = strfind (fs, filesep) ;
    if (isempty (slash))
        slash = 1 ;
    else
        slash = slash (end) + 1 ;
    end
    o = fs (slash:end) ;
    obj = [obj  ' ' o obj_extension] ;                                      %#ok
    s = sprintf ('mex %s -DDLONG -O %s -c %s.c', d, include, fs) ;
    kk = do_cmd (s, kk, details) ;
end

for k = 1:length(klu_zlsrc)
    ff = strrep (klu_zlsrc {k}, '/', filesep) ;
    slash = strfind (ff, filesep) ;
    if (isempty (slash))
        slash = 1 ;
    else
        slash = slash (end) + 1 ;
    end
    o = ff (slash:end) ;
    s = sprintf ('mex %s -DDLONG -O %s -c %s.c', d, include, ff) ;
    kk = do_cmd (s, kk, details) ;
    lobj = klu_lobj {k} ;
    obj = [obj  ' ' lobj obj_extension] ;                                   %#ok
    mvfile ([o obj_extension], [lobj obj_extension]) ;
    s = sprintf ('mex %s -DDLONG -DCOMPLEX -O %s -c %s.c', d, include, ff) ;
    kk = do_cmd (s, kk, details) ;
    zlobj = klu_zlobj {k} ;
    obj = [obj  ' ' zlobj obj_extension] ;                                  %#ok
    mvfile ([o obj_extension], [zlobj obj_extension]) ;
end

% compile the KLU mexFunction
s = sprintf ('mex %s -DDLONG -O %s -output klu klu_mex.c', d, include) ;
s = [s obj] ;                                                               %#ok
kk = do_cmd (s, kk, details) ;

% clean up
s = ['delete ' obj] ;
do_cmd (s, kk, details) ;

fprintf ('\nKLU successfully compiled\n') ;

%-------------------------------------------------------------------------------

function rmfile (file)
% rmfile:  delete a file, but only if it exists
if (length (dir (file)) > 0)                                                %#ok
    delete (file) ;
end

%-------------------------------------------------------------------------------

function cpfile (src, dst)
% cpfile:  copy the src file to the filename dst, overwriting dst if it exists
rmfile (dst)
if (length (dir (src)) == 0)    %#ok
    fprintf ('File does not exist: %s\n', src) ;
    error ('File does not exist') ;
end
try
    copyfile (src, dst) ;
catch ME
    % ignore errors of the form "cp: preserving permissions: ...
    % Operation not supported".  rethrow all other errors.
    if (isempty (strfind (ME.message, 'Operation not supported')))
        rethrow (ME) ;
    end
end

%-------------------------------------------------------------------------------

function mvfile (src, dst)
% mvfile:  move the src file to the filename dst, overwriting dst if it exists
cpfile (src, dst) ;
rmfile (src) ;

%-------------------------------------------------------------------------------
function kk = do_cmd (s, kk, details)
%DO_CMD: evaluate a command, and either print it or print a "."
if (details)
    fprintf ('%s\n', s) ;
else
    if (mod (kk, 60) == 0)
        fprintf ('\n') ;
    end
    kk = kk + 1 ;
    fprintf ('.') ;
end
eval (s) ;

