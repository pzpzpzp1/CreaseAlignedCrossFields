if ispc
    if strcmp(getenv('computername'),'DESKTOP-5UNMMEJ')
        tbb_include_path = 'D:/Documents/MATLAB/tbb-tbb_2019/include';
        tbb_libs = 'C:/Users/pzpzp/Documents/MATLAB/tbb-tbb_2019/build/vs2013/x64/Release-MT';
        mosek_include_path = 'C:\Program Files\Mosek\9.0\tools\platform\win64x86\h';
        mosek_libs = 'C:\Program Files\Mosek\9.0\tools\platform\win64x86\bin';
    else
        error('Fill in your Mosek and TBB directories!');
    end
    
    mex('-cxx', '-O', 'COMPFLAGS="\$COMPFLAGS -MT"',...
        ['-I' tbb_include_path],...
        ['-I' mosek_include_path],...
        ['-L' tbb_libs],...
        ['-L' mosek_libs],...
        '-lfusion64_9_0', '-lmosek64_9_0', 'MosekSoftCrossFields.cpp');
else
    error('Building on anything except pc is not tested. Some restructuring of the commented lines of code below should work though :)');
end

%{
if ismac
    %mex('-cxx', '-O', '-g', 'CXXFLAGS="\$CXXFLAGS -I/usr/local/include -I../../../external/mosek/8/tools/platform/osx64x86/h"', 'LDFLAGS="\$LDFLAGS -L../../../external/mosek/8/tools/platform/osx64x86/bin"', '-lmosek64', '-lfusion64', '-ltbb', 'MultiSdp.cpp');
    %!install_name_tool -change libmosek64.8.1.dylib @loader_path/../../../external/mosek/8/tools/platform/osx64x86/bin/libmosek64.8.1.dylib MultiSdp.mexmaci64
    %!install_name_tool -change libfusion64.8.1.dylib @loader_path/../../../external/mosek/8/tools/platform/osx64x86/bin/libfusion64.8.1.dylib MultiSdp.mexmaci64
elseif isunix
    %mex('-cxx', '-O', '-g', 'CXXFLAGS="\$CXXFLAGS -I/usr/local/include -I../../../external/mosek/8/tools/platform/linux64x86/h"', 'LDFLAGS="\$LDFLAGS -L../../../external/mosek/8/tools/platform/linux64x86/bin -Wl,-rpath-link,../../../external/mosek/8/tools/platform/linux64x86/bin ''-Wl,-rpath=$ORIGIN/../../../external/mosek/8/tools/platform/linux64x86/bin''"', '-lfusion64', '-lmosek64', '-ltbb', 'MultiSdp.cpp');
    %mex('-cxx', '-O', '-g', 'CXXFLAGS="\$CXXFLAGS -I/usr/local/include -I. -I/ascldap/users/pauzhan/Documents/eigen -I/ascldap/users/pauzhan/mosek/9.0/tools/platform/linux64x86/h -I/ascldap/users/pauzhan/mosek/9.0/tools/platform/linux64x86/bin -I/ascldap/users/pauzhan/tbb2019_20190605oss/include"',...
    %    'LDFLAGS="-L/ascldap/users/pauzhan/tbb2019_20190605oss/lib/intel64/gcc4.7 -L/ascldap/users/pauzhan/mosek/9.0/tools/platform/linux64x86/bin -Wl,-rpath-link,/ascldap/users/pauzhan/mosek/9.0/tools/platform/linux64x86/bin -Wl,-rpath=/ascldap/users/pauzhan/mosek/9.0/tools/platform/linux64x86/bin"', ...
    %    '-lfusion64', '-lmosek64', '-ltbb', 'MosekSoftCrossFields2.cpp');
end
%}