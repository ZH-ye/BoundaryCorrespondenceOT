% complie script

boost = 'path/to/boost';
stanmath='path/to/stanmath';
sundials = 'path/to/sundials/include';
Eigen = 'path/to/eigen';

flags = '-O2 -mfma -march=native';



%%%%%%%%%%%%%%%%%%%
LBFGS = './LBFGS';
ipath = {['-I',LBFGS], ...
        ['-I',sundials], ...
    ['-I',stanmath], ...
    ['-I',Eigen], ...
    ['-I',boost]};


mex('-v', ipath{:}, ...
    ['CXXFLAGS=$CXXFLAGS -std=c++14 ',flags],...
    'model.cpp', 'mexSinkhorn.cpp');