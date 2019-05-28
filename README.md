# Boundary Correspondence of Planar Domains for Isogeometric Analysis Based on Optimal Mass Transport[1]

<img style="display:inline;" src="./figures/example.png" height="200" title="example"> <img style="display:inline;" src="./figures/example2.png" height="200" title="example"> <img style="display:inline;" src="./figures/example3.png" height="200" title="example">



An automatic approach to compute a correspondence between the boundaries of a unit square and a planar domain based on optimal mass transport.  The result serves as the input of subsequent domain parameterization [2].

Tested on Windows 10 (matlab 2018b, MinGW) and Ubuntu 18.04 (matlab 2018b, gcc).

### Dependences

matlab toolboxs: 

- curve fitting
- mapping
- optimization (for parametrization [2])

`mex` files are provided on linux and Windows platform. If you want to compile the `mex` file yourself,  [stanmath](https://github.com/stan-dev/math) package is needed.



### References

[1]Zheng, Y., Pan, M., Chen, F., 2019. Boundary correspondence of planar domains for isogeometric analysis based on optimal mass transport. Computer-Aided Design 114, 28–36. 

[2] M. Pan, F. Chen, W. Tong, Low-rank parameterization of planar domains for isogeometric analysis, Computer Aided Geometric Design 63 (2018) 1–16.