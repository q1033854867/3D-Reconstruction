# 3D-Reconstruction
This is a MATLAB project of 3D-Reconstruction. By first running the prepross.m and then main.m, you can get the cloudpoints of demodata. Based on cloudpoints, you can reconstruct the object using other tools such as [MeshLab](http://www.meshlab.net/).

## Components
- demodata
  * A0...499.jpg
- automask
- centerline
- ...

## Platform
Tested under Windows 10 and MATLAB R2018b.<br>
You also can run this project in macOS or Linux. Before you run anything, you should compile all the c source file using `mex filename`. More details in [MATLAB mex documentation](https://ww2.mathworks.cn/help/matlab/ref/mex.html?lang=en)

## Acknowledgements
This work was supported by Independent Innovation Project in Wuhan Uni of Tech.