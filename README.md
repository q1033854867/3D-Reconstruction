# 3D-Reconstruction
This is a MATLAB project of 3D-Reconstruction. By first running the prepross.m and then main.m, you can get the cloudpoints of demodata. Based on cloudpoints, you can reconstruct the object using other tools such as [MeshLab](http://www.meshlab.net/).

## Components
- demodata
  * A0...499.jpg
- lib
- lib_c
- lib_aux
- prepeocess.m
- main.m
- rot_axis.jpg
- board.jpg

Total 500 demodata are in the `demodata`. And the measured object is an iron column. <br>
Files in `lib` and `lib_aux` are .m source file. Files in `lib_c` are c source file and its binary file which compiled in Windows 10. <br>
Two images `rot_axis.jpg` and `board.jpg` are axis of rotation and checherboard, respectively.<br>
`prepeocess.m` do some preprocessing. When you run it and get four important .mat file in root folder, you will be able to run `main.m` and get cloudpoints result finally. <br>
Note! The variable `config` in `main.m` control processing mode. Part of its default values are as follows:
```matlab
config.write_into_txt = false;
config.realtime_disp = true;
config.read_img_prev = true;
config.mask_dynamc_adj = true;
```
The first two lines mean we do not write cloudpoints into .txt file but realtime display the result. The third lines means we read all demo images into memory before processing loop. We do this because the MATLAB function `imread` is a little slow. The last lines means the mask of ROIs will dynamic adjust its size, but which is not yet perfect.

## Platform
Tested under Windows 10 and MATLAB R2018b.

You also can conveniently run this project in macOS or Linux. But before you run anything, you should compile all the c source file using `mex filename`. More details in [MATLAB mex documentation](https://ww2.mathworks.cn/help/matlab/ref/mex.html?lang=en)

In addition, we also have another Python version for Raspberry Pi, which aiming at realtime processing and are now in developing.

## Acknowledgements
This work was supported by Independent Innovation Project in Wuhan Uni of Tech.