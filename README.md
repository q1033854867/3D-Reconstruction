# 3D-Reconstruction
This is a MATLAB project of 3D-Reconstruction. By first running the `prepross.m` and then `main.m`, you can get the cloudpoints of demodata. Based on cloudpoints, you can reconstruct the object using other tools, such as [Geomagic Studio](https://www.3dsystems.com/) or [MeshLab](http://www.meshlab.net/).

## Usage
First, add the project fold with its subfold into your MATLAB search path.<br>
Then, in MATLAB command line, run:
```MATLAB
preprocess
main
```
If you need to recompile C source file: 
```MATLAB
cd lib_c
mex search_points.c
cd ..
```

## Components
- demodata/
- lib/
- lib_c/
- prepeocess.m
- main.m
- rot_axis.jpg
- board.jpg

Two images `rot_axis.jpg` and `board.jpg` are axis of rotation and checherboard, respectively.<br>
Files in `lib/` are MATLAB source file. Files in `lib_c/` are C source file and its binary file which compiled in Windows 10. <br>
The function `preprocess.m` do some preprocessing, such as find nodes of checkerboard, organize nodes into graph, generate mask of ROI automaticlly, get index of rotation axis in image, etc. After you run it and get four important .mat file in top directory, you will be able to run `main.m` and get cloudpoints result finally. <br>
Note! The global variable `config` in `main.m` controls processing mode. Parts of its default values are as follows:
```MATLAB
config.write_into_txt = false;
config.realtime_disp = true;
config.read_img_prev = true;
config.mask_dynamc_adj = true;
config.laser_algorithm = 'basic';
config.save_laser = false;
```
The first two lines mean we do not write cloudpoints into .txt file but realtime display the 3D cloudpoints. The third line means we read all images into memory before processing, which may take up 3.7634GBytes of memory space. We do this because the MATLAB function `imread()` is a little slow. If your memory is not enough, please set this item to false. The fourth line means the mask of ROIs will dynamic adjust its size, but the adjust algorithm is naive, which will be enhanced in fulture. The fifth line means we default use basic extraction algorithm to get the centerline of linear structured laser. Also, there are other alternative algorithms you can choose, see details in `lib/cloudpoints_imgidx.m`. 

## Platform
Tested under Windows 10 and MATLAB R2018b.

You also can conveniently run this project in macOS or Linux. But before you run anything, you should compile all the C source file using `mex filename`. More details in [MATLAB mex documentation](https://ww2.mathworks.cn/help/matlab/ref/mex.html?lang=en)

In addition, we also have another Python version for Raspberry Pi, which aiming at realtime processing and are now in developing.

## Acknowledgements
This work was supported by Independent Innovation Project in Wuhan Uni of Tech.