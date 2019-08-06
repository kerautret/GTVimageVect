# GTVimageVect
Source code of "Geometric Total Variation for Image Vectorization, Zooming and Pixel Art Depixelizing"



## Installation


   To use this source code, you need to install the following dependencies:
      - DGTal library, on Linux/Mac just follow these step:
        `git clone git@github.com:DGtal-team/DGtal.git`
        `cd DGtal; mkdir build; cd build`
        `cmake .. -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE:string="Release";`
        `make;`
        Not mandatory `make install`.