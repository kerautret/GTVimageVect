# GTVimageVect
Source code of "Geometric Total Variation for Image Vectorization, Zooming and Pixel Art Depixelizing"



## Installation


To use this source code, you need to install the following dependencies:
   - **DGTal library**, on Linux/Mac just follow these steps:
     - Install (if not already present), the boost dependancies (see https://www.boost.org).
     For instance on linux you can install the following package: sudo apt-get install `libboost-dev libboost-system-dev libboost-program-options-dev`(the `libboost-program-options-dev` is not mandatory for DGtal but is used by our project)
     - `git clone git@github.com:DGtal-team/DGtal.git`
     - `cd DGtal; mkdir build; cd build`
     - `cmake .. -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE:string="Release";`
     - `make;`
     - Not mandatory `make install`.
     - For any problem, don't hesitate to contact the DGtal team on the [GitHub repository](https://github.com/DGtal-team/DGtal)
     On windows or more details see instructions [here](https://dgtal-team.github.io/doc-nightly/moduleBuildDGtal.html).
   - **cairo**
      The installation details (Linux/MacOS/Windows) are on the website: https://www.cairographics.org/download/
   - **boost programm options** by default the application use the boost program options extension of 
