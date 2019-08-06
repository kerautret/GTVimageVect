# GTVimageVect
Source code of "Geometric Total Variation for Image Vectorization, Zooming and Pixel Art Depixelizing" Authored by  Bertrand Kerautret and Jacques-Oliver Lachaud

You can also access to [IPOL](http://www.ipol.im) online demonstration [here](https://ipolcore.ipol.im/demo/clientApp/demo.html?id=77777000076)

Build status:
 - Linux/MacOS [![Build Status](https://travis-ci.org/kerautret/GTVimageVect.svg?branch=master)](https://travis-ci.org/kerautret/GTVimageVect)
 - Windows [![Build status](https://ci.appveyor.com/api/projects/status/i1crefqj9j1e3lw2?svg=true)](https://ci.appveyor.com/project/kerautret/gtvimagevect)

## Installation


To use this source code, you need to install the following dependencies:
   - **DGTal library** (current version or at least [commit 0e13036](https://github.com/DGtal-team/DGtal/commit/0e13036afedee920373a2460afd02e2a21660baa)), on Linux/Mac just follow these steps:
     - Install (if not already present), the boost dependancies (see https://www.boost.org).
     For instance on linux you can install the following package: sudo apt-get install `libboost-dev libboost-system-dev libboost-program-options-dev`(the `libboost-program-options-dev` is not mandatory for DGtal but is used by our project)
     - `git clone git@github.com:DGtal-team/DGtal.git`
     - `cd DGtal; mkdir build; cd build`
     - `cmake .. -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE:string="Release";`
     - `make;` eventually:  `make install`.

     For any problem, don't hesitate to contact the DGtal team on the [GitHub repository](https://github.com/DGtal-team/DGtal).
     On windows or more details see instructions [here](https://dgtal-team.github.io/doc-nightly/moduleBuildDGtal.html).
   - **cairo**
      On *Linux*, the installation can be done from a package mananger (see [here](https://www.cairographics.org/download/) for more details.
      On *Windows*, you can follow the steps of the cairo website or run the script or upload the binary from the [cairo-windows](https://github.com/preshing/cairo-windows) repository. You can also consult the appveyor configuration file ([appveyor.yml](https://github.com/kerautret/GTVimageVect/blob/master/appveyor.yml) that uses it.
   - **boost programm options** by default the application use the package `libboost-program-options-dev`. You can use you default package manager to install it.


The installation of the code is then:
   - From the git command:
     - `git clone git@github.com:kerautret/GTVimageVect.git`
     - `cd GTVimageVect; mkdir build; cd build;`

     Then you start to build the code: (you can remove the DGtal path if you make a global installation of DGtal.
     - `cmake .. -DDGtal_DIR="/fullpath_to_yourParent_DGtal_dir/DGtal/build" -DCMAKE_BUILD_TYPE:string="Release"`
     - `make`
     
## Typical Use:
   The algorithm is run from executable `tv-triangulation-color`. For instance you can use it follows from the build directory:
   `tv-triangulation-color -i ../Input/dolphin.ppm -b 16 -D 16 -o output -C result.eps`
    Then you should obtain:
    <table>
    <tr><td><img width="200" src="https://user-images.githubusercontent.com/772865/62563570-931eb300-b883-11e9-8ee6-c6054d60040a.png"></td>
    <td><img width="200" src="https://user-images.githubusercontent.com/772865/62563720-e85ac480-b883-11e9-982c-01e3dedc316b.png"></td>
    </tr>
    <tr> <td> source</td> <td>result </td> </tr>
</table>
                                           


## Reproduction of Paper Figures

