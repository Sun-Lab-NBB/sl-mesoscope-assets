# sl-mesoscope-assets
Provides MATLAB assets used to acquire experiment data with ScanImage-controlled 2-Photon Random Access Mesoscope 
(2P-RAM).

![matlab](https://img.shields.io/badge/matlab-R2022b%2B-orange)
![license](https://img.shields.io/badge/license-GPLv3-blue)

___

## Detailed Description

This project provides additional assets that are not originally installed on the Mesoscope control PC sourced from 
[MBF Bioscience](https://www.mbfbioscience.com/). Primarily, the assets exposed by this project are used to detect and 
correct imaging plane motion and are originally developed in the 
[Pachitariu and Stringer lab](https://mouseland.github.io/). 

Additionally, the project provides the 'setupAcquisition' function, developed in the Sun lab to automate the preparation
for data acquisition and enable the [sl-experiment](https://github.com/Sun-Lab-NBB/sl-experiment) library to 
bidirectionally interface with the ScanImage software that controls the Mesoscope during runtime.

___

## Dependencies

### Main Dependency
- [MATLAB](https://www.mathworks.com/products/matlab.html) version R2022b+ with 
[ScanImage](https://www.mbfbioscience.com/products/scanimage/) version 2023.1.0 (Premium). This code will likely work 
with later MATLAB and ScanImage versions.

### Additional Dependencies
These additional toolbox dependencies must be installed via the MATLAB interface before using assets from this project:
- [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) to support fast online 
motion detection and correction.

### Hardware Dependencies
**Note!** Typically, the computer used to run this code is purchased from MBF Bioscience (together with ScanImage 
license). In this case, the vendor should satisfy all hardware dependencies and may be requested to also resolve all 
software dependencies.

- [Nvidia GPU](https://www.nvidia.com/en-us/). This library uses GPU hardware acceleration to support online motion 
detection and correction by processing Mesoscope-acquired images. The library was tested with **Nvidia GTX 1660**.
___

## Installation

1. Download this repository to the Mesoscope control PC. Use the latest stable release from 
   [GitHub](https://github.com/Sun-Lab-NBB/sl-mesoscope-assets/releases), as it always reflects the current state of 
   our data acquisition hardware and software.
2. Open MATLAB and navigate to the **Command Window**.
3. Use the `addpath("PATH_TO_THIS_REPOSITORY")` command, replacing the **PATH_TO_THIS_REPOSITORY** with the **absolute**
   path to the downloaded and unpacked repository on the local machine.

If installation worked as expected, the `setupAcquisition()` MATLAB function should now be available for calling from 
the Command Window.
___

## Usage

All assets from this library are intended to be used via the `setupAcquisition` function, which automatically configures
motion detection and correction. It also generates the motion estimation files as part of the acquisition preparation 
sequence.

The [sl-experiment](https://github.com/Sun-Lab-NBB/sl-experiment) advises the user to call the setup function as part of 
the broader Mesoscope preparation sequence. While it is possible to use provided assets outside the sl-experiment 
context, this is not the intended use pattern.

### Acquisition Setup Flow

In most cases, `setupAcquisition` function consists of three major steps:
1. **Motion estimation setup**. During this step, the function configures the acquisition according to the user-defined 
   parameters and establishes the single-plane or the z-stack of planes to be acquired at runtime. Then, it generates a 
   set of additional sub-planes within 20 microns above and below each target plane and acquires ~ 20 images at each 
   plane (traversing the entire stack on each acquisition pass). The produced z-stack is then used to set up the 
   MotionEstimator file used to detect and correct motion in the X, Y and Z axes.
2. **High-definition zstack acquisition**. After setting up the motion detection, the function increases the resolution
   of the target ROIs by 1.5 Scale Factor and repeats the zstack acquisition. This generates a high-definition 
   zstack.tiff file kept with the TIFF files acquired during experiment runtime.
3. **Data Acquisition**. Finally, the function configures the acquisition and motion detection parameters for the 
   runtime and enters the acquisition loop. While in the loop, the function starts or stops the acquisition depending on 
   the marker file(s) created by the sl-experiment library (see sl-experiment ReadMe for details).

**Note!** The function can also be used to resume an interrupted runtime by calling it with an optional `recovery` 
argument set to **true**. In this runtime mode, the function skips steps 1 and 2, going directly to step 3.
___

## API Documentation

Use the `doc` MATLAB command to generate the API documentation for all assets distributed under this project. 
After generating the documentation, use the `help` MATLAB command to view the documentation at any time.

**Note!** All users are highly encouraged to check the setupAcquisition documentation to familiarize themselves with 
function arguments and the range of supported runtime configuration options.

___

## Versioning

This project uses [semantic versioning](https://semver.org/). For the versions available, see the 
[tags on this repository](https://github.com/Sun-Lab-NBB/sl-mesoscope-assets/tags).

---

## Authors

- Ivan Kondratyev ([Inkaros](https://github.com/Inkaros))
- Marius Pachitariu
- Carsen Stringer
- Georg Jaindl

---

## License

This project is licensed under the GPL3 License: see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- All [Sun Lab](https://neuroai.github.io/sunlab/) members for providing the inspiration and comments during the
  development of this project.
- The members of the [Pachitariu and Stringer lab](https://mouseland.github.io/) that developed the motion detection
  and correction code featured in this project and the original acquisition setup function.

---