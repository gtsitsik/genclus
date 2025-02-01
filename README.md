[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
<div align="center">


# Multi-View Spectral Clustering for Graphs with Multiple View Structures

<img src="Images/genclus_light_mode.png#gh-light-mode-only" alt="Multi-Structure Multi-View Graph" width="500"/>
<img src="Images/genclus_dark_mode.png#gh-dark-mode-only" alt="Multi-Structure Multi-View Graph" width="500"/>

</div>



## Overview
GenClus is a graph clustering method aimed at modeling multi-view graphs where
each view may correspond to a different clustering of nodes. Its design is
built on principled foundations, generalizing various spectral clustering
methods while also facilitating the derivation of computationally efficient
numerical solutions.


This repo contains an implementation of GenClus along with the complete
experiment scripts of the associated paper

|Yorgos Tsitsikas and Evangelos E. Papalexakis. Multi-View Spectral Clustering for Graphs with Multiple View Structures. In _Proceedings of the 2025 SIAM International Conference on Data Mining (SDM)_, 2025|
| :--- |

The full version of this paper can be accessed at https://arxiv.org/abs/2501.11422

Also, BibTeX users can cite this work as follows:
```bibtex
@inproceedings{genclus,
  title={Multi-View Spectral Clustering for Graphs with Multiple View Structures},
  author={Tsitsikas, Yorgos and Papalexakis, Evangelos E.},
  booktitle={Proceedings of the 2025 SIAM International Conference on Data Mining (SDM)},
  year={2025},
}
```

## **Dependencies**

### **genclus.m**
- [MATLAB](https://www.mathworks.com/products/matlab.html) - **Required** (Tested on MATLAB 2019b)
- [MTIMESX](https://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support) - **Optional**


### **clustering_quality_on_artificial_data.m**
- [MATLAB](https://www.mathworks.com/products/matlab.html) - **Required** (Tested on MATLAB 2019b)
- [Multi-Graph Explorer](https://github.com/gtsitsik/multi-graph-explorer) - **Required**

### **real_world_case_study.m**
- [Openflights data](https://openflights.org/data):  
  - `airlines.dat`, `airports.dat`, `countries.dat`, `routes.dat` - **Required**
- [Continents data](https://www.kaggle.com/datasets/andradaolteanu/country-mapping-iso-continent-region):  
  - `continents2.csv` - **Required**
- [MATLAB](https://www.mathworks.com/products/matlab.html) - **Required** (Tested on MATLAB 2019b; Currently fails in MATLAB 2022b)
- [Multi-Graph Explorer](https://github.com/gtsitsik/multi-graph-explorer) - **Required**
- [Tensor Toolbox](https://www.tensortoolbox.org/) - **Required**  

### **execution_time_on_artificial_data.m**
- [MATLAB](https://www.mathworks.com/products/matlab.html) - **Required** (Tested on MATLAB 2019b)
- [Multi-Graph Explorer](https://github.com/gtsitsik/multi-graph-explorer) - **Required**


## License
[3-Clause BSD](./LICENSE)
