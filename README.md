# Structured Massive Access for Scalable Cell-Free Massive MIMO Systems

This is a code package is related to the following scientific article:

Shuaifei Chen, Jiayi Zhang, Emil Björnson, Jing Zhang, and Bo Ai, "[Structured Massive Access for Scalable Cell-Free Massive MIMO Systems](https://ieeexplore.ieee.org/abstract/document/9174860)," *IEEE Journal on Selected Areas in Communications*, vol. 39, no. 4, pp. 1086-1100, April 2021.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

How to meet the demand for increasing number of users, higher data rates, and stringent quality-of-service (QoS) in the beyond fifth-generation (B5G) networks? Cell-free massive multiple-input multiple-output (MIMO) is considered as a promising solution, in which many wireless access points cooperate to jointly serve the users by exploiting coherent signal processing. However, there are still many unsolved practical issues in cell-free massive MIMO systems, whereof scalable massive access implementation is one of the most vital. In this paper, we propose a new framework for structured massive access in cell-free massive MIMO systems, which comprises one initial access algorithm, a partial large-scale fading decoding (P-LSFD) strategy, two pilot assignment schemes, and one fractional power control policy. New closed-form spectral efficiency (SE) expressions with maximum ratio (MR) combining are derived. The simulation results show that our proposed framework provides high SE when using local partial minimum mean-square error (LP-MMSE) and MR combining. Specifically, the proposed initial access algorithm and pilot assignment schemes outperform their corresponding benchmarks, P-LSFD achieves scalability with a negligible performance loss compared to the conventional optimal large-scale fading decoding (LSFD), and scalable fractional power control provides a controllable trade-off between user fairness and the average SE.


## Content of Code Package

The package generates the simulation SE results which are used in Figure 6, Figure 7, Figure 8, and Figure 9. To be specific:

- `simulationSE`: Main function;
- `functionGenerateR`: Generate spatial correlation matrix and large-scale fading coefficients;
  - `functionRlocalscattering`: Generate the spatial correlation matrix for the local scattering model;
- `functionAPselection`: Perform initial access and generate AP selection results;
- `functionUEgroup`: Perform pilot assignment by using the proposed User-Group scheme;
- `functionPowerControl`: Perform uplink fractional power control; 
- `functionChannelEstimates`: Perform MMSE channel estimation;
- `functionComputeSE_MMSE_pLSFD`: Compute simulation SE results while using LP-MMSE combining and P-LSFD.

See each file for further documentation.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
