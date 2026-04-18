# SBTD: Block term decomposition of streaming tensors

We propose a novel tensor tracking method called SBTD for factorizing tensors derived from multidimensional data streams under the BTD format. Thanks to the alternating optimization framework, SBTD first applies a regularized least-squares solver to estimate the temporal factor of the underlying streaming tensor. Then, SBTD adopts an adaptive filter to track the nontemporal tensor factors over time by minimizing a weighted
least-squares cost function.

<img width="1007" height="344" alt="SBTD" src="https://github.com/user-attachments/assets/e2f2f992-fb18-4d92-bbdf-c6260b00f4af" />


## Demo
Please run 
+ `demo_comparison.m`: To illustrate the performance of SBTD in comparsion with [BTD-ALS](https://epubs.siam.org/doi/10.1137/070690730) and [onlineBTD](https://ieeexplore.ieee.org/document/9260061).
+ `demo_noise.m`: To illustrate the effect of Gaussian noise on the performance of SBTD
+ `demo_time_varying.m`: To illustrate the ability of SBTD in non-stationary environments

## Reference
This code is free and open source for research purposes. If you use this code, please acknowledge the following paper.

[1]  **L.T. Thanh**, K. Abed-Meraim, P. Ravier, O. Buttelli . "[*A Novel Tensor Tracking Algorithm For Block-Term Decomposition of Streaming Tensors*](https://ieeexplore.ieee.org/document/10208007)". **Proc.  IEEE SSP**, 2023. 
