# Disturbance Adaptive Rigid Tube Model Predictive Control on a quadrotor

## Installation of Required Packages

To run this code, you need to install the following packages:

1. **CasADi**: [https://web.casadi.org/](https://web.casadi.org/)  
2. **MPT3**: [https://www.mpt3.org/Main/Installation](https://www.mpt3.org/Main/Installation)  

> Note: Installing MPT3 should automatically install **YALMIP**, which is also required to run the code.
> On MacBooks with Apple Silicon, the code has only been tested on **MATLAB 2021b or earlier**.


## you should know
**drone:**for every flight we use 'InitializeDrone' to get struct 'drone'.
**Controller:** We use class UQRobustMPC to construct the tube MPC controller.
**Figure:** you can use figure_disturbances_flag = 1;
figure_states_flag ;
figure_inputs_flag ;
figure_path_flag = ;

# Quickly start
After making sure you have installed all the needed packages, you can directly run them in order:
1. **Case1_Conservative:** A basic LQR controller with input constraints truncated, flies in a dynamic wind field to get the conservative estimation of the disturbance of the system. You can change the parameter 'F_{wind}（the magnitude of the wind disturbance with respect to time t）'.
2. **Case1_compare:** fly to track a predefined Lissajous curve. by change  
``` Change the curve 
% Ax = 4; Ay = 8; Az = 2;    %  Amplitude
% p = 2; q = 1; r = 2;        % Frequency ratio

Ax = 4; Ay = 6; Az =8;     % Amptitude
p = 2; q = 1; r = 2;        % Frequency ratio
```
Tracking the different curves and outcomes will show different realizations of disturbances in the system.

3. **Case1_TMPC:** Comparison of LQR Controller and Rigid Tube Controller. 
4. **Case1_DATMPC:** comparison of LQR Controller and DATMPC. DATMPC will learn the new estimation of the system's disturbance and adaptively update the disturbance invariant set ('drone.W').



