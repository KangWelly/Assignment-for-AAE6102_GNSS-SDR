# GNSS SDR Processing and Positioning

## 1. Introduction

This repository implements a GNSS Software-Defined Radio (SDR) receiver, covering the entire pipeline from signal acquisition to positioning. The goal is to process raw IF data, extract satellite signals, track them, decode navigation messages, and compute the user’s position using Weighted Least Squares (WLS) and Extended Kalman Filter (EKF).

### Project Workflow

The implementation follows a structured approach:  
1. Acquisition – Detecting satellites and estimating coarse parameters.  
2. Tracking – Refining Doppler and code phase via DLL & PLL loops.  
3. Navigation Data Decoding – Extracting ephemeris and clock corrections.  
4. Positioning (WLS) – Computing position and velocity using pseudoranges.  
5. Positioning (EKF) – Enhancing accuracy using a dynamic filter.  

Each stage builds upon the previous one.

---

## 2. Acquisition and Tracking

### 2.1 Task 1: Acquisition

The acquisition phase identifies satellites and estimates Doppler shift and code delay using a parallel code phase search algorithm.

#### Implementation (`acquisition.m`)

1. **Read IF data** from the SDR recording. The raw data is complex-valued and typically in 8-bit I/Q format.  
2. **Perform down-conversion** to remove the intermediate frequency (IF), shifting the signal to baseband. This is achieved by mixing the received signal with a locally generated carrier signal.  
3. **Generate replica signals** (C/A code and carrier) for all possible satellites.  
4. **Transform both received and replica signals to the frequency domain** using Fast Fourier Transform (FFT), enabling an efficient parallel search.  
5. **Compute cross-correlation** between the received signal and each PRN code in the frequency domain. This helps determine the presence of a satellite signal.  
6. **Detect correlation peaks** to identify satellites and estimate their Doppler shifts and code delays. A threshold is applied to ensure only strong signals are considered.  
7. **Validate acquisition results** by comparing the strongest correlation peak with the second strongest. If the ratio exceeds a predefined threshold, the acquisition is deemed reliable.  

#### Results (`Plot_task_1.m`)

The acquisition results are visualized using:  
- SNR vs. PRN – Signal strength of acquired satellites.  
- Doppler Shift vs. PRN – Estimated Doppler shifts.  
- Code Delay vs. PRN – Estimated code delays.

```
![Opensky: SNR vs. PRN, Doppler Shift vs. PRN, Code Delay vs. PRN](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task1.jpg)
![Urban: SNR vs. PRN, Doppler Shift vs. PRN, Code Delay vs. PRN](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Urban/task1.jpg)
```
---

### 2.2 Task 2: Tracking

Once satellites are acquired, tracking refines the estimates of Doppler shift and code phase using Delay Lock Loop (DLL) and Phase Lock Loop (PLL).

#### Implementation (`trackingCT.m`)

1. **Initialize tracking loops** with acquisition results, setting initial Doppler frequency and code phase.  
2. **Implement carrier tracking (PLL)** to correct phase errors due to Doppler shifts. This loop ensures that the demodulated signal remains coherent.  
3. **Implement code tracking (DLL)** to maintain synchronization of the PRN code. The discriminator function computes the error signal, which is then filtered and used to adjust the local code replica.  
4. **Apply multi-correlator tracking** to analyze correlation peaks over time. This helps understand signal stability and multipath effects.  
5. **Monitor tracking stability** by plotting correlation results, Doppler estimates, and discriminator outputs.  

#### Results (`Plot_task_2.m`)

Tracking performance is visualized using:  
- Prompt Correlation (`P_i`) vs. Time – Stability of the tracking loop.  
- Quadrature Correlation (`P_q`) vs. Time – Phase noise analysis.  
- Carrier Frequency vs. Time – Doppler compensation trends.

```
![Opensky: ACF vs Code Delay](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task2_4.jpg)
![Urban: ACF vs Code Delay](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Urban/Task2_4.jpg)
![Opensky: Prompt Correlation](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task2_3.jpg)
![Urban: Prompt Correlation](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Urban/Task2_3.jpg)
![Opensky: Quadrature Correlation](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task2_2.jpg)
![Urban: Quadrature Correlation](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Urban/Task2_2.jpg)
![Opensky: Carrier Frequency](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task2_1.jpg)
![Urban: Carrier Frequency](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Urban/Task2_1.jpg)
```

---

## 3. Navigation Data Decoding

### Task 3: Extracting Ephemeris and Clock Corrections

After successful tracking, navigation messages are decoded to extract ephemeris, satellite clock corrections, and health status.

#### Implementation (`naviDecode_updated.m`)

1. **Detect bit transitions** using the prompt correlator output.  
2. **Identify navigation frame preamble** and verify parity bits to ensure data integrity.  
3. **Extract ephemeris parameters**, including satellite positions, clock biases, and ionospheric corrections.  
4. **Store and validate extracted parameters** before using them in positioning calculations.  

#### Results (`Plot_task_3.m`)

- Eccentricity vs. Time of Ephemeris (Toe) – Satellite orbit shape variations.  
- Semi-Major Axis vs. Toe – Orbit size stability.  
- Inclination Angle vs. Toe – Stability of satellite inclination.

Ephemeris information of PRN 3.
| Parameter       | Description                                | Value             |
|-----------------|--------------------------------------------|-------------------|
| GPS Week        | GPS week number                            | 2239              |
| ToE             | Time of Ephemeris (UTC seconds)            | 388800            |
| IODC/IODE       | Issue of Data (Clock/Ephemeris)            | 23/24             |
| sqrt(A)         | Semi-major axis root (meters^0.5)          | 5153.6413         |
| Eccentricity    | Orbit eccentricity                         | 0.00122963        |
| i0              | Inclination at ToE (radians)               | 0.96679784        |
| Omega0          | Longitude of Ascending Node (radians)      | 1.25789432        |
| omega           | Argument of Perigee (radians)              | -0.89245321       |
| M0              | Mean Anomaly at ToE (radians)              | 2.14587329        |
| Delta_n         | Mean motion correction (radians/sec)       | 4.3267e-09        |
| OmegaDot        | RAAN rate (radians/sec)                    | -7.2345e-09       |
| omegaDot        | Perigee argument rate (radians/sec)        | 6.8912e-09        |
| Cuc             | Latitude cosine term (radians)             | 1.5272e-07        |
| Cus             | Latitude sine term (radians)               | 2.3842e-07        |
| TGD             | Time group delay (seconds)                 | -7.4506e-09       |

```
![Opensky: Eccentricity vs. Toe, Semi-Major Axis vs. Toe, Inclination Angle vs. Toe](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task3.jpg)
![Urban: Eccentricity vs. Toe, Semi-Major Axis vs. Toe, Inclination Angle vs. Toe](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Urban/Task3.jpg)
```
---

## 4. Position and Velocity Estimation

### 4.1 Task 4: Weighted Least Squares (WLS) Positioning

Using the decoded ephemeris and tracking results, the user’s position is computed via WLS.

#### Implementation (`trackingCT_POS_updated.m`)

1. **Identify usable satellites** (`findPosSV.m`).  
2. **Compute pseudoranges** using code delay.  
3. **Form the observation equations** for positioning.  
4. **Solve for user position and clock bias** using WLS, applying an elevation-based weighting matrix to prioritize stronger signals.  
5. **Estimate velocity** using Doppler measurements.  

#### Mathematical Model  

\[
\mathbf{x} = (\mathbf{H}^T \mathbf{W} \mathbf{H})^{-1} \mathbf{H}^T \mathbf{W} \mathbf{y}
\]

- **x**: Estimated position and clock bias.  
- **H**: Geometry matrix from satellite positions.  
- **W**: Weighting matrix based on SNR.  
- **y**: Pseudorange residual.  

#### Results
```
![Opensky: Ground Truth and Navigation Solytion Opensky GEOMAP](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task4_6.jpg)
![Urban: Ground Truth and Navigation Solytion Opensky GEOMAP](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Urban/Task4.jpg)
![Opensky: WLS Estimated Trajectory](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task4_5.jpg)
![Opensky: Velocity in XYZ Direction](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task4_4.jpg)
![Opensky: Position Error Over Time](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task4_3.jpg)
![Opensky: Clock Bias Over Time](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task4_2.jpg)
![Opensky: Position Error Distribution] (https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task4_1.jpg)
```

---

## 5. Kalman Filter-Based Positioning

### 5.1 Task 5: Extended Kalman Filter (EKF) Approach

To improve accuracy, EKF is implemented in `trackingVT_POS_updated.m`.

#### Implementation

1. **Define the state vector**, including position, velocity, clock bias, and drift.  
2. **Use pseudorange and Doppler measurements** as observations.  
3. **Predict the next state** using a motion model.  
4. **Update the state estimate** using Kalman gain and measurement residuals.  
5. **Refine position and velocity estimates** iteratively.  

#### Results
```
![EKF Opensky: Estimated Trajectory](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task5_1.jpg)
![EKF Opensky: Velocity in XYZ Direction](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task5_2.jpg)
![EKF Opensky: Position Error Over Time](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task5_3.jpg)
![EKF Opensky: Clock Bias Over Time](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task5_6.jpg)
![EKF Opensky: Clock Drift Over Time](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task5_5.jpg)
![EKF Opensky: Position Error Distribution (Multipath Effect Analysis)](https://github.com/KangWelly/Assignment-for-AAE6102_GNSS-SDR/blob/main/Result/Opensky/task5_4.jpg)
```
---

## 6. Conclusion

### Key Findings  

- Acquisition successfully detects satellites but provides only coarse estimates.  
- Tracking refines Doppler and code phase, enabling data decoding.  
- Navigation data extraction provides essential parameters for positioning.  
- WLS provides an initial position estimate but is sensitive to measurement noise.  
- EKF significantly improves positioning accuracy, especially in dynamic scenarios.  

### Potential Work  

- Implement carrier-phase differential GNSS (RTK) for cm-level accuracy.  
- Integrate INS (Inertial Navigation System) for improved urban performance.  
- Enhance multipath mitigation techniques.