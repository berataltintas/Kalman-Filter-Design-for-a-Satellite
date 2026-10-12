# Kalman Filter Design for a Satellite

This repository contains a MATLAB-based simulation and a multi-**Multiplicative Extended Kalman Filter (MEKF)** approach for satellite attitude estimation.

## Contents

Three different estimation flows work together in the project:

- **Coarse mode (Sun + Magnetometer):** Initial attitude via TRIAD and MEKF update from vector measurements.
- **Eclipse mode (Gyroscope + Magnetometer):** MEKF update from magnetic measurements with gyro-bias estimation when the sun sensor is unavailable (eclipse/shadow).
- **Fine mode (Star Tracker + Gyroscope):** MEKF with high-accuracy star tracker measurements and gyroscope bias estimation.

The simulation includes orbit propagation, ground-truth attitude generation, sensor measurement generation, and comparative visualization of filter results.

## Project Flow

The `main.m` file performs the following steps:

1. Sets orbit/sensor/filter parameters with `init.m`.
2. Generates coarse and fine datasets with `simulate_mode.m`.
3. Computes the coarse initial attitude with `coarse_triad.m`.
4. Runs filters with `coarse_mekf.m`, `eclipse_mekf.m`, and `mekf.m`.
5. Produces quaternion- and Euler (RPY)-based error metrics and plots.

## File Summary

- `main.m`: Simulation entry point, filter calls, metrics, and plots.
- `init.m`: Orbit, sensor noise, bias, and covariance parameters.
- `simulate_mode.m`: Data generation for coarse/fine modes.
- `mekf.m`: MEKF for fine mode (star tracker + gyro bias).
- `coarse_mekf.m`: Direction-vector-based MEKF for coarse mode.
- `eclipse_mekf.m`: Gyro + magnetometer-based MEKF under eclipse conditions.
- `coarse_triad.m`: TRIAD-based estimation of coarse initial attitude.
- `simulate_*.m`: Sensor measurement generators (sun sensor, magnetometer, gyro, star tracker).
- `propagate_orbit.m`, `attitude_truth.m`, `check_eclipse.m`: Orbit, ground-truth attitude, and eclipse modeling helpers.

## Requirements

- MATLAB (the repository runs with MATLAB scripts/functions)
- No additional toolbox dependency is strictly required by the code; core MATLAB linear algebra functions are used.

## Running

In the MATLAB command window at the repository root:

```matlab
main
```

## Generated Outputs

At the end of the run:

- Attitude comparison plots (True vs Estimate) for coarse, eclipse, and fine modes,
- Orbit plot,
- Sensor measurement plots,
- Quaternion component comparisons,
- RMSE and RMS/max error summaries in the command window.

## License

This project is licensed under the [MIT License](LICENSE).
