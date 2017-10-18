# it-ift

An R package for likelihood estimation of time-homogeneous jump diffusion processes.

$$dX_t = \mu(\theta, X_t)dt + \sigma(\theta, X_t)dW_t + c(\theta,\eta_t)dN_t$$

The R package builds a [TMB](https://github.com/kaskr/adcomp) program, depending on several header files under /includes, that via the [exact saddlepoint approximation](url_to_come_to_arkivx_manuscript) applied Itô-Taylor approximated characteristic functions, returns the negative log-likelihood for a jump-diffusion process.
The package is currently in development.

# Usage
Create ```DATA``` and ```PARAMETER``` lists, following TMB conventions.

* The ```DATA$X``` vector contains the descrete observations of the jump diffusion
* The ```DATA$process```integer specifies the integer valued ```process```. See lists of processes and corresponding ```process``` values.
* The ```DATA$dt``` scalar specifies the equidistant timesteps.
* The ```PARAMETER$par``` object corresponds to the ```DATA$process``` value.

# Advanced notes
Information to come regarding:

* User specified jump diffusion process
* Tuning of number of Newton iterations and quadrature points.
