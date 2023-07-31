# Optical feedback simulation

This Python code reproduces the figures of the article "Optical feedback linear cavity-enhanced absorption spectroscopy" by J.Tian et al, DOI: 10.1364/OE.431934. 
The code simulates the free-running frequency of the laser. Then it computes the real frequency jumps. The last part of the code shows the evolution of the optical feedback rate beta, the locking range and some other parameters.

## Code structure
In the beginning, you can find a few constants useful for the simulation. Then you can set two boolean variables _non-res_ and _experimental_. With the first variable, you can choose to add the non-resonant part of the equation into the calculations. The second variable is here to switch between the article's parameters and the experimental ones we had.

The function _coupling_condition()_ returns the value of $\tau_a$ and $\tau_c$, the roundtrip delay time between the laser and the cavity front mirror and the roundtrip delay time inside the cavity. It returns also the frequency $\nu_q$ corresponding to our wavelength of 1500 nm, and lastly the arrays of frequency and angular frequency $\nu$ and $\omega$.

The function _main()_ calculates the free running frequency with the resonant part and the non-resonant part with the equation below. It returns also the locking range and $\Delta \nu_{cav}$.

$$\begin{align}
\omega_{\text{free}}=\omega+ K_1~\frac{\sin \left[\omega\left(\tau_{\text{a}}+\tau_{\text{c}}\right)+\theta\right]-r_{\text{m}}^2 \sin \left[\omega~\tau_{\text{c}}+\theta\right]}{1+F^2 \sin ^2\left(\omega~\tau_{\text{c}} / 2\right)} - K_2~\sin \left(\omega~\tau_{\text{a}}+\theta\right)
\end{align}$$

The function _free_running_freq()_ recreats the plot of the article.

The function _beta_evolution()_ plots the evolution of the locking range as a function of $\beta$ the optical feedback rate for 3 finesse values.

The function _finesse_evolution_ plots the evolution of the locking range and $\Delta \nu_{cav}$ as a function  of the finesse.

The function _upsample()_ adds more points to known data with a new span between each point. It returns the two _x_upsamp_ and _y_upsamp_.

The function _freq_jumps()_ computes the real frequency jumps as we scan the laser with a ramp for example. It returns two lists _list_jumps_x_ and _list_jumps_y_ which are the concatenation of a small part of the free-running frequency and the jumps.

The function _transmission()_ takes the results of _freq_jumps()_ and feeds the Airy function (below) with it. It then plots the spectrum of the transmission of the cavity with and without optical feedback.

$$\begin{align}
A(\omega) = \frac{(1-R)^2}{(1-R)^2 + 4R\sin^2(\omega L/c)}
\end{align}$$

These are the main plots of the code :

<img src="https://github.com/MaloBriend/optical_feedback_simulation/blob/main/free_running_frequency.png" width="400" />

<img src="https://github.com/MaloBriend/optical_feedback_simulation/blob/main/frequency_jumps.png" width="400" />
