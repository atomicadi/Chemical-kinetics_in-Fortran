# Conventional Transition-State Theory in Fortran
# FORT: Fortran Operator for Rate Theory, A CTST rate calculator
<img width="1440" height="1080" alt="Image" src="https://github.com/user-attachments/assets/e2052a93-ca0e-4fd0-8428-3042c930489b" />

The theory of reaction rates that was published almost simultaneously by Henry Eyring, and by M. G. Evans and M. Polanyi in 1935 is referred to now as conventional transtition-state theory (CTST). This theory involves the same assumptions and approximations that are made in the calculations of equilibrium constants using statistical mechanics.\
So after margeing the statistical mechanics and chemical equilibrium and to calculate the CTST rate constant, the partition functin of both the reactants and activated complex is needed. The total partition function of a molecule can be written as:
<p align="center">


$$
q = \sum_i g_i e^{-\frac{ε_i}{k_BT}} ...... (1)
$$


</p>

Where the summation is taken over all energy lavels. The energy ε<sub>i</sub> is the energy of the i<sup>th</sup> state relative to the zero-point energy, and g<sub>i</sub> is the degeneracy. Usually, is is assume that the various types of energy- electronic, vibrational, rotational, and translational- are independent of one another. So, the total energy corresponding to the i<sup>th</sup> energy state is thus expressed as the sum of the energies of the different types:
<p align="center">


$$
ε_i = e_i + ν_i + r_i + t_i  ...... (2)
$$


</p>

Where e<sub>i</sub> = electronic energy, ν<sub>i</sub> = vibrational energy, r<sub>i</sub> = rotational energy, and t<sub>i</sub> = translational energy. So, using equation (2), equation (1) can be written as:
<p align="center">


$$
q = \sum_i g_{e_i} e^{-\frac{e_i}{k_BT}}  g_{ν_i} e^{-\frac{ν_i}{k_BT}}  g_{r_i} e^{-\frac{r_i}{k_BT}}  g_{t_i} e^{-\frac{t_i}{k_BT}}  ...... (3)
$$


</p>

Omitting the degeneracy, equation (3) may be written as:
<p align="center">


$$
q = q_{e} q_{ν} q_{r} q_{t}  ...... (4)
$$


</p>

Where q<sub>e</sub> (electronic), q<sub>ν</sub> (vibrational), q<sub>r</sub> (rotational), and q<sub>t</sub> (translational) are separate partition functions, each referring to one type of energy. At ordinary temperatures the excited electronic levels of an atom or molecule are
usually too high to make a significant contribution to the partition function and equation (4) becomes,
<p align="center">


$$
q = q_{ν} q_{r} q_{t}  ...... (5)
$$


</p>

The formula of separate partition functions are:
<p align="center">


$$
q_{ν} =  \frac{1}{(1 - e^{-\frac{hν}{k_BT}})}  ...... (6)
$$

$$
q_{r} =  \frac{8π^2Ik_BT}{h}  ...... (7)
$$

</p>
Command line: ```./TST```
