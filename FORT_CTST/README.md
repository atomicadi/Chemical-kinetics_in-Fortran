# Conventional Transition-State Theory in Fortran
# FORT: Fortran Operator for Rate Theory, A CTST rate calculator
<img width="1440" height="1080" alt="Image" src="https://github.com/user-attachments/assets/e2052a93-ca0e-4fd0-8428-3042c930489b" />

The theory of reaction rates that was published almost simultaneously by Henry Eyring, and by M. G. Evans and M. Polanyi in 1935 is referred to now as conventional transtition-state theory (CTST). This theory involves the same assumptions and approximations that are made in the calculations of equilibrium constants using statistical mechanics.\
So after merging the statistical mechanics and chemical equilibrium and to calculate the CTST rate constant, the partition functins of both the reactants and activated complex are needed. The total partition function of a molecule can be written as:
<p align="center">


$$
q = \sum_i g_i e^{-\frac{ε_i}{k_BT}} ...... (1)
$$


</p>

Where, the summation is taken over all energy lavels. The energy ε<sub>i</sub> is the energy of the i<sup>th</sup> state relative to the zero-point energy, and g<sub>i</sub> is the degeneracy. Usually, it is assume that the various types of energy- electronic, vibrational, rotational, and translational- are independent of one another. So, the total energy corresponding to the i<sup>th</sup> energy state is thus expressed as the sum of the different types of the energies:
<p align="center">


$$
ε_i = e_i + ν_i + r_i + t_i  ...... (2)
$$


</p>

Where, e<sub>i</sub> = electronic energy, ν<sub>i</sub> = vibrational energy, r<sub>i</sub> = rotational energy, and t<sub>i</sub> = translational energy. So, using equation (2), equation (1) can be written as:
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

Where, q<sub>e</sub> (electronic), q<sub>ν</sub> (vibrational), q<sub>r</sub> (rotational), and q<sub>t</sub> (translational) are separate partition functions, each referring to one type of energy. At ordinary temperatures the excited electronic levels of an atom or molecule are
usually too high to make a significant contribution to the partition function and equation (4) becomes,
<p align="center">


$$
q = q_{ν} q_{r} q_{t}  ...... (5)
$$


</p>

The formula of separate partition functions are:
<p align="center">


$$
q_{ν} (per-normal-mode) =  \frac{1}{(1 - e^{-\frac{hν}{k_BT}})}  ...... (6)
$$

$$
q_{r} (linear-molecule) =  \frac{8π^2Ik_BT}{σh^2}  ...... (7)
$$

$$
q_{r} (nonlinear-molecule) =  \frac{8π^2(8π^3I_AI_BI_C)^{1/2}(k_BT)^{3/2}}{σh^3}  ...... (8)
$$

$$
q_{t} =  \frac{(2πmk_BT)^{3/2}}{h^2}  ...... (9)
$$
</p>

Where, h = Planck constant, ν = normal-mode vibrational frequency, k<sub>B</sub> = Boltzmann constant, T = absolute temperature, I = moment of inertia for the linear molecule, (I<sub>A</sub>,I<sub>B</sub>, and I<sub>C</sub>) = moment of inertia for a nonlinear molecule about three axes at right angles, σ = symmetry number, m = mass of molecule

Now, if a reaction is
<p align="center">


$$
A + B => a^{‡} => C + D 
$$


</p> 
then using all the formulas, finally the CSTS rate constant equation becomes,
<p align="center">


$$
k = \frac{k_BT}{h} \frac{q_{‡}}{q_A q_B} e^{-\frac{E_0}{RT}}  ...... (10)
$$


</p>

Where, q<sub>‡</sub> is the total partition function of activated complex, q<sub>A</sub> and q<sub>B</sub> are the total partition function of reactant A and B respectively, E<sub>0</sub> = Barrier hight, and R = Gas constant.


In this project the CTST rate constant is calculated by solving equation (10) through Fortran programming. Herein 3 files (mod_1_requirement_scratch.f90, mod2_partition.f90, and mother.f90) are generated in which mod_1_requirement_scratch.f90 is a module file that conatains the requirement values and calculations to solve the partition functions (eg Total mass, moment of inertia etc.). For solving the moment of inerta, the moment of inerta tensor is generated first using the equations,
<p align="center">


$$
I_{xx} = \sum_im_i(y_i^2 + z_i^2)  ......(11)
$$


$$
I_{yy} = \sum_im_i(z_i^2 + x_i^2)   ......(12)
$$


$$
I_{zz} = \sum_im_i(x_i^2 + y_i^2)  ......(13)
$$

$$
I_{xy} = -\sum_im_ix_iy_i = I_{yx}  ......(14)
$$

$$
I_{xz} = -\sum_im_ix_iz_i = I_{zx}  ......(15)
$$

$$
I_{yz} = -\sum_im_iy_iz_i = I_{zy}  ......(16)
$$

</p>


Finally using equation (11)-(16), the moment of inerta tensor is generated,\
$\begin{bmatrix}a&b\\c&d
\end{bmatrix}$
Command line: ```./TST```
