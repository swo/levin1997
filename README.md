# Levin et al. simulation

This is an implementation of the simulations described in Levin *et al*., "The population genetics of antibiotic resistance", *Clinical Infectious Diseases* (1997). The paper is [here](https://academic.oup.com/cid/article/24/Supplement_1/S9/283564/The-Population-Genetics-of-Antibiotic-Resistance).

# Derivation of pg. S10 equations

The new environmental resistant fraction is the result of replacing a fraction $f$ with the average host resistance $\overline{H}$ and enforcing selection in the remainder:
$$
\begin{aligned}
E' &= f \overline{H} + (1-f) \frac{(1-s)E}{(1-E) + (1-s)E} \\
   &= f \overline{H} + (1-f) \frac{(1-s)E}{1 - sE}
\end{aligned}
$$
Thus, the *change* in resistance is:
$$
\begin{aligned}
\Delta E &= E' - E \\
 &= f \overline{H} + (1-f) \frac{(1-s)E}{1-sE} - \left[ \frac{1-sE}{1-sE} (1-f) + f \right] E \\
 &= f (\overline{H} - E) + (1-f) \frac{sE(1-E)}{1-sE},
\end{aligned}
$$
which makes it look like there is a missing term in the equation.

Similarly, I would expect:
$$
H_i' = gE + (1-g) \frac{(1-s)H_i}{1 - sH_i}
$$
