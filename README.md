% Extension to Levin 1997 simulation

This is an implementation and extension of the simulations described in Levin *et al*., "The population genetics of antibiotic resistance", *Clinical Infectious Diseases* (1997). The paper is [here](https://academic.oup.com/cid/article/24/Supplement_1/S9/283564/The-Population-Genetics-of-Antibiotic-Resistance).

# The original simulation

Sensitive and resistant strains of a commensal organism live in $N$ hosts and a
shared environment, replicating with generation time $t_g$ for some number of
generations (which defines the duration of the simulation). Every host gets an
average of $T$ antibiotic treatments per unit time (e.g., per year), so that,
at each generation, each host is treated with probability $p = t_g T$.
(Clearly, $T$ cannot exceed the number of generations per year.)

During each generation, three things happen:

- Competition between the two strains occurs in each host and in the environment.
- Simultaneously, a fraction $g$ of each host's strains are replaced by a representative mix of environmental strains, and a fraction $f$ of the environment's strains are replaced by a representative mix of all the hosts' strains.
- Then, the resistant strain comes to complete domination in treated hosts.

During competition, the sensitive strain has fitness $w_S = 1$ and the
resistant strain has fitness $w_R = 1-s$. The initial fraction resistant in the
host and environment are set to a small number.

# Extensions

I made two extensions to the simulation.

First, I created diversity in the hosts. A fraction $f_c$ of the hosts are
consumers, and the remainder never get treated. The average treatment rate $T$
is now the consumption among consumers. Setting $f_c = 1$ recapitulates the
original model.

Second, I allow non-complete selection of resistant bacteria. Thus, there are
now four selection coefficients. The original simulation set the fitness of the
sensitive strain to zero during treatment.

Context      Sensitive Resistant
-----------  --------- ---------
No treatment $1$       $1-s_R$
Treatment    $1-s_S$   $1-s_R$

Table: Fitnesses of the two strains in the two contexts. Setting $s_S = 0$ recapitulates the original simulation.

# Results

First, I recapitulate the results from Figure 2 in the original manuscript,
which shows that increasing treatment rate $T$ leads to increasing resistance.
The slope of this relationship is highest when the fitness cost $s$ is
smallest. The curve is always concave down or linear.

Second, I show that the dependence of resistance on $f_c$ is strongest when $T$
is higher, and similarly the dependence of resistance on $T$ is higher when
$f_c$ is higher. This is no surprise: the average consumption over all consumers
is $f_c T$.

Third, I show that the consumption/resistance curve becomes sigmoidal when the
fitness costs of antibiotic treatment for the susceptible strain is smaller. In
other words, increasing consumption among consumers gives you increasing gains
in resistance when multiple treatments are needed to push a host to full
resistance. When just one dose will push a host to full resistance, a second
dose has no effect, which leads to the concave-down part of the curve. This
concavity is only observable when the fitness cost of resistance is small.
Otherwise, the resistant strain dies out too quickly for the successive
treatments to build on one another.

# Addendum: Technical errata

I think there is an error in the equations on page S10 in the original manuscript. There it says that:
$$
\begin{gathered}
\Delta E = f \left( \sum_{i=1}^N \frac{H_i}{N} - E \right) - \frac{sE(1-E)}{1-sE} \\
\Delta H_i = 1 - H_i \text{ with probability $p$} \\
\Delta H_i = g(E - H_i) - \frac{s H_i (1 - H_i)}{1 - s H_i} \text{ with probability $1-p$},
\end{gathered}
$$
where I have taken the liberty of writing $T/219$ as $p$. Here $E$ is the
fraction of population in the environment that is resistant; $H_i$ is the
fraction resistant in host $i$. I think there are missing terms in the first and
last equations.

For example, imagine the mathematically correct but implausible scenario where
$f>0$, $s=1$, and the average host resistance $\overline{H} \equiv \sum_i (H/N) = 0$.
Then the first equation reduces to $\Delta E = -(1+f)E$, i.e., the new
environmental resistance will be $E + \Delta E < 0$.

I derive a correction. The new environmental resistant fraction is the result
of replacing a fraction $f$ with the average host resistance $\overline{H}$ and
enforcing selection in the remainder:
$$
\begin{aligned}
E' &= f \overline{H} + (1-f) \frac{ w_R E }{w_S (1 - E) + w_R E} \\
   &= f \overline{H} + (1-f) \frac{(1-s)E}{(1-E) + (1-s)E} \\
   &= f \overline{H} + (1-f) \frac{(1-s)E}{1 - sE}
\end{aligned}
$$
Thus, the change in resistance is:
$$
\begin{aligned}
\Delta E &= E' - E \\
 &= f \overline{H} + (1-f) \frac{(1-s)E}{1-sE} - \left[ \frac{1-sE}{1-sE} (1-f) + f \right] E \\
 &= f (\overline{H} - E) + (1-f) \frac{sE(1-E)}{1-sE},
\end{aligned}
$$
which makes it look like there is a missing $(1-f)$ in the paper's equation. Similarly, I would expect:
$$
\begin{gathered}
H_i' = gE + (1-g) \frac{(1-s)H_i}{1 - sH_i} \\
\Delta H_i = g(E - H_i) + (1-g) \frac{s H_i (1-H_i)}{1 - sH_i}
\end{gathered}
$$
