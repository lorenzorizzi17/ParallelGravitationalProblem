---
marp: true
theme: gaia
backgroundColor: #221
color: #eee
paginate: true
math: mathjax
---

<style>
/* 1. CONFIGURAZIONE GENERALE */
section {
    font-size: 25px; 
    text-align: left;      /* Testo normale a sinistra */
    padding-top: 50px;     /* Spazio in alto */
}

/* 2. CONFIGURAZIONE TITOLI */
h1, h2, h3 {
    text-align: left;    /* Titoli centrati per tutte le slide */
    width: 100%;
}

/* 3. CONFIGURAZIONE SPECIFICA PER LA COPERTINA (LEAD) */
section.lead {
    text-align: center;    /* Centra tutto il testo */
    justify-content: center; /* Centra verticalmente */
}
/* Nella copertina i titoli devono essere grandi */
section.lead h1 { font-size: 60px; margin-bottom: 0.2em;}
section.lead h2 { font-size: 40px; margin-bottom: 1em; color: #aaa;}

/* 4. UTILITY: COLONNE E BOX */
.columns {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 1rem;
  align-items: center; /* Centra verticalmente le colonne tra loro */
}

.box-rosso {
    background-color: #800000;
    color: white;
    padding: 20px;
    border-radius: 10px;
    border-left: 5px solid red;
    box-shadow: 5px 5px 10px rgba(0,0,0,0.5);
}

.box-codice {
    background-color: #333;
    border: 4px solid #555;
    padding: 15px;
}

pre, code { font-size: 18px; }
ul, ol { margin-left: 20px; }
</style>

# N-body gravitational problem
## Computational Physics Project


**Lorenzo Rizzi**
Modern computing for physics @ Physics of Data
January 2026

<br>

<div class="columns">
<div>

![center height:280px](fig/copertina.png)

</div>
<div>

![center height:280px](fig/copertina2.png)

</div>
</div>

---

## The physical problem

<div class="columns">
<div>

Fix $N$. The gravitational interaction at time $t$ reads:

$$\vec{F}_i = G \sum_{j \neq i}^N \frac{m_i m_j}{|\vec{r}_{ij}|^3} \vec{r}_{ij}$$

Requiring $O(N^2)$ computations per iteration.

To avoid numerical instabilities, one usually introduces a **softening** $\epsilon$:

$$\frac{1}{(r^2 + \epsilon^2)^{3/2}}$$

We will use $G = 1, \epsilon = 10^{-1}$

</div>
<div>

<br>

![width:90%](fig/Nbody.png)

</div>
</div>

---

## Index

The code was realized in three distinct implementations:

1. **Sequential CPU**: one-threaded *reference line*.
2. **Parallel CPU**: using OpenMP enabling multi-threading.
3. **Parallel GPU**: hardware acceleration with CUDA-C.

<br>

**In this presentation:**

* Energy conservation analysis (using either Euler or Verlet scheme)
* Physically relevant properties (radial distribution)
* Performance comparaison between implementations and benchmarking

---

## Energy conservation: integration schemes

<br>
<br>

<div class="columns">
<div>

### 1. Euler (Explicit)

The simplest approach, but **not symplectic**. Energy drifts with time (explodes), spurious injection!

$$
\begin{aligned}
\vec{x}_{t+1} &= \vec{x}_t + \vec{v}_t \Delta t \\
\vec{v}_{t+1} &= \vec{v}_t + \vec{a}(\vec{x}_t) \Delta t
\end{aligned}
$$

<br>

**Bad for orbital mechanics**:
The potential is not well behaved.
Particles spiral out due to energy gain.

</div>
<div>

### 2. Velocity Verlet

**Symplectic** integrator. Reversible in time, conserves energy for long periods.

$$
\begin{aligned}
\vec{v}_{t+\frac{1}{2}} &= \vec{v}_t + \frac{1}{2}\vec{a}_t \Delta t \\
\vec{x}_{t+1} &= \vec{x}_t + \vec{v}_{t+\frac{1}{2}} \Delta t \\
\vec{a}_{t+1} &= F(\vec{x}_{t+1}) / m \\
\vec{v}_{t+1} &= \vec{v}_{t+\frac{1}{2}} + \frac{1}{2}\vec{a}_{t+1} \Delta t
\end{aligned}
$$

**Standard for N-Body**:  Stable orbits over time

</div>
</div>

---

## Energy conservation

![height:500px](fig/energyConservationAll.png)

---

## Energy conservation

![height:500px](fig/energyConservationVerlet.png)

---

## Sanity check: two-bodies orbital dynamics
<style scoped>
/* Trasforma l'immagine in un blocco e imposta i margini automatici ai lati */
img {
    display: block;
    margin: 0 auto;
}
</style>
<br>

<div class="columns">
<div>

A special initial condition, where $N = 2$:

$m_1 \gg m_2$

$\vec v_1 \approx 0,\>\>  \vec v_2 = v \hat u_z$

$\vec r_1 = 0, \>\> \vec r_2 = (L/2, 0, 0)$

A symplectic integrator (Verlet) should always oscillate aroung the analytical trajectory!

</div>
<div>


![height:220px](fig/2bodyTrajectoriesVerlet.png)



![height:220px](fig/2bodyTrajectoriesEuler.png)
</div>
<div>


---

## Virial theorem
<style scoped>
img {
    display: block;
    margin: 0 auto;
}
</style>
If the total enery of the system is negative, then it should converge to a stationary state where the virial theorem holds true. In particular, at **equilibrium**, one expects to see:
$$
2 \langle K \rangle = - \langle U \rangle
$$
Given the total energy constraint $\langle K \rangle + \langle U \rangle = E$, then, if we wait long enough, we should measure: $
K \approx - E, 
U \approx 2E
$

![center height:300px](fig/virial.png)

---

## Radial distribution function

<br>

<style scoped>
/* Trasforma l'immagine in un blocco e imposta i margini automatici ai lati */
img {
    display: block;
    margin: 0 auto;
}
</style>
<br>

<div class="columns">
<div>

When the initial condition is random (provided $E < 0$), all the bodies will collaps and form a dense cluster

With OVITO, we can extract the frame-by-frame radial distribution function $g(r)$:



Theoretical models (isothermal sphere) predicts a scaling behavior $\rho(r) \sim g(r) \sim r^{-2}$
</div>
<div>


![height:370px](fig/radialDistributionFunction.png)

</div>
<div>

---

## Code structure (CPU sequential)

The entire code is written in C++. A main class is defined `Simulation`:

```cpp
class Simulation {
    private:
    int m_N;        
    Real m_dt;          
    Real m_L;            
    // Fundamental data structures
    std::vector<Real> m_x, m_y, m_z;      
    std::vector<Real> m_vx, m_vy, m_vz;    
    std::vector<Real> m_fx, m_fy, m_fz;    
    std::vector<Real> m_mass;  

    public:
    // Constructor (random masses, velocities, positions); implememented in simulation.cpp
    Simulation(int n_particles, Real timeStep, Real length);
    Simulation(Real timeStep, Real length);
    ...

    // Crucial part, computes the forces (implemented in computeForces.cpp)
    Real computeForces(bool);

    // The numerical integrators (choose which to use in main; implemented in numericalIntegrator.cpp)
    void integrateEuler(int nSteps, int saveEvery, std::string saveEnergy, std::string saveTrajectory);
    void integrateVerlet(int nSteps, int saveEvery, std::string saveEnergy, std::string saveTrajectory);
};
```
The sequential logic flow is quite easy: while $t < T_{max}$, call `computeForce()` and run the numerical integration. Energy is computed along with the forces

---

## Code structure: OpenMP