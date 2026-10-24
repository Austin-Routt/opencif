# Ladder-Array Design Specification — v5 (revision 5.1, 2026-10-23)

A closed system of equations and an algorithm mapping

> **(shape function, three targets, fabrication floors) → (a drawn planar geometry)**

together with the validators that decide whether an implementation is correct.

---

## 0. Contract

**Produces.** A row-indexed planar geometry: one primary duct, two secondary ducts, a periodic sequence of transverse links; every width, every drawn path length, and the row count $N$ determined.

**Guarantees — provable inside this document.**

| quantity | status |
|---|---|
| resistance closure $r(w_s(i))=\rho(i)r(w_c(i))$ | **exact** to $\varepsilon_{\rm inv}$ |
| $w\!\leftrightarrow\!h$ symmetry of $r$ | **exact series identity** $B(1/\alpha)=\alpha^{-2}B(\alpha)$ (§3) |
| mass identity $P(i)+2S(i)=1+2\sigma_0$ | **exact for any $\sigma_0$** (§10) |
| terminal fraction: the schedule's $F_{\rm sched}=1-P(N)$, and the delivered $F=2/(\rho(N)+2)$ | **exact for any $\sigma_0$** (§10; rev. 5.1) |
| forward↔inverse round-trip | **exact, no special cases** (§11) |
| lamina $\delta(i)=\delta^\ast$ under S1 | **exact at any aspect ratio** (§7) |
| $F$ reproduced by streamline migration, **zero diffusion** | **exact** (§6.2 lemma) |
| depth factor $K_{\rm depth}(\beta)$ | **derived** from the same series as $r$ (§4) |
| canonical taper ($N$, $F$, $w_c(i)$, inversion) | **exact**, one monotone potential $\Psi$ (§13) |
| $\Lambda$-, GS-, $\mu$-, $Q_{\rm in}$-invariance | **exact** (§17) |

**Does not guarantee.** $\kappa$, $\beta_{\rm cal}$, $\phi_{\rm cal}$, $u_{\rm lo}$, $u_{\rm ff}$, $\eta$, $f_{\rm soft}$, $\gamma_{\rm GS}$, $\mathrm{AR}_{\max}$, $R_{\min}$, $\ell_g$, $c_{\rm ent}$, $L_{\rm eq}$. These are **supplied**. §20 lists them with the exact consequence of each being wrong. **No validator in §19 can check any of them.**

---

## 1. Objects and conventions

| symbol | meaning |
|---|---|
| $i\in\{0,\dots,N\}$ | row index; $N$ is an **output** |
| $h$ | fixed second dimension of every duct |
| $\mu,\ \rho_f,\ T$ | viscosity (scales $r$); density (Re); temperature (§20 Péclet only) |
| $w_c(i),\ w_s(i)$ | primary width; secondary **equivalent straight width** after row $i$ |
| $\mathcal W:[0,1]\to\mathbb R^+$ | shape function, $w_c(i)=\mathcal W(i/N)$ |
| $\beta(i)=h/w_c(i)$ | **local aspect ratio** |
| $f(i)\in(0,\tfrac12)$ | **per-side** extraction schedule |
| $\sigma(i)=Q_s(i)/Q_c(i)$ | march variable; $\rho(i)=1/\sigma(i)=R_s(i)/R_c(i)$ |
| $\delta(i)$ | flat-lamina width (§4); $D=\kappa\delta$ = derived scale |
| $\Psi(w)=w^2/c_2(h/w)$ | **taper potential** (§13) — strictly increasing |
| $F$ | terminal fraction |
| $\Lambda(i),\ \Phi(i),\ \hat L(i)$ | row pitch; fold factor; **measured** drawn traversal length |
| $\ell_{\rm stroke}$ | fold stroke length (sets $n_{\rm turns}$, §14.6) |
| $G,\ \ell_g,\ g_p$ | link width; link length; post width |
| $\phi=G/\Lambda$ | **wall porosity** — the open fraction of the "wall" the kernels assume is solid |
| $\mathrm{GS},\ W,\ F_{\min},\ R_{\min}$ | fold width; transverse budget; minimum feature; minimum inner-corner radius |
| $R_{\rm out}^{c},R_{\rm out}^{s}$ | run-out resistances, last link → terminals (**per side** for $s$) |
| $m,\ m_{\rm budget}$ | terminal-layer depth; the allowance asserted against it |
| $Q_c,Q_s,Q_{\rm gap}$ | flows (scale-free — §17) |

**Flow conventions, fixed once.**
$$\boxed{\;Q_{\rm gap}(i)=f(i)Q_c(i-1),\qquad Q_c(i)=Q_c(i-1)\bigl[1-2f(i)\bigr],\qquad Q_s(i)=Q_s(i-1)+Q_{\rm gap}(i)\;}$$

$f$ is **per side**; the primary loses $2f$ per row. $R_c(i),R_s(i)$ are the segments **after** row $i$, carrying $Q_c(i),Q_s(i)$.

> **Segment rule (settles turns, §14.5).** A segment is the span between two consecutive rungs. **Flows change only *at* rungs.** So any gap-free span — including a 180° turn — carries $Q_c(i),Q_s(i)$ **exactly** and takes the ratio $\rho(i)$ of its **upstream** rung. No interpolation. If a span contains rungs, it is not one segment.

---

## 2. Constants

### 2.1 Derived — exact, verified in §19

| symbol | value | where |
|---|---|---|
| $192/\pi^5$ | $0.6274106$ | §3 |
| $\sum_{n\,\rm odd}n^{-5}=\tfrac{31}{32}\zeta(5)$ | $1.0045238$ | §3 |
| $\sum_{n\,\rm odd}n^{-4}=\pi^4/96$ | $1.0146780$ | §4 |
| $\sum_{n\,\rm odd}n^{-2}=\pi^2/8$ | $1.2337006$ | §4 |
| $B(1)$ | $0.4217310$ | §3 (square) |
| $12/B(1)$ | $28.454154$ | §3 (square) |
| $24/B(1)$ | $\mathbf{56.908}$ | §3 — **square-duct $f\!\cdot\!\mathrm{Re}$; textbook 56.91** |
| $f\!\cdot\!\mathrm{Re}$, parallel plates | $\mathbf{96}$ | §3 — **textbook 96** |
| $\lim_{\alpha\to\infty}\alpha[1-B(\alpha)]$ | $0.6302489$ | §3 (slender) |
| $B(\alpha)$ range for $\alpha\ge1$ | $[0.4217310,\,1)$ | §3 — sorting eliminates cancellation |
| **reciprocity** $B(1/\alpha)$ | $\alpha^{-2}B(\alpha)$ (**exact**) | §3 — *why* sorting is the exact branch |
| $c_2(\infty)$ | $3$ (exact) | §4 (cubic limit) |
| $c_2(1)$ | $3.5568$ | §4 |
| $\lim_{\beta\to0}\beta\,c_2(\beta)$ | $\approx1.7$ | §4 — $c_2$ diverges as $1/\beta$ |
| $F_{\rm lam}(\tfrac12,\beta)$ | $\tfrac12$ (exact, **all** $\beta$) | §4 (hard bound) |
| $C(1)$ | $9.6078$ | §16-b (square) |
| $C(\infty)$ | $6$ (exact) | §16-b (slender) |

**$c_2(\beta)$ and $K_{\rm depth}(\beta)=\sqrt{3/c_2(\beta)}$ — derived:**

| $\beta$ | 0.10 | 0.45 | 0.53 | 0.64 | 0.80 | 1.00 | 1.25 | 1.40 | 1.87 | 2.00 | 2.69 | 4.00 | $\infty$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $c_2$ | 17.378 | 5.038 | 4.577 | 4.163 | 3.805 | 3.557 | 3.388 | 3.324 | 3.210 | 3.190 | 3.127 | 3.078 | **3.000** |
| $K_{\rm depth}$ | 0.416 | 0.772 | 0.810 | 0.849 | 0.888 | 0.918 | 0.941 | 0.950 | 0.967 | 0.970 | 0.979 | 0.987 | **1.000** |

### 2.2 Error laws (not constants)

| law | where |
|---|---|
| cubic seed error $\;|x_0/x-1|\simeq\sqrt{f/27}$ | §4 E-f — *an $f$-dependent law; do not table it as a constant* |
| terminal-layer decay $\;|e(i)|=A\,z_-^{\,N-i}$ | §15 — *$A$ is the boundary excitation, not 1* |

### 2.3 Calibrated — supplied; **unverifiable here** (§20)

| symbol | value | role |
|---|---|---|
| $\kappa$ | $4.9138$, 95 % CI $[4.750,5.078]$ | $D=\kappa\delta$. **Fitted with the cubic kernel.** §5. |
| $\beta_{\rm cal}$ | **unknown** | aspect ratio(s) at which $\kappa$ was fitted. **Controls a 27 % swing in $N$.** §5. |
| $\phi_{\rm cal}$ | **unrecorded** | wall porosity at which $\kappa$ was fitted. §20. |
| $u_{\rm lo},\ u_{\rm ff}$ | $0.035$, $0.012$ | edges of the $\kappa$-constant window |
| $\eta$ | $0.903$ | realized/model **terminal fraction**, $F_{\rm real}=\eta F$. **A yield correction — not a lamina correction.** |
| $f_{\rm soft}$ | $5\times10^{-4}$ | soft upper bound on $f$ |
| $\gamma_{\rm GS}$ | $[1.4,1.5]$ | admissible $\mathrm{GS}/G$ |
| $\mathrm{AR}_{\max}$ | $10$ | $w\ge h/\mathrm{AR}_{\max}$ |
| $R_{\min}$ | process | minimum printable inner-corner radius |
| $\ell_g,\ c_{\rm ent}$ | process | $R_{\rm gap}=r(G)[\ell_g+c_{\rm ent}D_h(G,h)]$ |
| $L_{\rm eq}$ | process | Stokes-corner equivalent length per turn (§14.6) |

### 2.4 Tolerances — explicit inputs; no literal in code

$$\varepsilon_r,\ \varepsilon_{\rm inv},\ \varepsilon_{\rm cub},\ \varepsilon_{\rm id},\ \varepsilon_{\rm geo},\ \varepsilon_{\rm net},\ \varepsilon_{\rm lam},\ \varepsilon_{\rm taper},\ \mathrm{Re}_{\max}$$

---

## 3. Resistance kernel

$$\boxed{\;r(w)=\frac{12\mu}{a\,b^{3}\,B(\alpha)},\qquad a=\max(w,h),\ b=\min(w,h),\ \alpha=\frac ab\ (\ge1)\;}$$
$$\boxed{\;B(\alpha)=1-\frac{192}{\pi^{5}\alpha}\sum_{n=1,3,5,\dots}^{n_{\max}}\frac{1}{n^{5}}\tanh\frac{n\pi\alpha}{2}\;}$$

**R-a.** $B$ depends on $\alpha$ alone. Tabulate once; **share the call site with §4 and §16-b.**

**R-b — the reciprocity identity (why sorting is not merely faster, but the *exact* branch).**

$$\boxed{\;B(1/\alpha)=\alpha^{-2}B(\alpha)\;}\qquad\text{(exact; verified to }7\times10^{-13}\text{ at }\alpha=1.4)$$

This *is* the $w\!\leftrightarrow\!h$ symmetry of $r$: substituting it into $r$ with the slots exchanged returns the same value identically. Two consequences:

1. **$B(1/\alpha)\sim\alpha^{-2}$ is small.** Evaluating the unsorted branch computes $1-(\text{something}\to1)$, losing $\approx2\log_{10}\alpha$ digits to cancellation, *and* amplifying the series-truncation error by $1/B(1/\alpha)\approx\alpha^{2}$.
2. **For $\alpha\ge1$, $B$ is monotone increasing with range $[0.4217310,\,1)$** — bounded away from zero. **Sorting eliminates the cancellation entirely.**

**Measured unsorted relative error, at $n_{\max}=$ 50 *odd* terms ($n\le99$):**

| $\alpha$ | 1.4 | 2.7 | 5.0 | 11.7 | 23.3 | 46.7 |
|---|---|---|---|---|---|---|
| unsorted | $\sim\!3\times10^{-9}$ | $2.0\times10^{-8}$ | $1.1\times10^{-7}$ | $1.3\times10^{-6}$ | $1.0\times10^{-5}$ | $8.1\times10^{-5}$ |
| sorted | $<10^{-11}$ and **falling** with $\alpha$ | | | | | |

> ⚠️ **State the term-count convention.** "50 terms" here means **50 odd terms, $n\le99$**. Under "$n\le50$" (25 odd terms) the magnitudes are ~15× larger. Set $n_{\max}$ **adaptively** (K4); never hard-code.

**R-c — limits and independent anchors.**
$$\alpha\to\infty:\ B\to1-\frac{0.6302489}{\alpha},\quad r\to\frac{12\mu}{ab^3};\qquad \alpha=1:\ B=0.4217310,\quad r=\frac{28.454154\,\mu}{h^{4}}$$

$$\boxed{\;f\!\cdot\!\mathrm{Re}=\frac{24\,D_h^{2}}{b^{2}B(\alpha)}\;\Longrightarrow\;\begin{cases}\alpha=1:& 24/B(1)=\mathbf{56.908}\quad(\text{textbook }56.91)\\[2pt]\alpha\to\infty:& 24(2b)^2/b^2=\mathbf{96}\quad(\text{textbook }96)\end{cases}}$$

These are **external anchors**: they validate $B$ against tabulated duct data without reference to anything in this document. (K10, K11.)

**R-d — inversion.** $r$ is strictly decreasing on $(0,\infty)$, hence bijective. $w=r^{-1}(y)$ is the unique root of $r(w)-y$, bracketed from the slender seed $w_0=(12\mu/(hy))^{1/3}$, solved to $\varepsilon_{\rm inv}$.

---

## 4. Extraction kernel

The exact flux fraction of the **flat** slab $0\le y\le\delta$ in a rectangular duct $w\times h$. Let $\beta=h/w$, $x=\delta/w$, $t_n=n\pi\beta/2$.

$$\boxed{\;f=F_{\rm lam}(x,\beta)=\frac{\displaystyle\sum_{n\,\rm odd}A_n(\beta)\bigl[1-\cos(n\pi x)\bigr]}{\displaystyle2\sum_{n\,\rm odd}A_n(\beta)},\qquad A_n(\beta)=\frac{2}{\pi n^{5}}\Bigl[t_n-\tanh t_n\Bigr]\;}$$

*(Equivalently $A_n=\tfrac1{n^4}[\beta-\tfrac{2}{n\pi}\tanh t_n]$ — but that form cancels; see E-e.)*

**E-a — positivity.** $t-\tanh t>0$ for $t>0$ ⟹ $A_n>0$ ∀$n,\beta$.

**E-b — cubic limit.** $\beta\to\infty\Rightarrow A_n\to\beta/n^4$ and $F_{\rm lam}\to g(x)\equiv3x^2-2x^3$. (Verified: $F_{\rm lam}(0.3,10^4)=0.216000=g(0.3)$; $F_{\rm lam}(0.2,10^6)=0.104000=g(0.2)$.)

**E-c — hard bound, every $\beta$.** $\cos(n\pi/2)=0$ for odd $n$ ⟹ $F_{\rm lam}(\tfrac12,\beta)=\tfrac12$ **exactly**. Hence $f<\tfrac12\iff x<\tfrac12\iff\delta<w_c/2$, with **no aspect-ratio caveat**.

**E-d — leading coefficient and the depth factor.**
$$\boxed{\;c_2(\beta)=\lim_{x\to0}\frac{F_{\rm lam}}{x^{2}}=\frac{\pi^{2}}{4}\cdot\frac{\sum n^{2}A_n(\beta)}{\sum A_n(\beta)},\qquad K_{\rm depth}(\beta)=\sqrt{\frac{3}{c_2(\beta)}}\;}$$
$c_2$ is strictly decreasing, $\to3$ as $\beta\to\infty$, and **diverges as $c_2\simeq1.7/\beta$ as $\beta\to0$** (the shallow limit, where $F_{\rm lam}\to x$).

**E-e — two numerical hazards, both mandatory.**
1. **Cancellation.** $\beta-\tfrac{2}{n\pi}\tanh t_n$ loses up to 3.5 digits at small $t_n$. Compute
$$t-\tanh t=\frac{t^{3}}{3}-\frac{2t^{5}}{15}+\frac{17t^{7}}{315}-\frac{62t^{9}}{2835}+\frac{1382t^{11}}{155925}-\cdots\quad(t<\tfrac12),$$
directly otherwise.
2. **Slow tail in $c_2$.** $\sum A_n$ converges like $\sum n^{-4}$; $\sum n^2A_n$ only like $\sum n^{-2}$. **Add the analytic tail** $\sum_{n\,\rm odd>n_{\max}}n^2A_n\simeq\beta/(2n_{\max})$, or $c_2$ is wrong in the 4th digit at $n_{\max}=10^3$. ($F_{\rm lam}$ itself needs only $n_{\max}\!\sim\!200$.)

**E-f — inversion.** $F_{\rm lam}$ is strictly increasing in $x$ ⟹ unique root on $[0,\tfrac12]$. **Seed $x_0=\sqrt{f/c_2(\beta)}$.** **Never seed with the linear $f$**: at $\beta=0.53$, $f=1.04\times10^{-4}$, the exact root is $4.79\times10^{-3}$ — the linear seed misses by **46×**. In the cubic branch, $x_0$ may be used as the *answer* when $\sqrt{f/27}<\varepsilon_{\rm cub}$.

**E-g — thin-lamina condition (a design gate, not a remark).** The near-wall quadratic $F_{\rm lam}\approx c_2x^2$ is valid for $x\ll\beta$, i.e.
$$\boxed{\;\frac{\delta(i)}{h}\ \le\ \varepsilon_{\rm lam}\;}$$
The plug law ($F_{\rm lam}\to x$) needs $\delta\gtrsim h$ and is **never reached** in practice ($\delta/h\sim3\times10^{-3}$–$3\times10^{-2}$). The real finite-$\beta$ effect is the prefactor $c_2(\beta)\ne3$: a **3–19 %** correction that **shrinks as the channel deepens**. *(Assert it — M14. The exact inversion needs no such bound; the seed, S1′, and the steric picture of §6.2 all do.)*

$$\boxed{\;\delta(i)=x\bigl(f(i);\beta(i)\bigr)\cdot w_c(i)\;}$$

> ### ⚠️ **What is exact, and what is a model.**
> $F_{\rm lam}$ is the **exact** flux fraction of a **flat** slab $y<\delta$ in a **fully developed, unperturbed** duct. Identifying that slab with the fluid actually skimmed through a link is a **modelling step**: the true dividing streamsurface is $z$-dependent (thicker near the top and bottom walls, where $u$ is slower), and the wall it hugs is $\phi\approx0.4$ **open** (§20). Both effects are $\beta$- and $\phi$-dependent and are **absorbed into $\kappa$**, not into $\eta$. **$\eta$ is a yield correction on the flow split; it has nothing to do with the lamina.**

---

## 5. Kernel–prefactor pairing (the rule that cannot be broken)

$\kappa$ is meaningless without the kernel it was fitted with. Declare a pair $\mathcal K\in\{g,\ F_{\rm lam}\}$ and $\kappa_{\mathcal K}$, and use $\mathcal K$ at **all three** call sites — inversion $x=\mathcal K^{-1}(f;\beta)$, schedule S1, scale $D=\kappa_{\mathcal K}\delta$. **P1 asserts one kernel object at all three.**

**Conversion.** $\ \boxed{\kappa_{F_{\rm lam}}=\kappa_{g}\big/K_{\rm depth}(\beta_{\rm cal})}$ — requires $\beta_{\rm cal}$, which is **unknown**.

**Transfer bias (mandatory report, never applied).**
$$\boxed{\;\frac{D_{\rm true}}{D_{\rm predicted}}=\sqrt{\frac{c_2(\beta_{\rm cal})}{c_2(\beta(i))}}\;}$$

The supplied $\kappa=4.9138$ was fitted on an array of $w_c\in\{100,125,150\}$ µm whose **depth is unreported**. Sweeping that one number, for a design tapering $75\to52$ µm at $h=140$ µm:

| assumed $h_{\rm cal}$ | mean $c_2(\beta_{\rm cal})$ | bias | $\kappa_{F_{\rm lam}}$ | $N$ | vs cubic ($N=2496$) |
|---|---|---|---|---|---|
| 80 µm | 4.177 | 1.141–1.156 | 5.798 | 3168 | **+26.9 %** |
| 100 µm | 3.817 | 1.090–1.105 | 5.543 | 2895 | +16.0 % |
| 120 µm | 3.605 | 1.060–1.074 | 5.387 | 2734 | +9.5 % |
| 140 µm | 3.470 | 1.040–1.053 | 5.285 | 2632 | +5.4 % |
| 150 µm | 3.421 | 1.032–1.046 | 5.247 | 2594 | **+3.9 %** |

> **One unreported number controls a 27 % swing in device length.** Until it exists: **use the cubic + $\kappa_g$ (self-consistent with the fit), and report the bias band as an uncertainty.** Do not silently apply a correction whose calibration you cannot verify.
>
> ⚠️ **Switching to $F_{\rm lam}$ removes only the *flux-fraction* part of $\kappa$'s $\beta$-dependence.** The dividing-surface-shape and wall-porosity parts (§4 box, §20) remain inside $\kappa$ and are *also* $\beta$-dependent. **Only calibrating at the design $\beta$ and $\phi$ removes all of it.**

---

## 6. Derived scale, window, and the transport lemma

### 6.1 Scale and window
$$\boxed{\;D(i)=\kappa_{\mathcal K}\,\delta(i),\qquad u(i)\equiv\frac{D(i)}{w_c(i)}=\kappa_{\mathcal K}\,x\bigl(f(i);\beta(i)\bigr)\;}$$
$\kappa$ is constant only for $u\in[u_{\rm lo},\sim\!0.10]$. Outside, §18 emits a **band**, never a scalar.

### 6.2 Transport lemma — why $\delta$ is the right scale, and where $\kappa$'s residual must live

Let $\varphi$ be the flux fraction of a fluid element measured **from one wall**. A link at that wall removes $fQ_c$; the opposite link removes another $fQ_c$ but does not change $\varphi Q_c$. Hence

$$\varphi' = \frac{\varphi-f}{1-2f}\qquad\Longleftrightarrow\qquad \Bigl(\varphi'-\tfrac12\Bigr)=\frac{\varphi-\tfrac12}{1-2f}$$

**The centreline $\varphi=\tfrac12$ is a repelling fixed point.** $(\varphi-\tfrac12)$ grows by exactly $1/P(i)$: fluid **migrates deterministically toward the wall**, and is extracted when $\varphi\le f$. Solving the threshold:

$$\boxed{\;\varphi_0\le\frac{1-P(N)}{2}=\frac F2\;\Longrightarrow\;\text{extracted; both walls give }2\times\frac F2=F\;}$$

**Verified numerically to 6 decimals** ($\varphi_{\rm crit}=0.216643$, $F/2=0.216643$).

Three consequences, all load-bearing:

1. **The device requires no cross-streamline transport.** Diffusion is a perturbation, not a mechanism. (Per rung: deterministic migration $\sim\delta/2=0.22$ µm vs diffusive $\sqrt{2D_{\rm diff}\Delta t}=0.033$ µm — **6.7× smaller**.)
2. **The cutoff is a hard steric wall constraint**, not a statistical position criterion: a sphere of radius $a$ cannot bring its centre closer than $a$ to the wall, so it is extracted iff $a\le\delta$. **This gives $\kappa_{\rm steric}=2$ exactly** ($D_c=2\delta$).
3. **Therefore the entire measured residual $\kappa/2=2.46$ must live in the three model gaps, and nowhere else:** the flat-lamina idealisation (§4 box), the wall porosity $\phi$ (§20), and finite-particle displacement. **The lemma pins the list.**

---

## 7. Schedule laws

**S1 — invariant lamina (canonical). Exact at every aspect ratio.**
$$\boxed{\;f(i)=\mathcal K\!\left(\frac{\delta^\ast}{w_c(i)};\ \beta(i)\right)\;}\Longrightarrow\;\delta(i)=\delta^\ast,\ D(i)=\kappa_{\mathcal K}\delta^\ast\ \text{exactly},\ \forall i,\ \forall\mathcal W.$$

**S1′ — leading-order shortcut** (required only for §13's closed form):
$$f(i)=c_2\bigl(\beta(i)\bigr)\left(\frac{\delta^\ast}{w_c(i)}\right)^{2}\ \xrightarrow{\ \beta\gg1\ }\ 3\left(\frac{\delta^\ast}{w_c(i)}\right)^{2}$$
**S1′ ≠ S1**: the difference is the $O(x^3)$ tail; relative $\delta$-drift $\simeq\sqrt{f/27}$ in the cubic branch ($\approx0.2\,\%$ at $f=10^{-4}$). **Declare which you used.** M8 tests S1 exactly.

**S2 — free schedule.** Any $f:\{1..N\}\to(0,\tfrac12)$; $\delta,D$ drift; §18 reports $D(i)$ as a profile and §12's input guard does not apply.

---

## 8. March

$$\boxed{\;\sigma(0)=\sigma_0,\qquad \sigma(i)=\frac{\sigma(i-1)+f(i)}{1-2f(i)},\qquad \rho(i)=\frac{1}{\sigma(i)}\;}$$

| mode | $\sigma_0$ | meaning |
|---|---|---|
| **FOLD** | $\sigma_0=0$ **exactly** | sealed inlet. $\rho(1)=\dfrac{1-2f(1)}{f(1)}$; $w_s(1)\ll w_{\min}$ typically ⟹ fold. |
| **DIRECT** | $\sigma_0=\dfrac{1-2f(1)}{\rho_1}-f(1)$, $\rho_1:=\dfrac{r(w_{\rm start})}{r(w_c(1))}$ | prescribed printable start width. **$\sigma_0>0$ is a fictitious inlet filtrate = the mass defect** (§10). Row 1 is off-design by construction; §15 quantifies it. |

**Assert $\sigma_0\ge0$** ($\sigma_0<0$ ⟹ $w_{\rm start}$ narrower than the closure requires — unphysical).

**Closed forms.** With $P(i)=\prod_{j=1}^i[1-2f(j)]$, $S(i)=\sigma_0+\sum_{j=1}^i f(j)P(j-1)$:
$$\boxed{\;\sigma(i)=\frac{S(i)}{P(i)},\quad\frac{Q_c(i)}{Q_c(0)}=P(i),\quad\frac{Q_s(i)}{Q_c(0)}=S(i),\quad\frac{Q_{\rm gap}(i)}{Q_c(0)}=f(i)P(i-1)\;}$$

---

## 9. Closure

$$\boxed{\;r\bigl(w_s(i)\bigr)=\rho(i)\,r\bigl(w_c(i)\bigr)\ \Longrightarrow\ w_s(i)=r^{-1}\!\Bigl(\rho(i)\,r\bigl(w_c(i)\bigr)\Bigr)\;}$$
$\Lambda(i)$ cancels identically. If folds carry a bend correction (§14.6), close against $R_s^{\rm eff}$.

---

## 10. Terminal fraction and the mass identity

**Theorem.** $\ \boxed{P(i)+2S(i)=1+2\sigma_0}\ $ for all $i$.
*Proof.* True at $i=0$. If true at $i-1$: $P(i)+2S(i)=P(i{-}1)[1-2f(i)]+2S(i{-}1)+2f(i)P(i{-}1)=P(i{-}1)+2S(i{-}1)$. $\square$

**Corollary — the schedule's terminal fraction, for any $\sigma_0$:**
$$\boxed{\;F_{\rm sched}=1-P(N)=2\sum_{j=1}^{N}f(j)P(j-1)\;}$$
(The fraction of $Q_c(0)$ that the schedule removes from the primary; a FOLD start delivers it. Revision 5 called it the physical $F$ for any $\sigma_0$, which a DIRECT start does not deliver: §21, "Fixed in revision 5.1".) **Independently confirmed by §6.2's migration lemma.**

**The other two forms carry the seed:** $\ 2S(N)=F_{\rm sched}+2\sigma_0$, $\ \dfrac{2}{\rho(N)+2}=\dfrac{F_{\rm sched}+2\sigma_0}{1+2\sigma_0}$. They coincide with $F_{\rm sched}$ **iff $\sigma_0=0$**. **The second is the delivered fraction:** a DIRECT start, sealed in the drawn network as every start is (§15), delivers $\boxed{F=2/(\rho(N)+2)=1-P(N)/(1+2\sigma_0)}$ (rev. 5.1).

**Outlet resistance ratio — generalized.** $S(N)=\sigma_0+F_{\rm sched}/2$, so, with the delivered $F$,
$$\boxed{\;\rho^{\rm out}=\frac{P(N)}{S(N)}=\frac{1-F_{\rm sched}}{\sigma_0+F_{\rm sched}/2}=\frac{2(1-F)}{F}\ \text{ for any }\sigma_0\;}$$

> ### ⚠️ **The run-out is load-bearing, and it is the most sensitive un-gated input in the document.**
> $R_c,R_s$ in $F=2R_c/(R_s+2R_c)$ are the **run-out** resistances (last link → terminals), reducing to $\rho(N)$ **iff** the two run-outs have equal length and preserve their widths. **A $3\times$ run-out mismatch collapses $F$ from $0.519$ to $\approx0.27$.** This is gated by **G13** and must be.

---

## 11. Inverse map (widths → schedule)

$\rho(i)=r(w_s(i))/r(w_c(i))$, $\sigma(i)=1/\rho(i)$, and
$$\boxed{\;f(i)=\frac{\sigma(i)-\sigma(i-1)}{1+2\,\sigma(i)}\;}\qquad\text{(canonical: no singularity, no special case, any }\sigma_0)$$

*Equivalently* $f(i)=\dfrac{1-\rho(i)/\rho(i-1)}{\rho(i)+2}$ with $\rho(0)=1/\sigma_0$ — but that form needs $\rho(0)=\infty$ when $\sigma_0=0$ and is a bug source. **Use the $\sigma$-form.**

**Admissibility.** $f(i)>0\iff\sigma$ strictly increasing $\iff\rho$ strictly decreasing. $f(i)<\tfrac12$ automatically (since $\sigma_{i-1}\ge0$).

**Scope.** This inverts the **design map** (§8–§9): *given these widths, what schedule does the closure imply?* It is **not** a physics solver — for that, §15.

---

## 12. Design inversions

$$\boxed{\;\delta^\ast=\frac{D_{\rm target}}{\kappa_{\mathcal K}}\ \Longrightarrow\ f^\ast=\mathcal K\!\left(\frac{\delta^\ast}{w_c(0)};\ \beta(0)\right)\;}$$
$$\boxed{\;\rho^{\rm out}=\frac{2(1-F_{\rm target})}{F_{\rm target}}\ \ (F_{\rm target}\text{ delivered, any }\sigma_0)\;}\qquad\boxed{\;w_{\min}=\max\!\Bigl(F_{\min},\ \frac{h}{\mathrm{AR}_{\max}}\Bigr)\;}$$

**Termination — on $P$, not on $w_s$.** $P$ is strictly decreasing, so
$$\boxed{\;\text{stop at the first }i\text{ with }P(i)\le(1-F_{\rm target})(1+2\sigma_0);\qquad N:=i\;}$$
is monotone-safe for **any** $\mathcal W$. Then $w_s(N)$ is an **output**; assert $w_s(N)\ge w_{\min}$.

**Plateau guard on inputs (S1 only).** $u(i)=D_{\rm target}/w_c(i)$ is smallest where $w_c$ is largest:
$$\boxed{\;\max_{t\in[0,1]}\mathcal W(t)\ \le\ \frac{D_{\rm target}}{u_{\rm lo}}\;}$$

---

## 13. Canonical taper — one potential, no fixed point

$$\boxed{\;\Psi(w)\ \equiv\ \frac{w^{2}}{c_2\bigl(h/w\bigr)}\;}$$

**Ψ is strictly increasing**, hence invertible:
$$\frac{\mathrm d\ln\Psi}{\mathrm d\ln w}=2+\frac{\mathrm d\ln c_2}{\mathrm d\ln\beta}\in[1,2]$$
because $c_2$ decreases in $\beta$ with $\mathrm d\ln c_2/\mathrm d\ln\beta\in(-1,0)$ ($c_2\to3$ as $\beta\to\infty$; $c_2\sim1.7/\beta$ as $\beta\to0$). *Verified: $1.905$ at $w=75$, $1.947$ at $w=52$.*

**Objective: $Q_{\rm gap}=$ const.** Under S1′,
$$Q_{\rm gap}=f\,Q_c=c_2(\beta)\Bigl(\frac{\delta^\ast}{w_c}\Bigr)^{2}Q_c=\frac{(\delta^\ast)^{2}\,Q_c}{\Psi(w_c)}\qquad\Longrightarrow\qquad Q_{\rm gap}\ \text{const}\iff Q_c\propto\Psi(w_c)$$

Continuity ($\mathrm dQ_c/\mathrm di=-2Q_{\rm gap}=$ const) then gives **everything, from one linear relation**:

$$\boxed{\;\frac{\mathrm d\Psi\bigl(w_c(i)\bigr)}{\mathrm di}=-2(\delta^\ast)^{2}\qquad\Longrightarrow\qquad \Psi\bigl(w_c(i)\bigr)=\Psi\bigl(w_c(0)\bigr)-2(\delta^\ast)^{2}\,i\;}$$

$$\boxed{\;w_c(i)=\Psi^{-1}\!\Bigl[\Psi(w_c(0))-2(\delta^\ast)^{2}i\Bigr],\qquad P(i)=\frac{\Psi(w_c(i))}{\Psi(w_c(0))}\ \ \text{(linear in }i)\;}$$

$$\boxed{\;N=\frac{\Psi(w_c(0))-\Psi(w_c(N))}{2(\delta^\ast)^{2}},\qquad F=1-\frac{\Psi(w_c(N))}{\Psi(w_c(0))},\qquad w_c(N)=\Psi^{-1}\!\Bigl[(1-F_{\rm target})\,\Psi(w_c(0))\Bigr]\;}$$

> ⚠️ **v4 wrote $F=1-\bigl(w_c(N)/w_c(0)\bigr)^2$. That is the cubic special case only.** In the $F_{\rm lam}$ branch it omits the factor $c_2(\beta_0)/c_2(\beta_N)=1.0265$ and gives $F=0.51929$ against the true $0.50654$ — **off by $0.0126$, i.e. 13× the M12 tolerance.** (The march is unaffected: it terminates on $P$, which is exact.) **The $\Psi$ formulation makes the error structurally impossible: the same potential appears in $w_c(i)$, $N$, $P$, $F$ and the inversion.**

**Cubic branch** ($c_2\equiv3$, so $\Psi=w^2/3$): everything reduces to
$$w_c(i)=\sqrt{w_c(0)^2+[w_c(N)^2-w_c(0)^2]\tfrac iN},\quad N=\frac{w_c(0)^2-w_c(N)^2}{6(\delta^\ast)^2},\quad F=1-\Bigl(\tfrac{w_c(N)}{w_c(0)}\Bigr)^2,\quad w_c(N)=w_c(0)\sqrt{1-F}$$

**Equivalence chain — with its scope conditions.**
$$Q_{\rm gap}\ \text{const}\iff\Delta P_{\rm cross}\ \text{const}\iff R_c(i)Q_c(i)=R_s(i)Q_s(i)\ \text{exactly}\iff\Psi(w_c)\ \text{linear in }i$$
holds **iff all links are identical**, because the underlying identity is
$$R_c(i)Q_c(i)-R_s(i)Q_s(i)=R_{\rm gap}(i)Q_{\rm gap}(i)-R_{\rm gap}(i{+}1)Q_{\rm gap}(i{+}1)$$
which collapses only for uniform $R_{\rm gap}$. **The invariant that always holds is $\Delta P_{\rm cross}=$ const.** Exact under S1′; under S1 an $O(x)$ residual remains in $Q_{\rm gap}$ flatness (N6).

**This taper does not extremize the shear functional** (§16-b). For any other $\mathcal W$: no closed form; use the damped fixed point (§18).

---

## 14. Realization map

$$\boxed{\;p=\mathrm{GS}+\max\bigl(F_{\min},2R_{\min}\bigr),\qquad R_{\rm bend}=\frac p2\ \text{(derived)},\qquad \Phi(i)=\frac{r(w_s(i))}{r(\mathrm{GS})},\qquad L(i)=\Phi(i)\Lambda(i)\;}$$

$\Phi$ is independent of $\Lambda$. Two adjacent strokes at pitch $p$ force a U-turn of centreline radius $p/2$ and inner-edge radius $(p-\mathrm{GS})/2$; that is what $R_{\min}$ constrains.

### 14.1 Dispatch
$$\text{row }i=\begin{cases}\textbf{STRAIGHT},& w_s(i)\ge w_{\min}\\ \textbf{FOLD},& w_s(i)<w_{\min}\end{cases}$$
STRAIGHT: width $w_s(i)$, length $\Lambda(i)$. FOLD: width $\mathrm{GS}$, drawn traversal $\hat L(i)$; requires $\Lambda(i)\ge2p$.

### 14.2 The three fold requirements
1. **Length.** $\hat L(i)=L(i)$ (or $L^{\rm eff}$, §14.6).
2. **Cross-section.** Width $\equiv\mathrm{GS}$ along the whole traversal.
3. **Topology.** Link $i$ has degree $\ge2$ in the drawn secondary graph, and no path $i\to i{+}1$ is shorter than $L(i)$.

### 14.3 Packing
A slot ($\Lambda$ along flow $\times$ $W$ across) admits strokes at pitch $p$, so its maximum path is $\Lambda W/p$:
$$\boxed{\;\Phi(i)\le\frac Wp\;}\qquad\boxed{\;\mathrm{arms}(i)=\Bigl\lceil\frac{\Phi(i)}{W/p}\Bigr\rceil-1\;}$$
Global **necessary** condition: $\sum_i\Phi(i)\Lambda(i)\le\frac Wp\sum_i\Lambda(i)$.

> ⚠️ **Necessary, not sufficient — and there is no routability guarantee.** These conditions bound *area*; they do not prove a planar, non-self-intersecting placement of the arms exists. **Global routability is the implementer's obligation.** The spec supplies necessary conditions and G2/G3/G5 to catch failures after the fact.

### 14.4 Links, posts, pitch, porosity
$$\Lambda(i)\ge\max\bigl(G+g_p,\ 2p\ \text{[if folded]}\bigr),\quad G\ge w_{\min},\quad g_p\ge F_{\min},\quad \frac{\mathrm{GS}}{G}\in\gamma_{\rm GS},\quad \phi=\frac{G}{\Lambda}$$
$$\boxed{\;R_{\rm gap}=r(G)\bigl[\ell_g+c_{\rm ent}D_h(G,h)\bigr],\qquad D_h=\frac{2Gh}{G+h}\;}$$
The bare $r(G)\ell_g$ is a **lower bound**: for $\ell_g\lesssim D_h$ the link is an orifice, not a duct. Declare $c_{\rm ent}$.

> ⚠️ **Report $\phi$.** The kernels of §3–§4 assume a *smooth, solid* wall. The real wall is $\phi\approx0.4$ **open**. The effective-slip offset is $O(0.3\ \mu\text{m})$ — the same order as $\delta$ itself. It is absorbed into $\kappa$ at $\phi_{\rm cal}$ (§20), and it is **not** controlled by gap dominance: porosity is geometric and is present at zero gap flow.

### 14.5 Turns
A gap-free 180° turn **is one segment** (§1) and takes $\rho(i)$ of its **upstream** rung — no interpolation. With primary arc $L_p$ and unequal secondary arcs $L_{\rm in}<L_{\rm out}$:
$$\boxed{\;r\bigl(w_s^{\rm in}\bigr)L_{\rm in}=r\bigl(w_s^{\rm out}\bigr)L_{\rm out}=\rho(i)\,r\bigl(w_c^{\rm turn}\bigr)L_p\;}$$
$L_{\rm in}<L_{\rm out}\Rightarrow w_s^{\rm in}<w_s^{\rm out}$. **§15's two-rail folding is valid only if this holds** (N8).

**Leg partition.** Contiguous blocks with $\sum_{i\in\text{leg}_k}\Lambda(i)\le\ell_{\rm leg}$.

### 14.6 Bend correction (Stokes-corner)
$$\boxed{\;R_s^{\rm eff}(i)=r(\mathrm{GS})\bigl[\hat L(i)+n_{\rm turns}(i)L_{\rm eq}\bigr]\ \Longrightarrow\ \hat L(i)=\Phi(i)\Lambda(i)-n_{\rm turns}(i)L_{\rm eq}\;}$$
with $n_{\rm turns}(i)\simeq\hat L(i)/\ell_{\rm stroke}$ (one iteration, or $\hat L\simeq\Phi\Lambda/[1+L_{\rm eq}/\ell_{\rm stroke}]$).

**Relative excess $=L_{\rm eq}/\ell_{\rm stroke}$, independent of $\hat L$** — typically 2–6 %, applying **only to folded rows** (a localized inlet-side error).

$L_{\rm eq}$ is a **length**, calibrated once at low Re. It is **not** a Dean correction: the secondary rail carries $Q_s$, giving $De\sim10^{-4}$–$0.15$ (§21). Do not import $De$-based correlations.

---

## 15. Network model

Fold the two secondary rails into one by symmetry (halve their resistances; the merged rail carries $2Q_s$). **Valid only if G10 holds** (N8).

| edge | resistance | carries |
|---|---|---|
| $(p_i,p_{i+1})$ | $R_c(i)=r(w_c(i))\Lambda(i)$ | $Q_c(i)$ |
| $(q_i,q_{i+1})$ | $\tfrac12R_s(i)$, $R_s(i)=r(\mathrm{GS})\hat L(i)$ or $r(w_s(i))\Lambda(i)$ | $2Q_s(i)$ |
| $(p_i,q_i)$, $i=1..N$ | $\tfrac12R_{\rm gap}$ | $2Q_{\rm gap}(i)$ |
| $(p_N,p_\ast)$, $(q_N,q_\ast)$ | $R_{\rm out}^{c}$, $\tfrac12R_{\rm out}^{s}$ | run-outs |

Inject $Q_{\rm in}$ at $p_0$; ground $p_\ast,q_\ast$. **$q_1$ has no upstream edge, so $Q_s(0)=0$ is structural.** Solve $\mathbf{GP}=\mathbf J$; extract $Q_{\rm gap}^{\rm em}(i)=[P(p_i)-P(q_i)]/R_{\rm gap}$, $f^{\rm em}(i)=Q_{\rm gap}^{\rm em}(i)/Q_c^{\rm em}(i-1)$, $F^{\rm em}=2Q_s^{\rm em}(N)/Q_{\rm in}$.

**Terminal layer — derived.** Let $V(i)=P(p_i)-P(q_i)$. The **bilateral** primary balance $Q_c(i{-}1)-Q_c(i)=2Q_{\rm gap}(i)$ gives
$$\boxed{\;V(i{+}1)-2V(i)+V(i{-}1)=q\,V(i),\qquad q=\frac{2R_c+R_s}{R_{\rm gap}}\;}$$
*(the **2** on $R_c$ is the same bilateral factor as everywhere else; a discriminating uniform-ladder experiment confirms it — two rails give $\lambda=1.390$, one rail $1.739$)*. Solutions $V\propto z_\pm^{\,i}$, $z_\pm=\tfrac12[(2+q)\pm\sqrt{(2+q)^2-4}]$, and $\lambda=1/\ln z_+$.

> ⚠️ **$q$ is NOT monotone.** $\rho$ falls, but $R_c$ **grows** as $w_c$ narrows, and near the outlet that wins: for the lineage device $q$ has a shallow interior minimum ($q_{\min}=0.501$ at $86\,\%$ of span) against $q(N)=0.510$, and $\mathrm dq/\mathrm di>0$ at $i=N$. **v4's justification ("$\rho$ falls monotonically") was incomplete.** Evaluate $\lambda$ as $\max$ over the terminal window $[N-m_{\rm budget},N]$ — conservative and correct.

**Layer depth — an amplitude, not just a rate.** The decay is $|e(i)|=A\,z_-^{\,N-i}$, so
$$\boxed{\;m=\Bigl\lceil\lambda\,\ln\!\bigl(A/\varepsilon_{\rm net}\bigr)\Bigr\rceil,\qquad A=\Bigl|\tfrac{f^{\rm em}(N)}{f(N)}-1\Bigr|\;}$$
$A$ is set by the **run-out mismatch**. With matched run-outs $A\approx0.05$ and $m$ is *smaller* than the unit-amplitude estimate; with a $3\times$ mismatch $A\approx5\times10^{2}$ and $m$ roughly doubles. **The exponential rate $\lambda$ is exact; the prefactor is a boundary condition.** Also **measure** $m$ (largest $i^\ast$ with bulk agreement; $m=N-i^\ast$) and assert $m\le m_{\rm budget}$ (N7).

---

## 16. Optional functionals and regime gates

### 16-a. Reynolds — three of them, all gated (M14)
$$\mathrm{Re}(i)=\frac{\rho_f\,Q_c(i)\,D_h(i)}{\mu\,h\,w_c(i)},\qquad D_h(i)=\frac{2w_c(i)h}{w_c(i)+h}$$
**Use $Q_c(i)$, not $Q_{\rm in}$** — maximum at $i=0$. Also compute:
- $\mathrm{Re}_{\rm gap}=\rho_fU_{\rm gap}D_h(G,h)/\mu$, $U_{\rm gap}=Q_{\rm gap}/(Gh)$ — **must be $\ll1$**: the skimming junction has to be Stokes, or §6.2's streamline picture fails.
- $\mathrm{Re}_\delta=\rho_f\,\dot\gamma\,\delta^2/\mu$ — **must be $\ll1$**.

*Lineage values: $\mathrm{Re}=4.7$ (inlet) → $2.3$ (outlet); $\mathrm{Re}_{\rm gap}=6.5\times10^{-4}$; $\mathrm{Re}_\delta=1.0\times10^{-3}$.*

### 16-b. Maximum wall shear rate
$$\boxed{\;\dot\gamma_{\max}(i)=C(\alpha_i)\frac{Q_c(i)}{a_ib_i^{2}},\qquad C(\alpha)=\frac{48}{\pi^{2}B(\alpha)}\sum_{n\,\rm odd}\frac{1}{n^{2}}\Bigl[1-\operatorname{sech}\frac{n\pi\alpha}{2}\Bigr]\;}$$
Same $B(\alpha)$ as §3 and §4 (assert one call site — S3). The maximum sits at the **long-wall midpoint** at every $\alpha$.

| $\alpha$ | 0.53 | 0.80 | 1.00 | 1.40 | 1.87 | 2.00 | 2.69 | 4.00 | $\infty$ |
|---|---|---|---|---|---|---|---|---|---|
| $C$ | 12.51 | 10.35 | **9.6078** | 8.80 | 8.25 | **8.133** | 7.65 | 7.10 | **6** |

$C\to6$ **only** as $\alpha\to\infty$. **Substituting 6 at $\alpha=2$ understates peak shear by 26.2 %** *and* discards a real drift along any taper. **Constant $Q_{\rm gap}\ne$ constant $\dot\gamma_{\max}$** at finite $\alpha$; for constant shear integrate $\mathrm d(w_c^2)/\mathrm di=-6(\delta^\ast)^2/[1-\tfrac{w_c}{2C}\tfrac{\mathrm dC}{\mathrm dw_c}]$ and accept $Q_{\rm gap}\propto1/C$ drift. **One or the other.**

### 16-c. Entrance length and quasi-developed taper (M15, M16)
The kernels assume **fully developed** flow. The correct laminar entrance length has a **Stokes floor** and does *not* vanish as $\mathrm{Re}\to0$:
$$\boxed{\;\frac{L_e}{D_h}\simeq\Bigl[0.619^{1.6}+(0.0567\,\mathrm{Re})^{1.6}\Bigr]^{1/1.6}\ \xrightarrow{\ \mathrm{Re}\to0\ }\ 0.62\;}$$
*(The naive $L_e=0.05\,\mathrm{Re}\,D_h$ **goes to zero** in this regime and is wrong here.)*

For the lineage device $L_e\approx70\ \mu$m against $\Lambda=49\ \mu$m: **the flow does not classically redevelop between rungs.** This is **benign, and the reason must be stated**: the per-rung perturbation is $O(f)=10^{-4}$ of the base flow, so the resistance error is $O(f^2)\sim10^{-8}$ — four orders below the closure tolerance. **Report $L_e/\Lambda$; do not gate on it. Gate instead on the perturbation, which is M13.**

**Quasi-developed taper.** The kernels are applied to a slowly varying width. Assert
$$\boxed{\;\Bigl|\frac{\mathrm dw_c}{\mathrm dx}\Bigr|=\Bigl|\frac{w_c(i)-w_c(i-1)}{\Lambda(i)}\Bigr|\le\varepsilon_{\rm taper}\;}$$
*(lineage: $1.7\times10^{-4}$ — four orders of margin.)*

### 16-d. Pressure budget — reported, not gated
$$\Delta P_{\rm tot}=\sum_iR_c(i)Q_c(i)=Q_{\rm in}\sum_ir(w_c(i))\Lambda(i)P(i)$$
Substrate-compliance limits are **out of scope** (they need a modulus and a deflection model, are fabrication-specific, and have nothing to check against here). Report; gate externally.

### 16-e. Sensitivity map — **mandatory report** (M17)
$$\frac{\partial\ln(\,\cdot\,)}{\partial\ln(\,\cdot\,)}\quad\text{at the lineage design point:}$$

| | $w_s(\rm out)$ | $w_c(\rm out)$ | $w_c$ | $h$ | $f$ | $\kappa$ |
|---|---|---|---|---|---|---|
| $F$ | $\mathbf{+1.21}$ | $-1.19$ | — | $\mathbf{-0.03}$ | — | — |
| $D$ | — | — | $+0.95$ | $\mathbf{+0.05}$ | $+0.50$ | $+1$ |
| $N$ | — | — | — | $+0.14$ | — | $\mathbf{+2}$ |

> **Read it once and it changes your process priorities: the *widths* carry the device; the *etch depth* barely does.** A 1 % mask error on the side-channel outlet moves $F$ by 1.2 %; a 1 % depth error moves it by 0.03 %. And $N\propto\kappa^2$, so $\kappa$'s ±3.4 % CI alone is ±6.8 % on row count.

---

## 17. Invariance table (every entry provable)

| rescale | $\rho,w_s$ | $f,\delta,D$ | $F$ | $N$ | $\Phi$ | $L=\Phi\Lambda$ | $R_s=r({\rm GS})L$ | $Q_{\rm gap}$ | $R_{\rm gap}/R_c$ | Re |
|---|---|---|---|---|---|---|---|---|---|---|
| $\mu\to\lambda\mu$ | — | — | — | — | — | — | $\lambda$ | — | — | $\lambda^{-1}$ |
| $Q_{\rm in}\to\lambda Q_{\rm in}$ | — | — | — | — | — | — | — | $\lambda$ | — | $\lambda$ |
| $\Lambda\to\lambda\Lambda$ | — | — | — | — | — | $\lambda$ | $\lambda$ | — | $\lambda^{-1}$ | — |
| $\mathrm{GS}\to\lambda\,\mathrm{GS}$ | — | — | — | — | $\frac{r({\rm GS})}{r(\lambda{\rm GS})}$ | $\frac{r({\rm GS})}{r(\lambda{\rm GS})}$ | **—** | — | — | — |
| $\ell_g\to\lambda\ell_g$ | — | — | — | — | — | — | — | — | $\approx\lambda$ | — |
| choose $f^\ast$ | ✓ | ✓ | **—** | ✓ | ✓ | ✓ | ✓ | ✓ | — | — |
| choose $F_{\rm target}$ | ✓ | **—** | ✓ | ✓ | — | — | — | — | — | — |

("—" = exactly invariant.)

> ⚠️ **v4's GS row said $L$ = "same". Wrong.** $L=\Phi\Lambda$ scales with $\Phi$; what is invariant is the **segment resistance** $R_s=r(\mathrm{GS})L$. The column is now explicit. (Verified: $r({\rm GS})L=r(2{\rm GS})L'$ to machine precision.)

**Decoupling lemmas.** $F$ is a functional of the outlet ratio alone; $D$ of $(w_c,f)$ alone. Orthogonal **because $N$ is free.** Pin $N$ and they re-couple.

---

## 18. Algorithm

```
INPUT
  h, mu, rho_f, T
  W : [0,1] -> R+                                  # primary shape function
  D_target, F_target
  F_min, R_min, AR_max, G, g_p, l_g, c_ent, GS, W_budget, l_stroke, L_eq
  n_legs, leg_length, R_out_c, R_out_s, m_budget
  KERNEL in {cubic, F_lam};  kappa_KERNEL          # PAIRED (Sec 5). Never mix.
  beta_cal, phi_cal                                # for the transfer-bias report only
  u_lo, u_ff, eta, f_soft, gamma_GS                # CALIBRATED (Sec 2.3)
  eps_r, eps_inv, eps_cub, eps_id, eps_geo, eps_net, eps_lam, eps_taper, Re_max
  mode in {FOLD, DIRECT};  if DIRECT: w_start

DERIVE
  n_max      : adaptive; |r(n)-r(2n)|/r(2n) < eps_r AT THE LARGEST alpha IN THE DESIGN   # K4
  B(alpha)   : shared series.  ASSERT alpha >= 1 at every call site.                     # K1
  r(w)       := 12 mu / (a b^3 B(a/b)),  a=max(w,h), b=min(w,h)
  r_inv(y)   := unique root of r(w)-y, seeded (12 mu/(h y))^(1/3)                        # K8
  A_n(beta)  := (2/(pi n^5)) * (t - tanh t), t = n pi beta/2,  TAYLOR branch for t<0.5   # E-e.1
  c2(beta)   := (pi^2/4) * [ sum n^2 A_n  +  beta/(2 n_max) ] / sum A_n                  # E-e.2 TAIL
  K(x;beta)  := 3x^2-2x^3          if KERNEL == cubic
                F_lam(x,beta)      if KERNEL == F_lam
  x(f;beta)  := unique root of K(x;beta)-f on [0,1/2], seeded sqrt(f/c2(beta))           # NEVER seed f
  Psi(w)     := w^2 / c2(h/w)                      # strictly increasing; Psi_inv by bracketing
  w_min      := max(F_min, h/AR_max)
  p          := GS + max(F_min, 2*R_min);   R_bend := p/2
  delta_star := D_target / kappa_KERNEL
  R_gap      := r(G) * ( l_g + c_ent * 2*G*h/(G+h) )
  phi        := G / Lambda

PRE-FLIGHT GATES  (fail, do not warn)
  ASSERT max_{t in [0,1]} W(t) <= D_target / u_lo                       # Sec 12 (MAX, not W(0))
  ASSERT GS/G in gamma_GS  and  GS >= w_min  and  G >= w_min
  ASSERT one kernel object at all three call sites                      # Sec 5 pairing -> P1
  ASSERT | R_out_s / R_out_c  /  rho_out  - 1 | < eps_geo               # RUN-OUT MATCH -> G13
  REPORT transfer bias sqrt( c2(beta_cal)/c2(beta(i)) )                 # Sec 5. DO NOT APPLY.
  REPORT phi                                                            # Sec 14.4

# ---------- MARCH  (damped fixed point on N; Sec 13 seeds it exactly for the canonical taper)
N := N_seed
repeat
    sigma := 0.0                                                        # FOLD: EXACT ZERO
    if mode == DIRECT:
        f1 := K(delta_star/W(1/N); h/W(1/N));  rho1 := r(w_start)/r(W(1/N))
        sigma := (1 - 2*f1)/rho1 - f1                                   # seed = the mass defect
        ASSERT sigma >= 0
    sigma_0 := sigma;  P := 1.0;  S := sigma;  rows := []
    for i = 1, 2, 3, ...
        w_c := W(i/N);  beta := h/w_c
        f   := K(delta_star/w_c; beta)                                  # S1 exact (or table for S2)
        ASSERT 0 < f < 0.5
        sigma := (sigma + f)/(1 - 2*f);   rho := 1/sigma
        if i == 1 and mode == FOLD:
            ASSERT |f*(rho + 2) - 1| < eps_id                           # ROW-1 CLOSURE   (2x tripwire)
        w_s := r_inv( rho * r(w_c) )
        ASSERT |r(w_s)/(rho*r(w_c)) - 1| < eps_inv                      # PE CLOSURE      (2x tripwire)
        delta := x(f;beta)*w_c;  D := kappa_KERNEL*delta;  u := D/w_c
        ASSERT delta / h <= eps_lam                                     # THIN LAMINA     (Sec 4 E-g)
        ASSERT |w_c - w_c_prev| / Lambda[i] <= eps_taper                # QUASI-DEVELOPED (Sec 16-c)
        S := S + f*P;   P := P*(1 - 2*f)
        rows.append(i, w_c, beta, f, rho, w_s, delta, D, u, P, S)
        if P <= (1 - F_target)*(1 + 2*sigma_0): N_new := i; break       # TERMINATE ON P (monotone-safe; rev. 5.1)
    if N_new == N: break
    N := floor(0.5*N + 0.5*N_new)                                       # damping: non-monotone W

ASSERT |P + 2*S - (1 + 2*sigma_0)| < eps_id                             # MASS IDENTITY (seed tripwire)
F := 1 - P/(1 + 2*sigma_0)                                              # DELIVERED, any sigma_0 (rev. 5.1)
ASSERT |F - 2/(rho+2)| < eps_id;  if sigma_0 == 0: ASSERT |F - 2*S| < eps_id
else:            REPORT mass defect 2*sigma_0 and F_sched := 1 - P
ASSERT w_s[N] >= w_min                                                  # w_s_out is an OUTPUT

# ---------- INVERSE ROUND-TRIP  (Sec 11, sigma-form) --------------------------------
sig[0] := sigma_0
for i = 1..N:  ASSERT | ((sig[i]-sig[i-1])/(1+2*sig[i])) / f[i] - 1 | < eps_id     # A->B->A, both modes

# ---------- REGIME & WINDOW GATES ---------------------------------------------------
if min(u) < u_ff:   FAIL "outside the calibrated window"
if min(u) < u_lo:   D := BAND[ kappa_lo*delta, kappa_hi*delta ]         # never a scalar
if max(f) > f_soft: WARN
ASSERT max(R_gap * f / R_c) < 0.05                                      # GAP DOMINANCE   (M13)
ASSERT max(Re) < Re_max ; ASSERT Re_gap << 1 ; ASSERT Re_delta << 1     # STOKES GATES    (M14)
REPORT L_e/Lambda  (Sec 16-c: >1 is EXPECTED and benign; the gate is M13, not L_e)

# ---------- RENDER  (Sec 14) --------------------------------------------------------
for i = 1..N:
    Phi := r(w_s[i]) / r(GS)
    if w_s[i] >= w_min:  draw STRAIGHT, width w_s[i], length Lambda[i]
    else:
        ASSERT Lambda[i] >= 2*p
        L_hat := Phi*Lambda[i] / (1 + L_eq/l_stroke)                    # Stokes-corner allowance
        arms  := ceil( Phi / (W_budget/p) ) - 1
        draw FOLD, width GS, traversal L_hat, borrowing `arms` slots
    draw LINK width G length l_g ; draw POST width g_p
at each leg boundary (gap-free turn = ONE segment):
    solve w_s_in, w_s_out_arc from Sec 14.5 using rho(i) of the UPSTREAM rung   # NOT an average

# ---------- OUTPUT ------------------------------------------------------------------
REPORT  N, F, eta*F, D (scalar or BAND), sum(Lambda), max(Phi), max(arms),
        max Re, Re_gap, Re_delta, L_e/Lambda, max gamma_dot via C(alpha) [NEVER 6],
        Delta_P_tot, phi, DIRECT mass defect 2*sigma_0,
        terminal-layer m (predicted with amplitude A, AND measured),
        SENSITIVITY MAP (Sec 16-e),
        transfer-bias band sqrt(c2(beta_cal)/c2(beta(i)))    <-- the dominant unquantified risk
```

---

## 19. Validators

Self-contained. No external data. Every one can fail.

### K — resistance kernel
| id | assertion |
|---|---|
| K1 | Every $B$ call receives $\alpha\ge1$ (sorting enforced at the call site). |
| K2 | $B(\alpha)>0$; monotone increasing on $[1,\infty)$ with $B(1)=0.4217310$. |
| K3 | $r$ strictly decreasing on a log grid spanning $[w_{\min}/10,\,10h]$. |
| K4 | $|r(n_{\max})/r(2n_{\max})-1|<\varepsilon_r$ **at the largest $\alpha$ in the design**. Term count = **odd terms**; state it. |
| K5 | $|r(w;h)/r(h;w)-1|<\varepsilon_r$ (both sorted). |
| K6 | $r\,h^4/\mu=28.454154$ at $w=h$. |
| K7 | $\alpha[1-B(\alpha)]\to0.6302489$; $r\,ab^3/(12\mu)\to1$ as $\alpha\to\infty$. |
| K8 | $|r(r^{-1}(y))/y-1|<\varepsilon_{\rm inv}$ over the design range. |
| **K9** | **Reciprocity:** $|B(1/\alpha)/[\alpha^{-2}B(\alpha)]-1|<\varepsilon_r$ for $\alpha\in\{1.4,2.7,10,47\}$. |
| **K10** | **External anchor:** $24/B(1)=56.908$ (square-duct $f\!\cdot\!\mathrm{Re}$; textbook 56.91). |
| **K11** | **External anchor:** $f\!\cdot\!\mathrm{Re}\to96$ in the parallel-plate limit. |

### E — extraction kernel
| id | assertion |
|---|---|
| E1 | $F_{\rm lam}\to3x^2-2x^3$ as $\beta\to\infty$: $|F_{\rm lam}(0.3,10^4)-0.216|<\varepsilon_r$ **and** $|F_{\rm lam}(0.2,10^6)-0.104|<\varepsilon_r$. |
| E2 | $F_{\rm lam}(\tfrac12,\beta)=\tfrac12$ to machine precision for $\beta\in\{0.1,1,10\}$. |
| E3 | $F_{\rm lam}$ strictly increasing in $x$ on $[0,\tfrac12]$; $A_n(\beta)>0$ ∀$n$. |
| E4 | $c_2(\infty)=3.0000$, $c_2(1)=3.5568$ to $\varepsilon_r$; $c_2$ strictly decreasing. |
| E5 | **Tail:** $|c_2(n_{\max})-c_2(2n_{\max})|/c_2<\varepsilon_r$ **with the analytic tail applied**. |
| E6 | **Cancellation:** $t-\tanh t$ uses the Taylor branch for $t<\tfrac12$; agrees with a high-precision reference at $\beta=0.02$. |
| E7 | $K_{\rm depth}(\beta)=\sqrt{3/c_2}\in(0,1]$, $\to1$ as $\beta\to\infty$. |
| E8 | Inversion seed is $\sqrt{f/c_2(\beta)}$; **assert the linear $f$ is never used as a seed**. |
| E9 | $A_n$, $B$, and $C$ derive from one shared duct series — assert one call site. |
| **E10** | **Shallow scaling:** $\beta\,c_2(\beta)\to\approx1.7$ as $\beta\to0$ (the $F_{\rm lam}\to x$ limit). |
| **E11** | **Thin lamina:** $\delta(i)/h\le\varepsilon_{\rm lam}$ ∀i — the validity condition of the seed, S1′, and the steric picture. |

### C — cubic branch (only when $\mathcal K=g$)
| id | assertion |
|---|---|
| C1 | $g'(x)=6x(1-x)\ge0$ on $[0,1]$. |
| C2 | $g(\tfrac12)=\tfrac12$ to machine precision. |
| C3 | Unique root for every $f\in(0,\tfrac12)$. |
| C4 | $\bigl||x_0/x-1|-\sqrt{f/27}\bigr|<10^{-3}\sqrt{f/27}$ for $f\in[10^{-5},10^{-2}]$ (**an error law, not a constant**); if $\sqrt{f/27}>\varepsilon_{\rm cub}$, the exact root is required. |

### P — pairing
| id | assertion |
|---|---|
| P1 | **The same kernel object at all three call sites** (inversion, schedule, $D$) — by identity, not by value. |
| P2 | $\kappa_{\mathcal K}$ is tagged with the kernel it was fitted with; a mismatch is a **hard fail**. |
| P3 | The transfer-bias band is computed and **reported** for every design. It is **never applied**. |

### M — march, regime, sensitivity
| id | assertion |
|---|---|
| M1 | FOLD: $\sigma(0)=0$ bit-exact. DIRECT: $\sigma_0\ge0$. |
| M2 | **Row-1 closure:** $|f(1)(\rho(1)+2)-1|<\varepsilon_{\rm id}$ (FOLD only — it is an identity at $\sigma_0=0$ and necessarily fails otherwise). |
| M3 | $0<f(i)<\tfrac12$; $\sigma$ strictly increasing (⟺ $\rho$ strictly decreasing). |
| M4 | **Closure:** $|r(w_s(i))/[\rho(i)r(w_c(i))]-1|<\varepsilon_{\rm inv}$ ∀i. |
| M5 | **Mass identity:** $|P(N)+2S(N)-(1+2\sigma_0)|<\varepsilon_{\rm id}$. |
| M6 | $F_{\rm sched}=1-P(N)$ and the delivered $F=1-P(N)/(1+2\sigma_0)$ with $|F-2/(\rho(N)+2)|<\varepsilon_{\rm id}$; and if $\sigma_0=0$, $|F-2S(N)|<\varepsilon_{\rm id}$ (rev. 5.1). |
| M7 | **Round-trip A→B→A** ($\sigma$-form): $|f_B(i)/f(i)-1|<\varepsilon_{\rm id}$ ∀i, **in both modes**. |
| M8 | Under S1: $|\delta(i)/\delta(1)-1|<\varepsilon_{\rm id}$ ∀i, at any aspect ratio. |
| M9 | **Invariance sweeps** (§17), $\lambda\in\{\tfrac12,2\}$: every "—" invariant; every scaling cell exact; **$R_s=r({\rm GS})L$ invariant under GS**. |
| M10 | **Damped fixed point** converges in $\le K$ iterations, for a non-monotone $\mathcal W$ test case. |
| M11 | **Grid refinement:** doubling $\mathcal W$'s resolution moves $F$, $D$, $\sum\Lambda$ by $<\varepsilon_{\rm id}$. |
| **M12** | **§13 taper, general form:** $|N-N_\Psi|\le1$ **and** $\bigl|F-\bigl[1-\Psi(w_c(N))/\Psi(w_c(0))\bigr]\bigr|<10^{-3}$. *(The old reference $1-(w_c(N)/w_c(0))^2$ is the **cubic special case** and fails by 13× in the $F_{\rm lam}$ branch.)* |
| M13 | **Gap-dominance gate:** $\max_iR_{\rm gap}f(i)/R_c(i)<0.05$. **Fail**, do not warn. |
| **M14** | **Stokes gates:** $\max\mathrm{Re}<\mathrm{Re}_{\max}$; $\mathrm{Re}_{\rm gap}\ll1$; $\mathrm{Re}_\delta\ll1$. |
| **M15** | **Entrance length:** report $L_e/\Lambda$ with the Stokes-floor formula. $L_e>\Lambda$ is **expected and benign** (the perturbation is $O(f)$, so the resistance error is $O(f^2)$). Do not gate on $L_e$; gate on M13. |
| **M16** | **Quasi-developed taper:** $|\mathrm dw_c/\mathrm dx|\le\varepsilon_{\rm taper}$ ∀i. |
| **M17** | **Sensitivity map** (§16-e) is computed and emitted with every design. |

### D — validity window
| id | assertion |
|---|---|
| D1 | $u(i)\ge u_{\rm ff}$ ∀i, else **FAIL**. |
| D2 | $u(i)\ge u_{\rm lo}$ ∀i, else $D$ is a **band**, never a scalar. |
| D3 | $f(i)\le f_{\rm soft}$ ∀i, else WARN. |
| D4 | $\delta(i)\le w_c(i)/2$ ∀i. |

### G — geometry
| id | assertion |
|---|---|
| G1 | **Resistance DRC** using the **measured** traversal: $|r({\rm GS})[\hat L(i)+n_{\rm turns}L_{\rm eq}]/(\rho(i)r(w_c(i))\Lambda(i))-1|<\varepsilon_{\rm geo}$; straight rows analogously. |
| G2 | **Series topology:** shortest-path$(\text{link}_i\to\text{link}_{i+1})\ge(1-\varepsilon_{\rm geo})L(i)$. |
| G3 | **No stub:** $\deg(\text{link}_i)\ge2$. |
| G4 | **Constant cross-section** along every traversal. |
| G5 | **Packing:** $\Phi(i)\le W/p$, or $\mathrm{arms}(i)$ slots allocated and pairwise non-overlapping. |
| G6 | **Global area:** $\sum_i\Phi\Lambda\le(W/p)\sum_i\Lambda$ (**necessary, not sufficient**). |
| G7 | **Min feature:** widths $\ge w_{\min}$; walls $\ge F_{\min}$; inner corners $\ge R_{\min}$. |
| G8 | **Fold pitch:** $\Lambda(i)\ge2p$, $p=\mathrm{GS}+\max(F_{\min},2R_{\min})$, $R_{\rm bend}=p/2$. |
| G9 | **Length budget:** $\sum_{\rm leg_k}\Lambda\le\ell_{\rm leg}$ ∀k; total $\le n_{\rm legs}\ell_{\rm leg}$. |
| G10 | **Turn rebalance** with the **upstream** $\rho(i)$. **Precondition for all of §15** (N8). |
| G11 | **Handover:** folded **iff** $w_s(i)<w_{\min}$. No other threshold. |
| G12 | **Bend accounting:** $n_{\rm turns}(i)$ counted from the drawn polyline, not assumed. |
| **G13** | **Run-out match:** $\bigl|\,(R_{\rm out}^{s}/R_{\rm out}^{c})\,/\,\rho(N)-1\bigr|<\varepsilon_{\rm geo}$. **The most sensitive un-gated input in v4: a $3\times$ mismatch collapses $F$ from $0.519$ to $\approx0.27$.** |
| **G14** | **Routability:** the drawn geometry is planar and non-self-intersecting. *(Not implied by G5/G6 — those bound area only. This must be checked on the drawing.)* |

### N — network
| id | assertion |
|---|---|
| N1 | Laplacian symmetric, non-singular after grounding; residual $<\varepsilon_{\rm net}$. |
| N2 | **Bulk:** $\max_{i\le N-m}|f^{\rm em}(i)/f(i)-1|<\varepsilon_{\rm net}$. |
| N3 | **Terminal:** $|F^{\rm em}/F-1|<\varepsilon_{\rm net}$, with $F$ the delivered fraction (rev. 5.1). |
| N4 | **$R_{\rm gap}$ insensitivity:** sweep $[10^{-1},10^{+1}]\times$; bulk $f^{\rm em}$ invariant. If not, M13 was violated. |
| N5 | **Structural BC:** $Q_s^{\rm em}(0)=0$ emerges un-imposed. |
| N6 | §13 taper: $\max_i|Q_{\rm gap}(i)/Q_{\rm gap}(1)-1|<10^{-3}$ (S1′) or $<10^{-2}$ (S1). |
| N7 | **Measured terminal layer** $m_{\rm meas}=N-i^\ast$ within a factor 2 of the predicted $m$, and $m_{\rm meas}\le m_{\rm budget}$. |
| N8 | **Symmetry precondition:** G10 passes. Otherwise §15's two-rail folding is invalid and N1–N7 are void. |
| **N9** | **Boundary amplitude:** $A=|f^{\rm em}(N)/f(N)-1|$ is measured and used in $m=\lceil\lambda\ln(A/\varepsilon_{\rm net})\rceil$. $\lambda$ is taken as $\max$ over $[N-m_{\rm budget},N]$ (**$q$ is non-monotone**). |

### S — shear (only if §16-b is used)
| id | assertion |
|---|---|
| S1 | $C(\alpha)\to6$ as $\alpha\to\infty$ ($C(10^6)=6.000004$). |
| S2 | $C(1)=9.6078$; $C(2)=8.133$ to $\varepsilon_r$. |
| S3 | $C$ uses the **same** $B(\alpha)$ as §3 and §4 — bitwise identical call site. |
| S4 | Constant-shear objective: $\max_i|\dot\gamma(i)/\dot\gamma(1)-1|<\varepsilon_{\rm id}$. **Do not assert on the §13 taper — it fails by design.** |

### R — reproducibility
| id | assertion |
|---|---|
| R1 | Determinism: hash(inputs) ↦ hash(outputs), stable across runs and platforms. |
| R2 | Every $\varepsilon$ of §2.4 is an explicit input. **No literal tolerance in the code.** |
| R3 | Every §2.3 constant is an explicit, tagged input, echoed in the output header — **including $\beta_{\rm cal}$, $\phi_{\rm cal}$, and the kernel $\kappa$ was fitted with.** |

---

## 20. What no validator here can check

These are **inputs**. Every equation above is exact *given* them; §19 passes for any values you supply.

| constant / effect | fixes | consequence if wrong |
|---|---|---|
| **$\kappa$** | $D=\kappa\delta$ | $D$ wrong by that factor; and $N\propto\kappa^2$ (§16-e). Geometry, $F$, $\Phi$ unaffected. |
| **$\beta_{\rm cal}$** | which kernel $\kappa$ implicitly carries | **Dominant risk.** A **27 % swing in $N$** and a **0–16 % bias in $D$** (§5), from the depth of an array that was never published. **Measure it, or re-fit $\kappa$ at the design depth.** |
| **$\phi_{\rm cal}$** | the wall porosity at which $\kappa$ was fitted | The kernels assume a **solid** wall; the real wall is $\phi\approx0.4$ **open**, with an effective-slip offset $O(\delta)$. §6.2 proves the residual $\kappa/2=2.46$ must live here, in the flat-lamina idealisation, or in finite-particle displacement — **and nowhere else**. Not controlled by gap dominance (porosity is geometric, present at zero gap flow). |
| **flat-lamina idealisation** | $\delta$ ↔ the actually-skimmed fluid | The true dividing streamsurface is $z$-dependent (thicker near top/bottom). Absorbed into $\kappa$, **$\beta$-dependently** — so it is *not* removed by switching to $F_{\rm lam}$ (§5). |
| $u_{\rm lo},\ u_{\rm ff}$ | the $\kappa$-constant window | $D$ wrong **silently** outside the true window. **The only failure mode the algorithm cannot detect at all.** |
| $\eta$ | $F_{\rm real}=\eta F$ | Realized yield misses; model $F$ unaffected. **A flow-split correction — it does not carry any lamina model gap.** |
| **diffusion (Péclet)** | out of the model entirely | For calibration beads the transverse spread over transit is $1$–$4\times\delta$ — *not* small. But §6.2 shows the cutoff is a **hard steric wall constraint**, and diffusion cannot change a sphere's radius: it **broadens the transition, it does not move the cutoff** (deterministic migration beats diffusion $6.7\times$ per rung). Becomes first-order for sub-micron particles or long transits. |
| $f_{\rm soft}$ | a warning threshold | nothing. |
| $\gamma_{\rm GS},\mathrm{AR}_{\max},F_{\min},R_{\min}$ | fabricability | the mask may be unbuildable. |
| $\ell_g,\ c_{\rm ent}$ | $R_{\rm gap}$ | §15 only, and $m$. Design insensitive while M13 passes. |
| $L_{\rm eq}$ | bend correction | folded rows 2–6 % under-resistant; a localized inlet-side closure violation. |
| **routability** | global planarity | G5/G6 bound area only. **A feasible packing is not proven to exist; the implementer must draw it** (G14). |

> **$\kappa$, $\beta_{\rm cal}$, $\phi_{\rm cal}$ and $u_{\rm lo}$ are the entire external dependency of the selectivity axis.** The resistance axis and the yield axis are closed.

---

## 21. Errata

### Fixed in v5

| # | defect (v4) | why wrong | fix |
|---|---|---|---|
| **1** | §13: $F=1-(w_c(N)/w_c(0))^2$, and **M12 built on it** | the cubic special case. $Q_c\propto\Psi=w_c^2/c_2$, not $w_c^2$. Omits $c_2(\beta_0)/c_2(\beta_N)=1.0265$ ⟹ $F=0.51929$ vs true $0.50654$ — **off by $0.0126$, 13× the M12 tolerance**. My own $N$ formula carried $c_2$; the $F$ formula silently dropped it. | **§13 rebuilt on one potential $\Psi$** (strictly increasing, $\mathrm d\ln\Psi/\mathrm d\ln w\in[1,2]$). $w_c(i)$, $P$, $N$, $F$ and the design inversion all read off it — **the error becomes structurally impossible**. M12 repointed. |
| **2** | §15: "$q$ is smallest at $i=N$, since $\rho$ falls monotonically" | incomplete: $R_c$ **grows** as $w_c$ narrows and wins near the outlet. $q$ has an interior minimum ($0.501$ at 86 % of span vs $q(N)=0.510$); $\mathrm dq/\mathrm di>0$ at $N$. | $\lambda$ = max over the terminal window; **and measure $m$** (N7/N9). |
| **3** | §15: $m=\lceil\lambda\ln(1/\varepsilon)\rceil$ | a **unit-amplitude** yardstick. The decay is $A z_-^{N-i}$; $A$ is set by run-out mismatch (matched $A\approx0.05$; $3\times$ mismatch $A\approx5\times10^2$). | $m=\lceil\lambda\ln(A/\varepsilon_{\rm net})\rceil$ with $A$ measured (N9). |
| **4** | §17: GS row, $L$ column = "same" | $L=\Phi\Lambda$ scales with $\Phi$. What is invariant is $R_s=r({\rm GS})L$. | column added; M9 asserts it. |
| **5** | §2.1 listed $\sqrt{f/27}$ among **constants** | it is an $f$-dependent **error law**. | moved to §2.2 "Error laws". |
| **6** | run-out condition stated (§10) but **never gated** | a $3\times$ mismatch collapses $F$ from $0.519$ to $\approx0.27$ — **the most sensitive un-gated input in the document**. | **G13**. |
| **7** | Re computed, never bounded; no entrance-length or taper-gradient condition | the Stokes picture of §6.2 requires $\mathrm{Re}_{\rm gap}\ll1$; the kernels require quasi-developed flow | **M14, M15, M16**; correct Stokes-floor $L_e$ (§16-c). |
| **8** | no thin-lamina condition | the seed, S1′, and the steric picture all need $\delta\ll h$ | **E11**, $\delta/h\le\varepsilon_{\rm lam}$. |
| **9** | no sensitivity map | cheap, and it reorders process priorities | **§16-e, M17**. |
| **10** | R-b table's term count unstated | "50 terms" is ambiguous (odd terms vs $n\le50$ differ ~15×) | stated: **50 odd terms, $n\le99$**. |
| **11** | "exact" kernel unqualified | $F_{\rm lam}$ is the exact flux of a **flat** slab in an **unperturbed** duct; the skimmed fluid is neither | §4 box; §20 rows. |

### Added — new results

| | |
|---|---|
| **Reciprocity identity** $B(1/\alpha)=\alpha^{-2}B(\alpha)$ | turns the $w\!\leftrightarrow\!h$ symmetry from an assertion into a **series identity**, and explains the unsorted digit loss quantitatively ($\approx2\log_{10}\alpha$ digits + $\alpha^2$ truncation amplification). **K9.** |
| **External anchors** $24/B(1)=56.908$ (textbook 56.91) and $\to96$ | validate $B$ against tabulated duct data **with no reference to this document**. **K10, K11.** |
| **Transport lemma (§6.2)** | streamline migration reproduces $F$ **exactly, with zero diffusion** ($\varphi_{\rm crit}=F/2$, verified to 6 decimals). Proves the device needs no cross-streamline transport; establishes $\kappa_{\rm steric}=2$; and **pins the entire $\kappa/2=2.46$ residual to exactly three named model gaps.** |
| **Taper potential $\Psi$** | one monotone object generating $w_c(i)$, $P(i)$, $N$, $F$, and the design inversion. |

### Rejected — and why

| claim | verdict |
|---|---|
| $L_e=0.05\,\mathrm{Re}\,D_h$ as an entrance-length gate | ❌ **Vanishes as $\mathrm{Re}\to0$** — wrong formula in this regime. The Stokes floor is $L_e\approx0.62\,D_h$, giving $L_e\approx70\ \mu$m $>\Lambda=49\ \mu$m. The correct response is not a gate but the **argument**: the per-rung perturbation is $O(f)$, so the resistance error is $O(f^2)\sim10^{-8}$. Report $L_e/\Lambda$; gate on M13. |
| "$\eta=0.903$ carries the flat-lamina model gap" | ❌ **Mis-assigned.** $\eta$ is a **yield** correction (a 3-D effect on the flow split). The lamina model gap lives in **$\kappa$** (the cutoff prefactor). They are independent axes — that is the whole content of the decoupling lemma (§17). |
| "wall segmentation ... controlled by gap dominance" | ❌ **Wrong mechanism.** Gap dominance bounds the *flow* through the links. Porosity is *geometric* and is present at **zero** gap flow. The effect is real ($\phi\approx0.4$, slip $O(\delta)$) and is now §20's $\phi_{\rm cal}$ — but nothing in M13 controls it. |
| "Péclet ... decisive for separation performance" | ❌ **Overstated.** §6.2 shows the cutoff is a **hard steric wall constraint**; diffusion cannot change a sphere's radius. Deterministic migration beats diffusion $6.7\times$ per rung. It **broadens the transition; it does not move the cutoff.** Named in §20 — not promoted to a mechanism. |
| "evaluate $\lambda$ at $\arg\min q$" | ⚠️ **Half right.** $q$ *is* non-monotone (adopted). But the $\arg\min$ sits at 86 % of span — far outside the terminal layer, which is governed by the *local* $q$ near $i=N$. Correct: $\max\lambda$ over the **terminal window**, and measure. *(Their claimed $m(10^{-3})$: 10→11 does not reproduce; I get 10 either way.)* |

### Fixed in revision 5.1 (2026-10-23)

| # | defect (v5) | why wrong | fix |
|---|---|---|---|
| **12** | §10: "the physical terminal fraction, for any $\sigma_0$: $F=1-P(N)$", with the stop on $P(N)\le1-F_{\rm target}$ (§12, §18), $\rho^{\rm out}$ from the same $F$ (§10, §12), and M6 | A DIRECT start draws its first side width at the printable start width with a **sealed** end, so no side flow enters it ($Q_s(0)=0$ is structural, §15), whereas the closure sizes every side channel for $\sigma(i)=S(i)/P(i)$ with $S(0)=\sigma_0>0$. With the pressures equalized on every segment from row 1 on, the sealed device's flows satisfy $Q_s(i)/Q_c(i)=\sigma(i)$ and $Q_c(i)+2Q_s(i)=Q_{\rm in}$, so they are the schedule's flows scaled by $1/(1+2\sigma_0)$, row 1 taking up the difference (on a drawn network with gap drops, the first few entrance rows). The delivered fraction is therefore $2\sigma(N)/(1+2\sigma(N))=2/(\rho(N)+2)=1-P(N)/(1+2\sigma_0)$, the form that v5 said "carries the seed", and not $1-P(N)$. The drawn network confirms it: on a Gifford Fig. 4 device with a 20 µm floor start ($\sigma_0=2.009\times10^{-2}$), the network of the drawn mask delivers $0.90647$ against a target of $0.90644$ when the march stops on the delivered value, where $1-P(N)$ would be $0.9027$, and the final side width is then $299.9$ µm against the paper's printed $300$ µm. The design code's own $\sigma_0$ compensation, which iterates its $F$ target until the network's $F^{\rm em}$ meets the request, treats the network's value as the delivered one in the same sense. | §10 separates $F_{\rm sched}=1-P(N)$ (the schedule's fraction, which the mass identity, M5 and §6.2's lemma concern, unchanged) from the delivered $F=2/(\rho(N)+2)$; §12 and §18 stop on $P(i)\le(1-F_{\rm target})(1+2\sigma_0)$; $\rho^{\rm out}=2(1-F)/F$ for any $\sigma_0$; M6 and N3 check the delivered $F$. A FOLD start ($\sigma_0=0$) is unchanged. Found by the independent re-audit of opencif 1.0.0's alignment matrix (its rows AA16 and S3; opencif's erratum record D108). |
