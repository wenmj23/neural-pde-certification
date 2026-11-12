# 中期报告定稿正文存档（2026-11-12）

> 用户已确认报告完成，尚无提交回执。本文件从最终 8 页 Word/PDF 对应正文提取，仅供研究接续；省去学生/导师姓名、邮箱及个人页眉，不复制排版。原始 Word/PDF 保留在用户文件中。本文件不是重新排版的提交版。11 月 12 日补充实验见单独研究报告，不自动视为已写入本文。

A Posteriori Error Certification for Neural PDE Solvers

## 1. Introduction and Research Objectives

### 1.1 Background and motivation

Partial differential equations (PDEs) describe many continuous physical systems. Physics-informed neural networks (PINNs) approximate a solution with a neural network trained using the governing equation and boundary data [1]. Neural solvers can produce useful approximations, but a small training loss only describes the objective evaluated during optimization. It does not, by itself, guarantee a small error at every point of the physical domain.

This project studies a posteriori certification: once a candidate solution has been fixed, a separate procedure uses the equation and the candidate to assess its error. The focus is a bound on the maximum absolute error over a continuous domain. This differs from improving the PINN training algorithm and from evaluating the error at a finite collection of test points.

### 1.2 Scope of the midterm work

The opening report proposed starting with a simple elliptic PDE before considering more general operators or neural architectures. The midterm work makes this first stage precise: the current main problem is the one-dimensional Poisson equation on [0, 1], with zero Dirichlet boundary conditions. The candidate is a frozen, sufficiently smooth neural-network function. The implementation and experiments focus on quadratic and cubic residual reconstruction; higher dimensions and neural operators remain possible future directions.

−u*″(x) = f(x),   0 < x < 1;       u*(0) = u*(1) = 0.                         (1)

Given a candidate v, the quantity of interest is E = ‖v − u*‖∞. The certifier seeks a computable interval [L, U] that contains E, subject to the stated mathematical and computational assumptions. For an absolute tolerance τ, U ≤ τ establishes acceptance, L > τ establishes exceedance, and an interval crossing τ is inconclusive. The ratio U/L describes enclosure tightness when L > 0; it is not the relative error of the neural solution.

### 1.3 Progress since the opening report

The opening report identified a solver–certifier–validation pipeline as the minimum outcome and left the exact error norm and model problem open. The present work has made two deliberate choices: use a one-dimensional Poisson problem with homogeneous Dirichlet boundary conditions and target the supremum norm of the solution error. These choices make the stability relation and the boundary assumptions explicit and permit numerical checks against problems with known solutions. They do not establish a certificate for the broader class of PDEs originally discussed as a possible extension.

Work to date spans more than a single set of numerical runs. It includes formulating the error equation; examining a global residual estimate; constructing quadratic and cubic residual reconstructions; implementing continuous range and remainder calculations; testing frozen sine and polynomial networks; comparing grids, tolerance decisions and selected timings; constructing artificial controls; and checking result persistence and reproducibility. These activities have different evidential roles: a proof step, a source-code check, a timed run and an empirical test should not be counted as interchangeable demonstrations of correctness.

## 2. Comparison and Selection of Research Approaches

### 2.1 Point checks and global residual bounds

Pointwise comparisons are easy to compute when a reference solution is available and are useful for debugging or illustration. However, even a dense finite grid may miss an error peak between adjacent points. Therefore the maximum observed on a grid is not treated as a continuous-domain certificate.

For the zero-boundary Poisson problem, the residual can be connected to the solution error through the Green function. If a validated global residual bound is available, the estimate E ≤ ‖r‖∞/8 follows. It gives a simple baseline, but a global maximum and an absolute-value estimate can lose spatial information and signed cancellation. In the experiments, a spatial residual baseline is also considered; its performance depends on the particular implementation and mesh.

### 2.2 Signed residual reconstruction

The selected approach reconstructs the signed residual on each mesh interval, then controls the difference between the true residual and this reconstruction. The reconstructed source is propagated through an auxiliary Poisson problem, while a nonnegative envelope accounts for the remainder. This separates a signed component, which can preserve cancellation, from a rigorously bounded remainder. The objective is a useful error enclosure rather than a claim that higher reconstruction order is universally superior.

### 2.3 Related certification and uncertainty methods

Residual and derivative verification for PINNs has been studied by Eiras et al. [2], whose work targets continuous-domain residual bounds; a residual bound is not automatically a bound on the solution error. Hillebrecht and Unger investigate rigorous a posteriori error bounds in PDE-defined PINNs [3], while Berrone et al. analyze residual-type estimators for variational PINNs in an energy norm [4]. Their settings and assumptions differ from the present one-dimensional maximum-error target. Liu et al. [5] also study residual-based bounds for neural approximations to differential equations, including linear ODEs. The present study focuses on a zero-boundary Poisson problem and a signed piecewise-polynomial reconstruction with a propagated remainder.

Conformal prediction is another useful perspective on PINN uncertainty. Yu, Ho and Wang [6] study calibrated prediction intervals under statistical assumptions. Such coverage is meaningful for the corresponding random prediction problem but does not by itself establish a deterministic worst-case error guarantee over every point in a continuous domain. It is therefore retained as a comparison, while residual-based maximum-error certification remains the main route.

## 3. Mathematical Formulation and Certificate Construction

### 3.1 Error equation and boundary assumptions

Write e = v − u* and r = −v″ − f. If the candidate satisfies the same boundary conditions as the exact solution, then the error solves

−e″ = r,       e(0) = e(1) = 0.                                                      (2)

The experiments commonly use v(x) = x(1 − x)N(x), with a fully connected tanh network N, so that the zero boundary conditions are imposed by construction. The neural weights are frozen before certification. The present mathematical target treats the stored floating-point weights as exact real constants; it does not certify every rounding error in a separate floating-point inference program.

Zero boundary values for the error are essential to this simplified formulation; they hold when the candidate matches the prescribed boundary values. If the candidate violates the boundary conditions, an additional boundary contribution must be controlled. A small interior residual alone would not exclude such an error. Nonzero and uncertain boundary treatments have been explored as research prototypes, but they are outside the default public command-line scope.

### 3.2 From residual to error

The zero-boundary solution of the one-dimensional error equation can be written using the Green kernel G(x,s) = min(x,s)[1 − max(x,s)]. For 0 ≤ x,s ≤ 1 this kernel is nonnegative, and integration gives ∫₀¹ G(x,s) ds = x(1−x)/2. Its maximum over x is 1/8. Consequently, an independently valid upper bound R ≥ sup|r| implies E ≤ R/8. This argument shows precisely where the equation, the boundary conditions and a continuous residual bound enter the estimate. An estimate based only on residual values at finitely many points would not justify replacing R by their maximum.

If the candidate has nonzero boundary mismatch, the representation also contains a linear boundary contribution. It is therefore unsafe to use the zero-boundary formula unchanged merely because the interior residual is small. The present main experiments avoid that issue with hard-enforced zero boundary data; boundary-error treatment is kept as a separate research task.

### 3.3 Reconstruction and error enclosure

Partition [0, 1] into mesh intervals and let p be a signed piecewise-polynomial reconstruction of r. Define z by −z″ = p with zero boundary conditions. Suppose a nonnegative envelope ρ bounds |r − p| on every interval, and a computable D_U bounds the maximum effect of that envelope. If the continuous maximum of |z| is enclosed by Z_L and Z_U, then the triangle inequality yields

max(0, Z_L − D_U)  ≤  ‖v − u*‖∞  ≤  Z_U + D_U.                   (3)

The computational task is consequently divided into a range enclosure for the reconstructed component and a validated upper bound for the propagated reconstruction remainder. In one dimension the Green kernel for the zero-boundary problem is nonnegative, which supports propagation of a nonnegative source envelope. The argument requires the envelope and the bounds used to construct it to be valid over entire mesh cells, not merely at selected points.

### 3.4 Reconstruction order and range computation

The prototype implements quadratic and cubic residual reconstructions. Their interpolation remainders depend on validated bounds for the third and fourth derivatives of the residual, respectively. With cell width h and suitable derivative bounds K₃ and K₄, the corresponding constant remainder estimates scale as K₃h³ and K₄h⁴. Shape and moment information can further reduce the propagated bound in the implemented construction. These formulas require sufficient smoothness and correct derivative enclosures.

Polynomial range enclosure uses Bernstein coefficients and subdivision; refining that range calculation is distinct from increasing the number of PDE mesh intervals. Both operations may improve an enclosure but also add work. The current comparison therefore records the mesh, reconstruction order and computational budget together, rather than interpreting one isolated number as an overall performance result.

## 4. Work Completed and Evaluation Setup

### 4.1 Implementation and experimental cases

The project has developed a public command-line workflow for selected sine-forcing cases and research prototypes for polynomial forcing, boundary variants and additional mesh strategies. The workflow produces readable reports and machine-readable records. A frozen network is a fixed input to the certifier; repeated meshes, tolerances or timings for the same network are not additional independently trained networks.

Experiments include sine forcing f(x) = (kπ)² sin(kπx), with known solution u*(x) = sin(kπx) for positive integers k; polynomial forcing with known analytic solutions, and deliberately constructed error functions for stress testing. The latter are explicitly distinguished from trained PINNs. Where feasible, results are compared with converged continuous interval references or exact analytic values. These comparisons help check specific cases but do not constitute an independent proof of every implementation component.

### 4.2 Experiments and record provenance

The experimental progression distinguishes different kinds of test. Previously trained shallow and deep sine networks exercise the public validation workflow. Polynomial-forcing networks broaden the source terms while retaining problems for which explicit solutions or converged numerical interval references can be constructed. Manufactured sine waves and polynomial errors provide controls with analytically known maxima; they are not represented as additional trained PINNs. The project also retains difficult or poorly trained networks so that a narrow enclosure is not confused with an accurate candidate.

Each experimental comparison is meaningful only together with its configuration. Relevant choices include the frozen network and source, reconstruction degree, mesh interval count, numerical precision, stopping criterion and timing scope. A result obtained by changing the mesh or tolerance on one network is another configuration, not another independently trained model. The project keeps selected configuration files, result records and post-exit integrity checks [7], although not all historical raw outputs and dependencies are currently available from the public repository.

### 4.3 Comparison criteria

The main criteria are: (i) whether the computed enclosure contains the selected reference interval or exact value; (ii) whether the enclosure is tight enough for the chosen U/L target; (iii) whether a tolerance is accepted, exceeded or remains undecided; and (iv) the measured computational cost under a stated timing protocol. The reference midpoint is used only as an approximate denominator in some comparative statistics. A narrow enclosure does not imply that the neural approximation itself is accurate.

## 5. Midterm Results

### 5.1 Three possible tolerance decisions

A representative shallow sine case illustrates why an inconclusive result must be preserved. With a 16-cell budget, a tolerance of 10⁻³ is above the reported upper bound and is accepted. With a tolerance of 10⁻⁴, a finer run can establish exceedance, whereas the 16-cell enclosure crosses the threshold and remains undecided. Values below are rounded for presentation.

| Setting | Error enclosure [L, U] | Decision |
| --- | --- | --- |
| τ = 10⁻³, 16 cells | [7.1293 × 10⁻⁵, 1.3631 × 10⁻⁴] | Within tolerance |
| τ = 10⁻⁴, 32 cells | [1.0079 × 10⁻⁴, 1.0583 × 10⁻⁴] | Exceeds tolerance |
| τ = 10⁻⁴, 16 cells | [7.1293 × 10⁻⁵, 1.3631 × 10⁻⁴] | Undetermined |

These are results for one previously studied candidate and selected budgets. A budget-limited “undetermined” outcome is not evidence that the method has proved failure or success. Earlier public examples were run in a Windows environment; the later research replays were executed in a Linux environment. The table uses quadratic reconstruction; its interval endpoints are rounded outwards for presentation.

### 5.2 Matched comparison on polynomial-forcing networks

A selected comparison used 12 trained polynomial-forcing networks, evaluating four concrete methods on matched 16-, 32- and 64-cell meshes. The reported 144 enclosures contained their corresponding converged continuous interval references. At 64 cells, the median of the upper bound divided by the midpoint of the reference interval is summarized below. The midpoint is not claimed to be the exact error.

| Concrete method | Median U / reference midpoint |
| --- | --- |
| Global residual baseline | 149.56 |
| Spatial residual baseline | 20.18 |
| Signed quadratic reconstruction | 1.019 |
| Signed cubic reconstruction | 1.0003 |

For the global baseline, a bound on the maximum absolute residual is computed from the quadratic residual data and divided by 8. The value 149.56 is an upper-bound-to-reference ratio, not the network error. Using one global maximum discards spatial variation and sign cancellation, so this bound can be conservative. In this finite set, signed reconstruction yields substantially tighter upper bounds than these particular residual baselines. For absolute tolerance 10⁻⁵, the reference classified three of the twelve networks as within tolerance, and the quadratic and cubic methods each proved acceptance for those three. The baseline methods did not prove acceptance in that setting; this fact alone does not prove that any network exceeded the tolerance.

Cost depends on what is held fixed. For fixed meshes on the same machine, the reported median cubic-to-quadratic time ratios were approximately 1.39 at 16 cells and 1.36 at 64 cells. A separate time-to-target experiment uses different stopping conditions, so its results cannot be combined with fixed-mesh timings to claim a universal speed ordering.

### 5.3 Sine frequency and weak approximations

Another experiment used 12 networks trained for sine frequencies k = 1, 3 and 5, with two network depths and two seeds. Five meshes and two reconstruction orders produced 120 enclosures; all contained their corresponding converged references. Two k = 5, seed 502 networks showed poor training behavior. For one reported 64-cell cubic case, the enclosure was approximately [4.40686, 4.40711], illustrating that a narrow certificate can accurately describe a poor approximation. The certifier and the neural solver have distinct roles.

### 5.4 Mesh choice and diagnostic checks

For the studied sine networks, selected equal-budget comparisons did not show a consistent advantage for the implemented adaptive mesh. By contrast, adaptive refinement gave narrower enclosures for the chosen artificial narrow-peak examples. The latter are constructed stress tests, not trained PINNs. These results motivate reporting positive and negative cases together and accounting for the cost of refinement.

Further numerical diagnostics compared residual derivatives against a separate high-precision point calculation. In a selected polynomial-network test, 864 checks were consistent with the computed intervals; the count includes repeated split settings and is not a count of independent networks. Deliberate omission of a source term also demonstrated that checking only third- and fourth-order derivatives can miss an error in the original residual when the source is constant or linear: the relevant higher derivatives are then zero. These are deliberately introduced test faults, not observed production failures. Point checks cannot rule out errors between the sampled locations.

### 5.5 Fixed-tolerance sequential evaluation

A subsequent online study used two common absolute tolerances for an existing set of twelve sine PINNs rather than selecting a different threshold from each reference error. Two policies started at 16 or 64 mesh intervals and then increased the mesh as needed, stopping when the current enclosure proved either within or above tolerance. Across both reconstruction orders and both starting policies, the final classification for tolerance 10⁻³ was seven within and five above; for tolerance 10⁻² it was eight within and four above. These classifications describe the same twelve networks under repeated settings, not additional trained networks.

A second experiment applied the same early-stopping idea to twenty-four previously studied polynomial-forcing networks, with two common tolerances. At 10⁻⁵, the selected cases classified seventeen within and seven above; at 10⁻⁶ they classified three within and twenty-one above. Both tested starting policies obtained the same final classification under the stated 64-cell maximum. A coarser initial mesh can avoid some later work when the decision is easy, but an undecided initial evaluation may add work. This is a budget and stopping-rule observation, not evidence that refining the certifier improves the already frozen neural solution.

For the sine follow-up at 10⁻³, each online setting received one warmup and three randomized repeated measurements. For each network the paired start-16/start-64 in-process time ratios were summarized by a median; the median across the twelve networks was approximately 0.72 for quadratic and 0.45 for cubic reconstruction. The median paired time ratio was below one for six of twelve networks in the quadratic comparison and nine of twelve in the cubic comparison. The timing includes candidate construction, mesh updates, certification and decision arithmetic, but excludes process startup and final output. These are single-machine results for known cases, not a claim that either starting mesh is generally optimal.

### 5.6 Independent numerical checks and an exact sampling counterexample

To check selected polynomial-source cases from another direction, a separate high-precision mpmath forward evaluation was compared with an explicit zero-boundary polynomial solution. Twelve previously studied networks were sampled on nested 64- and 256-interval grids. For both grids, the selected maximum point value was recomputed at two decimal precisions; all twenty-four selected point checks were stable under the chosen numerical criterion and lay below the corresponding saved reference upper bound. Each also lay below both the quadratic and cubic 64-cell certificate upper bounds. These are finite-point, finite-precision diagnostics, not an independent verification of the entire continuous-domain implementation.

The sampling limitation can be made explicit. For any positive integer m, let e(x) = 10⁻⁴ sin(2mπx). At every node x=i/m the error is zero, while at x=1/(4m) it reaches 10⁻⁴. Thus even a grid with zero observed error at all its nodes can miss a nonzero maximum between them. This constructed function is a mathematical control, not a trained network and not an observed failure of the certifier. Its purpose is to motivate why the continuous residual and remainder bounds matter.

Exact rational polynomial arithmetic was also used to check that the explicit polynomial reference formulas satisfy the selected forcing equations and zero boundary conditions. Separately, interval-versus-point comparisons and intentionally introduced source-omission controls exposed testing blind spots: for constant or linear forcing, its third and fourth derivatives vanish, so higher-derivative checks alone can agree even when the original residual formula omits the source. Such injected faults test the diagnostic method; they should not be described as faults found in production code.

## 6. Current Limitations and Reliability

The mathematical enclosure is conditional on the PDE, boundary data, smoothness, interpolation and Green-function argument, as well as correct interval arithmetic and software implementation. The existing numerical checks and replayed cases provide evidence for specific inputs; they are not a formal verification or an independent audit of the complete certifier. The present evidence does not establish correctness for all neural networks, PDEs, boundary conditions or numerical platforms.

Some comparisons reuse previously studied networks or were designed after earlier outcomes were observed; they should not be described as blind tests. Small same-machine timing samples are sensitive to environment and measurement scope. Published project materials contain selected scripts, inputs and records, while some older experimental details and dependencies remain absent. Thus not every historical study can currently be replayed from the repository alone.

Earlier runs also exhibited mismatches between reported totals and saved details. The original cause has not been established. Post-process record-count and hash checks, along with staging or fresh output directories, help detect certain inconsistencies but do not prove that the storage issue has been eliminated. The project must continue to preserve unsuccessful or inconclusive outcomes and distinguish data integrity from mathematical correctness.

## 7. Remaining Work and Tentative Plan

The immediate priority is to make the present one-dimensional result understandable and auditable. The remaining tasks are to review the error equation, boundary assumptions, derivative enclosures and remainder propagation, and to select a small number of representative experiments with traceable inputs and outputs.

Before extending the problem class, the project should compare quadratic and cubic reconstruction under clearly stated budgets, report inconclusive cases, and reproduce selected public examples on another platform if available. A meaningful evaluation should report both enclosure tightness and computation cost rather than optimizing one number in isolation.

After the midterm stage, possible extensions include explicit treatment of nonzero boundary error, broader forcing functions, additional network families or a different elliptic problem. Each extension requires a new statement of assumptions and a consistent derivation of the residual and error propagation. Applications to FNO or GNO outputs and higher-dimensional PDEs remain exploratory and will depend on the strength and cost of the one-dimensional baseline.

### 7.1 Proposed milestones after the midterm

First, check the derivation against the code, review boundary and derivative assumptions, and retain complete input–configuration–result packages for representative cases. Compare methods under predefined tolerances and stopping rules, reporting successful, exceeded and unresolved outcomes.

Next, investigate nonzero boundary contributions and broader smooth forcing functions, with corresponding derivations and checks. Extensions to neural operators, other equations or higher dimensions remain optional. The thesis will report certificate conditions, tested cases, computational cost and limitations together.

## 8. References

[1] M. Raissi, P. Perdikaris and G. E. Karniadakis, “Physics-informed neural networks: A deep learning framework for solving forward and inverse problems involving nonlinear partial differential equations,” J. Comput. Phys., vol. 378, pp. 686–707, 2019. doi:10.1016/j.jcp.2018.10.045.

[2] F. Eiras, A. Bibi, R. R. Bunel, K. D. Dvijotham, P. Torr and M. P. Kumar, “Efficient Error Certification for Physics-Informed Neural Networks,” Proceedings of Machine Learning Research, vol. 235, pp. 12318–12347, 2024. https://proceedings.mlr.press/v235/eiras24a.html.

[3] B. Hillebrecht and B. Unger, “Rigorous a Posteriori Error Bounds for PDE-Defined PINNs,” IEEE Transactions on Neural Networks and Learning Systems, vol. 36, no. 1, pp. 1583–1593, 2025. doi:10.1109/TNNLS.2023.3335837.

[4] S. Berrone, C. Canuto and M. Pintore, “Solving PDEs by Variational Physics-Informed Neural Networks: an a posteriori error analysis,” Annali dell’Università di Ferrara, vol. 68, pp. 575–595, 2022. doi:10.1007/s11565-022-00441-6.

[5] S. Liu, X. Huang and P. Protopapas, “Residual-based error bound for physics-informed neural networks,” in Proc. UAI, PMLR, vol. 216, pp. 1284–1293, 2023. https://proceedings.mlr.press/v216/liu23b.html.

[6] Y. Yu, C. H. Ho and Y. Wang, “A conformal prediction framework for uncertainty quantification in physics-informed neural networks,” Journal of Computational Physics, vol. 561, article 114979, 2026. doi:10.1016/j.jcp.2026.114979.

[7] Neural PDE Certification project repository, research reports, configurations and selected experiment records. https://github.com/wenmj23/neural-pde-certification. Accessed November 2026.

## 9. AI-Assisted Tools Declaration

Generative AI tools (ChatGPT) were used to assist with this project and with language polishing, structural organization, and typesetting during the preparation of this report. The author remains responsible for the final technical content, the verification and appropriate use of cited sources, and the submitted version.
