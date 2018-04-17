[![DOI](https://zenodo.org/badge/82066723.svg)](https://zenodo.org/badge/latestdoi/82066723)

An Efficient Two-Level Preconditioner for Multi-Frequency Wave Propagation Problems: Numerical Examples and MAPLE derivations
-----------------------------------------------------------------------------------------------------------------------------

This repository contains additional material to [BvG17]:

* complete the proof of Lemma 4.1 with two Maple scripts (files in subfolder `/maple`),
* visualization of Lemma 4.1: https://manuelmbaumann.github.io/optimal-tau/,
* numerical experiments as presented in Section 6 (subfolder `/num_exper`).

Numerical experiments:
----------------------
Our numerical experiments are done in Python3 and are stored in the subfolder `/num_exper`. The experiments can be reproduced with the scripts `exp1.sh` - `exp5.sh`, respectively. An example run from the command line looks like:

`python3 elast_wedge.py --ndims=2 --freq=[1.0,10.0] --Nom=5 --damping=0.7`

<img src="/num_exper/figs/circ_pic.png" width="425"/> <img src="/num_exper/figs/msconv-plot.png" width="425"/> 
*(extreme frequencies displayed in red)*

Dependencies:
-------------
For the proof of Lemma 4.1:
* Maple [v 18.02]

For the [visualization](https://delft.diskos.nl/bokeh/opt_tau_bokeh):
* [Bokeh](http://bokeh.pydata.org/en/latest/) [v 0.12.4]

For the numerics:
* [nutils](http://www.nutils.org/):  `pip install git+https://github.com/joostvanzwieten/nutils@955bc67d219496e26b037f47855709a222850d7c`
* NumPy [v 1.8.2], SciPy [v 0.14.0], matplotlib [v 1.4.2]

Hardware info:
--------------
All numerical examples presented in Section 6 of [BvG17] have been implemented in Python-3 running over GNU/Debian Linux, and executed on a computer with 4 CPUs Intel I5 with 32 GB of RAM.

References:
-----------
* [BvG17]: M. Baumann and M.B. van Gijzen (2017). [An Efficient Two-Level Preconditioner for Multi-Frequency Wave Propagation Problems.](https://d1rkab7tlqy5f1.cloudfront.net/EWI/Over%20de%20faculteit/Afdelingen/Applied%20Mathematics/Reports/2017/BvG_TR-17-03.pdf) Tech. report, DIAM Report 17-03, Delft University of Technology.
