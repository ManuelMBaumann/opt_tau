An Efficient Two-Level Preconditioner for Multi-Frequency Wave Propagation Problems: Numerical Examples and MAPLE derivations
-----------------------------------------------------------------------------------------------------------------------------

This repository contains additional material to [BvG17]:

* complete the proof of Lemma 4.1 with two Maple scripts (files in subfolder `/maple`),
* visualization of Lemma 4.1: https://delft.diskos.nl/bokeh/opt_tau_bokeh,
* numerical experiments as presented in Section 6 (subfolder `/num_exper`).

Numerical experiments:
----------------------
Our numerical experiments are done in Python3 and are stored in the subfolder `/num_exper`. The experiments can be reproduced with the scripts `exp1.sh` - `exp5.sh`, respectively. An example run from the command line looks like:

`python3 elast_wedge.py --ndims=2 --freq[1.0,9.0] --Nom=5`

<img src="/num_exper/figs/circ_pic.png" width="425"/> <img src="/num_exper/figs/msconv-plot.png" width="425"/> 

Dependencies:
-------------
For the proof of Lemma 4.1:
* Maple [v 18.02]
For the [visualization](https://delft.diskos.nl/bokeh/opt_tau_bokeh):
* [Bokeh](http://bokeh.pydata.org/en/latest/) [v 0.12.4]
For the numerics:
* [nutils](http://www.nutils.org/):  `pip install git+https://github.com/joostvanzwieten/nutils@955bc67d219496e26b037f47855709a222850d7c`
* NumPy [v 1.8.2], SciPy [v 0.14.0], matplotlib [v 1.4.2]

References:
-----------
* [BvG17]: M. Baumann and M.B. Van Gijzen (2017). [An Efficient Two-level Preconditioner for Multi-Frequency Wave Propagation Problems.](http://www.ewi.tudelft.nl/en/the-faculty/departments/applied-mathematics/reports/) Tech. report, DIAM Report 17-03, Delft University of Technology.
