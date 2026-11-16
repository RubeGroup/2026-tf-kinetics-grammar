# 2026-tf-kinetics-grammar

Gillespie (stochastic) simulations of a continuous-time Markov chain (CTMC) model of
transcription factor (TF) binding to DNA. This code generates panels b, c, d and f of
Figure 2 in a manuscript by S. Emami, G. Ben Zeev and H. T. Rube (2026, submitted).

## Contents

| Path | Description |
|---|---|
| `Markov_Chains/SingleStrandChain.py` | Builds the rate matrix `Q` of the single-stranded model (SSM) for a DNA of length `L`, plus maps between state names and indices. |
| `Markov_Chains/MarkovChain.py` | Minimal base class. |
| `Notebooks/Figure Elements.ipynb` | Numba-compiled Gillespie simulator, trajectory analysis (coarse-graining, holding times, off-binding sampling), and the cells that draw Figure 2. |

## Installation

Tested with Python 3.14 on macOS (Apple silicon). The figures were originally produced with
Python 3.11, and they reproduce identically.

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt jupyter
```

## Reproducing Figure 2

Launch Jupyter from the `Notebooks/` directory, so that the notebook can import
`Markov_Chains`, and run all cells in order:

```bash
cd Notebooks
jupyter notebook "Figure Elements.ipynb"
```

| Notebook section | Figure panel |
|---|---|
| "Figure 2b" | 2b: example trajectory over a 20-time-unit window |
| "Figure 2c" | 2c: holding-time distributions from a 10^5-time-unit trajectory |
| "Figure 2d" | 2d: bound tracks sampled at off-binding switch points |
| "Figure 2f" | 2f: off-binding duration histogram with length-biased theory overlay |

Model parameters (`L = 7`, `q01 = 0.1`, `q10 = 1`, `q11 = 2`, `q21 = 0.1`, `q12 = 2` at the fourth
offset (index 3) and `0.01` elsewhere) are set in the "CTMC used to represent the single stranded chain"
cell and match the paper's Materials and Methods. Random seeds are set inside
numba-compiled code (`set_numba_seed`), so results are reproducible. A full run takes
about 15 seconds on a laptop; the first run in a new environment is slower while
packages and numba functions compile.

## License

Released under the MIT License; see [`LICENSE`](LICENSE).

## Citation

If you use this code, please cite the manuscript above and this software; citation
metadata is in [`CITATION.cff`](CITATION.cff).
