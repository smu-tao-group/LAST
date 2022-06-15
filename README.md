# LAST: Latent Space Assisted Adaptive Sampling for Protein Trajectories
[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/smu-tao-group/ADMET_XGBoost.svg?style=square)](https://lgtm.com/projects/g/HTian1997/getarticle)
[![DOI](http://img.shields.io/badge/DOI-arXiv:2204.13040-B31B1B.svg)](https://arxiv.org/abs/2204.13040)


## News

- (2022.06) This paper has been accepted by [ICML 2022 Workshop AI4Science](http://www.ai4science.net/icml22/).
- (2022.04) This work has been awarded the Chemical Computing Group (CCG) Excellence Award for Graduate Students by ACS Computers In Chemistry (COMP) division in ACS Fall 2022.

## Installation

```bash
git clone https://github.com/smu-tao-group/LAST
cd LAST
conda env create -f environment.yml
conda activate last
```

## Usage

1. Prepare the input files (INPCRD and PRMTOP) to `inputs/`.
2. Go to `src/` and run `chmod +x ./run.sh`.
3. Run `./run.sh PDB_ID MAX_ITERATION PATIENCE`.

## License

GPL-3.0
