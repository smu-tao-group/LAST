# LAST: Latent Space Assisted Adaptive Sampling for Protein Trajectories
[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/smu-tao-group/ADMET_XGBoost.svg?style=square)](https://lgtm.com/projects/g/HTian1997/getarticle)
[![DOI:10.1021/acs.jcim.2c01213](http://img.shields.io/badge/DOI-10.1021/acs.jcim.2c01213-B31B1B.svg)](https://doi.org/10.1021/acs.jcim.2c01213)


## News

- (2022.12) This paper has been published in [JCIM](https://pubs.acs.org/doi/10.1021/acs.jcim.2c01213).
- (2022.11) This paper has been accepted by [MLCB 2022](https://sites.google.com/cs.washington.edu/mlcb2022/) as Spotlight (top 12%).
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
