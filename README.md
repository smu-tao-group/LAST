# LAST: Latent Space Assisted Adaptive Sampling for Protein Trajectories

## Install

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
