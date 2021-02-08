# Demonstrations

REGALS can be applied to SAXS datasets from a variety of experimental setups. Demos 1-4 included here reproduce the analysis in the REGALS paper ([Meisburger, Xu & Ando, 2021]). Each demo is contained in a live notebook and requires an input data file provided in the `data` folder. Both MATLAB and python versions are included. See below for a list of files and reference information for each dataset.

The easiest way to run the demos is to clone (or download) the REGALS repository and run the live notebooks in place. The jupyter notebooks can also be viewed directly within GitHub.

> [Meisburger, Xu, & Ando, 2021]: https://doi.org/10.1107/S2052252521000555
Meisburger, S.P., Xu, D. & Ando, N. (2021). _REGALS_: a general method to deconvolve X-ray scattering data from evolving mixtures. _IUCrJ_ **8**(2). https://doi.org/10.1107/S2052252521000555


## 1. AEX-SAXS

**Experiment:** SAXS with in-line anion exchange chromatography (AEX) collected on as-isolated _B. subtilis_ NrdE, the large subunit of ribonucleotide reductase (from [Parker et al. 2018]).

**Files:**
- `data/NrdE_mix_AEX.mat` - input data (MATLAB-formatted hdf5)
- `NrdE_mix_AEX.ipynb` - jupyter notebook
- `NrdE_mix_AEX.mlx` - MATLAB live notebook

_**Data is included here with permission from the authors. If you use it please cite the original publication:**_
> [Parker et al. 2018]: https://doi.org/10.1073/pnas.1800356115
Parker, M. J. _et al._ (2018). An endogenous dAMP ligand in Bacillus subtilis class Ib RNR promotes assembly of a noncanonical dimer for regulation by dATP. *Proc Natl Acad Sci USA* 115, E4594–E4603. https://doi.org/10.1073/pnas.1800356115

## 2. Ligand titration

**Experiment:** Phenylalanine hydroxylase (PheH) with varying amounts of allosteric ligand, phenylalanine (from [Meisburger et al. 2016]).

**Files:**
- `data/PheH_titration.mat` - input data (MATLAB-formatted hdf5)
- `PheH_titration.ipynb` - jupyter notebook
- `PheH_titration.mlx` - MATLAB live notebook

_**Data is included here with permission from the authors. If you use it please cite the original publication:**_
> [Meisburger et al. 2016]: https://doi.org/10.1021/jacs.6b01563
Meisburger, S. P. _et al._ (2016). Domain Movements upon Activation of Phenylalanine Hydroxylase Characterized by Crystallography and Chromatography-Coupled Small-Angle X-ray Scattering. *J Am Chem Soc* 138, 6506–6516. https://doi.org/10.1021/jacs.6b01563

## 3. Time-resolved mixing

**Experiment:** Stopped-flow mixing of ATP with nucleotide binding domain (NBD) from membrane transporter MsbA (from [Josts et al. 2020]).

**Preprocessing:** Downloaded scattering curves from the SASBDB (ID: [SASDGV5](https://www.sasbdb.org/data/SASDGV5/)) and reformatted.

**Files:**
- `data/MsbA_time_resolved.mat` - preprocessed input data (MATLAB-formatted hdf5)
- `MsbA_time_resolved.ipynb` - jupyter notebook
- `MsbA_time_resolved.mlx` - MATLAB live notebook


_**If you use this data please include the SASBDB ID and cite the original publication:**_
> [Josts et al. 2020]: https://doi.org/10.1016/j.str.2019.12.001
Josts, I. _et al._ (2020). Structural Kinetics of MsbA Investigated by Stopped-Flow Time-Resolved Small-Angle X-Ray Scattering. *Structure* 28, 348-354.e3. https://doi.org/10.1016/j.str.2019.12.001

## 4. Temperature-jump

**Experiment:** SAXS/WAXS measurement following rapid temperature jump (IR laser pulse) to 29.9 degC (from [Thompson et al. 2019]).

**Preprocessing:** downloaded raw data from the NIH Figshare Archive
 ([Fraser, Anfinrud & Thompson, 2019]), calculated difference profiles (laser on - laser off), scaled in WAXS regime, subtracted buffer blanks, and truncated to _q_ < 1 A<sup>-1</sup> (see [Thompson et al. 2019] and [Meisburger, Xu & Ando, 2021] for details).

**Files:**
- `data/CypA_Tjump.mat` - preprocessed input data (MATLAB-formatted hdf5)
- `CypA_Tjump.ipynb` - jupyter notebook
- `CypA_Tjump.mlx` - MATLAB live notebook

_**If you use this data please cite the Figshare entry and the original publication:**_

> [Thompson et al. 2019]: https://doi.org/10.1038/s41557-019-0329-3
Thompson, M. C. _et al._ (2019). Temperature-Jump Solution X-Ray Scattering Reveals Distinct Motions in a Dynamic Enzyme. *Nature Chemistry* 11, 1058–1066. https://doi.org/10.1038/s41557-019-0329-3

>[Fraser, Anfinrud & Thompson, 2019]: https://doi.org/10.35092/yhjc.9177143.v1
Fraser, J. Anfinrud, P. & Thompson, M. (2019). X-ray scattering curves (SAXS/WAXS) used for the analysis described in _Temperature-Jump Solution X-ray Scattering Reveals Distinct Motions in a Dynamic Enzyme_. The NIH Figshare Archive. Dataset. https://doi.org/10.35092/yhjc.9177143.v1
