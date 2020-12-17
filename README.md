# JUNOasimov

## Components

* `juno_fit.py`: main code for fitting
* `juno_osc.py`: class `Spectra` for calculations done for JUNO analysis (cross-section, reactor flux, oscillating the spectra, smearing with energy resolution, detector response etc.)
* `oscillation.py`: class for `Model` governing neutrino oscillation. Used for general purpose e.g. calculating Majorana mass
* `myplot.py`: framework for plotting. More [here](https://github.com/sagitta42/myplot)


## Example command

```console
python juno_fit.py Nev=10000 fit_method=chi2 data_method=asimov Nbins=200 save
```

## Parameters

* `Nev`: number of events in the sample
* `fit_method`: `chi2` or `lkl` (default `chi2`) (lkl currently removed)
* `data_method`: `sample` (sample events randomly based on PDF), `asimov` (scale PDF) or `gaus` (debug mode) (default `asimov`). Note: for `gaus` mode variables of the function `chisquare()` have to be changed from `chisquare(th12, th13, a, b)` to `chisquare(mu, sig, a, b)` and "JUNO mode" has to be replaced by "debug mode" (uncomment); same applies to `Minuit` initialization
* `Nbins`: number of bins (best value 250 for resolution with a single parameter `a`, 200 when `b` is included)

The additional argument `save` is meant for the `Myplot` class with prompts it to save the image rather than show it

## Free fit variables

The variables of the fit are defined on top of the file

```python
all_params = ['th12', 'th13', 'a', 'b']
free_params = ['th12', 'th13', 'a']
```

## Example output

```console
$ python juno_fit.py Nev=100000 fit_method=chi2 data_method=asimov Nbins=200 save
Oscillation parameters based on NuFIT
---------------
            mean     sigma
th12    0.590270  0.013614
th23    0.848230  0.024435
th13    0.150098  0.002269
dm2sol  0.000074  0.000002
dm2atm  0.002528  0.000031
---------------
### Get unoscillated reactor spectrum
### Get standard oscillated NO spectrum
### Smear it
### Get corresponding PDFs (norm to 1)
Mass ordering: NO
Original parameters:
         th12 : 0.5902703530244823
         a : 0.030151134457776358
         b : 0.05
### Construct NO data histogram for 100000 events based on PDF
Fit method: Chi^2
┌───┬──────┬───────────┬───────────┬────────────┬────────────┬─────────┬─────────┬───────┐
│   │ Name │   Value   │ Hesse Err │ Minos Err- │ Minos Err+ │ Limit-  │ Limit+  │ Fixed │
├───┼──────┼───────────┼───────────┼────────────┼────────────┼─────────┼─────────┼───────┤
│ 0 │ th12 │   0.590   │   0.006   │            │            │         │         │       │
│ 1 │ th13 │  0.1501   │  0.0015   │            │            │         │         │  yes  │
│ 2 │ a    │ 30.15e-3  │  0.30e-3  │            │            │         │         │       │
│ 3 │ b    │  50.0e-3  │  0.5e-3   │            │            │         │         │       │
└───┴──────┴───────────┴───────────┴────────────┴────────────┴─────────┴─────────┴───────┘
Minimizing...
Running time: 27.29336667060852 s
Result: <ValueView of Minuit at 25b6f28>
  th12: 0.5902703530244823
  th13: 0.15009831567151233
  a: 0.030151134457776358
  b: 0.05
Chi2/Ndof: 0.0
### Get unoscillated reactor spectrum
data
oscNO
Image: result_asimov_chi2_NO_Nev100000_smeared_Nbins200_th12-a-b.png
(saved)
-> Contours with 10 bins
('th12', 'a')
Running time: 68.17097401618958 s
Image: contour_asimov_chi2_NO_Nev100000_smeared_Nbins200_th12-a-b_th12VSa.png
(saved)
('th12', 'b')
Running time: 67.53335380554199 s
Image: contour_asimov_chi2_NO_Nev100000_smeared_Nbins200_th12-a-b_th12VSb.png
(saved)
('a', 'b')
Running time: 71.57248640060425 s
Image: contour_asimov_chi2_NO_Nev100000_smeared_Nbins200_th12-a-b_aVSb.png
(saved)
```
