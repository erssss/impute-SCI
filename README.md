
# SCI

This repository is for paper "Win-Win: On Simultaneous Clustering and Imputing over Incomplete Data".

## File Structure

* code: source code of algorithms.
* data: dataset source files of all eight public data collections used in experiments.

## Dataset

* Banknote: https://archive.ics.uci.edu/dataset/267/banknote+authentication
* LED: https://sci2s.ugr.es/keel/dataset.php?cod=63
* Ecoli: https://archive.ics.uci.edu/dataset/39/ecoli
* Crx: https://sci2s.ugr.es/keel/dataset.php?cod=59
* Dermatology: https://sci2s.ugr.es/keel/dataset.php?cod=60
* Horse: https://sci2s.ugr.es/keel/dataset.php?cod=180
* Soybean: https://archive.ics.uci.edu/dataset/90/soybean+large
* Solar Flare: https://archive.ics.uci.edu/dataset/89/solar+flare

## Dependencies

python 3.9

```
gurobipy==10.0.3
missingpy==0.2.0
numpy==1.23.5
pandas==2.1.3
scikit_learn==0.24.0
scipy==1.11.4
```

## Instruction

``` sh
cd code
python main.py
```