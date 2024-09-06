[:arrow_left: Back to portfolio](https://github.com/adrien-b-git/Portofolio)

# Geometrical Modelling

The purpose of the project is to create an easy to use interface allowing a user to draw curves from a set of points.

To compute these curves, I implemented several methods to choose from:
- Lagrange polynomials ;
- Béziers curves ;
- Hermite C1 splines;
- Hermite C2 splines;

The user can set the resolution of the curve, the tension, and display the curvature.

# How to install

```console
conda env create -f environment.yml
conda activate curve
python3 ./CurveEditor.py
```

