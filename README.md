# Stellar Parameters Determination by Spectra Comparison

## Introduction

This is a simple python module designed to facilitate the comparison of spectra by equivalent width comparison. The design if purposely modular, so the user can edit and hack together exactly what they need for their use case with minimal adaptation.

## Features

- IO functions to easily read and manipulate spectrum files (limit range, normalize, etc).
- Equivalent width calculation method (various methods for estimating continuum).
- Temperature estimation method from growth curves.
- Database lookup and EW comparison.
- Synthetic spectra processing (convolution methods, experimental profiles, rotational profiles).
