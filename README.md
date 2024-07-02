# Inferring Microbial Interactions from Time Series Data
## TODO
  - data preprocessing
  - build true $ \rho(D) $

Let a $2\times n$ matrix $M$ represent the relative abundance timeseries of two taxa $i : \{i_{1}..i_{n}\}$ and $j : \{j_{1}..j_{n}\}$  where the rows represent the taxa and the indices columns represent the chronologically ordered set of timeponts $T$ $\{1, 2, 3, … n\}$.A permutation is a both a one-to-one and onto function $\sigma: \{1, 2, 3, … n\} → \{1, 2, 3, … n\}$ that reorders the indices of $M$.
$$
\begin{align*}
\sigma(1) &= 3 \\
\sigma(2) &= 2 \\
\sigma(3) &= 4 \\
\sigma(4) &= 1
\end{align*}
$$
The empirically observed Spearman coefficient is calculated as
$$\rho(\Delta i_{t},j)$$
where
$$\Delta i_{t} = i_{t+1} - i_{t}$$
The permuted Spearman coefficient is calculated as
$$\rho(\Delta i_{\sigma(t)},j)$$
where
$$\Delta i_{t} = i_{\sigma(t)+1} - i_{\sigma(t)}$$# bangladesh_time_series
