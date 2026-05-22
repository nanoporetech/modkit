# Direct RNA Isoform DMR details

The RNA-specific DMR functions in Modkit use a Dirichlet-multinomial model.
The model assumes that the underlying methylation-state proportions are drawn from a Dirichlet distribution:
\\[
  p_i \sim \text{Dirichlet}(\alpha_{1}, \alpha_{2}, \dots, \alpha_{K})
\\]

Where \\( K \\) is the number of methylation states (m6A and Inosine, for example). 
The `score` is then calculated as the likelihood ratio statistic:
\\[
  \text{LRT} = 2\ (\text{log}\ L_{alt} - L_{null})
\\]

Where the null model is one where all isoforms (or both conditions) share the same proportions and the alternative model is one where each isoform or condition has it's own proportions. 

The degrees of freedom is then calculated as
\\[
  df = (M - 1)(K - 1)
\\]
where \\( M \\) is the number of isoforms (or \\( 2 \\) in the case of `compare-tx-sites`).
Finally, the p-value is calculated as
\\[
  p = P(\chi^2_{df}\\ge\ LRT)
\\]


