# Examples of how thresholds affect base modification calls

The following examples are meant to demonstrate how individual base modification calls will be made during
`pileup` or `call-mods`. Important to remember that the probability of the modification needs to be >= the
pass threshold value.

## Two-way base modification calls
```text
probability of canonical cytosine: p_C
probability of 5mC: p_m
threshold for all cytosine base modifications: threshold_C 

{p_C: 0.25, p_m: 0.75}, threshold_C: 0.7 => call: 5mC (m), most likely.
{p_C: 0.25, p_m: 0.75}, threshold_C: 0.8 => call: Fail, both probabilities are below threshold.
```

## Three-way base modification calls
```text
probability of canonical cytosine: p_C
probability of 5mC: p_m
probability of 5hmC: p_h
threshold for all cytosine base modifications: threshold_C 

{p_C: 0.05, p_m: 0.7, p_h: 0.25}, threshold_C: 0.7 => call: 5mC (m), most likely.
{p_C: 0.15, p_m: 0.6, p_h: 0.25}, threshold_C: 0.7 => Fail, all below threshold.
```

## Three-way base modification calls with modification-specific thresholds
```text
probability of canonical cytosine: p_C
probability of 5mC: p_m
probability of 5hmC: p_h
threshold for all canonical cytosine: mod_threshold_C 
threshold for all 5mC: mod_threshold_m 
threshold for all 5hmC: mod_threshold_h 

command line: --filter-threshold C:0.7 --mod-threshold m:0.8 --mod-threshold h:0.9

filter_threshold C: 0.7
mod_threshold_m: 0.8
mod_threshold_h: 0.9
{p_C: 0.05, p_m: 0.85, p_h: 0.1} => call: 5mC (m) most likely

filter_threshold C: 0.7
mod_threshold_m: 0.8
mod_threshold_h: 0.9
{p_C: 0.75, p_m: 0.05, p_h: 0.2} => call: C (canonical) most likely
{p_C: 0.05, p_m: 0.75, p_h: 0.2} => Fail, all below respective thresholds
{p_C: 0.1, p_m: 0.05, p_h: 0.85} => Fail, all below respective thresholds
```



