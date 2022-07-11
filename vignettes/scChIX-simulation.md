## Simulating scChIX input data

First we set up some parameters to simulate scChIX data. We will set
half the bins to be “mutually exclusive” while the other half will be
overlapping.

For each of the three cell types, we generate three count matrices.
First represents histone mark 1, second represents histone mark 2, and
the third represents histone mark 1+2 double-incubation.

    #> simATAC is:
    #> ...updating parameters...
    #> ...setting default parameters...
    #> ...setting up SingleCellExperiment object...
    #> Your data has different number of bins compared to the provided genome positions. Please give a file of bin information consistent with your input data with three columns and header of "chr start end" as the bin.coordinate.file parameter. If you don't give a file containing the information of bins, simATAC considers the bin.coordinate.file parameter as "None" and names the bins {Bin1 to BinX} with X number of bins. In this case, you wont be able to get the coordinate information of bins. Please make sure the "species" parameter of the simATACCount object is set correctly.
    #> ...simulating library size...
    #> ...simulating non-zero cell proportion...
    #> ...simulating bin mean...
    #> ...generating final counts...
    #> ...Done...
    #> [1] 0.4140963
    #>             Bin  BinMean Bin.Orig
    #> Bin127  Bin4499 4.069776   Bin127
    #> Bin126  Bin4782 3.867244   Bin126
    #> Bin781  Bin4023 2.236004   Bin781
    #> Bin2382 Bin1391 2.112053  Bin2382
    #> Bin4674 Bin1918 1.952029  Bin4674
    #> Bin1983 Bin1116 1.804295  Bin1983
    #> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    #>   ..@ i       : int [1:1552861] 27 36 38 43 53 60 63 76 83 97 ...
    #>   ..@ p       : int [1:751] 0 2459 5705 6915 8289 9792 11767 13086 14551 15927 ...
    #>   ..@ Dim     : int [1:2] 5000 750
    #>   ..@ Dimnames:List of 2
    #>   .. ..$ : chr [1:5000] "Bin1" "Bin2" "Bin3" "Bin4" ...
    #>   .. ..$ : chr [1:750] "Cell1" "Cell2" "Cell3" "Cell4" ...
    #>   ..@ x       : num [1:1552861] 1 3 1 1 3 1 1 1 1 1 ...
    #>   ..@ factors : list()
    #> NULL
    #> 5 x 5 sparse Matrix of class "dgCMatrix"
    #>         Cell1_A Cell2_A Cell3_A Cell4_A Cell5_A
    #> Bin1422       .       .       .       .       .
    #> Bin1017       .       .       .       .       .
    #> Bin4775       .       1       .       .       .
    #> Bin2177       .       .       .       .       .
    #> Bin1533       .       .       .       .       .
    #> simATAC is:
    #> ...updating parameters...
    #> ...setting default parameters...
    #> ...setting up SingleCellExperiment object...
    #> Your data has different number of bins compared to the provided genome positions. Please give a file of bin information consistent with your input data with three columns and header of "chr start end" as the bin.coordinate.file parameter. If you don't give a file containing the information of bins, simATAC considers the bin.coordinate.file parameter as "None" and names the bins {Bin1 to BinX} with X number of bins. In this case, you wont be able to get the coordinate information of bins. Please make sure the "species" parameter of the simATACCount object is set correctly.
    #> ...simulating library size...
    #> ...simulating non-zero cell proportion...
    #> ...simulating bin mean...
    #> ...generating final counts...
    #> ...Done...

<img src="scChIX-simulation_files/figure-markdown_github/generate_counts-1.png" width="80%" />

    #> [1] 0.4261821
    #>             Bin  BinMean Bin.Orig
    #> Bin127   Bin195 4.069776   Bin127
    #> Bin126  Bin3092 3.867244   Bin126
    #> Bin2382 Bin1210 2.208168  Bin2382
    #> Bin781  Bin3399 2.091723   Bin781
    #> Bin4674 Bin3662 1.816923  Bin4674
    #> Bin2362 Bin3444 1.810604  Bin2362
    #> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    #>   ..@ i       : int [1:1598183] 0 10 23 26 36 56 122 125 126 149 ...
    #>   ..@ p       : int [1:751] 0 1575 2846 5039 7016 7830 9893 12104 14257 15534 ...
    #>   ..@ Dim     : int [1:2] 5000 750
    #>   ..@ Dimnames:List of 2
    #>   .. ..$ : chr [1:5000] "Bin1" "Bin2" "Bin3" "Bin4" ...
    #>   .. ..$ : chr [1:750] "Cell1" "Cell2" "Cell3" "Cell4" ...
    #>   ..@ x       : num [1:1598183] 1 1 1 1 1 1 1 22 15 1 ...
    #>   ..@ factors : list()
    #> NULL
    #> 5 x 5 sparse Matrix of class "dgCMatrix"
    #>         Cell1_B Cell2_B Cell3_B Cell4_B Cell5_B
    #> Bin2463       1       .       .       .       .
    #> Bin2511       .       .       .       .       .
    #> Bin2227       .       .       .       .       .
    #> Bin526        .       .       .       .       .
    #> Bin4291       .       1       .       .       .
    #> simATAC is:
    #> ...updating parameters...
    #> ...setting default parameters...
    #> ...setting up SingleCellExperiment object...
    #> Your data has different number of bins compared to the provided genome positions. Please give a file of bin information consistent with your input data with three columns and header of "chr start end" as the bin.coordinate.file parameter. If you don't give a file containing the information of bins, simATAC considers the bin.coordinate.file parameter as "None" and names the bins {Bin1 to BinX} with X number of bins. In this case, you wont be able to get the coordinate information of bins. Please make sure the "species" parameter of the simATACCount object is set correctly.
    #> ...simulating library size...
    #> ...simulating non-zero cell proportion...
    #> ...simulating bin mean...
    #> ...generating final counts...
    #> ...Done...

<img src="scChIX-simulation_files/figure-markdown_github/generate_counts-2.png" width="80%" /><img src="scChIX-simulation_files/figure-markdown_github/generate_counts-3.png" width="80%" />

    #> [1] 0.4240136
    #>             Bin  BinMean Bin.Orig
    #> Bin127  Bin2559 4.069776   Bin127
    #> Bin126  Bin4901 3.794844   Bin126
    #> Bin2382 Bin2098 2.236004  Bin2382
    #> Bin4674 Bin3961 1.893565  Bin4674
    #> Bin3303 Bin4328 1.887121  Bin3303
    #> Bin781  Bin2801 1.855056   Bin781
    #> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    #>   ..@ i       : int [1:1590051] 9 11 19 23 36 49 83 102 117 120 ...
    #>   ..@ p       : int [1:751] 0 2670 3640 5598 8178 11791 14587 16569 19987 22007 ...
    #>   ..@ Dim     : int [1:2] 5000 750
    #>   ..@ Dimnames:List of 2
    #>   .. ..$ : chr [1:5000] "Bin1" "Bin2" "Bin3" "Bin4" ...
    #>   .. ..$ : chr [1:750] "Cell1" "Cell2" "Cell3" "Cell4" ...
    #>   ..@ x       : num [1:1590051] 1 1 1 2 2 1 1 1 1 1 ...
    #>   ..@ factors : list()
    #> NULL
    #> 5 x 5 sparse Matrix of class "dgCMatrix"
    #>         Cell1_C Cell2_C Cell3_C Cell4_C Cell5_C
    #> Bin923        .       .       .       .       1
    #> Bin2409       .       .       .       .       .
    #> Bin1034       .       .       .       .       .
    #> Bin2339       .       .       1       .       .
    #> Bin3351       .       .       .       .       .

We do a quick check the UMAPs of these count matrices that they are
three different cell types for each histone mark.

    #> 5 x 5 sparse Matrix of class "dgCMatrix"
    #>      Cell1_A Cell2_A Cell3_A Cell4_A Cell5_A
    #> Bin1       .       .       .       .       .
    #> Bin2       3       5       .       1       1
    #> Bin3       1       .       1       .       .
    #> Bin4       1       1       .       .       .
    #> Bin5       5      10       2       .       3
    #> 5 x 5 sparse Matrix of class "dgCMatrix"
    #>      Cell251_A Cell252_A Cell253_A Cell254_A Cell255_A
    #> Bin1         2         .         3         .         2
    #> Bin2         .         .         .         .         1
    #> Bin3         .         1         1         .         .
    #> Bin4         .         .         1         .         .
    #> Bin5         .         .         .         .         .
    #> 5 x 5 sparse Matrix of class "dgCMatrix"
    #>      Cell501_AxCell501_A Cell502_AxCell502_A Cell503_AxCell503_A Cell504_AxCell504_A Cell505_AxCell505_A
    #> Bin1                   2                   6                   .                   3                   3
    #> Bin2                   .                   6                   1                   .                   1
    #> Bin3                   .                   .                   .                   .                   .
    #> Bin4                   .                   .                   .                   .                   .
    #> Bin5                   3                   3                   2                   1                   4
    #> $mark1

<img src="scChIX-simulation_files/figure-markdown_github/check-1.png" width="80%" />

    #> 
    #> $mark2

<img src="scChIX-simulation_files/figure-markdown_github/check-2.png" width="80%" />

    #> 
    #> $`mark1-mark2`

<img src="scChIX-simulation_files/figure-markdown_github/check-3.png" width="80%" />

## Running scChIX on input data

Write the count matrices to a directory (see
`inst/extdata/countmat_var_filt.mark1.rds` for example) and then run the
snakemake workflow `snakemake_workflow/run_snakemake.simulation_data.sh`
changing the paths to the correct directories.

The snakemake workflow takes several hours to complete, so we will just
load and analyze the results in this notebook.

## Downstream analysis of scChIX to check simulated data

Load the scChIX outputs for the simulated data for three overlapping
scenarios: `frac.mutual.excl=0.01, 0.5, 0.99`

We can plot the empirical 95% confidence intervals from our estimates of
the degree of overlaps.

<img src="scChIX-simulation_files/figure-markdown_github/ci-1.png" width="80%" /><img src="scChIX-simulation_files/figure-markdown_github/ci-2.png" width="80%" /><img src="scChIX-simulation_files/figure-markdown_github/ci-3.png" width="80%" /><img src="scChIX-simulation_files/figure-markdown_github/ci-4.png" width="80%" />

We can show that the two UMAPs can be linked together using the
deconvolved double-incubated cells as anchors.

<img src="scChIX-simulation_files/figure-markdown_github/umap-1.png" width="80%" />

Finally we can plot the distribution of overlap estimates across the
bins and compare how these distributions look compared to ground truth.

    #> $`0.01`

<img src="scChIX-simulation_files/figure-markdown_github/hist-1.png" width="80%" />

    #> 
    #> $`0.5`

<img src="scChIX-simulation_files/figure-markdown_github/hist-2.png" width="80%" />

    #> 
    #> $`0.99`

<img src="scChIX-simulation_files/figure-markdown_github/hist-3.png" width="80%" />
<img src="scChIX-simulation_files/figure-markdown_github/hist-4.png" width="80%" />

    #> $`0.01`

<img src="scChIX-simulation_files/figure-markdown_github/hist-5.png" width="80%" />

    #> 
    #> $`0.5`

<img src="scChIX-simulation_files/figure-markdown_github/hist-6.png" width="80%" />

    #> 
    #> $`0.99`

<img src="scChIX-simulation_files/figure-markdown_github/hist-7.png" width="80%" />
<img src="scChIX-simulation_files/figure-markdown_github/hist-8.png" width="80%" />
