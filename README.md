selection_stats
===============

Tests of neutrality based on the Site Frequency Spectrum for VCF files.

This Java program implements many known (cross-population) tests of natural selection, based on allele freuencies and the site frequency spectrum (SFS).

*Overview*
The program takes VCF files (target and control populations, and optionally an outgroup population) and computes statistic values in a sliding window (default 50kb) starting at the first variant, and in each consecutive 2kb interval until the last variant. The program runs on an individual chromosome (i.e. should run 22 times for the complete human autosomal genome).

To run, use:
$> java -jar s_stats.jar -case case_chr1.vcf -control cont_chr1.vcf -outg outg_chr1.vcf

For help, use:
$> java -jar s_stats.jar -h

*Output, by column:*
1. Chromosome
2. Window start
3. Window end
4. Num. sites in case population
5. Num. non-fixed sites in case population (used in all but Fst/PBS)
6. Num. sites in control population
7. Num. non-fixed sites in control population (used in all but Fst/PBS)
8. Num. sites in out-group
9. Num. non-fixed sites in out-group (used in all but Fst/PBS)
10. Watterson's theta, case population
11. Watterson's theta, control population
12. Tajima's theta, case population
13. Tajima's theta, control population
14. Negative log case over control Tajima's theta
15. Tajima's D, case population
16. Tajima's D, control population
17. Sum frequencies theta, case population
18. Sum frequencies theta, control population
19. Negative log case over control Sum Frequencies theta
20. Fst, case vs. control
21. Fst, control vs. out-group
22. Fst, case vs. out-group
23. PBS
