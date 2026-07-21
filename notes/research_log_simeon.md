# Simeon's Research Log (Teachers Project)
Log initiated on July 15, 2026.
## July 15, 2026
### Manuscript
Reviewed version 1.00 of manuscript. This is a major revision.

Here is a list of necessary updates (non-exhaustive):
1. Add $L$ locations, indexed by $\ell = 1,\ldots,L$. Index variables with $\ell$, where necessary.
2. Change notation: the teacher's wage function needs to be defined more generally since her compensation will depend on the ability compensation of the students, which is not identical across locations, in general. The second argument has to be a location-specific index (or moment) of the students' ability (on site)
3. Section 2.2.2: review the discussion of sorting forces since the teacher's preference for higher-ability students only matters when she is paid her marginal product (rather than based on an exogenous wage function with her human capital as the only argument). The second paragraph needs to be reviewed too, to make sure the section is internally coherent.
4. Footnote 6: when the wage profile is exogenous, the occuptional choice boundary **does** depend on $\tilde{H}^T$!
5. Section 3.3, Proposition 2: the choice boundary depends on $\tilde{H}^T$ when the teacher's wage function is exogenous; the proposition is only correct when teachers are paid their marginal product.
6. Section 4.1, Tables 1 & 2:
    + Report single decimal point for all percentages
    + There are discrepancies between Table 1 and Table 2: the total number of females is not the same and the percentages are inconsistent with the raw numbers (in Table 1, 47,417 females = 51.0% of the total; in Table 2, 43,859 females = 51.4% of the total)
7. Section 4.3: why do we need the O*Net skill requirements and the decomposition of job-specific abilities into their $CA$ and $AA$ components? Is the idea that we somehow quantify the extent of mismatch in the data and in the model across the three steady-state calibrations? I doesn't look like we actually do anything with these abilities in the calibration in Section 5.
8. Section 5:
    + Do we need to revisit the exclusion of unemployed workers? It's a bit strange since we assign workers with less than 15 hours per week to home production. A more "natural" would be to assign them to home production, no?
    + Normalization of labor market barriers:
        + teaching is a natural choice since we are assuming there is no discrimination between men and women
        + earlier this year, we kicked around some alternative ideas but I don't quite recall what they were (I should check the slides from Ananth's presentation in March or April)
    + Normalization of occupational productivities: home production for men (review what we did for Ananth's spring 2026 presentation; we may have kicked around some additional options too)
    + $\lambda_f$ is introduced in Section 5.1: this is, most likely, a relic from the early version where teachers are paid their marginal product and the the calibration was built around this assumption. It appears out of nowhere and is not defined when it is first used on page 26; it is used again later in Appendix B.
    + Review Figures 7 and 8: I believe we have better plots in the slide deck (spring 2026)
9. New sections:
    + #5.3: untargeted moments, sensitivity checks, Jacobian to illustrate identification
    + #6: counterfactual experiments
    + #7: conclusion (previously #6)
### Previous Literature
Need to re-read a couple of papers to get a better handle on the calibration strategy:
1. Guvenen et al. (2020): "Multidimensional Skill Mismatch"
2. Hsieh et al (2019): "The Allocation of Talent (...)"