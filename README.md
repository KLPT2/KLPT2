There are a few parallel files that you can run. They are very messy at the moment, so apologies!

You can run example.sage, or genEG.sage.
These two files require cost.py and klpt_panny.py in the same folder.

The first file example.sage runs a small example to show that it works. The example begins with a reduced matrix as defined in Definition 3.11 in the paper.

The second file is genEG.sage which is used to generate small examples for you to use.

The final file is G2KLPT.sage which contains all the functions required for the G2KLPT algorithm to complete.
Key functions in this file include:
 - RandomPolarisation( O, sbound=20 ) which produces a random polarisation
 - Compute_ac_LLL( O, g ) which computes values a and c of the transformation matrix such that the resulting polarisation matrix has a prime top left entry. This is the first step to get a reduced matrix.
 - Compute_bd_KLPT( O, a, c, L=2 ) which computes values b and d of the transformation matrix such that the transformation matrix has determinant a power of L. This is used in conjunction with Compute_ac_LLL to get a transformation matrix that will transform a given polarisation matrix into a reduced matrix.
 - Compute_ac( O, g, L=2 ) which computes values a and c of the transformation matrix such that the top left entry is a power of L. Note that the input g should be a reduced matrix as defined in Definition 3.11. Then using the Compute_bd_KLPT function, this returns a transformation matrix which is close to the final form.
 - FindAlpha( g, O ) which computes the alpha which is used in Section 3.3.
 - ChoosePolarisationPrimePower( g, O, L=2 ) which combines Compute_ac_LLL, Compute_bd_KLPT, Compute_ac, Compute_bd_KLPT to return a transformation matrix which is close to the final form.
 - ChoosePolarisationPrimePower_Reduced( g, O, L=2 ) which combines Compute_ac, Compute_bd_KLPT to return a transformation matrix which is close to the final form. But it accepts as input a g which is a reduced matrix.
