This is a proof-of-concept for the correctness of the KLPT2 algorithm.

Description of files:
 - G2KLPT.sage contains most of the functions required to run the KLPT2 algorithm. This is meant to act as a library for future computations.
 - cost.py and klpt_panny.py are files that have been pulled from https://github.com/friends-of-quaternions/deuring/tree/main. Functions in these files are called in this proof-of-concept.
 - genEG.sage generates examples to be run and is mainly used for debugging.
 - example.sage contains an example that can be ran to demonstrate the correctness of the KLPT2 algorithm.

There are options to run 3 examples in example.sage FastExample(), SlowExample(), NewExample():
 - FastExample() begins with a reduced matrix as defined in Definition 3.11 in the paper. Solutions obtained from previous computations have been hardcoded into the example to demonstrate the correctness of the algorithm.
 - SlowExample() begins with a reduced matrix as defined in Definition 3.11 in the paper (exactly the same as FastExample() so far), but instead of hardcoding the intermediate values as in FastExample(), SlowExample() will re-compute these values. Due to the randomised nature of the intermediate steps, errors and incomplete computations are to be expected. Apologies about that!
 - NewExample() aims to generate a totally new example of a reduced matrix and show that KLPT2 works. Again, errors and incomplete computations are to be expected, and apologies!

To use the functions in this repo, note the key functions required for the KLPT2 algorithm:
 - RandomPolarisation( O, sbound=20 ) which produces a random polarisation
 - Compute_ac_LLL( O, g ) which computes values a and c of the transformation matrix such that the resulting polarisation matrix has a prime top left entry. This is the first step to get a reduced matrix.
 - Compute_bd_KLPT( O, a, c, L=2 ) which computes values b and d of the transformation matrix such that the transformation matrix has determinant a power of L. This is used in conjunction with Compute_ac_LLL to get a transformation matrix that will transform a given polarisation matrix into a reduced matrix.
 - Compute_ac( O, g, L=2 ) which computes values a and c of the transformation matrix such that the top left entry is a power of L. Note that the input g should be a reduced matrix as defined in Definition 3.11. Then using the Compute_bd_KLPT function, this returns a transformation matrix which is close to the final form.
 - FindAlpha( g, O ) which computes the alpha which is used in Section 3.3.
 - ChoosePolarisationPrimePower( g, O, L=2 ) which combines Compute_ac_LLL, Compute_bd_KLPT, Compute_ac, Compute_bd_KLPT to return a transformation matrix which is close to the final form.
 - ChoosePolarisationPrimePower_Reduced( g, O, L=2 ) which combines Compute_ac, Compute_bd_KLPT to return a transformation matrix which is close to the final form. But it accepts as input a g which is a reduced matrix.
