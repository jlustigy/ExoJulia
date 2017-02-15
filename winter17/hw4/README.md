Homework 4:

1. Expreiment with the routine ttv_nplanet.jl, which is a wrapper for compute_ttv.jl which calls N(N-1) pairs of planets to compute the TTV of each planet due to all of the others.
2. Carry out an initial fit to the transit times of the two planets, ttv_planet1.txt & ttv_planet2.txt, including only 2 planets in your model. Assume a 30-second timing precision for each. Note 2: you cannot set e=0 in compute_ttv.jl
3. Add in a third (non-transiting) planet to your model. Make a grid in period, and optimize the fit over periods from 500-10,000 days & over ~10 phases for each period. Plot the maximum likelihood versus period. Note 3: you'll need to use a wrapper for curve_fit which allows the 3rd period to be fixed, ttv_wrapper_fixp3.jl
4. Starting at the best-fit, run a Markov chain. Obtain the masses of the planets, period of the outer planet, eccentricity vectors of the planets, and uncertainties on all.
5. Which planets are these? What are the implications of this result?
