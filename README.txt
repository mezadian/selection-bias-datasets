This README describes the data for the manuscript:

"Correcting for selection bias in mutation accumulation experiments: bias is amplified when rare lineages may be lost"

by
Mojgan Ezadian (Westtern University, London, Canada, mezadian@uwo.ca)
Lindi M. Wahl (Western University, London, Canada, lwahl@uwo.ca)

All data are stored in csv format.  Note that there are no new original data associated with this study; these data are provided for convenience but can all be reproduced exactly using the simulation code and mathematical techniques described in the manuscript.
For reference, all filenames give the appropriate figure number and panel from the manuscript.

Columns in each file as described below.  For each data column, we provide the following in order:

Variable name, Type of computation, Type of distribution, number of generations, Mathematical symbol, Units, Allowable values (any relevant parameter values that distinguish that column, in parentheses)
In any field, "--" indicates not applicable.

Fig1A:
column 1: --, --, --, --, x, unitless, [0,1]
column 2: Probability generating function, numerical, fission distribution, 1, f^1(x,1), unitless, [0,1]
column 3: Probability generating function, symbolic,  fission distribution, 1, f^1(x,1), unitless, positive number
column 4: Probability generating function, numerical, fission distribution, 3, f^3(x,1), unitless, positive number
column 5: Probability generating function, symbolic,  fission distribution, 3, f^3(x,1), unitless, positive number
column 6: Probability generating function, numerical, fission distribution, 5, f^5(x,1), unitless, positive number
column 7: Probability generating function, symbolic,  fission distribution, 5, f^5(x,1), unitless, positive number
column 8: Probability generating function, numerical, fission distribution, 7, f^7(x,1), unitless, positive number
column 9: Probability generating function, symbolic,  fission distribution, 7, f^7(x,1), unitless, positive number
Fig1B:
column 1:  --,--, --, --, x, unitless, [0,1]
column 2:  Probability generating function, numerical, Poisson distribution, 1, f^1(x,1), unitless, positive number
column 3:  Probability generating function, symbolic,  Poisson distribution, 1, f^1(x,1), unitless, positive number
column 4:  Probability generating function, numerical, Poisson distribution, 3, f^3(x,1), unitless, positive number
column 5:  Probability generating function, symbolic,  Poisson distribution, 3, f^3(x,1), unitless, positive number
column 6:  Probability generating function, numerical, Poisson distribution, 5, f^5(x,1), unitless, positive number
column 7:  Probability generating function, symbolic,  Poisson distribution, 5, f^5(x,1), unitless, positive number
column 8:  Probability generating function, numerical, Poisson distribution, 7, f^7(x,1), unitless, positive number
column 9:  Probability generating function, symbolic,  Poisson distribution, 7, f^7(x,1), unitless, positive number
column 10: Probability generating function, numerical, Poisson distribution, 15, f^15(x,1), unitless, positive number
column 11: Probability generating function, symbolic, Poisson distribution, 15, f^15(x,1), unitless, positive number

Fig2:
column 1:  selective effect, --, --, --, s, unitless, [-1,0.2]
column 2:  relative bias, symbolic, fission distribution, 3, b(s), unitless, positive number
column 3:  relative bias, numerical FFT, fission distribution, 3, b(s), unitless, positive number
column 4:  relative bias, direct-bias estimation, fission distribution, 3, b(s), unitless, positive number
column 5:  relative bias, symbolic, fission distribution, 4, b(s), unitless, positive number
column 6:  relative bias, numerical FFT, fission distribution, 4, b(s), unitless, positive number
column 7:  relative bias, direct-bias estimation, fission distribution, 4, b(s), unitless, positive number
column 8:  relative bias, symbolic, fission distribution, 5, b(s), unitless, positive number
column 9:  relative bias, numerical FFT, fission distribution, 5, b(s), unitless, positive number
column 10: relative bias, direct-bias estimation, fission distribution, 5, b(s), unitless, positive number
column 11: relative bias, symbolic, fission distribution, 6, b(s), unitless, positive number
column 12: relative bias, numerical FFT, fission distribution, 6, b(s), unitless, positive number
column 13: relative bias, direct-bias estimation, fission distribution, 6, b(s), unitless, positive number

Fig3A:
column 1: selective effect, --, --, --, s, unitless, [-1,0.25]
column 2: relative bias, numerical, pure-birth model, 5,  b(s), unitless, positive number
column 3: relative bias, numerical, pure-birth model, 10, b(s), unitless, positive number
column 4: relative bias, numerical, pure-birth model, 15, b(s), unitless, positive number
column 5: relative bias, numerical, pure-birth model, 20, b(s), unitless, positive number
column 6: relative bias, numerical, pure-birth model, 27, b(s), unitless, positive number

Fig3B:
column 1: selective effect, --, --, --, s, unitless, [-1,0.25]
column 2: relative bias, direct-bias estimation, fission distribution, 5,  b(s), unitless, positive number
column 3: relative bias, direct-bias estimation, fission distribution, 10, b(s), unitless, positive number
column 4: relative bias, direct-bias estimation, fission distribution, 15, b(s), unitless, positive number
column 5: relative bias, direct-bias estimation, fission distribution, 20, b(s), unitless, positive number
column 6: relative bias, direct-bias estimation, fission distribution, 27, b(s), unitless, positive number

Fig3C:
column 1: selective effect, --, --, --, s, unitless, [-1,0.25]
column 2: relative bias, direct-bias estimation, Poisson distribution, 5,  b(s), unitless, positive number
column 3: relative bias, direct-bias estimation, Poisson distribution, 10, b(s), unitless, positive number
column 4: relative bias, direct-bias estimation, Poisson distribution, 15, b(s), unitless, positive number
column 5: relative bias, direct-bias estimation, Poisson distribution, 20, b(s), unitless, positive number

Fig4A:
column 1: selective effect, --, --, --, s, unitless, [-1,0.11]
column 2: relative bias, numerical, pure-birth model for d=0.1, 20,  b(s), unitless, positive number
column 3: relative bias, direct-bias estimation, fission distribution for d=0.1, 20, b(s), unitless, positive number
column 4: relative bias, direct-bias estimation, Poisson distribution for d=0.1, 20, b(s), unitless, positive number
 
Fig4B:
column 1: selective effect, --, --, --, s, unitless, [-1,0.25]
column 2: relative bias, numerical, pure-birth model for d=0.2, 20,  b(s), unitless, positive number
column 3: relative bias, direct-bias estimation, fission distribution for d=0.2, 20, b(s), unitless, positive number
column 4: relative bias, direct-bias estimation, Poisson distribution for d=0.2, 20, b(s), unitless, positive number
 

Fig4C:
column 1: selective effect, --, --, --, s, unitless, [-1,0.42]
column 2: relative bias, numerical, pure-birth model for d=0.3, 20,  b(s), unitless, positive number
column 3: relative bias, direct-bias estimation, fission distribution for d=0.3, 20, b(s), unitless, positive number
column 4: relative bias, direct-bias estimation, Poisson distribution for d=0.3, 20, b(s), unitless, positive number

Fig4D:
column 1: selective effect, --, --, --, s, unitless, [-1,0.66]
column 2: relative bias, numerical, pure-birth model for d=0.4, 20,  b(s), unitless, positive number
column 3: relative bias, direct-bias estimation, fission distribution for d=0.4, 20, b(s), unitless, positive number
column 4: relative bias, direct-bias estimation, Poisson distribution for d=0.4, 20, b(s), unitless, positive number



Fig7:
column 1: generations, --, --, 5 10 15 20 27, tau, unitless, positive integer
column 2: observed beneficial fraction, simulated DFE, fission distribution, 5 10 15 20 27, f_b, unitless, positive integer
column 3: corrected beneficial fraction, simulated DFE, fission distribution, 5 10 15 20 27, f_b, unitless, positive integer
column 4: true beneficial fraction, simulated DFE, fission distribution, 5 10 15 20 27, f_b, unitless, positive integer





