#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
//#include"gammln.c"
//#include"bnldev.c"
//#include"ran1.c"
//#include"gamdev.c"
#define MAXTYPES 1501   // maximum number of fitness classes to keep track of
#define MAXGENS  30     // maximum number of generations of growth
#define MAXMUTES 10001  // maximum number of entries in the DFE

/* This simulates an experimental mutation accumulation protocol */

/* running from the command line, a single argument (positive integer)
   can be included to set the seed for the random number generator */

/* In this version, the underlying fitness distibution is the sum of
   two exponential functions that decay away from s=0  */

void main(int argc, char** argv)

{

    FILE* fpout, * fps, * fpfit;
    char filename[100];
    long Nvisible;   // minimum number of indls in a visible colony
    int i, igen, ngens, iline, nlines, itransfer, ntransfers, extinctionflag = 0;
    int initialanc, threshold, muthreshold, nNtypes, mtwt;
    long long Nend, nlucky, ilucky, Ntmp, new, nmu;
    float survivalprob;
    double Wsave, r, dfe[3][MAXMUTES], fmax, fben, alphab, alphad, delW, death;
    int muput, jsave = 0;
    double mutotal;
    double twoplus2s;
    int nsinglets;   // total number of single-mutation entries in DFE
    /*  Ntypes stores the population size of each fitness class
     The zeroth class has relative fitness in the bin 0 to delW,
     the first class has relative fitness in the bin delW to 2*delW, etc.
     There are two columns in Ntypes.
     In the zeroth column we count members of the colony that differ
     from the founding indl of that colony (immediate "ancestor")
     by exactly one mutation (or do not differ at all).
     In the first column we keep any members of the colony that
     differ from the founding indl of that colony by 2 or more mutations.
     nextNtypes is a working copy of Ntypes that we'll need later.  */
    long** Ntypes = (long**)malloc(2 * sizeof(long*));
    for (i = 0;i <= 1;i++) Ntypes[i] = (long*)malloc(MAXTYPES * sizeof(long));
    long** nextNtypes = (long**)malloc(2 * sizeof(long*));
    for (i = 0;i <= 1;i++) nextNtypes[i] = (long*)malloc(MAXTYPES * sizeof(long));

    int ianc;  // fitness class of the ancestor of the currently growing colony
    long N[MAXGENS];  // total population size of the colony
    int j, iold, imu, indmu;
    float ran1(long*), gamdev(int, long*), bnldev(float, int, long*);
    long seed = -30, seedsave;
    /* We will also sample the underlying DFE nsdist times and save the
       the frequency of each fitness class in sdist.  When the number of
       mutations is larger than muthreshold, we will use this distribution
       to distribute them among classes in the next generation, rather than
       choosing each mutation's fitness value stochastically  */
    double* sdist = (double*)malloc(MAXTYPES * sizeof(double));

    // set the random number generator seed according to the command line
    // argument, or else it stays at -1
    if (argc > 1) seed = -((long)(atof(argv[1])));
    seedsave = seed;
    death = 0.2;    // death probability of wildtype individuals
    Nvisible = 1e5;  // colonies chosen below this threshold trigger a warning
    threshold = 500;  // lineage above this size grows deterministically
    muthreshold = 1000;  // more than 1000 mutations from one lineage, then distribute s values deterministically
    ngens = 20;   // generations of growth = number of doublings btwn transfers
    fprintf(stdout, "%d generations between transfers\n", ngens);

    nlines = 1000;  // number of replicate lines to simulate
    ntransfers = 30;   // transfers per line
    mutotal = 0.00043;   // need genome-wide mutation rate, 4.3e-4 
    
    fben = 0.25;    // fraction beneficial in the underlying DFE
    alphab = 30;    // decay rate of DFE for beneficial mutations
    alphad = 20;    // decay rate of DFE for deleterious mutations
    fmax = (double)(1.0 / (1.0 - death)) - 1.0;  // highest relative fitness in the fitness bins will be 1+fmax
    //  This is different from max value of DFE because beneficial mutations
    //  could possibly accrue over the course of the experiment.
    nNtypes = MAXTYPES - 1;
    delW = (double)((1.0 + fmax) / nNtypes);   // width of fitness bins
    initialanc = (int)(1.0 / delW);    // initial ancestor has relative fitness 1
                                       // so this gives its bin index
    nsinglets = 0;  // total number of entries in dfe (so far)

    // check that we can open the output file properly before doing the work
    sprintf_s(filename,100, "madata/exp_%dtrans_%dlines_%1.2ffben_%2.1fab_%2.1fad_%1.5fmu_%dgens_%ldseed.out", ntransfers, nlines, fben, alphab, alphad, mutotal, ngens, -seed);
    int pout = fopen_s(&fpout, filename, "w");
    if (pout != 0) { fprintf(stdout, "Unable to open output file\n"); exit(1); }
    sprintf_s(filename,100, "madata/exp_%dtrans_%dlines_%1.2ffben_%2.1fab_%2.1fad_%1.5fmu_%dgens_%ldseed_fitness.out", ntransfers, nlines, fben, alphab, alphad, mutotal, ngens, -seed);
    int pfit = fopen_s(&fpfit, filename, "w");
    if (pfit != 0) { fprintf(stdout, "Unable to open output file\n"); exit(1); }

    /* set up sdistn so that when there are many mutants,
       we can distribute them deterministically  */
    double sdistsum = 0, stmp;
    for (i = 0;i < nNtypes;i++) {
        stmp = -1 + delW * i;
        if (stmp <= 0) sdist[i] = (1 - fben) * alphad * exp(alphad * stmp);
        else sdist[i] = fben * alphab * exp(-alphab * stmp);
        sdistsum += sdist[i];
    }
    for (i = 0;i < nNtypes;i++) sdist[i] /= sdistsum;
    // save the underlying DFE to an output file
    int ps = fopen_s(&fps, "madata/sdist.out", "w");
    for (i = 0;i < nNtypes;i++) fprintf(fps, "%f\n", sdist[i]);
    fclose(fps);

    // preliminary set-up is finished.  Main loops begin below.  Backing the
    // indentation out to get as much horizontal space as possible
    // loop for each line
    for (iline = 1;iline <= nlines;iline++) {
        // initially start with one indl with fitness 1
        for (i = 0;i < nNtypes;i++) {
            Ntypes[0][i] = 0; Ntypes[1][i] = 0;
        }
        ianc = initialanc;
        Ntypes[0][ianc] = 1;
        Wsave = ianc * delW;   // save the ancestor fitness here
        // loop to do ngens of growth and then sampling, repeated ntransfers times
        for (itransfer = 1;itransfer <= ntransfers;itransfer++) {
            extinctionflag = 0;
            for (igen = 1;igen <= ngens;igen++) {
                // calculate the current size of the colony
                N[igen] = 0;
                for (i = 0;i < nNtypes;i++) N[igen] += Ntypes[0][i] + Ntypes[1][i];
                // if the total size of current colony becomes zero 
                // this may happen with the birth-death process
                if (N[igen] == 0) {
                    extinctionflag = 1;
                    igen = ngens + 1;    // move to end of growth loop
                }
                /* Poisson reproduction, birth-death process:
                Each fitness class gets copied into the next genn (new = Ntypes).
                We allow mutations to happen to individuals as well.
                The pupolation have a Poisson-generated number of offspring with mean Ntypes*2*(1-death)*(1+s).
                If a fitness class has over 500, it grows deterministically (Ntypes*2*(1-death)*(1+s)).
                Mutation can independently occur with each birth.
                If the number of mutants > 1000, they are added deterministically. */

                if (extinctionflag == 0) {
                    for (iold = 0;iold < nNtypes;iold++) nextNtypes[0][iold] = nextNtypes[1][iold] = 0;
                    for (iold = 0;iold < nNtypes;iold++) {
                        for (j = 0;j <= 1;j++) {
                            if (Ntypes[j][iold] > 0) {
                                if (Ntypes[j][iold] > threshold) {
                                    new = 2 * Ntypes[j][iold] * delW * iold * (1 - death);
                                    if (new < 0) {  // this happens if we exceed MAX_LONG
                                        fprintf(stderr, "neg new %d %d %d\n\n", iline, itransfer, igen);
                                        exit(1);
                                    }
                                }
                                else {
                                    // twoplus2s = 2*(1+s) = 2+2*s
                                    twoplus2s = 2.0 + 2.0 * (delW * iold - 1.0);
                            
                                    new = (long)poidev((float)(Ntypes[j][iold] * twoplus2s * (1 - death)), &seed);
                                }

                                // Now that we know how many new individuals to add, let them mutate
                                if (new > 0) nmu = (long)poidev((float)mutotal * (long)new, &seed);
                                else nmu = 0;
                                //  The next 2 lines tell us to put any new indls into column 1
                                // of nextNtypes, unless they are new mutations of the ancestor,
                                // in which case they will be stored in column 0
                                muput = 1;
                                if (iold == ianc) muput = 0;
                                if (nmu > muthreshold) {  //determ mutns
                                    for (i = 0;i < nNtypes;i++) {  // i gives the fitness class of DFE
                                        // newfitness = old fitness * r, which in other words is:
                                        //       newW = (iold*delW)*(delW*i);
                                        // Then newindex = newfitness/delW, which in other words is:
                                            //          indmu = (int)(newW/delW);
                                        // We combine the two steps above as follows:
                                        indmu = (int)(iold * delW * i);
                                        if (indmu >= 0 && indmu < nNtypes) nextNtypes[muput][indmu] += (long)(sdist[i] * nmu);
                                    }
                                }
                                else {
                                    for (imu = 1; imu <= nmu;imu++) {  //stochastic mutn
                                        if (ran1(&seed) < fben)  r = 1.0 - log((double)ran1(&seed)) / alphab;
                                        else r = 1.0 + log((double)ran1(&seed)) / alphad;
                                        //newW = iold*delW*r;  indmu = newW/delW;
                                        indmu = (int)(iold * r);
                                        if (indmu >= 0 && indmu < nNtypes) nextNtypes[muput][indmu] += 1;
                                    }
                                }
                                //for new-nmu births, mutation did not occur
                                nextNtypes[j][iold] += new - nmu;
                            } // if non-zero Ntypes
                        } // loop on i and j, 1 step and multistep mutations
                    } // loop on iold

                    // new becomes the old
                    for (i = 0;i < nNtypes;i++) {
                        Ntypes[0][i] = nextNtypes[0][i];
                        Ntypes[1][i] = nextNtypes[1][i];
                    }
                } // end if extinctionflag == 0
            } // loop on igen, growth is finished

             // If the total size of current colony becomes zero
            // just start again with a single individual at ianc
            if (extinctionflag == 1) {
                Ntypes[0][ianc] = 1; // put back the ancestor and start again
                Wsave = ianc * delW;
                fprintf(fpfit, " %f ", Wsave);

            } // end If loop: when the colony size becomes zero

            // else the population does not go extinct: now we randomly sample a single individual.
           /* This is the individual that lands closest to the marked spot on the next plate.
            We will only simulate the growth process for this colony. */
            else {

                // Nend gives total size of colony at end of growth
                Nend = 0;
                for (i = 0;i < nNtypes;i++) Nend += Ntypes[0][i] + Ntypes[1][i];
                if (ngens > 20 && Nend < Nvisible)
                    fprintf(stderr, "Warning: we chose a colony with only %d indls\n", Nend);
                nlucky = (long)(ran1(&seed) * Nend);
                Ntmp = 0; jsave = 1;
                ilucky = ianc;
                for (i = 0;i < nNtypes;i++) {
                    for (j = 0;j <= 1;j++) {
                        Ntmp += Ntypes[j][i];
                        if (Ntmp > nlucky) {
                            ilucky = i;
                            jsave = j;
                            i = nNtypes + 1; j = 2;  // jump out of both loops
                        }
                    }
                }

                /* add this mutation to the list if previous ancestor
               and this lucky differ by  a single mutation  */
                if ((ilucky != ianc) && (jsave == 0)) {
                    nsinglets++;
                    if (nsinglets > MAXMUTES) fprintf(stderr, "Error: MAXMUTES exceeded\n");
                    dfe[0][nsinglets] = ilucky * delW / Wsave;
                    dfe[1][nsinglets] = (double)itransfer;
                    dfe[2][nsinglets] = ilucky * delW;
                }
                // now lucky becomes the new ancestor for the next growth
                ianc = (int)ilucky;
                for (i = 0;i < nNtypes;i++) {
                    Ntypes[0][i] = Ntypes[1][i] = 0;
                }
                Ntypes[0][ianc] = 1; // this type is the new ancestor
                Wsave = ianc * delW;
                fprintf(fpfit, " %f ", Wsave);
            }// loop on else: when the population does not go extinct

        } // loop on itransfers, do the next growth/sampling now
        fprintf(fpfit, "\n");
        fprintf(stdout, "%d ", iline);

    } // at this point, this line is done, return to do a new line

     // at this point, all the lines are done.  Save results.
    fprintf(stdout, "\n");
    for (i = 1;i <= nsinglets;i++)
        fprintf(fpout, "%f %f %f\n", dfe[0][i], dfe[1][i], dfe[2][i]);
    fclose(fpout);
}
