/*
**
**
**
**      evolve.c
**
**
**
**
**      Peter Turney
**
**      February 3, 2000
**
**
**
**      Copyright (c) 2000
**
**      National Research Council of Canada
**
**
**
**
**      experimental evolutionary system
**
**      experiments reported in:
**
**      A Simple Model of Unbounded Evolutionary Versatility
**      as a Largest-Scale Trend in Organismal Evolution
**
**
**
**
**
**      - genome:
**
**          [ 10 bits encode mutation rate ] [ 0 to N bits encode phenotype ]
**
**      - initially 0 phenotype bits
**
**      - phenotype bit is added / deleted with probability equal to
**        mutation rate
**
**      - fitness = number of matches to random target
**
**      - if target is static, mutation rate will eventually approach 0,
**        since mutation is more likely to disrupt matches than to find
**        new matches
**
**      - if target is dynamic, mutation rate will be non-zero
**
**      - as genome length increases, absolute number of bit changes per
**        era increases (one indication of increasing evolvability)
**
**      - if 10% target change per era, fitness drops 10% but recovers;
**        absolute magnitude of recovery increases (10% of N as N grows)
**
**
*/

/*
**      include memory management
*/

#include <stdlib.h>
#include <malloc.h>

/*
**      for output
*/

#include <stdio.h>

/*
**      include math library
*/

#include <math.h>

/*
**      include random number generator
**
**      - generates a random number between 1 and 2147483647 
**        (2**31 - 1)
**
**      - built-in C function rand only generates between 1 
**        and 32767 (2**15 - 1)
*/

#include "longrand.h"

/*
**      macros for random numbers
**
**      - fracrand() returns random real number between 0 and 1
**
**      - bitgenc() selects either '0' or '1' with equal probability
**
**      - randomain(hi, lo) selects an integer between lo and hi (inclusive)
*/

#define    fracrand()        ((double) longrand() / 2147483647L)
#define    bitgenc()         ((char) ((fracrand() > 0.5) ? '1' : '0'))
#define    randomain(hi, lo) ((int) floor(fracrand() * ((hi - lo) + 0.999999)) + lo)

/*
**      mutation macro
*/

#define    mutate(x)         ((char) ((x == '1') ? '0' : '1'))

/*
**      main procedure
*/

int main(int argc, char *argv[])
{
  /*
  **    local variable declarations
  */

  FILE *config_file;
  int pop_size, num_children, seed;
  int len_mut_code, report_interval;
  int tournament_size, era_length;
  double era_change_rate;
  int check;
  char **population;
  int i, j;
  char *individual;
  char *target, *new_target;
  int *pop_lengths;
  int *pop_fitnesses;
  int *pop_ages;
  int child_num;
  int max_phenotype_length;
  int mum_pos, dad_pos, min_len, cross_point;
  int child_length;
  int child_fitness;
  char *mum, *dad, *child;
  int max_mut_code, mut_code;
  double mutation_rate;
  int least_fit_fitness, least_fit_pos, least_fit_age;
  int first_best_pos, first_best_fit;
  int second_best_pos, second_best_fit;
  int random_pos, random_fit;
  int fitness, length;
  double average_fitness, average_length;
  double average_age, average_mutation_rate;
  double total_start_era_fitness, average_start_era_fitness;
  double total_end_era_fitness, average_end_era_fitness;
  int num_start_era_fitness, num_end_era_fitness;
  double start_era_fitness, end_era_fitness;
  double fitness_incr_start_end_era;
  double fitness_decr_end_start_era;
  double total_fit_incr_start_end;
  double total_fit_decr_end_start;
  double avg_fit_incr_start_end;
  double avg_fit_decr_end_start;
  int num_fit_incr_start_end;
  int num_fit_decr_end_start;

  /*
  **    read the configuration file
  **
  **    - the name of the configuration file is passed as a command
  **      line argument
  **
  **    - the format of the file is:
  **
  **      2000        Population Size
  **      1000000     Number of Children
  **      34524       Random Number Seed
  **      10          Bits for Encoding Mutation Rate
  **      2000        Reporting Interval
  **      400         Tournament Size
  **      500         Length of an Era
  **      0.2         Rate of Change for a New Era
  **
  **    - the order is important; the numbers are read and the strings
  **      are ignored
  */

  if (argc != 2) {
    printf("\n\n");
    printf("Usage:  > %s my_config.dat \n\n", argv[0]);
    return 0;
  }
  config_file = fopen(argv[1], "r");
  if (config_file == NULL) {
    printf("\n\n");
    printf("Error: Could not open file %s. \n\n", argv[1]);
    return 0;
  }

  /* Population Size */

  check = fscanf(config_file, "%d%*[^\n\r]%*[\n\r]", &pop_size);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Population Size. \n\n");
    return 0;
  }
  printf("\n\n");
  printf("Population Size:                    %d \n", pop_size);

  /* Number of Children */

  check = fscanf(config_file, "%d%*[^\n\r]%*[\n\r]", &num_children);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Number of Children. \n\n");
    return 0;
  }
  printf("Number of Children:                 %d \n", num_children);

  /* Random Number Seed */

  check = fscanf(config_file, "%d%*[^\n\r]%*[\n\r]", &seed);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Random Number Seed. \n\n");
    return 0;
  }
  printf("Random Number Seed:                 %d \n", seed);

  /* Bits for Encoding Mutation Rate */

  check = fscanf(config_file, "%d%*[^\n\r]%*[\n\r]", &len_mut_code);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Bits for Encoding Mutation Rate. \n\n");
    return 0;
  }
  printf("Bits for Encoding Mutation Rate:    %d \n", len_mut_code);

  /* Reporting Interval */

  check = fscanf(config_file, "%d%*[^\n\r]%*[\n\r]", &report_interval);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Reporting Interval. \n\n");
    return 0;
  }
  printf("Reporting Interval:                 %d \n", report_interval);

  /* Tournament Size */

  check = fscanf(config_file, "%d%*[^\n\r]%*[\n\r]", &tournament_size);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Tournament Size. \n\n");
    return 0;
  }
  printf("Tournament Size:                    %d \n", tournament_size);

  /* Length of an Era */

  check = fscanf(config_file, "%d%*[^\n\r]%*[\n\r]", &era_length);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Length of an Era. \n\n");
    return 0;
  }
  printf("Length of an Era:                   %d \n", era_length);

  /* Rate of Change for a New Era */

  check = fscanf(config_file, "%lf%*[^\n\r]%*[\n\r]", &era_change_rate);
  if (check != 1) {
    printf("\n\n");
    printf("Error: Error reading Rate of Change for a New Era. \n\n");
    return 0;
  }
  printf("Rate of Change for a New Era:       %f \n", era_change_rate);
  printf("\n\n");
  fclose(config_file);

  /*
  **    the mutation rate varies from 0 to 1
  **
  **    the mutation rate is encoded using len_mut_code bits
  **
  **    the bits encode a number from 0 to 2**len_mut_code - 1
  */

  max_mut_code = (int) pow(2, len_mut_code) - 1;

  /*
  **    initialize the random number generator
  */

  slongrand(seed);

  /*
  **    allocate memory for the initial population
  */

  population = (char **) malloc(pop_size * sizeof(char *));
  for (i = 0; i < pop_size; i++) {
    population[i] = (char *) malloc(len_mut_code * sizeof(char));
  }

  /*
  **    initialize the population
  */

  for (i = 0; i < pop_size; i++) {
    individual = population[i];
    for (j = 0; j < len_mut_code; j++) {
      individual[j] = bitgenc();
    }
  }

  /*
  **    initialize the lengths
  */

  pop_lengths = (int *) malloc(pop_size * sizeof(int));
  for (i = 0; i < pop_size; i++) {
    pop_lengths[i] = len_mut_code;
  }

  /*
  **    initialize the fitnesses
  */

  pop_fitnesses = (int *) malloc(pop_size * sizeof(int));
  for (i = 0; i < pop_size; i++) {
    pop_fitnesses[i] = 0;
  }

  /*
  **    initialize the ages
  */

  pop_ages = (int *) malloc(pop_size * sizeof(int));
  for (i = 0; i < pop_size; i++) {
    pop_ages[i] = 0;
  }

  /*
  **    the length of the longest phenotype so far
  */

  max_phenotype_length = len_mut_code;

  /*
  **    initialize target
  */

  target = NULL;

  /*
  **    the fitness of the least fit member of the population
  */

  least_fit_fitness = 0;

  /*
  **    initialize variables for statistics for fitness at start of era
  **
  **    - include the initial population
  */

  total_start_era_fitness    = 0;
  average_start_era_fitness  = 0;
  num_start_era_fitness      = 1;
  start_era_fitness          = 0;

  /*
  **    initialize variables for statistics for fitness at end of era
  */

  total_end_era_fitness      = 0;
  average_end_era_fitness    = 0;
  num_end_era_fitness        = 0;
  end_era_fitness            = 0;

  /* 
  **    initialize variables for fitness increment from start to end
  */

  fitness_incr_start_end_era = 0;
  total_fit_incr_start_end   = 0;
  avg_fit_incr_start_end     = 0;
  num_fit_incr_start_end     = 0;

  /*
  **    initialize variables for fitness decrement from end to start
  */

  fitness_decr_end_start_era = 0;
  total_fit_decr_end_start   = 0;
  avg_fit_decr_end_start     = 0;
  num_fit_decr_end_start     = 0;

  /*
  **    MAIN LOOP
  **
  **    create and evaluate each child
  **
  **    - this is a steady-state GA, not a generational GA
  */

  for (child_num = 1; child_num <= num_children; child_num++) {

    /*
    **    randomly select two parents
    **
    **    - use tournament selection
    **
    **    - randomly select tournament_size individuals (sampling
    **      with replacement) and select the best two individuals
    **      as parents
    */

    first_best_pos  = -1;
    second_best_pos = -1;
    first_best_fit  = -1;
    second_best_fit = -1;
    for (i = 0; i < tournament_size; i++) {
      random_pos = randomain((pop_size - 1), 0);
      random_fit = pop_fitnesses[random_pos];
      if (random_fit > first_best_fit) {
        first_best_pos = random_pos;
        first_best_fit = random_fit;
      }
      else if ((random_fit > second_best_fit) && (random_pos != first_best_pos)) {
        second_best_pos = random_pos;
        second_best_fit = random_fit;
      }
    }
    if (fracrand() < 0.5) {
      mum_pos = first_best_pos;
      dad_pos = second_best_pos;
    }
    else {
      mum_pos = second_best_pos;
      dad_pos = first_best_pos;
    }

    /*
    **    pick a cross-over point
    */

    if (pop_lengths[mum_pos] > pop_lengths[dad_pos]) min_len = pop_lengths[dad_pos];
    else min_len = pop_lengths[mum_pos];
    cross_point = randomain((min_len - 2), 1);

    /*
    **    cross-over
    **
    **    make the child
    **
    **    - take the left half of the mum's genome and the right
    **      half of the dad's genome
    **
    **    - thus the length of the child's genome is determined by
    **      the length of the dad's genome
    **
    **    - when allocating memory, add one bit to the child's length, 
    **      in case we later decide to add to the child's genome
    */

    mum = population[mum_pos];
    dad = population[dad_pos];
    child_length = pop_lengths[dad_pos];
    child = (char *) malloc((child_length + 1) * sizeof(char));
    for (i = 0; i < child_length; i++) {
      if (i <= cross_point) child[i] = mum[i];
      else child[i] = dad[i];
    }

    /*
    **    decode the child's mutation rate
    */

    mut_code = 0;
    for (i = 0; i < len_mut_code; i++) {
      mut_code *= 2;
      if (child[i] == '1') mut_code++;
    }
    mutation_rate = (double) mut_code / max_mut_code;

    /*
    **    apply mutation to genome
    */

    for (i = 0; i < child_length; i++) {
      if (fracrand() < mutation_rate) child[i] = mutate(child[i]);
    }

    /*
    **    decide whether to shrink or grow the genome by one bit
    */

    if (fracrand() < mutation_rate) {
      if (fracrand() < 0.5) {

        /*
        **    shrink, unless length is already minimal
        */

        if (child_length > len_mut_code) child_length--;
      }
      else {

        /*
        **    grow, and randomly generate new bit
        */

        child[child_length] = bitgenc();
        child_length++;
      }
    }

    /*
    **    evaluate the fitness of the child by the number
    **    of matches with the target
    */

    if (child_length == len_mut_code) {

      /*
      **    if the part of the child's genome that encodes
      **    the phenotype is empty, then the fitness (the
      **    number of matches with the target) is zero
      */

      child_fitness = 0;
    }
    else if (child_length > max_phenotype_length) {

      /*
      **    if this is the longest child we have seen so far,
      **    then we need to increase the target length
      */

      new_target = (char *) malloc((child_length - len_mut_code) * sizeof(char));
      for (i = 0; i < (max_phenotype_length - len_mut_code); i++) {
        new_target[i] = target[i];
      }
      for (i = (max_phenotype_length - len_mut_code); 
           i < (child_length - len_mut_code); i++) {
        new_target[i] = bitgenc();
      }
      free(target);
      target = new_target;
      max_phenotype_length = child_length;

      /*
      **    now that we have increased the target, we can
      **    evaluate the child
      */

      child_fitness = 0;
      for (i = len_mut_code; i < child_length; i++) {
        if (child[i] == target[(i - len_mut_code)]) child_fitness++;
      }
    }
    else {

      /*
      **    the target does not need to be increased
      */

      child_fitness = 0;
      for (i = len_mut_code; i < child_length; i++) {
        if (child[i] == target[(i - len_mut_code)]) child_fitness++;
      }
    }

    /*
    **    if the child is fitter than the least fit member of
    **    the population, then replace the least fit member of
    **    the population with the child
    **
    **    - if there is a tie for least fit, then replace the
    **      oldest of the least fit
    */

    if (child_fitness >= least_fit_fitness) {

      /*
      **    find the position of the oldest, least fit member of
      **    the population
      */

      least_fit_pos = -1;
      least_fit_age = -1;
      for (i = 0; i < pop_size; i++) {
        if (pop_fitnesses[i] == least_fit_fitness) {
          if (least_fit_pos == -1) {
            least_fit_pos = i;
            least_fit_age = pop_ages[i];
          }
          else if (pop_ages[i] < least_fit_age) {
            least_fit_pos = i;
            least_fit_age = pop_ages[i];
          }
        }
      }

      /*
      **    replace the oldest, least fit member
      */

      free(population[least_fit_pos]);
      population[least_fit_pos]     = child;
      pop_fitnesses[least_fit_pos]  = child_fitness;
      pop_ages[least_fit_pos]       = child_num;
      pop_lengths[least_fit_pos]    = child_length;

      /*
      **    update least_fit_fitness
      */

      least_fit_fitness = child_fitness;
      for (i = 0; i < pop_size; i++) {
        if (pop_fitnesses[i] < least_fit_fitness) {
          least_fit_fitness = pop_fitnesses[i];
        }
      }
    }

    /*
    **    otherwise, free the memory used by the child
    */

    else {
      free(child);
    }

    /*
    **    if it is the start of a new era, then change the
    **    target
    */

    if ((child_num % era_length) == 0) {

      /*
      **    calculate fitness statistics for end of era
      */

      average_fitness = 0;
      for (i = 0; i < pop_size; i++) {
        average_fitness += pop_fitnesses[i];
      }
      average_fitness /= pop_size;
      end_era_fitness = average_fitness;
      total_end_era_fitness += end_era_fitness;
      num_end_era_fitness++;
      average_end_era_fitness = total_end_era_fitness / num_end_era_fitness;

      /*
      **    calculate fitness increase
      */

      fitness_incr_start_end_era = end_era_fitness - start_era_fitness;
      total_fit_incr_start_end += fitness_incr_start_end_era;
      num_fit_incr_start_end++;
      avg_fit_incr_start_end = total_fit_incr_start_end / num_fit_incr_start_end;

      /*
      **    change the target
      */

      for (i = 0; i < (max_phenotype_length - len_mut_code); i++) {
        if (fracrand() < era_change_rate) target[i] = mutate(target[i]);
      }

      /*
      **    re-evaluate the fitness of the population
      */

      fitness = 0;
      for (i = 0; i < pop_size; i++) {
        individual = population[i];
        fitness    = 0;
        length     = pop_lengths[i];
        for (j = len_mut_code; j < length; j++) {
          if (individual[j] == target[(j - len_mut_code)]) fitness++;
        }
        pop_fitnesses[i] = fitness;
      }

      /*
      **    update least_fit_fitness
      */

      least_fit_fitness = fitness;
      for (i = 0; i < pop_size; i++) {
        if (pop_fitnesses[i] < least_fit_fitness) {
          least_fit_fitness = pop_fitnesses[i];
        }
      }

      /*
      **    calculate fitness statistics for start of era
      */

      average_fitness = 0;
      for (i = 0; i < pop_size; i++) {
        average_fitness += pop_fitnesses[i];
      }
      average_fitness /= pop_size;
      start_era_fitness = average_fitness;
      total_start_era_fitness += start_era_fitness;
      num_start_era_fitness++;
      average_start_era_fitness = total_start_era_fitness / num_start_era_fitness;

      /*
      **    calculate fitness decrease
      */

      fitness_decr_end_start_era = end_era_fitness - start_era_fitness;
      total_fit_decr_end_start += fitness_decr_end_start_era;
      num_fit_decr_end_start++;
      avg_fit_decr_end_start = total_fit_decr_end_start / num_fit_decr_end_start;
    }

    /*
    **    show progress
    */

    if ((child_num == 1) || ((child_num % report_interval) == 0)) {

      /*
      **    calculate population statistics
      */

      average_fitness        = 0;
      average_length         = 0;
      average_age            = 0;
      average_mutation_rate  = 0;
      for (i = 0; i < pop_size; i++) {
        average_fitness += pop_fitnesses[i];
        average_length  += pop_lengths[i];
        average_age     += pop_ages[i];
        individual = population[i];
        mut_code = 0;
        for (j = 0; j < len_mut_code; j++) {
          mut_code *= 2;
          if (individual[j] == '1') mut_code++;
        }
        mutation_rate = (double) mut_code / max_mut_code;
        average_mutation_rate += mutation_rate;
      }
      average_fitness        /= pop_size;
      average_length         /= pop_size;
      average_age            /= pop_size;
      average_mutation_rate  /= pop_size;

      /*
      **    print stats
      */

      printf("Child Number:                %d \n", child_num);
      printf("Average Fitness:             %f \n", average_fitness);
      printf("Average Genome Length:       %f \n", average_length);
      printf("Average Mutation Rate:       %f \n", average_mutation_rate);
      printf("Average Birth Date:          %f \n", average_age);
      printf("Average Start Era Fitness:   %f \n", average_start_era_fitness);
      printf("Average End Era Fitness:     %f \n", average_end_era_fitness);
      printf("Average Incr. Start-End:     %f \n", avg_fit_incr_start_end);
      printf("Average Decr. End-Start:     %f \n", avg_fit_decr_end_start);

      /*
      **    print separator
      */

      printf("\n");
      
      /*
      **    halt if average mutation rate is zero
      */

      if (average_mutation_rate == 0) {
        printf("\n\n");
        printf("Average Mutation Rate is Zero. \n\n");
        printf("HALTED \n\n");
        return 0;
      }

      /*
      **    reset start/end era stats
      */

      total_start_era_fitness    = 0;
      average_start_era_fitness  = 0;
      num_start_era_fitness      = 0;

      total_end_era_fitness      = 0;
      average_end_era_fitness    = 0;
      num_end_era_fitness        = 0;

      total_fit_incr_start_end   = 0;
      avg_fit_incr_start_end     = 0;
      num_fit_incr_start_end     = 0;

      total_fit_decr_end_start   = 0;
      avg_fit_decr_end_start     = 0;
      num_fit_decr_end_start     = 0;
    }
  }

  /*
  **    end of MAIN LOOP
  */

  /*
  **    free memory
  */

  for (i = 0; i < pop_size; i++) {
    free(population[i]);
  }
  free(population);
  free(pop_fitnesses);
  free(pop_lengths);
  free(pop_ages);
  free(target);

  /*
  **    end of main
  */

  return 0;
}



