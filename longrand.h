/*
**  longrand() -- generate 2**31-2 random numbers
**
**  public domain by Ray Gardner
** 
**  based on "Random Number Generators: Good Ones Are Hard to Find",
**  S.K. Park and K.W. Miller, Communications of the ACM 31:10 (Oct 1988),
**  and "Two Fast Implementations of the 'Minimal Standard' Random
**  Number Generator", David G. Carta, Comm. ACM 33, 1 (Jan 1990), p. 87-88
**
**  linear congruential generator f(z) = 16807 z mod (2 ** 31 - 1)
**
**  uses L. Schrage's method to avoid overflow problems
*/

/*
**  used internally in longrand()
*/

long nextlongrand(long seed);

/* 
**  generate a random number between 1 and 2147483647 (2**31 - 1)
*/

long longrand(void);

/* 
**  seed the random number generator
*/

void slongrand(unsigned long seed);




