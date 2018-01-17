#include "furry/common/rand.h"
#include <cassert>

//using std::assert;
using std::rand;

namespace furry
{

/* 

   deal: Choose a random k-subset of {0,1,...n-1}

   Author: Frederic Devernay, from public domain Fortran code.

   Synopsis:
   deal(int n, int k, int *a)

   Description:
   Chose a random subset of k elements from {0,1,...n-1}, and return these indices in
   the k first elements of the a[] array.

   This function uses two different algorithms, depending on the respective values of k and n.

   If k > n/2, use a simple algorithm which picks up each sample with a probability k/n.
   If k <= n/2, use the ranksb algorithm.

   ranksb is a translation of the Fortran subroutine RANKSB which appears
   on pp 38-9, Chapter 4 of Nijenhuis & Wilk (1975) Combinatorial Algorithms, Academic Press.

   Compilation/test:
   gcc -DMAIN deal.c -o deal
   ./deal

   The original public domain Fortran source code is available at:
   http://www.cs.sunysb.edu/~algorith/implement/wilf/distrib/processed/
   The Fortran source code was translated using f2c, and then hand-tuned for C.

   Some comments by Gregory P. Jaxon, from
   https://groups.google.com/group/alt.math.recreational/msg/26f549d66c5de8c8

*/

/* There are four different versions of ranksb in this file, which should all give the same result */
static void ranksb_(int n, int k, int *a);     /* f2c-translated version */
static void ranksb1(int n, int k, int *a);     /* hand-tuned f2c version, first stage */
static void ranksb2(int n, int k, int *a);     /* hand-tuned f2c version, second stage */
static void ranksb3(int n, int k, int *a);     /* translated from Gregory P. Jaxon's version */

static void deal_k_near_n(int n, int k, int *a); /* adapted from Gregory P. Jaxon's GRADEAL_K_NEAR_N */

#define ranksb ranksb2 /* ranksb2 should be the best choice */

/* PRINTA macros are used for debugging and comparing the different versions */
#define PRINTA_(I) (void)0
//#define PRINTA_(I) printf("a[%d]=%d\n",I,a[I])
//#define PRINTA_(I) printf("a[%d]=%d, line %d\n",I,a[I], __LINE__)

#define PRINTA(I) (void)0
//#define PRINTA(I) printf("a[%d]=%d\n",I,a[I]+1)
//#define PRINTA(I) printf("a[%d]=%d, line %d\n",I,a[I], __LINE__)

/* return a random number between 0 and n (0 and n included), with a uniform distribution
   reference: http://www.bourguet.org/v2/clang/random/ */
int
alea(int n)
{ 
    /* The following version is often used:
       return (n+1) * rand() / (RAND_MAX + 1.0);
       But it has several issues:
       * risk of overflow if RAND_MAX == MAX_INT (but using a floating-point
         addition fixes this)
       * part of the random bits are lost if RAND_MAX is bigger than the
         biggest integer that can be represented in a double (this is not
         really a problem with 64-bits double, but could be one with 32-bits
         float or 16-bits half)
       * but most importantly, the distribution is biased if RAND_MAX+1 is
         not a multiple of n+1 (the favored values are distributed regularly
         within the interval)
    */
    int partSize, maxUsefull, draw;
    if (n==0)
        return 0;
    assert (0 < n && n <= RAND_MAX);
    partSize =
        n == RAND_MAX ? 1 : 1 + (RAND_MAX-n)/(n+1);
    maxUsefull = partSize * n + (partSize-1);
    //draw;
    do {
      draw = std::rand();
    } while (draw > maxUsefull);
    return draw/partSize;
}
 
/* f2c-translated version (reference) */
void ranksb_(int n, int k, int *a)
{
    /* Local variables */
    int c__, i__, l, m, p, r__, s, x, m0, ds;

    /* Function Body */
    c__ = k;
    for (i__ = 1; i__ <= k; ++i__) {
	a[i__-1] = (i__ - 1) * n / k;
        PRINTA_(i__-1);
    }
L10:
    x = alea(n-1) + 1; //(int)(n * ((double)rand()/(RAND_MAX+1))) + 1;
    l = (x * k - 1) / n + 1;
    if (x <= a[l-1]) {
	goto L10;
    }
    ++a[l-1];
    PRINTA_(l-1);
    --c__;
    if (c__ != 0) {
	goto L10;
    }
    p = 0;
    s = k;
    for (i__ = 1; i__ <= k; ++i__) {
	m = a[i__-1];
	a[i__-1] = 0;
        PRINTA_(i__-1);
	if (m == (i__ - 1) * n / k) {
	    goto L20;
	}
	++p;
	a[p-1] = m;
        PRINTA_(p-1);
L20:
	;
    }
L30:
    l = (a[p-1] * k - 1) / n + 1;
    ds = a[p-1] - (l - 1) * n / k;
    a[p-1] = 0;
    PRINTA_(p-1);
    a[s-1] = l;
    PRINTA_(s-1);
    s -= ds;
    --p;
    if (p > 0) {
	goto L30;
    }
    l = k;
L40:
    if (a[l-1] == 0) {
	goto L50;
    }
    r__ = l;
    m0 = (a[l-1] - 1) * n / k + 1;
    m = a[l-1] * n / k - m0 + 1;
L50:
    x = m0 + alea(m -1);//m0 + (int)(m * ((double)rand()/(RAND_MAX+1)));
    i__ = l;
L60:
    ++i__;
    if (i__ <= r__) {
	goto L80;
    }
L70:
    a[i__ - 1-1] = x;
    PRINTA_(i__-1-1);
    --m;
    --l;
    if (l == 0) {
        /* decrement all elements of a to get C-style indexes */
        for(i__=0; i__<k; i__++)
            --a[i__];
	return;
    }
    goto L40;
L80:
    if (x < a[i__-1]) {
	goto L70;
    }
    ++x;
    a[i__ - 1 - 1] = a[i__-1];
    PRINTA_(i__-1-1);
    goto L60;
} /* ranksb_ */

/* hand-tuned f2c-translated version */
void ranksb1(int n, int k, int *a)
{
    /* Local variables */
    int c, i, l, m, p, r__, s, x, m0, ds;

    /* Function Body */

    /* Partition [0 : n-1] into k intervals:
       Store the least element of the I'th interval in a[i] */
    for (i = 0; i < k; ++i) {
	a[i] = i * n / k;
        PRINTA_(i);
    }

    /* Using a uniformly distributed random variable in the
       range 0 <= x < n, make k selections of x such that
       x lands in the remaining portion of some interval.
       At each successful selection, reduce the remaining
       portion of the interval by 1. */
    /* If k is close to n (say, bigger than n/2), the
       while loop may take many iterations. For this reason,
       it is better to use another algorithm in these
       situations (see deal_k_near_n() below). */
    for(c = k; c > 0; c--) {
        do {
            x = alea(n-1) + 1; //n * ((double)rand()/(RAND_MAX+1)) + 1;
            l = (x * k - 1) / n + 1;
        } while (x <= a[l - 1]);
        ++a[l - 1];
        PRINTA_(l-1);
    }

    /* Collect the least elements of any interval which
       absorbed a selection in the previous step into the
       low-order indices of a. */
    p = -1;
    for (i = 0; i < k; ++i) {
	m = a[i];
	a[i] = 0;
        PRINTA_(i);
	if (m != i * n / k) {
            /* A non-empty partition */
            ++p;
            a[p] = m;
            PRINTA_(p);
        }
    }
    /* Allocate space for each non-empty partition starting
       from the high-order indices.  At the last position
       in each partition's segment, store the interval index
       of the partitions's successor. */
    s = k;
    for(; p >=0; p--) {
        l = (a[p] * k - 1) / n + 1;
        ds = a[p] - (l - 1) * n / k;
        a[p] = 0;
        PRINTA_(p);
        a[s - 1] = l;
        PRINTA_(s-1);
        s -= ds;
    }

    for(l = k; l > 0; l--) {
        /* ranksb each of the sub-problems */
        x = a[l - 1];
        if (x != 0) {
            r__ = l;
            m0 = (x - 1) * n / k + 1;
            m = x * n / k - m0 + 1;
            /* The order of arithmetic operations is important!
               The same rounding errors must be produced in each
               computation of a boundary. */
        }

        /* base_l is the least element of the current (l'th)
           interval.  size_l is the count of the number of
           unselected members of the interval. */
        x = m0 + alea(m-1);//m0 + (int)(m * ((double)rand()/(RAND_MAX+1)));
        i = l;
        ++i;
        while (i <= r__ && x >= a[i - 1]) {
            a[i - 1 - 1] = a[i - 1];
            PRINTA_(i-1-1);
            ++x;
            ++i;
        }
        a[i - 1 - 1] = x;
        PRINTA_(i-1-1);
        --m;
    }
    /* decrement all elements of a to get C-style indexes */
    for(i=0; i<k; i++)
        --a[i];
} /* ranksb_ */

/* hand-tuned f2c-translated version, second (best) version */
/*  k << n, order(k) space & time */
void ranksb2(int n_, int k_, int *a)
{
    /* Local variables */
    int64_t i, l, m, p, r, s, x, m0, ds;
    int64_t n = n_, k = k_;

    /* Function Body */
    if (k == 0)
        return;
    if (k == 1) {
        a[0] = alea(n-1);//n * ((double)rand() / (RAND_MAX+1));
        return;
    }


    /* Partition [0 : n-1] into k intervals:
       Store the least element of the I'th interval in a[i] */
    for (i = 0; i < k; ++i) {
	a[i] = i * n / k;
        PRINTA_(i);
    }

    /* Using a uniformly distributed random variable in the
       range 0 <= x < n, make k selections of x such that
       x lands in the remaining portion of some interval.
       At each successful selection, reduce the remaining
       portion of the interval by 1. */
    /* If k is close to n (say, bigger than n/2), the
       while loop may take many iterations. For this reason,
       it is better to use another algorithm in these
       situations (see deal_k_near_n() below). */
    for(i = 0; i < k;  ++i) {
        do {
            x = alea(n-1) + 1;//n * ((double)rand()/(RAND_MAX+1)) + 1;
            l = (x * k - 1) / n;
        } while (x <= a[l]);
        ++a[l];
        PRINTA_(l);
    }

    /* Collect the least elements of any interval which
       absorbed a selection in the previous step into the
       low-order indices of a. */
    p = -1;
    for (i = 0; i < k; ++i) {
	m = a[i];
	a[i] = 0;
        PRINTA_(i);
	if (m != i * n / k) {
            /* A non-empty partition */
            ++p;
            a[p] = m;
            PRINTA_(p);
        }
    }
    /* Allocate space for each non-empty partition starting
       from the high-order indices.  At the last position
       in each partition's segment, store the interval index
       of the partitions's successor. */
    s = k-1;
    for(; p >=0; p--) {
        l = (a[p] * k - 1) / n;
        ds = a[p] - l * n / k;
        a[p] = 0;
        PRINTA_(p);
        a[s] = l + 1;
        PRINTA_(s);
        s -= ds;
    }

    for(l = k-1; l >= 0; l--) {
        /* ranksb each of the sub-problems */
        x = a[l];
        if (x != 0) {
            /* Start a new bin */
            r = l;
            m0 = (x - 1) * n / k;
            m = x * n / k - m0;
            /* The order of arithmetic operations is important!
               The same rounding errors must be produced in each
               computation of a boundary. */
        }

        /* m0 is the least element of the current (l'th)
           interval.  m is the count of the number of
           unselected members of the interval. */
        x = m0 + alea(m-1);//m0 + (int)(m * ((double)rand()/(RAND_MAX+1)));

        /* Bubble Merge the (x-base_l)'th unselected member
           of the current interval into the current interval's
           segment (a [l..r]). */

        i = l;
        while (i < r && x >= a[i+1]) {
            a[i] = a[i+1];
            PRINTA(i);
            ++x;
            ++i;
        }
        a[i] = x;
        PRINTA(i);
        --m;
    }
} /* ranksb_ */

/* This version is translated from GRADEAL, from
   https://groups.google.com/group/alt.math.recreational/msg/26f549d66c5de8c8
   a few arithmetic opreations had to be changed to be conformant with the original Fortran version
*/
void ranksb3(int n, int k, int *a)
{
    /*  k << n, order(k) space & time */
    if (k == 0)
        return;
    if (k == 1) {
        a[0] = alea(n-1);//n * ((double)rand() / (RAND_MAX+1));
        return;
    }

    /* a[<GRADEUP>a<-k?n] - Works only on the low order
       k elements of a, everything else is left untouched. */
    int x, l, r, i;
    int base_l, size_l;

    /* Partition [0 : n-1] into k intervals:
       Store the least element of the I'th interval in a[i] */
    for(i = 0; i<k; i++) {
        a[i] = i*n/k;
        PRINTA_(i);
    }

    /* Using a uniformly distributed random variable in the
       range 0 <= x < n, make k selections of x such that
       x lands in the remaining portion of some interval.
       At each successful selection, reduce the remaining
       portion of the interval by 1. */
    /* If k is close to n (say, bigger than n/2), the
       while loop may take many iterations. For this reason,
       it is better to use another algorithm in these
       situations (see deal_k_near_n() below). */
    for(i = k-1; i >= 0; i--) {
        do {
            x = alea(n-1) + 1;//(int)(n * ((double)rand()/(RAND_MAX+1))) + 1; // CHANGED from (int)(n * ((double)rand()/(RAND_MAX+1)))
            l = (x * k - 1 )/ n; // CHANGED from x * k / n
        } while(a[l] >= x);
        a[l]++;
        PRINTA_(l);
    }

    /* Collect the least elements of any interval which
       absorbed a selection in the previous step into the
       low-order indices of a. */
    l = -1;
    for(i=0; i<k; i++) {
        x = a[i];
        a[i] = 0;
        PRINTA_(i);
        if( x != i*n/k ) {
            /* A non-empty partition */
            l++;
            a[l] = x;
            PRINTA_(l);
        }
    }

    /* Allocate space for each non-empty partition starting
       from the high-order indices.  At the last position
       in each partition's segment, store the interval index
       of the partitions's successor. */
    r = k-1;
    for(i=l; i>= 0; i--) {
        l= (a[i]*k-1)/n; // CHANGED from a[i]*k/n
        size_l = a[i] - l*n/k;
        a[i] = 0;
        PRINTA_(i);
        a[r] = l+1;
        PRINTA_(r);
        r -= size_l;
    }

    for(l=k-1; l>=0; l--) {
        /* ranksb each of the sub-problems */
        x = a[l];
        if( x != 0 ) {
            /* Start a new bin */
            r = l;
            base_l = (x-1)*n/k;
            size_l = x*n/k - base_l;
            /* The order of arithmetic operations is important!
               The same rounding errors must be produced in each
               computation of a boundary. */
        }
        /* base_l is the least element of the current (l'th)
           interval.  size_l is the count of the number of
           unselected members of the interval. */

        x = base_l + alea(size_l-1);//base_l + (int)(size_l * ((double)rand()/(RAND_MAX+1)));

        /* Bubble Merge the (x-base_l)'th unselected member
           of the current interval into the current interval's
           segment (a [l..r]). */

        i = l;
        while( i < r && x >= a[i+1] ) {
            a[i] = a[i+1];
            PRINTA(i);
            x++;
            i++;
        }
        a[i] = x;
        PRINTA(i);
        size_l--;
    } /* for(l... */
}

/* k near n, o(n) time, o(k) space */
void deal_k_near_n(int n, int k, int *a)
{
    /* Warning: modifies k and n */
    /* Algorithm: go though all candidates from n-1 to 0, and pick each one with probability k/n */
    while((n > k) && (k > n/2)) {
        /* each number has probability k/n of being picked up */
        if (k > alea(n-1)) { //n*((double)rand()/(RAND_MAX+1))) {
            /* pick this one up */
            k--;
            n--;
            a[k]= n;
        }
        else {
            /* don't pick this one */
            n--;
        }
    }
    if (n == k) {
        /* we've got k numbers to sample from a set of k, easy... */
        for(n=n-1; n>=0; n--) {
            a[n] = n;
        }
        k = 0;
    }
    if (k > 0) {
        assert(k <= n/2);
        ranksb(n, k, a);                /* reduced to ranksb */
    }
}

/* k near n, o(n) time, o(k) space */
void deal(int n, int k, int *a)
{
    assert(k <= n);
    if (k <= n/2)
        ranksb(n, k, a);
    else
        deal_k_near_n(n, k, a);
}

std::vector<int> deal(int n, int m)
{
  std::vector<int> samples(m);
  deal(n, m, &samples[0]);
  return samples;
}

void SampleUnique(int n, int m, std::vector<int> *samples) {
  samples->resize(m);
  SampleUnique(n, m, samples->data());
}

void SampleUnique(int n, int m, int *samples) {
  deal(n, m, samples);
}

} // furry
